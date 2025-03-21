import ctypes
import itertools
import multiprocessing
import pathlib
from importlib.resources import files

import pandas as pd

ERROR = 1


class NgsFile(ctypes.Structure):
    _fields_ = [
        ("file", ctypes.POINTER(ctypes.c_void_p)),
        ("header", ctypes.POINTER(ctypes.c_void_p)),
        ("index", ctypes.POINTER(ctypes.c_void_p)),
    ]


class HaplotypeDetector:
    """
    This class implements methods to parse reads that span two snps and to detected
    instances where they support more than 2 haplotypes

    Attributes
    ----------
    threads : int
        number of threads to use when computing the haplotypes

    Methods
    ----------
    detectMultipleHaplotypes(pileup,
        alignmentFile,
        groups=files("ngsTroubleFinder.tools").joinpath("exonic_consecutive_snps_no_HVR.hg38.tsv"),
        referenceGenome=None,
        outFile=None,
    )
        Starts the computation to detect multiple haplotypes
    """

    def __init__(self, threads=4, minBaseQuality=30, minReadQuality=40):
        """
        Parameters
        ----------
        threads : int
            number of threads to use in the parallel computation

        """
        self.threads = threads
        self.minBaseQuality = minBaseQuality
        self.minReadQuality = minReadQuality

        self.libfile = pathlib.Path(__file__).parent / "_pyPaCBAM.so"

    def _filterPileupByGroups(self, pileup, groups, minCoverage=20):
        """
        Get usable snps by getting variants that are in a 150 base window,
        have at least a minimum coverage and are on autosomes

        Parameters
        ----------
        pileup : pd.DataFrame
            pileup dataframe in PyPaCBAM format
        groups : pd.DataFrame
            dataframe with [chr, pos, id, group] column where group indicates all the snps in a
            150 window
        minCoverage : int
            minimum coverage to consider a variatn

        Returns
        -------
        pd.DataFrame
            filtered pileup

        """
        pileup = pileup[pileup["rsid"].isin(groups["rsid"])]
        pileup = pileup.join(groups[["rsid", "group"]].set_index("rsid"), on="rsid")

        pileup = pileup[pileup["cov"] >= minCoverage]
        return pileup

    def _analyzeHaplotype(self, haplotypes, minimalSupportReads=3):
        """
        Analyze the summary of the reads to detect the number of haplotypes

        Parameters
        ----------
        haplotypes : dictionary
            dictionary with the counts of the various haplotypes
        minimalSupportReads : int
            minimal number of reads to say that an haplotype is observed

        Returns
        -------
        True if 3 or more haplotype are observed, False otherwise

        """

        foundHaplotypes = 0
        for k in haplotypes:
            if k != "Other" and haplotypes[k] >= minimalSupportReads:
                foundHaplotypes += 1

        if foundHaplotypes >= 3:
            return True
        else:
            return False

    @staticmethod
    def chunkise(a, n):
        """
        chunk a list in n similar sized chunks
        Parameters
        ----------
        a : list to be chunked
        n : number of chunks

        Returns
        -------
        generator
            generator for the chunks
        """
        k, m = divmod(len(a), n)
        return (a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n))

    def _processVariantPair(self, ngsFile, variant1, variant2):
        detectedHaplotypes = [0, 0, 0, 0, 0]

        # prepare pointer to save the results
        detectedHaplotypesPtr = (ctypes.c_int * 5)(*detectedHaplotypes)
        haplotypeError = self._pileupEngine.inferHaplotype(
            ngsFile,
            ctypes.c_char_p(variant1.chr.encode()),
            variant1.pos - 1,
            variant2.pos - 1,
            ctypes.c_char(variant1.ref.encode()),
            ctypes.c_char(variant1.alt.encode()),
            ctypes.c_char(variant2.ref.encode()),
            ctypes.c_char(variant2.alt.encode()),
            self.minBaseQuality,
            self.minReadQuality,
            detectedHaplotypesPtr,
        )

        haplotypes = {
            "RefRef": detectedHaplotypesPtr[0],
            "RefAlt": detectedHaplotypesPtr[1],
            "AltRef": detectedHaplotypesPtr[2],
            "AltAlt": detectedHaplotypesPtr[3],
            "Other": detectedHaplotypesPtr[4],
        }
        return haplotypes, haplotypeError

    def _analyzeChunk(self, fileName, chunk, chunkPosition, q, referenceGenome=None):
        """
        wrapper function that analyze a chunk and return the results to the main thread using
        a queue
        Parameters
        ----------
        ngsFile : str
            path to a bam/cram file
        chunk : list
            group of snp to analyze
        chunkPosition : int
        q : multiprocessing.Queue
            queue used to return the results to the main thread
        referenceGenome : str | None
            path to the reference genome (used only with crams)

        """

        self._pileupEngine = ctypes.CDLL(str(self.libfile))

        self._pileupEngine.openHtsFile.argtypes = (
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.POINTER(NgsFile),
        )

        self._pileupEngine.closeHtsFile.argtypes = (ctypes.POINTER(NgsFile),)

        self._pileupEngine.inferHaplotype.argtypes = (
            ctypes.POINTER(NgsFile),  # file
            ctypes.c_char_p,  # chrom
            ctypes.c_int,  # position
            ctypes.c_int,  # quality
            ctypes.c_char,  # ref base variant 1
            ctypes.c_char,  # alt base variant 1
            ctypes.c_char,  # ref base variant 2
            ctypes.c_char,  # alt base variant 2
            ctypes.c_int,  # baseQuality
            ctypes.c_int,  # readQuality
            ctypes.POINTER(ctypes.c_int),  # results
        )

        analyzedGroups = []

        ngsFile = NgsFile(None, None, None)

        if referenceGenome is not None:
            referenceGenome = ctypes.c_char_p(referenceGenome.encode())

        openError = self._pileupEngine.openHtsFile(
            ctypes.c_char_p(fileName.encode()),
            ctypes.c_char_p(fileName.encode()),
            referenceGenome,
            ngsFile,
        )
        if openError:
            q.put(ERROR)
            quit()

        for g in chunk:
            v1 = g.iloc[0]
            v2 = g.iloc[1]

            haplotypes, haplotypeError = self._processVariantPair(ngsFile, v1, v2)
            if haplotypeError:
                q.put(ERROR)
                quit()
            multipleHaplotypeSupport = self._analyzeHaplotype(haplotypes)

            analyzedGroups.append(
                [
                    v1.chr,
                    v1.pos,
                    v2.pos,
                    v1.ref,
                    v1.alt,
                    v2.ref,
                    v2.alt,
                    haplotypes["RefRef"],
                    haplotypes["RefAlt"],
                    haplotypes["AltRef"],
                    haplotypes["AltAlt"],
                    haplotypes["Other"],
                    multipleHaplotypeSupport,
                ]
            )

        self._pileupEngine.closeHtsFile(ngsFile)
        q.put((chunkPosition, analyzedGroups))

    def detectMultipleHaplotypes(
        self,
        pileup,
        alignmentFile,
        groups=files("ngsTroubleFinder.tools").joinpath(
            "exonic_consecutive_snps_no_HVR.hg38.tsv"
        ),
        referenceGenome=None,
        outFile=None,
    ):
        """
        Main function to start the analysis of the haplotypes
        It starts multiple threads and assigns to each threads a chunk of
        snps, collect the results and returns them to the caller
        Parameters
        ----------
        pileup : pd.DataFrame
            pileup dataframe in PyPaCBAM format
        alignmentFile : str
            bam or cram file path
        groups : str
            path to a tsv with the groups (already provided in the repository)
        referenceGenome : str | None
            path to the reference genome (only used for crams)
        outFile : str | None
            output path for the final statistic

        Returns
        -------
        pd.DataFrame
            dataframe containing the summary data for the multiple haplotypes

        """
        q = multiprocessing.Queue()

        groups = pd.read_csv(groups, sep="\t")
        pileup = self._filterPileupByGroups(pileup, groups)

        multipleSNPs = [x[1] for x in pileup.groupby("group") if len(x[1]) == 2]
        multipleSNPs = self.chunkise(multipleSNPs, self.threads)
        processes = []

        for position, chunk in enumerate(multipleSNPs):
            p = multiprocessing.Process(
                target=self._analyzeChunk,
                args=(alignmentFile, chunk, position, q, referenceGenome),
            )
            p.start()
            processes.append(p)

        analyzedGroups = []
        receivedMessages = 0
        errors = 0
        while receivedMessages < self.threads:
            m = q.get()
            if m == ERROR:
                print("Error while processing haplotypes")
                errors = 1
            else:
                analyzedGroups.append(m)
            receivedMessages += 1

        for p in processes:
            p.join()
            p.close()

        if errors:
            print("Quitting")
            quit()

        analyzedGroups = sorted(analyzedGroups, key=lambda x: x[0])
        analyzedGroups = [x[1] for x in analyzedGroups]
        analyzedGroups = list(itertools.chain(*analyzedGroups))

        haplotypeDataframe = pd.DataFrame(
            analyzedGroups,
            columns=[
                "chr",
                "SNP1",
                "SNP2",
                "SNP1Ref",
                "SNP1Alt",
                "SNP2Ref",
                "SNP2Alt",
                "RefRef",
                "RefAlt",
                "AltRef",
                "AltAlt",
                "Other",
                "multipleHaplotypes",
            ],
        )

        if outFile is not None:
            haplotypeDataframe.to_csv(outFile, sep="\t", index=False)

        return haplotypeDataframe
