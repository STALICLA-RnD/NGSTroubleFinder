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


class PyPaCBAM:
    """
    Attributes:
        minBaseQuality: min quality to consider a base
        minReadQuality: min quality to consider a read
        uncalledThreshold: threshold to consider a variant uncalled
        heteroMinValue: heterozygosity range
        threads: number of threads to use

    """

    def __init__(
        self,
        threads=4,
        minBaseQuality=30,
        minReadQuality=1,
        uncalledThreshold=0.02,
        heteroMinValue=0.2,
    ):
        self.minBaseQuality = minBaseQuality
        self.minReadQuality = minReadQuality
        libfile = pathlib.Path(__file__).parent / "_pyPaCBAM.so"
        self._pileupEngine = ctypes.CDLL(str(libfile))

        self._pileupEngine.openHtsFile.argtypes = (
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.POINTER(NgsFile),
        )

        self._pileupEngine.closeHtsFile.argtypes = (ctypes.POINTER(NgsFile),)

        self._pileupEngine.pileup.argtypes = (
            ctypes.POINTER(NgsFile),  # file
            ctypes.c_char_p,  # chrom
            ctypes.c_int,  # position
            ctypes.c_int,  # quality
            ctypes.c_int,  # read quality
            ctypes.POINTER(ctypes.c_int),  # results
        )

        if uncalledThreshold > heteroMinValue:
            raise ValueError(
                "Uncalled threshold must be lower than the heterozygosity value"
            )

        self.uncalledThreshold = uncalledThreshold
        self.heteroMinValue = heteroMinValue

        self.threads = threads

    def _positionPileup(self, ngsfile, chrom, pos):
        """
        Compute the pileup for a single variant
        Args:
            ngsfile: Pysam.AlignmentFile
            file to pileup
            chrom: str
            chromosome of the location
            pos: int
            position to pileup

        Returns: Dictionary{str:int}
            count of A, C, G, T in each position

        """
        pileup = [0, 0, 0, 0]

        pileupPtr = (ctypes.c_int * 4)(*pileup)
        pileupError = self._pileupEngine.pileup(
            ngsfile,
            ctypes.c_char_p(chrom.encode()),
            pos - 1,
            self.minBaseQuality,
            self.minReadQuality,
            pileupPtr,
        )

        pileup = {
            "A": pileupPtr[0],
            "C": pileupPtr[1],
            "G": pileupPtr[2],
            "T": pileupPtr[3],
        }
        return pileup, pileupError

    def _computeAF(self, pileup, refBase, altBase):
        """
        Compute the allelic fraction of a position based on the pileup
        This definition of allelic fraction considers only the
        bases matching the REF or the ALT alleles. It doesn't consider
        other bases
        E.G. ref = A, alt = G {A:10, C:2, G=10, T=0} AF = 10/(10=10) = 0.5
        Args:
            pileup: Dictionary{str:int}
            pileup dictionary counting the read supporting each bases
            refBase: str
                A, C, G, T based on the variant
            altBase:
                A, C, G, T based on the variant

        Returns: float
            allelic fraction observed for a variant

        """
        if (pileup[refBase] + pileup[altBase]) == 0:
            return 0
        else:
            return pileup[altBase] / (pileup[refBase] + pileup[altBase])

    def _positionGenotype(self, af):
        """
        Genotype a position based on the allelic fraction
        This method leaves a buffer zone of uncalled variant
        to detect contamination
        Args:
            af: float
            observed allelic fraction of a locus

        Returns: str
            "0/0" ref homozygous
            "0/1" heterozygous
            "1/1" alt homozygous
            "./." uncalled

        """
        if af <= self.uncalledThreshold:
            return "0/0"
        elif af >= 1 - self.uncalledThreshold:
            return "1/1"
        elif self.heteroMinValue <= af <= 1 - self.heteroMinValue:
            return "0/1"
        else:
            return "./."

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

    def genotypeChunk(self, fileName, chunk, chunkPosition, q, reference=None):
        """
        wrapper function that genotype a chunk using an independent process
        Args:
            ngsFile: str
            path to the file (every process must have an independent file or there will be a
            race condition on seek)
            chunk: List(str)
                list of vcf lines
            chunkPosition: int
                position of the chunk in the vcf (hack to sort the vcf back to the original order)
                it limits the intra-process comunications by sending only chunks and not lines
                it also improves sorting performance (maybe it's not needed)
            q: multiprocessing.Queue
                queue used to return the results to the main process
            reference: str
                path to the reference genome

        """
        results = []
        # ngsFile = pysam.AlignmentFile(ngsFile, reference_filename=reference)
        ngsFile = NgsFile(None, None, None)

        if reference is not None:
            reference = ctypes.c_char_p(reference.encode())

        openError = self._pileupEngine.openHtsFile(
            ctypes.c_char_p(fileName.encode()),
            ctypes.c_char_p((fileName + ".bai").encode()),
            reference,
            ngsFile,
        )

        if openError:
            q.put(ERROR)
            quit()

        for snp in chunk:
            chrom = snp[0]
            pos = int(snp[1])
            rsid = snp[2]
            ref = snp[3]
            alt = snp[4]
            pileup, pileupError = self._positionPileup(ngsFile, chrom, pos)
            cov = sum([pileup[k] for k in pileup])
            af = self._computeAF(pileup, ref, alt)
            genotype = self._positionGenotype(af)

            if pileupError:
                q.put(ERROR)
                quit()

            results.append(
                (
                    chrom,
                    int(pos),
                    rsid,
                    ref,
                    alt,
                    int(pileup["A"]),
                    int(pileup["C"]),
                    int(pileup["G"]),
                    int(pileup["T"]),
                    float(af),
                    int(cov),
                    genotype,
                )
            )

        self._pileupEngine.closeHtsFile(ngsFile)
        q.put((chunkPosition, results))

    def parseVcf(self, vcf):
        """
        Parse a vcf by keeping only the first 5 column
        Args:
            vcf: str
                path to the vcf

        Returns: generator(List(str))
            chunked generator of lists containing a list of snps to genotype
        """
        with open(vcf) as vcf:
            snps = vcf.readlines()
            snps = [x.strip().split("\t")[0:5] for x in snps if not x.startswith("#")]
            snps = self.chunkise(snps, self.threads)
            return snps

    def genotype(
        self,
        ngsFile,
        vcf=files("ngsTroubleFinder.tools").joinpath("exonic_no_HVR_sites.hg38.vcf"),
        reference=None,
        outFile=None,
    ):
        """
        function to genotype a vcf using multiprocessing
        The vcf get "parsed" (only relevant informations)
        chunked and fed into different processes
        The results can be saved into a file and returned to the caller
        Args:
            ngsFile: str
                path to a bam or cram file
            vcf: str
                path to a vcf (the vcf is already included in the model)
            reference: str
                path to the reference genome only used for crams
            outFile: str
                path to the output pileup file

        Returns: pd.Dataframe
            pandas dataframe in PaCBAM format
            (for each position it has the number of read supporting each base, allelic fraction
            and genotype)
        """
        q = multiprocessing.Queue()

        snps = self.parseVcf(vcf)

        processes = []

        for position, chunk in enumerate(snps):
            p = multiprocessing.Process(
                target=self.genotypeChunk,
                args=(ngsFile, chunk, position, q, reference),
            )
            p.start()
            processes.append(p)

        completePileup = []
        receivedMessages = 0
        errors = 0
        while receivedMessages < self.threads:
            m = q.get()
            if m == ERROR:
                print("Error while processing pileup")
                errors = 1
            else:
                completePileup.append(m)
            receivedMessages += 1

        for p in processes:
            p.join()
            p.close()

        if errors:
            print("Quitting")
            quit()

        completePileup = sorted(completePileup, key=lambda x: x[0])
        completePileup = [x[1] for x in completePileup]
        completePileup = list(itertools.chain(*completePileup))

        haplotypeDataframe = pd.DataFrame(
            completePileup,
            columns=[
                "chr",
                "pos",
                "rsid",
                "ref",
                "alt",
                "A",
                "C",
                "G",
                "T",
                "af",
                "cov",
                "genotype",
            ],
        )

        if outFile is not None:
            haplotypeDataframe.to_csv(outFile, sep="\t", index=False)

        return haplotypeDataframe
