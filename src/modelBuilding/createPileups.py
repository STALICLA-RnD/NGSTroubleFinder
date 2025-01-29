import glob
import multiprocessing
import random
import subprocess
import sys
from pathlib import PurePath

import pandas as pd
from ngsTroubleFinder.haplotypeDetector import HaplotypeDetector
from ngsTroubleFinder.pyPaCBAM import PyPaCBAM


def createContamination(contaminant, percentage, seed, basePath):
    contaminantName = PurePath(
        contaminant.rstrip("".join(PurePath(contaminant).suffixes))
    ).stem

    subsample = f"{basePath}/{contaminantName}_{percentage}_{seed}.bam"

    subprocess.run(
        f"samtools view {contaminant} -bs {seed}.{percentage} -o {subsample} -T {verifyBamIDPath}/GRCh38.genome.fa",
        shell=True,
    )

    subprocess.run(
        f"samtools index {subsample}",
        shell=True,
    )

    return subsample


def runVerifyBamID(verifyBamIDPath, bam, outFile, dataset="", pileup=False):
    if dataset == "exome":
        svdprefix = (
            f"{verifyBamIDPath}/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat"
        )
    else:
        svdprefix = f"{verifyBamIDPath}/resource/1000g.phase3.100k.b38.vcf.gz.dat"

    if not pileup:
        subprocess.run(
            f"{verifyBamIDPath}/VerifyBamID \
                --BamFile {bam} \
                --SVDPrefix {svdprefix} \
                --Output {outFile} \
                --Reference {verifyBamIDPath}/GRCh38.genome.fa \
                --NumThread 2 --OutputPileup",
            shell=True,
        )
    else:
        subprocess.run(
            f"{verifyBamIDPath}/VerifyBamID \
                --PileupFile {bam} \
                --SVDPrefix {svdprefix} \
                --Output {outFile} \
                --Reference {verifyBamIDPath}/GRCh38.genome.fa \
                --NumThread 2 --OutputPileup",
            shell=True,
        )


def cleanup(f):
    subprocess.run(
        f"rm {f}",
        shell=True,
    )


def mergePileup(original, subsample):
    columnsToMerge = ["A", "C", "G", "T", "cov"]
    for c in columnsToMerge:
        original[c] = original[c] + subsample[c]

    for i, row in original.iterrows():
        d = row[row["ref"]] + row[row["alt"]]
        if d > 0:
            af = row[row["alt"]] / d
        else:
            af = 0
        original.at[i, "af"] = af

    return original


def mergeHaplotypes(original, subsample):
    original = original.merge(
        subsample,
        how="outer",
        on=["chr", "SNP1", "SNP2", "SNP1Ref", "SNP1Alt", "SNP2Ref", "SNP2Alt"],
        suffixes=("", "_contaminated"),
    )
    original = original.fillna(0)

    original["RefRef"] = (original["RefRef"] + original["RefRef_contaminated"]).astype(
        int
    )
    original["RefAlt"] = (original["RefAlt"] + original["RefAlt_contaminated"]).astype(
        int
    )
    original["AltRef"] = (original["AltRef"] + original["AltRef_contaminated"]).astype(
        int
    )
    original["AltAlt"] = (original["AltAlt"] + original["AltAlt_contaminated"]).astype(
        int
    )
    original["Other"] = (original["Other"] + original["Other_contaminated"]).astype(int)

    original["multipleHaplotypes"] = (
        original[["RefRef", "RefAlt", "AltRef", "AltAlt"]] >= 3
    ).sum(axis=1) >= 3

    original = original[
        [
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
        ]
    ]

    return original


def mergeVBIDPileup(original, subsample, outpath):
    concatenated = original.merge(
        subsample,
        how="outer",
        on=["chrom", "pos", "base"],
        suffixes=("", "_contaminated"),
    )

    concatenated["depth"] = concatenated["depth"].fillna(0).astype(int)
    concatenated["l"] = concatenated["l"].fillna("")
    concatenated["q"] = concatenated["q"].fillna("")

    concatenated["depth_contaminated"] = (
        concatenated["depth_contaminated"].fillna(0).astype(int)
    )
    concatenated["l_contaminated"] = concatenated["l_contaminated"].fillna("")
    concatenated["q_contaminated"] = concatenated["q_contaminated"].fillna("")

    concatenated["depth"] = concatenated["depth"] + concatenated["depth_contaminated"]
    concatenated["l"] = concatenated["l"] + concatenated["l_contaminated"]
    concatenated["q"] = concatenated["q"] + concatenated["q_contaminated"]

    concatenated = concatenated[["chrom", "pos", "base", "depth", "l", "q"]]

    # flushing is needed otherwise sometimes the file won't be found
    outpath = open(outpath, "w")
    concatenated.to_csv(outpath, index=False, header=False, sep="\t")
    outpath.flush()
    outpath.close()


def createBam(chunk):
    # pileupEngine = PyPaCBAM(threads=2, minReadQuality=40)
    # haplotypeDetector = HaplotypeDetector(threads=2)
    for c in chunk:
        contaminant = c[0]
        contaminated = c[1]
        perc = c[2]
        index = c[3]
        basePath = c[4]

        contaminatedName = PurePath(
            contaminated.rstrip("".join(PurePath(contaminated).suffixes))
        ).stem

        contaminantName = PurePath(
            contaminant.rstrip("".join(PurePath(contaminant).suffixes))
        ).stem

        sampleName = f"{contaminantName}_{perc}_{index}_{contaminatedName}"

        print(f"Processing {sampleName}")

        contaminantBam = createContamination(contaminant, perc, index, basePath)

        pileup = pileupEngine.genotype(
            contaminantBam, reference=f"{verifyBamIDPath}/GRCh38.genome.fa"
        )

        contaminatedPileup = pd.read_table(
            f"{basePath}/original/pileups/{contaminatedName}.snps"
        )

        pileup = mergePileup(contaminatedPileup, pileup)
        pileup.to_csv(
            f"{basePath}/mixed/pileups/{sampleName}.snps",
            index=False,
            sep="\t",
        )

        haplotypes = haplotypeDetector.detectMultipleHaplotypes(
            pileup=pileup,
            alignmentFile=contaminantBam,
            referenceGenome=f"{verifyBamIDPath}/GRCh38.genome.fa",
        )

        contaminatedHaplotype = haplotypeDetector.detectMultipleHaplotypes(
            pileup=pileup,
            alignmentFile=contaminated,
            referenceGenome=f"{verifyBamIDPath}/GRCh38.genome.fa",
        )

        haplotypes.to_csv(
            f"{basePath}/mixed/haplotypes/{sampleName}.tsv",
            index=False,
            sep="\t",
        )

        haplotypes = mergeHaplotypes(haplotypes, contaminatedHaplotype)

        haplotypes.to_csv(
            f"{basePath}/mixed/haplotypes/{sampleName}.tsv",
            index=False,
            sep="\t",
        )

        runVerifyBamID(
            verifyBamIDPath,
            contaminantBam,
            f"{basePath}/{sampleName}_subsample",
            dataset="",
        )

        subsamplePileup = pd.read_table(
            f"{basePath}/{sampleName}_subsample.Pileup",
            header=None,
            names=["chrom", "pos", "base", "depth", "l", "q"],
        )
        subsamplePileup = subsamplePileup.fillna("NA")

        contaminatedPileup = pd.read_table(
            f"{basePath}/original/verifyBamID/{contaminatedName}.Pileup",
            header=None,
            names=["chrom", "pos", "base", "depth", "l", "q"],
        )

        mergeVBIDPileup(
            contaminatedPileup, subsamplePileup, f"{basePath}/{sampleName}.Pileup"
        )

        print("START SECOND VBI")
        runVerifyBamID(
            verifyBamIDPath,
            f"{basePath}/{sampleName}.Pileup",
            f"{basePath}/mixed/verifyBamID/{sampleName}",
            dataset="",
            pileup=True,
        )
        print("ENDING SECOND VBI")

        cleanup(f"{contaminantBam}")
        cleanup(f"{contaminantBam}.bai")
        cleanup(f"{basePath}/{sampleName}_subsample.Pileup")
        cleanup(f"{basePath}/{sampleName}_subsample.Ancestry")
        cleanup(f"{basePath}/{sampleName}_subsample.selfSM")


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


def analyzeBam(bam, path, reference=None):
    bamName = PurePath(bam.rstrip("".join(PurePath(bam).suffixes))).stem
    pileupEngine = PyPaCBAM(threads=4, minReadQuality=40)
    haplotypeDetector = HaplotypeDetector(threads=4)
    pileup = pileupEngine.genotype(
        bam, reference=reference, outFile=f"{path}/pileups/{bamName}.snps"
    )

    haplotypeDetector.detectMultipleHaplotypes(
        pileup=pileup,
        alignmentFile=bam,
        referenceGenome=reference,
        outFile=f"{path}/haplotypes/{bamName}.tsv",
    )

    runVerifyBamID(
        verifyBamIDPath,
        bam,
        f"{path}/verifyBamID/{bamName}",
        dataset="",
    )


def main(originalPath, basePath, samplesToCreate=50):
    bams = sorted(glob.glob(f"{originalPath}/bams/*.bam"))
    print(bams)

    for bam in bams:
        analyzeBam(bam, originalPath, reference=f"{verifyBamIDPath}/GRCh38.genome.fa")

    random.seed(0)

    contaminationPercentage = ["005", "01", "015", "020", "025", "03", "04", "05", "1"]
    samples = []
    for i, perc in enumerate(contaminationPercentage):
        for j in range(samplesToCreate):
            contaminant, contaminated = random.sample(bams, 2)
            assert contaminated != contaminant
            samples.append(
                (contaminant, contaminated, perc, i * samplesToCreate + j, basePath)
            )

    random.shuffle(samples)
    chunks = chunkise(samples, 6)
    processes = []
    for c in chunks:
        p = multiprocessing.Process(
            target=createBam,
            args=(c,),
        )
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
        p.close()


if __name__ == "__main__":
    inputBams = sys.argv[1]  # /ngsTroubleFinder/RNA/test/original/
    basePath = sys.argv[2]  # /ngsTroubleFinder/RNA/test/
    verifyBamIDPath = sys.argv[3]  # /verifyBamID/
    samples = int(sys.argv[4])  # 15
    random.seed(0)
    main(inputBams, basePath, samplesToCreate=samples)
