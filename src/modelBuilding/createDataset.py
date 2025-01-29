import glob
import os
import sys

import pandas as pd


def createDataset(path, contaminated=False):
    percentageMapping = {
        "005": 0.005,
        "01": 0.01,
        "015": 0.015,
        "020": 0.020,
        "025": 0.025,
        "03": 0.03,
        "04": 0.04,
        "05": 0.05,
        "1": 0.1,
    }
    trainingDataframe = {
        "Sample": [],
        "Contamination": [],
        "SnpInTails": [],
        "HomozygousSnps": [],
        "MultipleHaplotypes": [],
    }

    for f in sorted(glob.glob(f"{path}/pileups/*.snps")):
        print(f)
        if contaminated:
            try:
                contaminationPercentage = percentageMapping[
                    os.path.basename(f).split("_")[1]
                ]
            except KeyError:
                contaminationPercentage = percentageMapping[
                    os.path.basename(f).split("_")[2]
                ]
        else:
            contaminationPercentage = 0
        print(contaminationPercentage)

        pileup = pd.read_csv(f, sep="\t")
        pileup = pileup[pileup["cov"] >= 20]
        pileup = pileup[pileup["chr"] != "chrX"]
        pileup = pileup[pileup["chr"] != "chrY"]
        snpInTails = pileup[
            ((pileup["af"] > 0.05) & (pileup["af"] < 0.15))
            | ((pileup["af"] > 0.85) & (pileup["af"] < 0.95))
        ]

        homoSNPs = pileup[((pileup["af"] <= 0.02) | (pileup["af"] > 0.98))]

        fractionHomo = len(homoSNPs) / len(pileup)

        fractionInTails = len(snpInTails) / len(pileup)

        haplotype = f.replace("pileups", "haplotypes").replace("snps", "tsv")
        haplotype = pd.read_csv(haplotype, sep="\t")

        haplotype = haplotype.set_index(["chr", "SNP1", "SNP2"])

        multipleHaplotypes = len(haplotype[haplotype["multipleHaplotypes"]]) / len(
            haplotype
        )

        trainingDataframe["Sample"].append(os.path.basename(f))
        trainingDataframe["Contamination"].append(contaminationPercentage)
        trainingDataframe["SnpInTails"].append(fractionInTails)
        trainingDataframe["MultipleHaplotypes"].append(multipleHaplotypes)
        trainingDataframe["HomozygousSnps"].append(fractionHomo)

    return pd.DataFrame(trainingDataframe)


def main():
    pathContaminated = sys.argv[1]  # /ngsTroubleFinder/RNA/train/mixed
    pathOriginal = sys.argv[2]  # /ngsTroubleFinder/RNA/train/original
    outPath = sys.argv[3]  # /ngsTroubleFinder/RNA/

    print(pathContaminated)
    print(pathOriginal)

    pathTestContaminated = pathContaminated.replace("train", "test")
    pathTestOriginal = pathOriginal.replace("train", "test")

    contaminated = createDataset(pathContaminated, contaminated=True)
    original = createDataset(pathOriginal, contaminated=False)
    dataset = pd.concat([contaminated, original])
    dataset.to_csv(
        f"{outPath}/trainDataset.tsv",
        sep="\t",
        index=False,
    )

    test_contaminated = createDataset(pathTestContaminated, contaminated=True)
    test_original = createDataset(pathTestOriginal, contaminated=False)
    test_dataset = pd.concat([test_contaminated, test_original])
    test_dataset.to_csv(f"{outPath}/testDataset.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
