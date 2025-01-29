import glob
import os
import sys

import numpy as np
import pandas as pd
from skl2onnx import to_onnx
from sklearn.linear_model import LinearRegression
from sklearn.metrics import max_error
from sklearn.metrics import mean_absolute_error


def parseVerifyBamID(path, original):
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
    expected = []
    observed = []

    for f in glob.glob(original):
        expectedContamination = 0
        observedContamination = pd.read_table(f)["FREEMIX"].iloc[0]
        expected.append(expectedContamination)
        observed.append(observedContamination)

    for f in glob.glob(path):
        try:
            expectedContamination = percentageMapping[
                os.path.basename(f).split("/")[-1].split("_")[1]
            ]
        except KeyError:
            expectedContamination = percentageMapping[
                os.path.basename(f).split("/")[-1].split("_")[2]
            ]
        observedContamination = pd.read_table(f)["FREEMIX"].iloc[0]
        expected.append(expectedContamination)
        observed.append(observedContamination)

    expected = np.array(expected)
    observed = np.array(observed)

    return (
        mean_absolute_error(y_true=expected, y_pred=observed),
        max_error(y_true=expected, y_pred=observed),
    )


def testModel(model, testDataset):
    X = testDataset[["SnpInTails", "MultipleHaplotypes", "HomozygousSnps"]].to_numpy(
        np.float32
    )
    y = testDataset["Contamination"]
    yHat = model.predict(X)
    return (
        mean_absolute_error(y_true=y, y_pred=yHat),
        max_error(y_true=y, y_pred=yHat),
    )


def fitLinearRegression(dataset, outFile):
    X = dataset[["SnpInTails", "MultipleHaplotypes", "HomozygousSnps"]].to_numpy(
        np.float32
    )

    y = dataset["Contamination"].astype(np.float32)
    clf = LinearRegression()
    clf.fit(X, y)
    yHat = clf.predict(X)
    onx = to_onnx(clf, X[:1])

    with open(outFile, "wb") as f:
        f.write(onx.SerializeToString())

    return (
        clf,
        mean_absolute_error(y_true=y, y_pred=yHat),
        max_error(y_true=y, y_pred=yHat),
    )


def main(basePath):
    datasetDNA = pd.read_table(f"{basePath}/DNA/trainDataset.tsv")
    testDatasetDNA = pd.read_table(f"{basePath}/DNA/testDataset.tsv")

    datasetRNA = pd.read_table(f"{basePath}/RNA/trainDataset.tsv")
    testDatasetRNA = pd.read_table(f"{basePath}/RNA/testDataset.tsv")

    absVerifyBamIDErrorRNA, maxVerifyBaMIDErrorRNA = parseVerifyBamID(
        f"{basePath}/RNA/test/mixed/verifyBamID/*.selfSM",
        f"{basePath}/RNA/test/original/verifyBamID/*.selfSM",
    )

    absVerifyBamIDErrorDNA, maxVerifyBaMIDErrorDNA = parseVerifyBamID(
        f"{basePath}/DNA/test/mixed/verifyBamID/*.selfSM",
        f"{basePath}/DNA/test/original/verifyBamID/*.selfSM",
    )

    print("RNA")
    modelRNA, absErrorRNATrain, maxErrorRNATrain = fitLinearRegression(
        datasetRNA, "../ngsTroubleFinder/tools/RNA.onnx"
    )

    absErrorRNATest, maxErrorRNATest = testModel(modelRNA, testDatasetRNA)

    print("DNA")
    modelDNA, absErrorDNATrain, maxErrorDNATrain = fitLinearRegression(
        datasetDNA, "../ngsTroubleFinder/tools/DNA.onnx"
    )

    absErrorDNATest, maxErrorDNATest = testModel(modelDNA, testDatasetDNA)

    print(
        pd.DataFrame(
            {
                "Dataset": ["RNA", "DNA", "VerifyBamIDRNA", "VerifyBamIDDNA"],
                "absErrorTrain": [absErrorRNATrain, absErrorDNATrain, np.nan, np.nan],
                "maxErrorTrain": [maxErrorRNATrain, maxErrorDNATrain, np.nan, np.nan],
                "absErrorTest": [
                    absErrorRNATest,
                    absErrorDNATest,
                    absVerifyBamIDErrorRNA,
                    absVerifyBamIDErrorDNA,
                ],
                "maxErrorTest": [
                    maxErrorRNATest,
                    maxErrorDNATest,
                    maxVerifyBaMIDErrorRNA,
                    maxVerifyBaMIDErrorDNA,
                ],
            }
        )
    )


if __name__ == "__main__":
    basePath = sys.argv[1]  # /ngsTroubleFinder/
    main(basePath)
