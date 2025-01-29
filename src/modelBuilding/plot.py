import glob
import os.path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import max_error
from sklearn.metrics import mean_absolute_error


def getSingleValues(results, model):
    train = results[(results["model"] == model) & (results["dataset"] == "Train")][
        "performance"
    ].iloc[0]

    test = results[(results["model"] == model) & (results["dataset"] == "Test")][
        "performance"
    ].iloc[0]

    return train, test


def plotSingleValues(dataset):
    plt.hlines(
        y=dummyTrain,
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        color="r",
    )
    plt.hlines(
        y=dummyTest,
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        linestyles="dashed",
        color="r",
    )
    plt.hlines(
        y=linearTrain,
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        color="g",
    )
    plt.hlines(
        y=linearTest,
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        linestyles="dashed",
        color="g",
    )
    plt.hlines(
        y=bayesianTrain,
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        color="b",
    )
    plt.hlines(
        y=bayesianTest,
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        linestyles="dashed",
        color="b",
    )
    plt.hlines(
        y=verifyBamIdTrain[0],
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        color="black",
    )
    plt.hlines(
        y=verifyBamIdTest[0],
        xmin=min(dataset["estimators"]),
        xmax=max(dataset["estimators"]),
        linestyles="dashed",
        color="black",
    )


def plotDataset(dataset):
    for l in sorted(set(dataset["learningRate"])):
        df = dataset[(dataset["learningRate"] == l)]
        sns.lineplot(x="estimators", y="performance", data=df, style="dataset")
        plotSingleValues(dataset)
        plt.title(f"{df['model'].iloc[0]} Learning Rate: {l}")
        plt.ylim((0, 0.03))
        plt.savefig(
            f"/home/svalentini/coderHome/rnaseqqc/src/modelBuilding/ML/DNA/{df['model'].iloc[0]}_{str(l)}.png"
        )
        plt.clf()
        plt.close()
        # plt.show()


def plotBest(dataset, verifyBamIDTest, verifyBamIdTrain):
    bests = dataset.groupby(["model", "dataset"]).min().reset_index()
    verifyBam = pd.DataFrame(
        {
            "model": ["verifyBamID", "verifyBamID"],
            "dataset": ["Test", "Train"],
            "estimators": [0, 0],
            "learningRate": [0, 0],
            "performance": [verifyBamIDTest[0], verifyBamIdTrain[0]],
            "maxError": [verifyBamIDTest[1], verifyBamIdTrain[1]],
        }
    )
    bests = pd.concat([bests, verifyBam])
    print(bests)
    fig, ax = plt.subplots(2, 1)
    sns.scatterplot(x="model", y="performance", style="dataset", data=bests, ax=ax[0])
    sns.scatterplot(x="model", y="maxError", style="dataset", data=bests, ax=ax[1])
    ax[0].set_title("DNA: Mean absolute error across classifier")
    ax[1].set_title("DNA: Max absolute error across classifier")
    plt.show()


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
        if abs(expectedContamination - observedContamination) > 0.08:
            print(f)
            print(expectedContamination)
            print(observedContamination)
            print("_______")

    error = mean_absolute_error(expected, observed)
    maxError = max_error(expected, observed)
    return error, maxError


verifyBamIdTest = parseVerifyBamID(
    "/mnt/nfs/ngsTroubleFinder/DNA/test/mixed/verifyBamID/*.selfSM",
    "/mnt/nfs/ngsTroubleFinder/DNA/test/original/verifyBamID/*.selfSM",
)

verifyBamIdTrain = parseVerifyBamID(
    "/mnt/nfs/ngsTroubleFinder/DNA/train/mixed/verifyBamID/*.selfSM",
    "/mnt/nfs/ngsTroubleFinder/DNA/train/original/verifyBamID/*.selfSM",
)

results = pd.read_table("./gridSearchDNA.tsv")

plotBest(results, verifyBamIdTest, verifyBamIdTrain)


dummyTrain, dummyTest = getSingleValues(results, "Dummy")
linearTrain, linearTest = getSingleValues(results, "LinearRegression")
bayesianTrain, bayesianTest = getSingleValues(results, "BayesianRidge")


xboost = results[(results["model"] == "XBoost")]
plotDataset(xboost)

randomForest = results[(results["model"] == "RandomForest")]
plotDataset(randomForest)

adaBoosting = results[(results["model"] == "AdaBoosting")]
plotDataset(adaBoosting)

svm = results[(results["model"] == "SVM")]
plotDataset(svm)
