import argparse
import itertools
import os
from importlib.resources import files

import onnxruntime as rt
import pandas as pd

from .plotter import Plotter
from .relatednessProcessor import RelatednessProcessor
from .sample import Sample
from .sampleProcessor import SampleProcessor


class TroubleFinder:
    """
    This class instantiate the tool using a config file and starts the analysis
    Attributes:
        transcriptomic: str
            path to the transcriptomic (integer read count) file
        referenceGenome: str
            path to the reference genome
        metadata: str
            path to the metadata file
        output: str
            path to the output folder
        outFolders: dict
            internal attribute used to simplify the output creation

    """

    def __init__(self, args):
        self.transcriptomic = args.transcriptomic
        self.referenceGenome = args.reference
        self.metadata = args.metadata
        self.threads = args.processes
        self.output = args.output
        self.outFolders = {
            "pileups": f"{self.output}/pileups/",
            "haplotypes": f"{self.output}/haplotypes/",
            "relatedness": f"{self.output}/relatedness/",
            "images": f"{self.output}/images/",
        }

    def _parseMetadata(self):
        """
        Parse the metadata file for Trouble Finder
        It expects a tsv file with columns:
        Sample_Name, Bam_Path, Sequencing (RNA, DNA), Sex (Male, Female)
        Sample_Name must match with the one in transcriptomic (if present)
        Bam_Path is used to compute the pileup for the snp
        Sequencing is used to select the ML model to predict contamination
        Sex is used to check the metadata against the inference
        Returns: List(Sample)
            List of Samples to process

        """
        sample_list = []
        try:
            metadata = pd.read_csv(self.metadata, sep="\t")
        except FileNotFoundError:
            print(f"{self.metadata} file not found")
            quit()

        if not all(
            [
                item in metadata.columns
                for item in ["Sample_Name", "Bam_Path", "Sequencing", "Sex"]
            ]
        ):
            print(
                "Malformed metadata file: Metadata must contain Sample_Name, Bam_Path, Sequencing and Sex"
            )
            quit()

        for sample in metadata.itertuples():
            try:
                sample = Sample(
                    sample.Sample_Name, sample.Bam_Path, sample.Sequencing, sample.Sex
                )
                sample_list.append(sample)
            except ValueError as e:
                print(f"Invalid value in metadata: {sample}")
                print(e)

        return sample_list

    def _parseTranscriptomicFile(self):
        """
        Parse the transcriptomic file and extract the male biomarkers and female biomarkers
        It's supposed to work with an INTEGER count file from salmon (output of our rnaseq pipeline)
        Returns: Dictionary{(femaleCount, maleCount:pd.Series}
            Dictionary with the female biomarker counts and the male biomarker counts per sample


        """
        if self.transcriptomic is not None:
            maleGenes = [
                "ENSG00000129824",
                "ENSG00000198692",
                "ENSG00000067048",
                "ENSG00000012817",
            ]

            femaleGenes = ["ENSG00000229807"]

            try:
                transcriptomicProfile = pd.read_csv(
                    self.transcriptomic, sep="\t", index_col="gene_id"
                )

                transcriptomicProfile.index = transcriptomicProfile.index.str.split(
                    "."
                ).str[0]

                femaleGenesCount = transcriptomicProfile[
                    transcriptomicProfile.index.isin(femaleGenes)
                ].sum()

                maleGenesCount = transcriptomicProfile[
                    transcriptomicProfile.index.isin(maleGenes)
                ].sum()
                return {"femaleCount": femaleGenesCount, "maleCount": maleGenesCount}
            except FileNotFoundError:
                print(f"{self.transcriptomic} file not found. Please check your input")
                quit()

        else:
            return None

    def _initFileSystem(self):
        """
        Creates the output folder needed
        """
        for k in self.outFolders:
            os.makedirs(self.outFolders[k], exist_ok=True)

    def _loadModel(self, modelPath):
        """
        Load a ML predictor using pickle
        Standard ML model are provided in the package
        Args:
            modelPath: str
                path to the ONNX model

        Returns:
            ONNX Inference Session
            Session to run the inference using ONNX

        """
        predictorModel = rt.InferenceSession(modelPath)
        return predictorModel

    def _process_Batch(self, sample_list):
        """
        Performs quality control on a batch of samples
        It will save the result in the Sample class
        Args:
            sample_list: List(Sample)
                Sample list to process
        """
        transcriptomicProfile = self._parseTranscriptomicFile()

        dnaModel = self._loadModel(files("ngsTroubleFinder.tools").joinpath("DNA.onnx"))
        rnaModel = self._loadModel(files("ngsTroubleFinder.tools").joinpath("RNA.onnx"))

        for sample in sample_list:
            if sample.sequencing == "DNA":
                predictorModel = dnaModel
            elif sample.sequencing == "RNA":
                predictorModel = rnaModel
            else:
                print("Malformed sample sequencing must be DNA or RNA")
                quit()

            sample = SampleProcessor(
                sample,
                model=predictorModel,
                referenceGenome=self.referenceGenome,
                outFolders=self.outFolders,
                threads=self.threads,
            )
            sample.processSample(transcriptomicProfile=transcriptomicProfile)

    def _computeRelatedness(self, sample_list):
        """
        Compute the relatedness among all the samples
        Args:
            sample_list: List(Sample)
            samples on which we want to compute the relatedness

        Returns:
            pd.DataFrame
                DataFrame with relatedness and ibs2(the number of sites where both samples have the same genotype)
                for each pair of samples

        """
        relatednessTable = []
        processor = RelatednessProcessor()

        for sample1, sample2 in itertools.combinations(sample_list, 2):
            relatedness, ibs2 = processor.processRelatedness(sample1, sample2)
            relatednessTable.append([sample1.name, sample2.name, relatedness, ibs2])

        relatednessTable = pd.DataFrame(
            relatednessTable, columns=["Sample1", "Sample2", "relatedness", "ibs2"]
        )
        return relatednessTable

    def _qualityChecks(self, sample):
        """
        Performs some check on a sample to see if some parameter are over an arbitrary threshold
        SEX: check if all the detected sex are the same
        CONTAMINATION: fails if the ML model predicts more than 1% contamination
        NOISY_ALLELIC_FRACTIONS: fails if the fraction of variant with an allelic fraction
        between 0.05-0.15 and 0.95-0.85 is more than 1% (it can indicate contamination or bad trimming or high error rate or cancer)
        HAPLOTYPES: fraction of regions where multiple haplotypes are found (they can happen in contaminated samples, CNV or in areas that are difficult to align)
        Args:
            sample: Sample
                Sample to check

        Returns: str
            string with all the fail QC concatenated by ';', PASS if everything is ok

        """
        qc = []

        if sample.transcriptomicSex == "Unknown":
            if sample.originalSex != sample.genomicSex:
                qc += ["SEX MISMATCH"]
        else:
            if not (
                sample.originalSex == sample.genomicSex == sample.transcriptomicSex
            ):
                qc += ["SEX MISMATCH"]

        if sample.contamination > 0.01:
            qc += ["CONTAMINATION"]

        if qc == []:
            return "PASS"
        else:
            return ";".join(qc)

    def _writeReports(self, sample_list, relatednessTable):
        """
        Output the QC report in form of tsv files
        Args:
            sample_list: list of analyzed samples
            relatednessTable: relatedness table containing all the relatedness information

        """
        for sample in sample_list:
            qcStatus = self._qualityChecks(sample)
            sample.qcStatus = qcStatus

        relatednessTable.to_csv(
            f"{self.outFolders['relatedness']}Relatedness.tsv", sep="\t", index=False
        )

        report = pd.DataFrame.from_records(vars(sample) for sample in sample_list)
        report = report[
            [
                "name",
                "bam",
                "sequencing",
                "originalSex",
                "contamination",
                "heterozygosityRatio",
                "transcriptomicSex",
                "maleExpression",
                "femaleExpression",
                "genomicSex",
                "yCoverage",
                "xHeterozygosity",
                "snpsInTails",
                "multipleHaplotypesFraction",
                "homozygousSnps",
                "familyMembers",
                "familyMemberSamples",
                "replicates",
                "replicateSamples",
                "qcStatus",
            ]
        ]

        report["familyMemberSamples"] = report["familyMemberSamples"].str.join(",")
        report["replicateSamples"] = report["replicateSamples"].str.join(",")

        report.to_csv(f"{self.output}/qcReport.tsv", sep="\t", index=False)
        return report

    def batchQualityControl(self):
        """
        Main function of the module
        It perform the complete QC analysis
        It parses the metadata, init the file system, computes the pileup, haplotypes, relatedness
        it then plots the report and writes all the metric per sample as output
        """
        print("Start processing batch")
        sample_list = self._parseMetadata()
        self._initFileSystem()
        self._process_Batch(sample_list)

        print("Start processing kinship")
        relatednessTable = self._computeRelatedness(sample_list)

        print("Writing reports")
        qcReport = self._writeReports(sample_list, relatednessTable)
        plotter = Plotter(self.outFolders)
        plotter.plotQCBatch(
            qcReport, sample_list, relatednessTable, f"{self.output}/report.html"
        )


def main():
    def parseArgs():
        """
        Parse the command line and returns the arguments

        Returns
        -------

        argparse.namespace
            Namespace with parsed arguments

        """
        parser = argparse.ArgumentParser(
            prog="NGS Trouble Finder",
            description=""" NGS Trouble Finder
                            A tool for detection and quantification of contamination and kinship across NGS data
                            """,
        )
        parser.add_argument(
            "-m",
            "--metadata",
            action="store",
            type=str,
            required=True,
            help="""Path to the cohort metadata.
                    The file must be a tsv with Sample_Name, Bam_Path, Sequencing and Sex fields.
                    The sequencing must be DNA (WGS) or RNA (rna-seq).
                    Sex must be Male or Female.
                    """,
        )
        parser.add_argument(
            "-o",
            "--output",
            action="store",
            type=str,
            required=True,
            help="Path to the output folder",
        )

        parser.add_argument(
            "-t",
            "--transcriptomic",
            action="store",
            type=str,
            required=False,
            help="""Path of the salmon count file.
                    Sample Name must match the one provided in the metadata.
            """,
        )

        parser.add_argument(
            "-p",
            "--processes",
            action="store",
            type=int,
            required=False,
            default=4,
            help="Number of processes used by the pileup engines [Default: 4]",
        )
        parser.add_argument(
            "-r",
            "--reference",
            action="store",
            type=str,
            help="Path to the reference genome fasta. Used only for CRAMS",
        )
        args = parser.parse_args()

        return args

    args = parseArgs()
    TroubleFinder(args).batchQualityControl()


if __name__ == "__main__":
    main()
