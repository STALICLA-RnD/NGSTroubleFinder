import os

import numpy as np
import pandas as pd

from .haplotypeDetector import HaplotypeDetector
from .pyPaCBAM import PyPaCBAM


class SampleProcessor:
    """
    Class to compute the sample property.
    Should it become part of the sample class?
    Attributes:
        sample: Sample
            Sample to characterize
        referenceGenome: str
            path to the reference genome
        model: str
            path to the pickle model used to infer contamination
            models for rna and dna are provided in the package
            a model must provide a predict function
        outFolders: Dictionary{str:str}
            dictionary containing the various output paths


    """

    def __init__(self, sample, model, referenceGenome, outFolders, threads=4):
        self.sample = sample
        self.referenceGenome = referenceGenome
        self.model = model
        self.outFolders = outFolders
        self.threads = threads

    def _pileup(self):
        """
        Compute the pileup for a sample
        Returns: pd.Dataframe
            dataframe with the pileup in the PyPaCBAM format


        """
        if os.path.isfile(f"{self.outFolders['pileups']}/{self.sample.name}.snps"):
            return pd.read_csv(
                f"{self.outFolders['pileups']}/{self.sample.name}.snps", sep="\t"
            )

        pileupEngine = PyPaCBAM(threads=self.threads)
        return pileupEngine.genotype(
            self.sample.bam,
            reference=self.referenceGenome,
            outFile=f"{self.outFolders['pileups']}/{self.sample.name}.snps",
        )

    def _haplotypeDetector(self, pileup):
        """
        Detect haplotype with a more precise second pass pileup
        in areas with close enough snps
        If the pileup is ALREADY present just read it from file to save time
        Args:
            pileup: pd.Dataframe
                pileup in PyPaCBAM format

        Returns: pd.Dataframe
            dataframe with the multiple haplotype detected

        """
        if os.path.isfile(f"{self.outFolders['haplotypes']}/{self.sample.name}.tsv"):
            return pd.read_csv(
                f"{self.outFolders['haplotypes']}/{self.sample.name}.tsv", sep="\t"
            )

        haplotypeDetector = HaplotypeDetector(threads=self.threads)
        return haplotypeDetector.detectMultipleHaplotypes(
            pileup,
            self.sample.bam,
            referenceGenome=self.referenceGenome,
            outFile=f"{self.outFolders['haplotypes']}/{self.sample.name}.tsv",
        )

    def _detectSexFromTranscriptomic(self, transcriptomicProfile):
        """

        Args:
            transcriptomicProfile: Dictionary{(femaleCount, maleCount:pd.Series}
            Dictionary with the female biomarker counts and the male biomarker counts per sample

        Returns: str
            Male if the sample is identified as male
            Female if the sample is identified as female
            Anomaly if the sample has an unexpected transcriptomic profile
            Unknown if the sample doesn't have an associated expression
        """
        try:
            maleCount = transcriptomicProfile["maleCount"][self.sample.name]
            femaleCount = transcriptomicProfile["femaleCount"][self.sample.name]
            self.sample.maleExpression = maleCount
            self.sample.femaleExpression = femaleCount

            """Thresholds are the 99 percentile from GTeX version 8 in whole blood
            Horrible but it works"""
            if maleCount > 1300 and femaleCount < 142:
                return "Male"
            elif femaleCount > 104 and maleCount < 236:
                return "Female"
            else:
                return "Anomaly"
        except KeyError:
            return "Unknown"

    def _detectSexFromGenotype(self, pileup):
        """
        Detect sex from the heterozygosity ratio of the X chromosome and the coverage on Y
        Thresholds are empirical.
        Args:
            pileup: pd.Dataframe
                Pileup in PyPaCBAM format
        Returns:
            Male if the sample is identified as male
            Female if the sample is identified as female
            Anomaly if the sample has unexpected genetic profile
        """

        # Usually 20 as coverage threshold is used but males have only one X chromosome
        pileup = self._filterPileupByCoverage(
            pileup, coverageThreshold=10, sexualChromosomes=True
        )

        yCov = pileup.loc[(pileup["chr"] == "chrY")][["cov"]].sum().iloc[0]

        homoAltX = pileup.loc[(pileup["chr"] == "chrX") & (pileup["af"] > 0.98)]

        heteroDf = pileup.loc[
            (pileup["chr"] == "chrX") & (pileup["af"] <= 0.8) & (pileup["af"] >= 0.2)
        ]

        heteroHomoRatio = len(heteroDf) / len(homoAltX)

        thrHeteroHomoRatio = 1

        self.sample.yCoverage = yCov
        self.sample.xHeterozygosity = heteroHomoRatio

        if heteroHomoRatio >= (thrHeteroHomoRatio) and yCov <= 10:
            return "Female"
        elif heteroHomoRatio < (thrHeteroHomoRatio) and yCov > 10:
            return "Male"
        else:
            return "Anomaly"

    def _filterPileupByCoverage(
        self, pileup, coverageThreshold=20, sexualChromosomes=False
    ):
        """
        Filter a pileup by keeping only the variants with at least a coverage
        and removing sexual chromosomes
        Args:
            pileup: pd.DataFrame
                pileup in pyPaCBAM format
            coverageThreshold: int
                coverage threshold. Only loci with at least this coverage will be kept
            sexualChromosomes: boolean
                If false remove the variants on sexual chromosomes
                If True keep the variants on sexual chromosomes

        Returns: pd.Dataframe
            piluep in pyPaCBAM format
        """
        pileup = pileup[pileup["cov"] >= coverageThreshold]

        if not sexualChromosomes:
            pileup = pileup[pileup["chr"] != "chrX"]
            pileup = pileup[pileup["chr"] != "chrY"]

        return pileup

    def _getSNPsInTails(self, pileup):
        """
        Get the SNPs in the noisy tails (we are expecting a very small number
        of variants here because of how NGS works)
        (RNA can have a little bit more with respect to WGS due to allele specific expression events)
        Args:
            pileup: pd.DataFrame
                pileup in PyPaCBAM format

        Returns: float
            Fraction of variants in the "noisy" tails

        """
        pileup = self._filterPileupByCoverage(pileup)

        # Empirical thresholds usually we should observe a very small amount of reads
        # in the 0.05-0.15 area of the AF distribution
        snpsInTails = pileup[
            ((pileup["af"] > 0.05) & (pileup["af"] < 0.15))
            | ((pileup["af"] > 0.85) & (pileup["af"] < 0.95))
        ]
        return len(snpsInTails) / len(pileup)

    def _getHomozygousSnps(self, pileup):
        """
        Get the SNPs in the homozygous tails
        Args:
            pileup: pd.DataFrame
                pileup in PyPaCBAM format

        Returns: float
            Fraction of variants in the homozygous tails

        """
        pileup = self._filterPileupByCoverage(pileup)

        # Empirical thresholds. We are quite strict to avoid counting noise or contamination
        homozygousSNPS = pileup[(pileup["af"] <= 0.02) | (pileup["af"] > 0.98)]
        return len(homozygousSNPS) / len(pileup)

    def _getMultipleHaplotypeFraction(self, haplotype):
        """
        Computes the fractions of regions with two close snps
        that support multipe haplotypes
        Having multiple haplotypes support can happen:
            1. In case of contamination
            2. In presence of CNV
            3. In regions that are difficult to align

        So a non zero amount is expected

        Args:
            haplotype: pd.Dataframe
                Dataframe with the haplotype information

        Returns: float
            the fraction of regions supporting multiple haplotypes

        """
        multipleHaplotypesFraction = len(
            haplotype[haplotype["multipleHaplotypes"]]
        ) / len(haplotype)
        return multipleHaplotypesFraction

    def _predictContamination(self, snpsInTails, multipleHaplotypes, homozygousSNPS):
        """
        Predict the contamination using a ML model that uses the amount of snps in tails
        and the regions with multiple haplotype detected.
        The model have been trained on artificial contaminated data with features SnpInTails and MultipleHaplotypes
        Args:
            snpsInTails: float
                fraction of snps in the tails
            multipleHaplotypes:
                fraction of regions supporting multiple haplotypes

        Returns: float
            contamination estimate [0-1]

        """
        features = np.array(
            [[snpsInTails, multipleHaplotypes, homozygousSNPS]], np.float32
        )

        predictedContamination = self.model.run(None, {"X": features})[0][0][0]
        return max(predictedContamination, 0)

    def _estimateContamination(self, pileup, haplotype):
        """
        Estimate the contamination of a sample
        Args:
            pileup: pd.DataFrame
                pileup in pyPaCBAM format
            haplotype: pd.Dataframe
                dataframe with the haplotype estimation

        Returns: float
            contamination estimate [0-1]

        """
        snpsInTailFraction = self._getSNPsInTails(pileup)
        self.sample.snpsInTails = snpsInTailFraction

        multipleHaplotypesFraction = self._getMultipleHaplotypeFraction(haplotype)
        self.sample.multipleHaplotypesFraction = multipleHaplotypesFraction

        homozygousSNPS = self._getHomozygousSnps(pileup)
        self.sample.homozygousSnps = homozygousSNPS

        contaminationEstimate = self._predictContamination(
            snpsInTailFraction, multipleHaplotypesFraction, homozygousSNPS
        )
        return contaminationEstimate

    def _computeHeterozygosityRatio(self, pileup):
        """
        Compute the heterozygosityRatio of a sample
        Heterozygosity Ratio is the ratio between heterozygous SNPs over homozygous alternative SNPS
        Args:
            pileup: pd.Dataframe
                pileup in pyPaCBAM format
        Returns: float
            Heterozygosity ratio

        """
        pileup = self._filterPileupByCoverage(pileup)
        return len(pileup[pileup["genotype"] == "0/1"]) / len(
            pileup[pileup["genotype"] == "1/1"]
        )

    def processSample(self, transcriptomicProfile=None):
        """
        Wrapper function to compute all the required statistics
        Args:
            transcriptomicProfile: Dictionary{(femaleCount, maleCount:pd.Series}

        """
        print(f"Processing pileup for sample {self.sample.name}")
        pileup = self._pileup()
        self.sample.pileup = pileup

        print(f"Processing haplotype for sample {self.sample.name}")
        haplotypes = self._haplotypeDetector(pileup)
        self.sample.haplotypes = haplotypes

        genomicSex = self._detectSexFromGenotype(pileup)
        self.sample.genomicSex = genomicSex

        print(f"Processing contamination for sample {self.sample.name}")
        contamination = self._estimateContamination(pileup, haplotypes)
        self.sample.contamination = contamination

        heterozygosityRatio = self._computeHeterozygosityRatio(pileup)
        self.sample.heterozygosityRatio = heterozygosityRatio

        if transcriptomicProfile is not None:
            transcriptomicSex = self._detectSexFromTranscriptomic(transcriptomicProfile)
            self.sample.transcriptomicSex = transcriptomicSex
        else:
            self.sample.transcriptomicSex = "Unknown"

        print(f"Finished processing sample {self.sample.name}")
