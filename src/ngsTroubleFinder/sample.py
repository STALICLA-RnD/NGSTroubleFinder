import dataclasses
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class Sample:
    """
    Dataclass to wrap all the required information
    """

    # Name of the sample (must match the one in transcriptomic)
    name: str
    # Path to the bam file
    bam: str
    # DNA or RNA. Used to select the ML model
    sequencing: str
    # Sex from the metadata
    originalSex: str
    # Estimate of contamination
    contamination: float = np.nan
    # Heterozygosity Ratio estimate
    heterozygosityRatio: float = np.nan
    # Sex inferred from Transcriptomic
    transcriptomicSex: str | None = None
    # Expression of male biomarkers
    maleExpression: float = np.nan
    # Expression of female biomarkers
    femaleExpression: float = np.nan
    # Sex inferred from genomic information
    genomicSex: str | None = None
    # Y chromosome coverage
    yCoverage: float = np.nan
    # X heterozygosity
    xHeterozygosity: float = np.nan
    # fraction of snps in the tail of the af distribution
    snpsInTails: float = np.nan
    # fraction of regions with multiple haplotypes
    multipleHaplotypesFraction: float = np.nan
    # fraction of snps in the tail of the af distribution
    homozygousSnps: float = np.nan
    # number of inferred family members
    familyMembers: int = 0
    # family members names
    familyMemberSamples: list["str"] = dataclasses.field(default_factory=list)
    # inferred replicates/twins
    replicates: int = 0
    # inferred replicates/twins name
    replicateSamples: list["str"] = dataclasses.field(default_factory=list)
    # pileup of the sample
    pileup: pd.DataFrame | None = None
    # Regions with inferred haplotypes
    haplotypes: pd.DataFrame | None = None
    # QC status PASS or ';' concatenated string with all the fails
    qcStatus: str = ""

    def __post_init__(self):
        if not self.sequencing in ["DNA", "RNA"]:
            raise ValueError("Sequencing must be DNA or RNA")
        if not self.originalSex in ["Male", "Female", "Unknown"]:
            raise ValueError("Sex must be Male, Female or Unknown")

    def __repr__(self):
        return "{}: {}".format(self.__class__.__name__, vars(self))
