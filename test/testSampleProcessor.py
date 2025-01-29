import unittest
from io import StringIO
from unittest.mock import MagicMock
from unittest.mock import patch

import pandas as pd
from ngsTroubleFinder.sample import Sample
from ngsTroubleFinder.sampleProcessor import SampleProcessor


class TestMain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model = MagicMock()
        model.run.return_value = [[[0.5]]]

        sample1 = Sample(
            name="test1", bam="test.bam", sequencing="RNA", originalSex="Female"
        )

        output = "/tmp/"
        outFolders = {
            "sex": f"{output}/sex/",
            "pileups": f"{output}/pileups/",
            "haplotypes": f"{output}/haplotypes/",
            "relatedness": f"{output}/relatedness/",
            "images": f"{output}/images/",
        }

        cls.processor = SampleProcessor(
            sample1, model=model, referenceGenome="GRCh38", outFolders=outFolders
        )

        cls.pileup = pd.read_csv(
            StringIO(
                """chr	pos	rsid	ref	alt	A	C	G	T	af	cov	genotype
chr20	68351	rs757428359	A	G	13	0	0	0	0.000000	13	0/0
chr20	68363	rs200192457	A	T	129	0	0	0	0.000000	129	0/0
chr20	68373	rs745889706	T	C	0	0	70	70	0.500000	140	0/1
chr20	68373	rs745889706	T	C	0	10	0	90	0.100000	100	0/0
chr20	68373	rs745889706	T	C	0	90	0	10	0.900000	100	0/0
chr20	68373	rs745889706	T	C	0	50	0	50	0.500000	100	0/0
chr20	68373	rs745889706	T	C	0	130	0	0	1.000000	130	1/1
chrX	68373	rs745889706	T	C	0	130	0	0	1.000000	130	1/1
chr20	68375	rs754912258	A	G	5	0	5	0	0.480769	10	0/1"""
            ),
            sep="\t",
        )

        cls.haplotype = pd.read_csv(
            StringIO(
                """
    chr	SNP1	SNP2	SNP1Ref	SNP1Alt	SNP2Ref	SNP2Alt	RefRef	RefAlt	AltRef	AltAlt	Anomaly	multipleHaplotypes
chr1	633824	633887	T	C	G	T	22	0	0	0	0	False
chr1	919598	919695	A	T	C	G	0	0	0	0	19	False
chr1	3856982	3857024	A	C	C	G	3	16	11	0	0	True"""
            ),
            sep="\t",
        )

    def test_pileup(self):
        with patch("ngsTroubleFinder.pyPaCBAM.PyPaCBAM.genotype") as mock_pileupEngine:
            mock_pileupEngine.return_value = self.pileup

            pileup = self.processor._pileup()
            self.assertTupleEqual(pileup.shape, (9, 12))
            mock_pileupEngine.assert_called_once()

    def test_haplotypeDetector(self):
        with patch(
            "ngsTroubleFinder.haplotypeDetector.HaplotypeDetector.detectMultipleHaplotypes"
        ) as mock_haplotypeDetector:
            mock_haplotypeDetector.return_value = self.haplotype
            haplotypes = self.processor._haplotypeDetector(self.pileup)
            self.assertTupleEqual(haplotypes.shape, (3, 13))
            mock_haplotypeDetector.assert_called_once()

    def test_detectSexFromTranscriptomic(self):
        transcriptomicProfile = {
            "maleCount": pd.Series([100], ["test1"]),
            "femaleCount": pd.Series([100], ["test1"]),
        }
        sex = self.processor._detectSexFromTranscriptomic(transcriptomicProfile)
        self.assertEqual(sex, "Anomaly")

        transcriptomicProfile = {
            "maleCount": pd.Series([0], ["test1"]),
            "femaleCount": pd.Series([1000], ["test1"]),
        }
        sex = self.processor._detectSexFromTranscriptomic(transcriptomicProfile)
        self.assertEqual(sex, "Female")

        transcriptomicProfile = {
            "maleCount": pd.Series([1500], ["test1"]),
            "femaleCount": pd.Series([10], ["test1"]),
        }
        sex = self.processor._detectSexFromTranscriptomic(transcriptomicProfile)
        self.assertEqual(sex, "Male")

    def test_detectSexFromGenotype(self):
        sex = self.processor._detectSexFromGenotype(self.pileup)
        self.assertEqual(sex, "Anomaly")

    def test_filterPileupByCoverage(self):
        filteredPileup = self.processor._filterPileupByCoverage(self.pileup)
        self.assertTupleEqual(filteredPileup.shape, (6, 12))

    def test_getSNPsInTails(self):
        snpsInTails = self.processor._getSNPsInTails(self.pileup)
        self.assertAlmostEqual(snpsInTails, 1 / 3)

    def test_getMultipleHaplotypeFraction(self):
        haplotypeFraction = self.processor._getMultipleHaplotypeFraction(self.haplotype)
        self.assertAlmostEqual(haplotypeFraction, 1 / 3)

    def test_getHomozygousSnps(self):
        homozygousSnps = self.processor._getHomozygousSnps(self.pileup)
        self.assertAlmostEqual(homozygousSnps, 1 / 3)

    def test_predictContamination(self):
        contamination = self.processor._predictContamination(0.5, 0.3, 0.6)
        self.assertAlmostEqual(contamination, 0.5)

    def test_estimanteContamination(self):
        contamination = self.processor._estimateContamination(
            self.pileup, self.haplotype
        )
        self.assertAlmostEqual(self.processor.sample.snpsInTails, 1 / 3)
        self.assertAlmostEqual(self.processor.sample.multipleHaplotypesFraction, 1 / 3)
        self.assertAlmostEqual(self.processor.sample.homozygousSnps, 1 / 3)
        self.assertAlmostEqual(contamination, 0.5)

    def test_computeHeterozygosityRatio(self):
        heterozygosity = self.processor._computeHeterozygosityRatio(self.pileup)
        self.assertAlmostEqual(heterozygosity, 1)

    def test_processSample(self):
        pileupMock = MagicMock()
        self.processor._pileup = pileupMock
        pileupMock.return_value = self.pileup

        haplotypeMock = MagicMock()
        self.processor._haplotypeDetector = haplotypeMock
        haplotypeMock.return_value = self.haplotype

        self.processor.processSample()

        self.assertEqual(self.processor.sample.genomicSex, "Anomaly")
        self.assertAlmostEqual(self.processor.sample.contamination, 0.5)
        self.assertAlmostEqual(self.processor.sample.heterozygosityRatio, 1.0)
        self.assertEqual(self.processor.sample.transcriptomicSex, "Unknown")


if __name__ == "__main__":
    unittest.main()
