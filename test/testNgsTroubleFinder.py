import unittest
from collections import namedtuple
from importlib.resources import files
from io import StringIO
from unittest import mock
from unittest import skip
from unittest.mock import DEFAULT
from unittest.mock import MagicMock
from unittest.mock import patch

import pandas as pd
from ngsTroubleFinder.ngsTroubleFinder import TroubleFinder
from ngsTroubleFinder.sample import Sample
from pandas.testing import assert_frame_equal


class TestMain(unittest.TestCase):
    def setUp(self):
        self.transcriptomic = StringIO(
            """gene_id	gene_symbol	Test01	Test02	Test03
ENSG00000000001	test1	10	10	10
ENSG00000129824	test2	0	400	10
ENSG00000198692	test3	10	450	10
ENSG00000067048	test4	1	300	10
ENSG00000012817	test5	5	100	10
ENSG00000229807	test6	4000	10	10
"""
        )
        self.metadata = StringIO(
            """Sample_Name	Bam_Path	Sequencing	Sex
Test01	test.bam	RNA	Female
Test02	test.bam	DNA	Male"""
        )

        args = namedtuple(
            "args", ["transcriptomic", "reference", "metadata", "output", "processes"]
        )
        args = args(self.transcriptomic, "test.fasta", self.metadata, "/tmp/test/", 4)
        self.ngsTroubleFinder = TroubleFinder(args)

    @classmethod
    def setUpClass(self):
        self.transcriptomic = StringIO(
            """gene_id	gene_symbol	Test01	Test02	Test03
ENSG00000000001	test1	10	10	10
ENSG00000129824	test2	0	400	10
ENSG00000198692	test3	10	450	10
ENSG00000067048	test4	1	300	10
ENSG00000012817	test5	5	100	10
ENSG00000229807	test6	4000	10	10
"""
        )
        self.metadata = StringIO(
            """Sample_Name	Bam_Path	Sequencing	Sex
Test01	test.bam	RNA	Female
Test02	test.bam	DNA	Male"""
        )

        args = namedtuple(
            "args", ["transcriptomic", "reference", "metadata", "output", "processes"]
        )
        args = args(self.transcriptomic, "test.fasta", self.metadata, "/tmp/test/", 4)
        self.ngsTroubleFinder = TroubleFinder(args)

        self.pileup = pd.read_csv(
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

        self.haplotype = pd.read_csv(
            StringIO(
                """
    chr	SNP1	SNP2	SNP1Ref	SNP1Alt	SNP2Ref	SNP2Alt	RefRef	RefAlt	AltRef	AltAlt	Anomaly	multipleHaplotypes
chr1	633824	633887	T	C	G	T	22	0	0	0	0	False
chr1	919598	919695	A	T	C	G	0	0	0	0	19	False
chr1	3856982	3857024	A	C	C	G	3	16	11	0	0	True"""
            ),
            sep="\t",
        )

    def test_parseMetadata(self):
        samples = self.ngsTroubleFinder._parseMetadata()
        self.assertEqual(len(samples), 2)

    def test_parseTranscriptomicFile(self):
        transcriptomic = self.ngsTroubleFinder._parseTranscriptomicFile()
        femaleCount = transcriptomic["femaleCount"]
        maleCount = transcriptomic["maleCount"]
        self.assertEqual(femaleCount["Test01"], 4000)
        self.assertEqual(femaleCount["Test02"], 10)
        self.assertEqual(femaleCount["Test03"], 10)

        self.assertEqual(maleCount["Test01"], 16)
        self.assertEqual(maleCount["Test02"], 1250)
        self.assertEqual(maleCount["Test03"], 40)

    def test_processBatch(self):
        sample1 = Sample("Test01", "test.bam", "RNA", "Female")
        sample2 = Sample("Test02", "test.bam", "DNA", "Male")

        sample_list = [sample1, sample2]

        mock_model = MagicMock()
        mock_model._loadModel.return_value = ["DNA", "RNA"]

        self.ngsTroubleFinder._loadModel = mock_model

        with patch(
            "ngsTroubleFinder.sampleProcessor.SampleProcessor.processSample"
        ) as mock_processSample:
            # mock_transcriptomic.return_value = self.transcriptomic
            self.ngsTroubleFinder._process_Batch(sample_list)
            self.assertEqual(mock_processSample.call_count, 2)

    def test_computeRelatedness(self):
        sample1 = Sample("Test01", "test.bam", "RNA", "Female")
        sample2 = Sample("Test02", "test.bam", "DNA", "Male")

        sample_list = [sample1, sample2]

        with patch(
            "ngsTroubleFinder.relatednessProcessor.RelatednessProcessor.processRelatedness"
        ) as mock_relatedness:
            expected = pd.read_csv(
                StringIO(
                    """Sample1\tSample2\trelatedness\tibs2
Test01\tTest02\t0.5\t100"""
                ),
                sep="\t",
            )

            mock_relatedness.return_value = (0.5, 100)

            relatedness = self.ngsTroubleFinder._computeRelatedness(sample_list)

            assert_frame_equal(expected, relatedness)

    def test_qualityChecks(self):
        expected = "SEX;CONTAMINATION"
        sample = Sample("Test", "test.bam", "RNA", "Female")
        sample.transcriptomicSex = "Male"
        sample.genomicSex = "Female"
        sample.contamination = 0.02
        sample.snpsInTails = 0.005
        sample.multipleHaplotypesFraction = 0.007

        qc = self.ngsTroubleFinder._qualityChecks(sample)
        self.assertEqual(qc, expected)

        expected = "PASS"
        sample = Sample("Test", "test.bam", "RNA", "Female")
        sample.transcriptomicSex = "Female"
        sample.genomicSex = "Female"
        sample.contamination = 0.002
        sample.snpsInTails = 0.0005
        sample.multipleHaplotypesFraction = 0.0007

        qc = self.ngsTroubleFinder._qualityChecks(sample)
        self.assertEqual(qc, expected)

    def test_writeReports(self):
        sample1 = Sample(
            "Test01",
            "test.bam",
            "DNA",
            "Female",
            transcriptomicSex="Female",
            genomicSex="Female",
            contamination=0.002,
            snpsInTails=0.0005,
            multipleHaplotypesFraction=0.0007,
        )
        sample2 = Sample(
            "Test02",
            "test.bam",
            "DNA",
            "Female",
            transcriptomicSex="Female",
            genomicSex="Female",
            contamination=0.002,
            snpsInTails=0.0005,
            multipleHaplotypesFraction=0.0007,
        )
        relatedness = pd.read_csv(
            StringIO(
                """Sample1\tSample2\trelatedness\tibs0
Test01\tTest02\t0.5\t100"""
            ),
            sep="\t",
        )

        sample_list = [sample1, sample2]

        with patch("pandas.DataFrame.to_csv") as mock_write:
            self.ngsTroubleFinder._writeReports(sample_list, relatedness)
            self.assertEqual(mock_write.call_count, 2)


if __name__ == "__main__":
    unittest.main()
