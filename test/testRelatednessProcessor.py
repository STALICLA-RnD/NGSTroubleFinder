import unittest
from io import StringIO
from unittest import mock

import pandas as pd
from ngsTroubleFinder.relatednessProcessor import RelatednessProcessor
from ngsTroubleFinder.sample import Sample


class TestMain(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.processor = RelatednessProcessor()
        pileup1 = StringIO(
            """chr	pos	rsid	ref	alt	A	C	G	T	af	cov	genotype\n
                20	68351	rs757428359	A	G	130	0	0	0	0.000000	130	0/0\n
                20	68363	rs200192457	A	T	129	0	0	0	0.000000	129	0/0\n
                20	68373	rs745889706	T	C	0	0	0	130	0.000000	130	0/0\n
                20	68375	rs754912258	A	G	54	0	50	0	0.480769	104	0/1\n
                20	68396	rs138777928	C	T	0	141	0	0	0.000000	141	0/0\n
                20	68397	rs748102612	G	A	0	0	141	0	0.000000	141	0/0\n
                20	68406	rs771803424	A	G	140	0	0	0	0.000000	140	0/0\n
                20	76654	rs564320474	G	T	0	0	31	0	0.000000	31	0/0\n
                20	76658	rs745496891	C	A	0	49	0	0	0.000000	49	0/0"""
        )
        self.pileup1 = pd.read_csv(pileup1, sep="\t")

        pileup2 = StringIO(
            """chr	pos	rsid	ref	alt	A	C	G	T	af	cov	genotype\n
                20	68351	rs757428359	A	G	130	0	0	0	0.000000	130	1/1\n
                20	68363	rs200192457	A	T	129	0	0	0	0.000000	129	1/1\n
                20	68373	rs745889706	T	C	0	0	0	130	0.000000	130	1/1\n
                20	68375	rs754912258	A	G	54	0	50	0	0.480769	104	0/1\n
                20	68396	rs138777928	C	T	0	141	0	0	0.000000	141	0/0\n
                20	68397	rs748102612	G	A	0	0	141	0	0.000000	141	0/0\n
                20	68406	rs771803424	A	G	140	0	0	0	0.000000	140	0/0\n
                20	76654	rs564320474	G	T	0	0	31	0	0.000000	31	0/0\n
                20	76658	rs745496891	C	A	0	49	0	0	0.000000	49	0/0"""
        )
        self.pileup2 = pd.read_csv(pileup2, sep="\t")

        pileup3 = StringIO(
            """chr	pos	rsid	ref	alt	A	C	G	T	af	cov	genotype\n
                20	68351	rs757428359	A	G	130	0	0	0	0.000000	130	0/1\n
                20	68363	rs200192457	A	T	129	0	0	0	0.000000	129	0/1\n
                20	68373	rs745889706	T	C	0	0	0	130	0.000000	130	0/1\n
                20	68375	rs754912258	A	G	54	0	50	0	0.480769	104	0/0\n
                20	68396	rs138777928	C	T	0	141	0	0	0.000000	141	0/0\n
                20	68397	rs748102612	G	A	0	0	141	0	0.000000	141	0/0\n
                20	68406	rs771803424	A	G	140	0	0	0	0.000000	140	0/0\n
                20	76654	rs564320474	G	T	0	0	31	0	0.000000	31	0/0\n
                20	76655	rs564320474	G	T	0	0	31	0	0.000000	31	0/0\n
                20	76656	rs564320474	G	T	0	0	31	0	0.000000	31	0/1\n
                20	76658	rs745496891	C	A	0	49	0	0	0.000000	49	0/0"""
        )
        self.pileup3 = pd.read_csv(pileup3, sep="\t")

        pileup4 = StringIO(
            """chr	pos	rsid	ref	alt	A	C	G	T	af	cov	genotype\n
                20	68351	rs757428359	A	G	130	0	0	0	0.000000	130	0/1\n
                20	68363	rs200192457	A	T	129	0	0	0	0.000000	129	0/0\n
                20	68373	rs745889706	T	C	0	0	0	130	0.000000	130	0/1\n
                20	68375	rs754912258	A	G	54	0	50	0	0.480769	104	0/0\n
                20	68396	rs138777928	C	T	0	141	0	0	0.000000	141	0/1\n
                20	68397	rs748102612	G	A	0	0	141	0	0.000000	141	0/0\n
                20	68406	rs771803424	A	G	140	0	0	0	0.000000	140	0/0\n
                20	76654	rs564320474	G	T	0	0	31	0	0.000000	31	0/0\n
                20	76655	rs745496891	C	A	0	49	0	0	0.000000	49	0/0\n
                20	76657	rs745496891	C	A	0	49	0	0	0.000000	49	0/1\n
                20	76658	rs745496891	C	A	0	49	0	0	0.000000	49	0/1"""
        )
        self.pileup4 = pd.read_csv(pileup4, sep="\t")

    def test_createSNPSets(self):
        expected = {
            "ref": {(20, 68406), (20, 68396), (20, 76654), (20, 76658), (20, 68397)},
            "het": {(20, 68375)},
            "alt": {(20, 68373), (20, 68351), (20, 68363)},
            "total": {
                (20, 68406),
                (20, 68396),
                (20, 76654),
                (20, 76658),
                (20, 68397),
                (20, 68375),
                (20, 68373),
                (20, 68351),
                (20, 68363),
            },
        }

        snpSet = self.processor._createSNPSets(self.pileup2)

        self.assertDictEqual(expected, snpSet)

    @mock.patch("matplotlib.backends.backend_pdf.PdfPages.savefig")
    def test_plotAllelicFractions(self, mock_run):
        sample1 = Sample(
            name="test1",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Female",
            pileup=self.pileup1,
        )
        sample2 = Sample(
            name="test2",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            pileup=self.pileup2,
        )
        sample3 = Sample(
            name="test3",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            pileup=self.pileup3,
        )
        sample4 = Sample(
            name="test4",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            pileup=self.pileup4,
        )

        relatedness, ibs2 = self.processor._computeRelatedness(sample1, sample1)
        self.assertEqual(relatedness, 1.0)
        self.assertEqual(ibs2, 9)

        relatedness, ibs2 = self.processor._computeRelatedness(sample1, sample2)
        self.assertAlmostEqual(relatedness, -5.0)
        self.assertEqual(ibs2, 6)

        relatedness, ibs2 = self.processor._computeRelatedness(sample3, sample4)
        self.assertAlmostEqual(relatedness, 2 / 3)
        self.assertEqual(ibs2, 7)

    def test_identifyRelation(self):
        sample1 = Sample(
            name="test1",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Female",
            pileup=self.pileup1,
        )
        sample2 = Sample(
            name="test2",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            pileup=self.pileup2,
        )

        self.processor._identifyRelation(sample1, sample1, 1)
        self.assertEqual(sample1.replicates, 2)

        self.processor._identifyRelation(sample1, sample2, -6.0)
        self.assertEqual(sample1.familyMembers, 0)
        self.assertEqual(sample2.familyMembers, 0)

        self.processor._identifyRelation(sample1, sample2, 0.5)
        self.assertEqual(sample1.familyMembers, 1)
        self.assertEqual(sample2.familyMembers, 1)


if __name__ == "__main__":
    unittest.main()
