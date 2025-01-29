import unittest
from collections import namedtuple
from io import StringIO
from unittest.mock import MagicMock

import pandas as pd
from ngsTroubleFinder.haplotypeDetector import HaplotypeDetector


class TestHaplotypeDetector(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.detector = HaplotypeDetector()

        testPileup = (
            "chr	pos	rsid	ref	alt	A	C	G	T	af	cov\n"
            "chr20	64073787	.	A	C	21	0	0	0	0.000000	21\n"
            "chrX	135342337	.	C	G	0	21	0	0	0.000000	21\n"
            "chr12	297126	.	G	A	5	0	15	0	0.250000	21\n"
            "chr12	21501897	.	G	T	0	0	0	21	1.000000	21\n"
            "chr1	3869625	.	G	A	0	0	0	0	0.000000	0\n"
            "chr19	17826620	.	A	G	10	0	11	0	0.523810	21\n"
            "chr19	14381453	.	T	C	0	2	0	19	0.095238	21\n"
            "chrX	13502337	.	A	G	0	0	21	0	1.000000	21\n"
            "chr1	827212	rs3115850	A	G	5	0	5	0	0.523810	10\n"
            "chr1	827221	rs3115849	A	G	5	0	5	0	0.523810	10\n"
            "chr1	826372	rs1057213	A	G	10	0	11	0	0.523810	21\n"
            "chr1	826420	rs1064272	T	C	0	2	0	19	0.095238	21\n"
        )

        testGroups = (
            "chr	pos	rsid	group\n"
            "chr1	826372.0	rs1057213	6\n"
            "chr1	826420.0	rs1064272	6\n"
            "chr1	827209.0	rs3115848	7\n"
            "chr1	827212.0	rs3131950	7\n"
            "chr1	827221.0	rs3131949	7\n"
            "chr1	827252.0	rs3131948	7\n"
        )

        self.testPileup = pd.read_csv(StringIO(testPileup), sep="\t")
        self.testGroups = pd.read_csv(StringIO(testGroups), sep="\t")

    def test_filterPileupByGroups(self):
        filtered_pileup = self.detector._filterPileupByGroups(
            self.testPileup, self.testGroups
        )
        self.assertEqual(len(filtered_pileup), 2)

    def test_analyzeHaplotype(self):
        haplotypes = {"RefRef": 0, "RefAlt": 0, "AltRef": 3, "AltAlt": 3, "Other": 4}
        haplotypes = self.detector._analyzeHaplotype(haplotypes)
        self.assertFalse(haplotypes)
        haplotypes = {"RefRef": 3, "RefAlt": 0, "AltRef": 3, "AltAlt": 3, "Other": 4}
        haplotypes = self.detector._analyzeHaplotype(haplotypes)
        self.assertTrue(haplotypes)

    def test_chunkise(self):
        toChunk = [1, 2, 3, 4, 5]
        nchunk = 2

        expected = [[1, 2, 3], [4, 5]]

        result = [chunk for chunk in self.detector.chunkise(toChunk, nchunk)]

        self.assertListEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
