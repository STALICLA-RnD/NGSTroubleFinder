import unittest
from io import StringIO
from unittest import mock
from unittest.mock import DEFAULT
from unittest.mock import patch

import pandas as pd
from ngsTroubleFinder.plotter import Plotter
from ngsTroubleFinder.sample import Sample


class TestMain(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.plotter = Plotter({"images": "/tmp", "report": "/tmp/report.html"})

    @mock.patch("matplotlib.pyplot.savefig")
    def test_plotTranscriptomicSex(self, mock_run):
        sample1 = Sample(
            name="test1",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Female",
            maleExpression=0,
            femaleExpression=100,
        )
        sample2 = Sample(
            name="test2",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            maleExpression=100,
            femaleExpression=0,
        )

        sample_list = [sample1, sample2]

        self.plotter.plotTranscriptomicSex(sample_list)
        mock_run.assert_called()

    @mock.patch("matplotlib.pyplot.savefig")
    def test_plotGenomicSex(self, mock_run):
        sample1 = Sample(
            name="test1",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Female",
            xHeterozygosity=2,
            yCoverage=1,
        )
        sample2 = Sample(
            name="test2",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            xHeterozygosity=0.01,
            yCoverage=100,
        )

        sample_list = [sample1, sample2]

        self.plotter.plotGenomicSex(sample_list)
        mock_run.assert_called()

    @mock.patch("matplotlib.pyplot.savefig")
    def test_plotRelatedness(self, mock_run):
        relatednessTable = pd.DataFrame(
            {
                "Sample1": ["Sample1", "Sample2"],
                "Sample2": ["Sample2", "Sample3"],
                "relatedness": [0.3, 0.1],
                "ibs2": [40, 50],
            }
        )

        self.plotter.plotRelatedness(relatednessTable)
        mock_run.assert_called()

    @mock.patch("matplotlib.backends.backend_pdf.PdfPages.savefig")
    def test_plotAllelicFractions(self, mock_run):
        pileup = StringIO(
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
        pileup = pd.read_csv(pileup, sep="\t")

        sample1 = Sample(
            name="test1",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Female",
            pileup=pileup,
        )
        sample2 = Sample(
            name="test2",
            bam="test.bam",
            sequencing="RNA",
            originalSex="Male",
            pileup=pileup,
        )

        sample_list = [sample1, sample2]

        self.plotter.plotAllelicFractions(sample_list)
        mock_run.assert_called()

    def test_plotQCBatch(self):
        sample1 = Sample(
            name="test1", bam="test.bam", sequencing="RNA", originalSex="Female"
        )
        sample2 = Sample(
            name="test2", bam="test.bam", sequencing="RNA", originalSex="Male"
        )

        sample_list = [sample1, sample2]

        relatednessTable = pd.DataFrame(
            {
                "Sample1": ["Sample1", "Sample2"],
                "Sample2": ["Sample2", "Sample3"],
                "relatedness": [0.3, 0.1],
                "ibs2": [40, 50],
            }
        )

        report = pd.DataFrame.from_records(vars(sample) for sample in sample_list)

        with patch.multiple(
            self.plotter,
            plotTranscriptomicSex=DEFAULT,
            plotGenomicSex=DEFAULT,
            plotAllelicFractions=DEFAULT,
            plotRelatedness=DEFAULT,
        ) as mock:
            self.plotter.plotQCBatch(
                report, sample_list, relatednessTable, "/tmp/report.html"
            )
            mock["plotTranscriptomicSex"].assert_called_once()
            mock["plotGenomicSex"].assert_called_once()
            mock["plotAllelicFractions"].assert_called_once()
            mock["plotRelatedness"].assert_called_once()


if __name__ == "__main__":
    unittest.main()
