import unittest
from io import StringIO
from unittest.mock import MagicMock
from unittest.mock import patch

from ngsTroubleFinder.pyPaCBAM import PyPaCBAM


class TestHaplotypeDetector(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pileupEngine = PyPaCBAM()

    def test_positionPileup(self):
        expected = {"A": 2, "C": 2, "G": 2, "T": 1}
        pyPaCBAM = MagicMock()
        pyPaCBAM.return_value = 0
        side_effect = [2, 2, 2, 1]

        def mock_pileup(ngsfile, chrom, pos, quality, readQuality, pileupPtr):
            pileupPtr[0] = 2
            pileupPtr[1] = 2
            pileupPtr[2] = 2
            pileupPtr[3] = 1
            return 0

        pyPaCBAM.side_effect = mock_pileup
        self.pileupEngine._pileupEngine.pileup = pyPaCBAM

        pileup, pileupError = self.pileupEngine._positionPileup("test", "chr1", 1000)

        self.assertDictEqual(pileup, expected)
        self.assertEqual(pileupError, 0)

    def test_computeAF(self):
        basePileup = {"A": 2, "C": 2, "G": 2, "T": 1}
        af = self.pileupEngine._computeAF(basePileup, "A", "C")
        self.assertAlmostEqual(af, 0.5)

        basePileup = {"A": 0, "C": 0, "G": 0, "T": 0}
        af = self.pileupEngine._computeAF(basePileup, "A", "C")
        self.assertEqual(af, 0)

    def test_positionGenotype(self):
        af = 0
        genotype = self.pileupEngine._positionGenotype(af)
        self.assertEqual(genotype, "0/0")

        af = 0.5
        genotype = self.pileupEngine._positionGenotype(af)
        self.assertEqual(genotype, "0/1")

        af = 1.0
        genotype = self.pileupEngine._positionGenotype(af)
        self.assertEqual(genotype, "1/1")

        af = 0.05
        genotype = self.pileupEngine._positionGenotype(af)
        self.assertEqual(genotype, "./.")

    @patch("builtins.open")
    def test_parseVcf(self, mock_open):
        # We are expecting the chunks
        # Maybe the function should be renamed as parse and chunk
        # or the functionality should be splitted
        expectedChunks = [
            [
                ["chr1", "14464", "rs546169444", "A", "T"],
                ["chr1", "58814", "rs114420996", "G", "A"],
            ],
            [
                ["chr1", "63671", "rs80011619", "G", "A"],
                ["chr1", "597549", "rs12240002", "T", "C"],
            ],
            [[""]],
            [["chr1", "597799", "rs111501994", "A", "G"]],
        ]

        vcf = StringIO(
            """##fileformat=VCFv4.0
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
        chr1	14464	rs546169444	A	T	.	.	ASP;CAF=0.9042,0.09585;COMMON=1;G5;GENEINFO=DDX11L1:100287102|WASH7P:653635;KGPhase3;R3;RS=546169444;RSPOS=14464;SAO=0;SLO;SSR=0;TOPMED=0.83032301223241590,0.16967698776758409;VC=SNV;VLD;VP=0x050100040005150026000100;WGT=1;dbSNPBuildID=142
        chr1	58814	rs114420996	G	A	.	.	ASP;CAF=0.891,0.109;COMMON=1;G5;KGPhase1;KGPhase3;RS=114420996;RSPOS=58814;SAO=0;SLO;SSR=0;TOPMED=0.89643061926605504,0.10356938073394495;VC=SNV;VLD;VP=0x050100000005150036000100;WGT=1;dbSNPBuildID=132
        chr1	63671	rs80011619	G	A	.	.	ASP;CAF=0.8125,0.1875;COMMON=1;G5;GNO;KGPhase1;KGPhase3;RS=80011619;RSPOS=63671;SAO=0;SLO;SSR=0;TOPMED=0.80337347094801223,0.19662652905198776;VC=SNV;VLD;VP=0x050100000005150136000100;WGT=1;dbSNPBuildID=131
        chr1	597549	rs12240002	T	C	.	.	ASP;CAF=0.9197,0.08027;COMMON=1;G5;GENEINFO=LOC105378947:105378947;INT;KGPhase3;RS=12240002;RSPOS=597549;SAO=0;SLO;SSR=0;TOPMED=0.83668609836901121,0.16331390163098878;VC=SNV;VLD;VP=0x050100080005150026000100;WGT=1;dbSNPBuildID=120\n
        chr1	597799	rs111501994	A	G	.	.	ASP;CAF=0.97,0.02995;COMMON=1;G5;GENEINFO=LOC105378947:105378947;GNO;INT;KGPhase1;KGPhase3;RS=111501994;RSPOS=597799;SAO=0;SLO;SSR=0;TOPMED=0.88215150356778797,0.11784849643221202;VC=SNV;VLD;VP=0x050100080005150136000100;WGT=1;dbSNPBuildID=132"""
        )
        mock_open.return_value = vcf
        chunk = self.pileupEngine.parseVcf("test")
        chunk = [x for x in chunk]

        self.assertListEqual(chunk, expectedChunks)


if __name__ == "__main__":
    unittest.main()
