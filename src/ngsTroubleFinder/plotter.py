import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


class Plotter:
    """
    Class that implements the plotting logic for QC
    It uses STALICLA's color. Color should be modified once Anezka has prepared the new palette
    Attributes:
        outFolders: dictionary with the output folders
        staliclaColors: color palette for the datapoint
        colorText: color palette for the text
        colorLines: color for the line
        colorSpine: color for the spine

    """

    def __init__(
        self,
        outFolders,
        colorPalette=("#327CF3", "#FC8A00", "#53DD6C"),
        colorText="#090A0B",
        colorLines="#CFD3D6",
        colorSpine="#969CA3",
    ):
        self.outFolders = outFolders

        self.staliclaColors = colorPalette
        self.colorText = colorText
        self.colorLines = colorLines
        self.colorSpine = colorSpine

    def _resetImage(self):
        """
        Reset a matplotlib image after saving
        """
        plt.clf()
        plt.close()

    def plotTranscriptomicSex(self, sample_list):
        """
        Plots and saves (png, svg) the sexual transcriptomic profile of the samples
        X axis female expression (should be zero for males)
        Y axis male expression (should be zero for females)
        Args:
            sample_list: list of samples with the expression computed

        """
        maleExpression = [s.maleExpression + 1 for s in sample_list]
        femaleExpression = [s.femaleExpression + 1 for s in sample_list]
        originalSex = [s.originalSex for s in sample_list]
        sampleNames = [s.name for s in sample_list]
        transcriptomicSex = [s.transcriptomicSex for s in sample_list]

        toPlot = pd.DataFrame(
            {
                "Male Expression": maleExpression,
                "Female Expression": femaleExpression,
                "Original Sex": originalSex,
                "Sample Name": sampleNames,
                "Transcriptomic Sex": transcriptomicSex,
            }
        )

        sns.scatterplot(
            x="Female Expression",
            y="Male Expression",
            hue="Original Sex",
            data=toPlot,
            palette=sns.color_palette(self.staliclaColors, 2),
            hue_order=["Male", "Female"],
        )
        plt.title("Transcriptomic Sex")

        fig = plt.gcf()
        ax = plt.gca()
        fig.set_size_inches(10, 6)

        plt.axhline((236 + 1), color=self.staliclaColors[0])
        plt.axhline((1300 + 1), color=self.staliclaColors[0])

        plt.axvline((142 + 1), color=self.staliclaColors[1])

        ax.set(xscale="log", yscale="log")

        plt.savefig(
            f"{self.outFolders['images']}/transriptomicSex.png",
            bbox_inches="tight",
            dpi=300,
        )
        plt.savefig(
            f"{self.outFolders['images']}/transriptomicSex.svg",
            bbox_inches="tight",
            dpi=300,
        )
        self._resetImage()

        fig = px.scatter(
            toPlot,
            x="Female Expression",
            y="Male Expression",
            color="Original Sex",
            hover_name="Sample Name",
            hover_data="Transcriptomic Sex",
            color_discrete_map={
                "Male": self.staliclaColors[0],
                "Female": self.staliclaColors[1],
            },
        )
        fig.add_hline(236, fillcolor=self.staliclaColors[0])
        fig.add_hline(1300, fillcolor=self.staliclaColors[0])

        fig.add_vline(142, fillcolor=self.staliclaColors[1])

        return fig.to_html(full_html=False, include_plotlyjs=False)

    def plotGenomicSex(self, sample_list):
        """
        Plots and saves (png, svg) the genomic sexual profile of the samples
        X axis X chromosome heterozygosity (should be 0 for males)
        Y axis Y chromosome coverage (should be 0 for females)
        Args:
            sample_list: list of samples with the yCoverage and xHeterozygosity computed

        """

        yCoverage = [s.yCoverage for s in sample_list]
        xHeterozygosity = [s.xHeterozygosity for s in sample_list]
        originalSex = [s.originalSex for s in sample_list]
        sampleNames = [s.name for s in sample_list]
        genomicSex = [s.genomicSex for s in sample_list]

        toPlot = pd.DataFrame(
            {
                "Y Chromosome Coverage": yCoverage,
                "X Chromosome Heterozygosity": xHeterozygosity,
                "Original Sex": originalSex,
                "Sample Name": sampleNames,
                "Genomic Sex": genomicSex,
            }
        )

        sns.scatterplot(
            x="X Chromosome Heterozygosity",
            y="Y Chromosome Coverage",
            hue="Original Sex",
            data=toPlot,
            palette=sns.color_palette(self.staliclaColors, 2),
            hue_order=["Male", "Female"],
        )

        plt.title("Genomic Sex")

        fig = plt.gcf()
        fig.set_size_inches(10, 6)

        plt.axhline(10, color=self.staliclaColors[0])
        plt.axvline(1, color=self.staliclaColors[1])

        plt.savefig(
            f"{self.outFolders['images']}/genomicSex.png", bbox_inches="tight", dpi=300
        )
        plt.savefig(
            f"{self.outFolders['images']}/genomicSex.svg", bbox_inches="tight", dpi=300
        )
        self._resetImage()

        fig = px.scatter(
            toPlot,
            x="X Chromosome Heterozygosity",
            y="Y Chromosome Coverage",
            color="Original Sex",
            hover_name="Sample Name",
            hover_data="Genomic Sex",
            color_discrete_map={
                "Male": self.staliclaColors[0],
                "Female": self.staliclaColors[1],
            },
        )

        fig.add_hline(10, fillcolor=self.staliclaColors[0])
        fig.add_vline(1, fillcolor=self.staliclaColors[1])

        return fig.to_html(full_html=False, include_plotlyjs=False)

    def plotRelatedness(self, relatednessTable):
        """
        Plots and saves (png, svg) the relatedness of the samples
        X axis relatedness
        Y axis ibs2 (the number of sites where the samples have the same genotype)
        Args:
            sample_list: list of samples with the yCoverage and xHeterozygosity computed
        """
        sns.scatterplot(
            x="relatedness",
            y="ibs2",
            data=relatednessTable,
            color=self.staliclaColors[0],
        )
        plt.xlim(-0.3, 1)
        plt.xticks(np.arange(-0.3, 1.1, 0.1), rotation=45)

        plt.axvline(0.95, color=self.staliclaColors[0])

        plt.title("Relatedness")
        fig = plt.gcf()
        fig.set_size_inches(10, 6)

        plt.savefig(
            f"{self.outFolders['images']}/relatedness.png", bbox_inches="tight", dpi=300
        )
        plt.savefig(
            f"{self.outFolders['images']}/relatedness.svg", bbox_inches="tight", dpi=300
        )
        self._resetImage()

        fig = px.scatter(
            relatednessTable,
            x="relatedness",
            y="ibs2",
            hover_data=["Sample1", "Sample2"],
            color_discrete_sequence=self.staliclaColors,
        )
        fig.add_vline(0.95, fillcolor=self.staliclaColors[0])

        return fig.to_html(full_html=False, include_plotlyjs=False)

    def plotAllelicFractions(self, sample_list, coverageThreshold=20):
        """
        Plots and saves (pdf multiple pages) the per sample allelic fraction profile
        X axis histogram with allelic fractions in 20 bins
        Args:
            sample_list: list of samples with the pileup computed
            coverageThreshold: the minimum coverage to consider a locus
        """
        plots = []
        for sample in sample_list:
            pileup = sample.pileup
            pileup = pileup[pileup["cov"] >= coverageThreshold]
            pileup = pileup[pileup["chr"] != "chrX"]
            pileup = pileup[pileup["chr"] != "chrY"]

            fig = px.histogram(
                pileup, x="af", nbins=21, color_discrete_sequence=self.staliclaColors
            )
            fig.update_traces(xbins=dict(start=0.0, end=1.0, size=0.050000001))

            plots.append(
                f"""
                <div class="bg-base-200 collapse collapse-arrow">
                <input type="checkbox" class="peer" />
                <div
                    class="collapse-title text-primary-content peer-checked:bg-secondary peer-checked:text-secondary-content">
                <h2 class="text-3xl font-bold">{sample.name}</h2>
                </div>
                <div
                class="collapse-content text-primary-content peer-checked:bg-secondary peer-checked:text-secondary-content">
                    <div class="hero bg-base-200">
                    <div class="hero-content flex-col lg:flex-row">
                        {fig.to_html(
                        full_html=False, include_plotlyjs=False,)}
                    </div>
                </div>
            </div>
            </div>
                """
            )

        with PdfPages(
            f"{self.outFolders['images']}/allelicFractions.pdf",
        ) as pdf:
            for sample in sample_list:
                pileup = sample.pileup
                pileup = pileup[pileup["cov"] >= coverageThreshold]
                pileup = pileup[pileup["chr"] != "chrX"]
                pileup = pileup[pileup["chr"] != "chrY"]

                sns.histplot(x="af", data=pileup, bins=20, color=self.staliclaColors[0])
                plt.title(f"{sample.name}")

                fig = plt.gcf()
                fig.set_size_inches(10, 6)

                pdf.savefig()
                self._resetImage()

        return plots

    def renderQCTable(self, qcTable):
        qcTable = qcTable[
            [
                "name",
                "sequencing",
                "qcStatus",
                "contamination",
                "originalSex",
                "transcriptomicSex",
                "genomicSex",
                "familyMembers",
                "replicates",
            ]
        ].rename(
            columns={
                "name": "Sample Name",
                "sequencing": "Sequencing",
                "qcStatus": "QC Status",
                "contamination": "Inferred Contamination",
                "originalSex": "Original Sex",
                "transcriptomicSex": "Transcriptomic Sex",
                "genomicSex": "Genomic Sex",
                "familyMembers": "#Family Members",
                "replicates": "#Replicates/Twins",
            },
        )
        return qcTable.to_html(classes='qcTable" id = "qcTable', index=False)

    def renderKinshipTable(self, qcTable):
        kinshipTable = qcTable[
            [
                "name",
                "familyMemberSamples",
                "replicateSamples",
            ]
        ].rename(
            columns={
                "name": "Sample Name",
                "familyMemberSamples": "Inferred Family Members",
                "replicateSamples": "Inferred Replicates/Twins",
            },
        )
        return kinshipTable.to_html(
            classes='kinshipTable" id = "kinshipTable', index=False
        )

    def renderRelatednessTable(self, relatednessTable):
        relatednessTable = relatednessTable[
            ["Sample1", "Sample2", "relatedness", "ibs2"]
        ].rename(
            columns={
                "Sample1": "Sample",
                "Sample2": "Sample",
                "relatedness": "Relatedness",
                "ibs2": "ibs2",
            }
        )
        return relatednessTable.to_html(
            classes='relatednessTable" id = "relatednessTable', index=False
        )

    def plotQCBatch(self, qc, sample_list, relatedness, reportFile):
        """
        Function to plot all the QC metrics
        Args:
            sample_list: list of samples
            relatedness: relatedness table

        """
        qcTable = self.renderQCTable(qc)
        transcriptomicSexPlot = self.plotTranscriptomicSex(sample_list)
        genomicSexPlot = self.plotGenomicSex(sample_list)
        kinshipPlot = self.plotRelatedness(relatedness)
        kinshipTable = self.renderKinshipTable(qc)
        relatednessTable = self.renderRelatednessTable(relatedness)
        allelicFractionPlots = self.plotAllelicFractions(sample_list)

        html_string = f"""
        <!doctype html>
        <html data-theme="nord">
        <head>
            <meta charset="UTF-8">
            <title>NGSTroubleFinder QC Report</title>
        </head>

        <body>
        <link href="https://cdn.jsdelivr.net/npm/simple-datatables@latest/dist/style.css" rel="stylesheet" type="text/css">
        <link href="https://cdn.jsdelivr.net/npm/daisyui@4.12.14/dist/full.min.css" rel="stylesheet" type="text/css" />
        <script src="https://cdn.jsdelivr.net/npm/simple-datatables@latest" type="text/javascript"></script>
        <script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
        <script src="https://cdn.tailwindcss.com"></script>

        <div class="divider"></div>
        <h1 class="text-5xl font-bold" style="text-align:center">NGSTroubleFinder Quality Control Report</h1>
        <div class="divider"></div>

        <div role="tablist" class="tabs tabs-boxed">

        <input
        type="radio"
        name="my_tabs_1"
        role="tab" class="tab"
        aria-label="Quality Control Summary"
        checked="checked" />
        <div role="tabpanel" class="tab-content">
            <div class="divider"></div>
            <h2 class="text-3xl font-bold">Quality Control Summary</h2>
            <div class="divider"></div>
            <div class="hero bg-base-200">
                <div class="hero-content text-center">
                    {qcTable}
                    <script>const qcTable = new simpleDatatables.DataTable("#qcTable");</script>
                </div>
            </div>
        </div>

        <input
        type="radio"
        name="my_tabs_1"
        role="tab"
        class="tab"
        aria-label="Sex Inference"
        />
            <div role="tabpanel" class="tab-content">
            <div class="divider"></div>
                <div class="flex w-full">
                <div class="card bg-base-300 rounded-box grid h-20 flex-grow place-items-center">
            <h2 class="text-3xl font-bold">Transcriptomic Sex</h2>
                <div class="hero bg-base-200">
                    <div class="hero-content flex-col lg:flex-row">
                    {transcriptomicSexPlot}
                    </div>
                </div>
                </div>
                <div class="card bg-base-300 rounded-box grid h-20 flex-grow place-items-center">
                <h2 class="text-3xl font-bold">Genomic Sex</h2>
                <div class="hero bg-base-200">
                    <div class="hero-content flex-col lg:flex-row">
                    {genomicSexPlot}
                    </div>
                </div>

                </div>
                </div>

                        </div>

            <input type="radio"
            name="my_tabs_1"
            role="tab"
            class="tab"
            aria-label="Kinship Analysis" />
            <div role="tabpanel" class="tab-content">
            <div class="divider"></div>
            <h2 class="text-3xl font-bold">Kinship Analysis</h2>
            <div class="divider"></div>
                <div class="hero bg-base-200">
                    <div class="hero-content flex-col lg:flex-row">
                    {kinshipPlot}
                    </div>
                </div>
                Each point represent a pair of samples. Common Index Values
                <ul>
                <li>~0.25 - grandparent/grandkid</li>
                <li>~0.5 parent/kid or siblings </li>
                <li>~1 twins or replicas of the same individual.</li>
                </ul>

            <div class="bg-base-200 collapse collapse-arrow">
            <input type="checkbox" class="peer" />
            <div
            class="collapse-title text-primary-content peer-checked:bg-secondary peer-checked:text-secondary-content">
            <h2 class="text-3xl font-bold">Kinship Table</h2>
            </div>
            <div
            class="collapse-content text-primary-content peer-checked:bg-secondary peer-checked:text-secondary-content">
                <div class="hero bg-base-200">
                    <div class="hero-content flex-col lg:flex-row">
                    {kinshipTable}
                        <script>const kinshipTable = new simpleDatatables.DataTable("#kinshipTable");</script>
                    </div>
                </div>
            </div>
            </div>

            <div class="bg-base-200 collapse collapse-arrow">
            <input type="checkbox" class="peer" />
            <div
            class="collapse-title text-primary-content peer-checked:bg-secondary peer-checked:text-secondary-content">
            <h2 class="text-3xl font-bold">Relatedness Table</h2>
            </div>
            <div
            class="collapse-content text-primary-content peer-checked:bg-secondary peer-checked:text-secondary-content">
                <div class="hero bg-base-200">
                    <div class="hero-content flex-col lg:flex-row">
                    {relatednessTable}
                    <script>const relatednessTable = new simpleDatatables.DataTable("#relatednessTable");</script>
                    </div>
                </div>
            </div>
            </div>

        </div>

        <input type="radio"
        name="my_tabs_1"
        role="tab"
        class="tab"
        aria-label="Advanced" />
        <div role="tabpanel" class="tab-content">
        <div class="divider"></div>
        <h2 class="text-3xl font-bold">Allelic Fraction Distributions</h2>
        <div class="divider"></div>
        {''.join(allelicFractionPlots)}

        </div>
        </div>

        </body>
    </html>"""

        with open(reportFile, "w") as out:
            out.write(html_string)
