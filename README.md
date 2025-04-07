# NGSTroubleFinder

## Introduction
NGSTroubleFinder provides a simple and effective way to detect common issues in
sequencing datasets.

It analyzes:
  - distribution of the allelic fractions of common SNPs
  - detection of regions that support multiple haplotypes
  - family relationship inside a cohort
  - sex of the samples using transcriptomic biomarkers and variants to infer the sex of the sample

## Table of Contents

- [NGSTroubleFinder](#project-name)
  - [Table of Contents](#table-of-contents)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
  - [Usage](#usage)
  - [Output](#output)
  - [Example](#Example)
  - [File Format](#file-format)
  - [Development](#development)
    - [Running Locally](#running-locally)
    - [Testing](#testing)
  - [Train Dataset](#Train-Dataset)
  - [Contact](#contact)

## Getting Started

### Prerequisites

The project have been tested with
- python==3.10.10 (versions >=3.12 know to cause problems with dependencies in certain OS)
- matplotlib==3.8.2
- numpy>=1.26.0
- pandas==2.2.0
- onnxruntime==1.20.1
- scipy==1.12.0
- plotly==5.24.1
- seaborn==0.13.2

> [!CAUTION]
> The package also requires libhts-dev for the C module.

### Installation

Run
```
sudo apt update
sudo apt install -y --no-install-recommends python3-dev libhts-dev
python3 -m pip install .
```

Alternatively you can use Conda inside the project where the repository was cloned

```
conda create -n "my_env" python=3.10 bioconda::htslib
conda activate my_env
python3 -m pip install .
```

You can find the executable in ```<path_to_my_env>/bin/ngsTroubleFinder```


You can find a docker image [here](https://hub.docker.com/r/staliclarnd/ngstroublefinder)
```
docker pull staliclarnd/ngstroublefinder:1.0.3
```

## Usage

You can run the installed tool `ngsTroubleFinder` to run an analysis
As input you should provide:
- a metadata file (mandatory)
- an output folder path (mandatory)
- a salmon file with INTEGER counts (optional, used to infer sex in RNA)
- a reference genome (mandatory for CRAMs processing, not needed when processing bams)

```
ngsTroubleFinder -m metadata.tsv -o /analysis/output -t transcriptomicProfile.tsv -r GRCh38.fasta
```

Alternatively you can use the provided docker container in the same way.
```
docker run --rm -it staliclarnd/ngstroublefinder:1.0.3 \
  ngsTroubleFinder \
	-m metadata.tsv \
	-t transcriptomicProfile.tsv \
	-o /analysis/output/
```

## Output
- report.html: Interactive output containing all the information and plots from the tool
- qcReport.tsv: Machine readable output with all the collected information per sample
- images: folder containing all the plots in png and svg format
- pileups: folder cointaing all samples' pileups
- haplotypes: folder containing all the samples' haplotype informations
- relatedness: folder containig the relatdness table statistics


## Example
Download some data and the reference genome (may require some time)
```
mkdir /tmp/test/
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/exome_alignment/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram -O /tmp/test/NA12878.cram
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/exome_alignment/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai -O /tmp/test/NA12878.cram.crai

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12891/exome_alignment/NA12891.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram -O /tmp/test/NA12891.cram
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12891/exome_alignment/NA12891.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai -O /tmp/test/NA12891.cram.crai

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O /tmp/test/reference.fa
```

Create the file ```/tmp/test/metadata.tsv```
```
Sample_Name	Bam_Path	Sequencing	Sex
NA12878	/tmp/test/NA12878.cram	DNA	Female
NA12891	/tmp/test/NA12891.cram	DNA	Male
```

Run NGSTroubleFinder
```
ngsTroubleFinder -m /tmp/test/metadata.tsv -o /tmp/test/analysis/ -r /tmp/test/reference.fa
```

### File Formats

Metadata Format:
```
Sample_Name	Bam_Path	Sequencing	Sex
Sample1	/path/to/bam/Sample1.bam	DNA	Male
Sample2	/path/to/bam/Sample2.bam	DNA	Female
Sample1-RNA	/path/to/bam/Sample1_RNA.bam	RNA	Male
Sample2-RNA	/path/to/bam/Sample2_RNA.bam	RNA	Female
```

- Sample_Name: name of the sample
  If you are using the transcriptome option the name MUST match the names in the salmon file
- Bam_Path: path to the associated bam file
  It must be indexed and it doesn't support s3 paths. (You can mount the bucket with s3fs)
- Sequencing: RNA or DNA. Used to select the ML model to infer contamination
- Sex: Sex of the sample Male, Female, Unknown


What's PyPaCBAM format:
```
chr	pos	rsid	ref	alt	A	C	G	T	af	cov	genotype
20	68351	rs757428359	A	G	130	0	0	0	0.000000	130	0/0
20	68363	rs200192457	A	T	129	0	0	0	0.000000	129	0/0
20	68373	rs745889706	T	C	0	0	0	130	0.000000	130	0/0
20	68375	rs754912258	A	G	54	0	50	0	0.480769	104	0/1
20	68396	rs138777928	C	T	0	141	0	0	0.000000	141	0/0
20	68397	rs748102612	G	A	0	0	141	0	0.000000	141	0/0
20	68406	rs771803424	A	G	140	0	0	0	0.000000	140	0/0
20	76654	rs564320474	G	T	0	0	31	0	0.000000	31	0/0
20	76658	rs745496891	C	A	0	49	0	0	0.000000	49	0/0
...
```

 - chr: chromosome
 - pos: position
 - rsid: rsid in dbSNP
 - ref: reference allele
 - alt: alternative allele
 - A: number of reads supporting A
 - C: number of reads supporting C
 - G: number of reads supporting G
 - T: number of reads supporting T
 - af: reads supporting alt / (reads supporting ref + reads supporting alt)
 - cov: coverage (number of reads or A + C + G + T)
 - genotype: Genotype (0/0, 0/1, 1/1, ./.)

### Pileup Engine
Although it's not part of the intended usage of the tool the pileup module can be used as stand alone module to compute the pileups of NGS files using python.
After installing the tool you can invoke the pileup engine in python on your vcf.
```
from ngsTroubleFinder.pyPaCBAM import PyPaCBAM

pileupEngine = PyPaCBAM(threads=4)

pileupEngine.genotype("test.bam", vcf="positions.vcf", outFile="test.snps")
```
Here we are computing a pileup of `test.bam` on a personalized `positions.vcf` outputting the pileup into `test.snps`

> [!CAUTION]
> Genotyping is performed using a heuristic approach, which may not be optimal, especially at low coverage. These results are for internal use only and should be used with caution, as their accuracy is not guaranteed.


## Development
### Running Locally

After installing you can run the tool by invoking it as `ngsTroubleFinder`

### Testing

Tests are located in the test folder. You can run the test by running
`python3 test/<test_file.py>`

## Train Dataset
The train dataset used to train the model is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15166523.svg)](https://doi.org/10.5281/zenodo.15166523)


## Contact

If you have any issue contact: samuel.valentini@stalicla.com
