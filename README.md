# GenePy

_what is [genepy](https://en.wikipedia.org/wiki/G%C3%A9n%C3%A9pi)?_

A set of awesome functions & tools for Computational Genomists

![long genome](documentation/genome.jpg)

## Content

- **utils**: where a bunch of helper functions and usefull general scripts are stoed
  - **plots**: a set of plotting tools based on [matplotlib]() and [bokeh]() to make volcano plots / CNV maps etc..
  - **helper**: and additional helper functions to save data, do merging of dataframes...
- **terra**: contains a set of functions that uses [dalmatian]() to interact with the [GCP]() powered genomics HPC platform: [Terra](). 
- **sequencing**: contains a set of function to works with bed/bam/fastqs...
- **rna**: contains function to work with RNAseq (and related) data.
  - **pyDESeq2**: it is a python integration of [deseq2]() (the differential expression analyser) with [rpy2]()
- **mutations**: a set of functions to work with maf files, vcf files etc..
- **google**: functions and packages linked to google's apis
  - **google_sheet**: function to upload a df as a google sheet
  - **gcp**: sets of functions to interact with google storage (relies on gsutil)
- **epigenetics**: where we have things related to epigenomics
  - **chipseq**: has functions to read, merge, denoise, ChIP seq data, it contains a lot of functions required for the AML paper.

### Helper tools

_tools that you do not need to use directly as they have binding functions in GenePy._ 

- **epigenetics/rose:**: where an updated version of the rose algorithm is stored (as a git submodule) 
- **cell_line_mapping**: a set of functions to map cell line ids to other cell line ids based on an up to date google spreadsheet. 


## Install

### with pip (WIP)

`pip install GenePy`
### dev mode (better for now)

```bash
git clone git://github.com/BroadInstitute/GenePy.git
cd GenePy
git submodule update --init
```

then you can import files in python with e.g:
```python
from GenePy import TerraFunction as terra
```

if GenePy is not in your path, first do:

```python
import sys
sys.path.append(RELATIVE_PATH_TO_GenePy)
```

now you can install the necessary python packages:

```bash
pip install requirements.txt
pip install rpy2-bioconductor-extensions gseapy macs2 deeptools
```

or if not using the requirements.txt (computation results might change):

```bash
pip install numpy pandas
```

```bash
pip install bokeh dalmatian firecloud_dalmatian google_api_python_client gsheets gspread ipdb ipython matplotlib Pillow pybedtools pyBigWig pysam pytest requests rpy2 scikit_learn scipy seaborn setuptools taigapy taigapy typing venn rpy2-bioconductor-extensions gseapy macs2 deeptools
```

then install the following tools:
- [htslib/samtools](http://www.htslib.org/)
- [bwa](https://github.com/lh3/bwa)
just used once:
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

finaly you can install R packages (GSEABase, erccdashboard, GSVA, DESeq2):

```bash
R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};BiocManager::install(c("GSEABase", "erccdashboard", "GSVA", "DESeq2"));'
```
## About

please do contribute, we do not have time to fix all issues or work on feature requests

Jeremie Kalfon jkalfon@broadinstitute.org jkobject@gmail.com https://jkobject.com



Apache license 2.0.
