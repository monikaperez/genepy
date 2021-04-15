from setuptools import setup
import sys
import os
if sys.version_info.major < 3 or sys.version_info.minor < 2:
    raise ValueError("genepy is only compatible with Python 3.3 and above")
if sys.version_info.minor < 5:
    import warnings
    warnings.warn("genepy may not function properly on Python < 3.5")

os.system('git submodule init && git submodule sync')

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='Broad-genepy',
    version='1.0',
    description='A useful module for any CompBio',
    long_description=long_description,
    author='Jeremie Kalfon',
    author_email='jkobject@gmail.com',
    url="https://github.com/BroadInstitute/genepy",
    packages=['genepy/cell_line_mapping-master/python/cell_line_mapper',
              'genepy/epigenetics', 'genepy/mutations', 'genepy/google', 'genepy/sequencing/',
              'genepy/terra', 'genepy/rna', 'genepy/utils'],
    package_data={'genepy': ['data/*']},
    python_requires='>=3.8',
    install_requires=[
        'rpy2-bioconductor-extensions',
        'gseapy',
        'macs2',
        'deeptools',
        ## from requirements.txt
        "bokeh",
        "colorcet",
        "firecloud_dalmatian",
        "gseapy",
        "gsheets",
        "gspread",
        "matplotlib",
        "numpy",
        "oauth2client",
        "pandas",
        "Pillow",
        "pybedtools",
        "pyBigWig",
        "pysam",
        "pytest",
        "requests",
        "rpy2",
        "itermplot",
        "scikit_learn",
        "scipy",
        "seaborn",
        "statsmodels",
        "taigapy",
        "venn",
        ],  # external packages as dependencies
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)


print("You might want to install Bowtie2, samtools, bwa and R to be able to use all functions of this package:\n\
  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml\n\
  http://www.htslib.org/\n\
  https://github.com/lh3/bwa\n")

print("trying to install R packages")
os.system(
    "R -e \"if(!requireNamespace(\"BiocManager\", quietly = TRUE)){install.packages(\"BiocManager\", repos=\"http://cran.us.r-project.org\")};BiocManager::install(c(\"GSEABase\", \"erccdashboard\", \"GSVA\", \"DESeq2\"));\"")
print('if it did not work. please install R or check your R installation')
print("once R is installed you need to install erccdashboard, GSEABase GSVA, DESeq2 to have access to all the functions")

print("Finished!")
