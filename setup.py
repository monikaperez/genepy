from setuptools import setup
import sys
import os
if sys.version_info.major < 3 or sys.version_info.minor < 2:
    raise ValueError("genepy is only compatible with Python 3.3 and below")
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
    python_requires='>=3.5',
    install_requires=[
        'rpy2-bioconductor-extensions',
        'gseapy',
        'macs2',
        'deeptools',
        ## from requirements.txt
        "bokeh",
        "dalmatian",
        "firecloud_dalmatian",
        "google_api_python_client",
        "gsheets",
        "gspread",
        "ipython",
        "matplotlib",
        "numpy",
        "pandas",
        "Pillow",
        "pybedtools",
        "pyBigWig",
        "pysam",
        "pytest",
        "requests",
        "rpy2",
        "scikit_learn",
        "scipy",
        "seaborn",
        "setuptools",
        "taigapy",
        "typing",
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

print("once R is installed you need to have installed erccdashboard, GSEABase GSVA, DESeq2 to have access to aall the functions")

print("Finished!")
