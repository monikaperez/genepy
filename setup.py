from setuptools import setup

with open("README", 'r') as f:
    long_description = f.read()

setup(
    name='JKBio',
    version='1.0',
    description='A useful module for any CompBio',
    author='Jeremie Kalfon',
    author_email='jkobject@gmail.com',
    packages=['JKBio', 'cell-line-mapping', 'pyDESeq2'],  # same as name
    install_requires=[
            'pysam',
            'numpy',
	        'pandas',
	        'venn',
	        'sklearn',
	        'seaborn',
	        'scikit-learn',
	        'rpy2',
	        'rpy2-bioconductor-extensions',
	        'pysam',
	        'jupyter',
	        'gseapy',
	        'bokeh',
	        'igv'],  # external packages as dependencies
)
