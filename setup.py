from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='JKBio',
    version='1.0',
    description='A useful module for any CompBio',
    author='Jeremie Kalfon',
    author_email='jkobject@gmail.com',
    url="https://github.com/jkobject/JKBio",
    packages=['cell_line_mapping-master/python'],  # same as name
    python_requires='>=3.5',
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
