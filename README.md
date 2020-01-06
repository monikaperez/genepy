# JKBio

A list of awesome functions for every Computational Biologist


## Content

### Terra

_in ./TerraFunction.py

a file containing a set of functions that uses [dalmatian]() to interact with the [GCP]() powered genomics HPC platform: [Terra](). 
They contain a list of additional functions to do more than what is available in dalmatian

### pyDEseq2

_in ./helper/pyDESeq2.py

it is a python integration of [deseq2]() (the differential expression analyser) with [rpy2]()


### Datanalytic

a set of function to do better data analytics (supplement of numpy, pandas etc..)

### Helper

a set of plotting tools based on [matplotlib]() and [bokeh]() and additional helper functions to save data, do merging of dataframes, etc.. 
It mostly contains a list of early function that will then be respectively integrated.

### taigr

a version of taiga that do not requires RCurl (and can save you when you have a faulty RCurl-Curl link)


### data

should not contain anything when pulled but is used by any of the functions in the folder, to save some required data files


### cell_line_mapping

a set of functions to map cell line ids to other cell line ids based on an up to date google spreadsheet. 


### epigenomics/ChipHelper

ChipHelper is especially targeted to functions related to the analysis of epigenomics data. It has functions to read, merge, denoise, ChIP seq data, it contains a lot of functions required for the AML paper.


you can import files in python with e.g. `from JKBio import TerraFunction as terra`

## Reason

As I am working in different domains of computational genomics, I need to have a set of reusable function that will help me during my work.
It can be functions to work with different tools that I have to use. Functions to do some plots. etc..

I will be trying to keep each of these functions functional and documented. Feel free to pull and start use anything that might be useful to you.
If you see anything suspicious or not working. A pull request would definitely get reviewed within a day.

I hope to be able to give back to the community and maybe save a couple of hours to couple of researchers.

Best.

#/!\ Under construction /!\

jkalfon@broadinstitute.org

jkobject@gmail.com

https://jkobject.com

Apache license 2.0.