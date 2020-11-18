### Utils

a set of plotting tools based on [matplotlib]() and [bokeh]() and additional helper functions to save data, do merging of dataframes, etc.. 
It mostly contains a list of early function that will then be respectively integrated.

#### Available functions:


_in ./Helper.py_

- fileToList: convert a txt with a list of values to a python list
- listToFile: converts a list of values to a txt
- dictToFile: converts a dict to a json
- fileToDict: converts a json to a dict
- batchMove: move a lot of file in batch (can provide different locations)
- batchRename: rename a bunch of files in batch
- createFoldersFor: makes the required folders for a given filepath


- filterProteinCoding: removes all non protein coding genes from a list (you need taiga access)
- convertGenes: converts genes from a naming to another (you need taiga access)
- getSpikeInControlScales: extracts the spike in control values from a set of bam files
- GSEAonExperiments: perform GSEA to compare a bunch of conditions at once
- runERCC: creates an ERCC dashboard and extract the RNA spike ins from it (need rpy2 and ipython and R's ERCCdashboard installed)

- scatter: makes a hoverable/zoomable bokeh scatter plot
- CNV_Map: makes a hoverable Copy Number plot using bokeh
- volcano: makes a searchable volcano plot for a Differential expression experiment using bokeh
- plotCorrelationMatrix: makes a hoverable bokeh correlation matrix
- venn: makes a venn diagram from a list of sets
- mergeImages: merge multiple pngs/pdfs together into one
- addTextToImage adds a text in an image to a specific location


- grouped: to use in a forloop to group values in batch
- overlap: given two tuples, returns the overlap
- union: given two tuples, returns the union
- nans: gets nans from a panda df
- randomString: generate a random string for namings
- parrun: runs list of commands in parallel
- askif: ask the user a questions and returns the y/n answer
- inttodate: converts an int to a string date.
- datetoint: converts a date to an int

- getBamDate: get the date of creation of a bam file
- changeToBucket: move values to a bucket
