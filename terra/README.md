# terra

a file containing a set of functions that uses [dalmatian](github.com/broadinstitute/dalmatian) to interact with the [GCP](https://cloud.google.com/storage/docs/gsutil) powered genomics HPC platform: [Terra](www.terra.bio). 
They contain a list of additional functions to do more than what is available in dalmatian

The goal is to improve reproducibility and productionalization of pipelines working with Terra.

#### Available functions:

- createManySubmissions: allows you to create many terra jobs in parallel
- waitForSubmission: an await function on Terra jobs
- removeSamples: a function that removes samples on a workspace and takes care of more edge cases (linked sample sets and pair sets..).
- uploadFromFolder: uploads fastq samples from a folder into a Terra workspace with the right namings etc..
- updateAllSampleSet: updates a sample set with all samples
- addToSampleSet: updates a sample set with some new samples
- addToPairSet: updates a pair set with some new pairs
- saveOmicsOutput:
- changeGSlocation:
- renametsvs:
- findBackErasedDuplicaBamteFromTerraBucket:
- shareTerraBams:
- shareCCLEbams:
- saveConfigs:
- cleanWorkspace:
- changeToBucket:
- delete_job:
- removeFromFailedWorkflows:
- listHeavyFiles:
- findFilesInWorkspaces:
- updateWorkflows:
- uploadWorkflows:
## highly recommanded

- firecloud-dalmatian (python)
- gsutil
- nextflow (better than terra)