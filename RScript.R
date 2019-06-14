MofaRun <- function(valueList){
	library(reticulate)
	MOFAobject <- createMOFAobject(valueList)
	DataOptions <- getDefaultDataOptions()
	ModelOptions <- getDefaultModelOptions(MOFAobject)
	TrainOptions <- getDefaultTrainOptions()
	ModelOptions$numFactors <- 200
	TrainOptions$DropFactorThreshold <- 0.02
	MOFAobject <- prepareMOFA(
	  MOFAobject, 
	  DataOptions = DataOptions,
	  ModelOptions = ModelOptions,
	  TrainOptions = TrainOptions
	)
	MOFAobject <- runMOFA(MOFAobject)
	return(MOFAobject)
}