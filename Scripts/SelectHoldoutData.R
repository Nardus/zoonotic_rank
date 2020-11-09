#! Rscript
# Part 7 of Zoonosis prediction pipeline
#	 - Split data into complete holdout and training sets

library(argparse)
library(dplyr)
library(rprojroot)


# Input args:
parser <- ArgumentParser(description = paste('Split data into (complete) holdout and training sets.',
																						 'Note that only viruses for which whole genomes are',
																						 'available are considered candidates for the training',
																						 'and holdout datasets. Viruses with partial genomes',
																						 'are split into a second holdout dataset to allow',
																						 'testing the accrucay of predicting such viruses.',
																						 'Outputs are SplitData_Holdout.rds,', 
																						 'SplitData_PartialGenomes.rds and SplitData_Training.rds'))
parser$add_argument('randomSeed', type = 'integer', 
										help = 'a random seed to make this selection reproducible')
parser$add_argument('--holdoutProportion', metavar = 'N', type = 'double', default = 0.15, 
										help = 'proportion of data to select as holdout (default: 0.15)')

INPUT_ARGS <- parser$parse_args()

if (INPUT_ARGS$holdoutProportion < 0 | INPUT_ARGS$holdoutProportion >= 1)
	stop('holdoutProportion must be in [0, 1)')

set.seed(INPUT_ARGS$randomSeed)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
TopDir <- find_rstudio_root_file()
setwd(TopDir)

source(file.path('Utils', 'selection_utils.R'))

FinalData <- readRDS(file.path('CalculatedData', 'FinalData_Cleaned.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Extract viruses lacking complete genomes ---------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
partialData <- FinalData %>% 
	filter(WholeGenome == 0) %>% 
	select(-WholeGenome)

completeData <- FinalData %>% 
	filter(WholeGenome == 1) %>% 
	select(-WholeGenome)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Do the selection ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# To ensure this dataset is representative, this selection is stratified:
# 	- Selection is at the virus spp level, ensuring viruses with multiple strains do not end up
# 	  in both the holdout and train datasets (instead, all strains should go to one of the two
# 	  datasets)
# 	- Proportions of zoonotic / non-zoonotic viruses represent their frequency in full dataset

## Sample holdout virus species
holdoutNames <- completeData %>% 
	distinct(LatestSppName, InfectsHumans) %>% 
	group_by(InfectsHumans) %>% 
	sample_frac(size = INPUT_ARGS$holdoutProportion) %>% 
	.$LatestSppName

trainNames <- completeData %>% 
	distinct(LatestSppName) %>% 
	filter(! LatestSppName %in% holdoutNames) %>% 
	.$LatestSppName


## Select the actual datasets
holdoutData <- completeData %>% 
	filter(LatestSppName %in% holdoutNames)

trainData <- completeData %>% 
	filter(LatestSppName %in% trainNames)


stopifnot(!any(holdoutData$LatestSppName %in% trainData$LatestSppName))
stopifnot(!any(trainData$LatestSppName %in% holdoutData$LatestSppName))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat('SelectHoldoutData.R:', length(holdoutNames), 'species selected as holdout set, leaving',
		length(trainNames), 'for training\n')

saveRDS(holdoutData, file.path('CalculatedData', 'SplitData_Holdout.rds'))
saveRDS(partialData, file.path('CalculatedData', 'SplitData_PartialGenomes.rds'))
saveRDS(trainData, file.path('CalculatedData', 'SplitData_Training.rds'))
