## =================================================================================================
## Predict zoonotic status from virus family
## =================================================================================================

library(argparse)

parser <- ArgumentParser(description = 'Predict zoonotic status from taxonomy using a simple heuristic')

parser$add_argument('random_seed', type = 'integer', 
										help = 'a random seed (used by both R and xgboost)')

parser$add_argument('runid', type = 'character', 
										help = 'an identifier for this run, used as a base to name output files')


# Options to control training
trainGroup <- parser$add_argument_group('Training options')
trainGroup$add_argument('--trainProportion', metavar = 'p', type = 'double', default = 0.70, 
												help = 'proportion of virus species to select for each training set - all remaining virus species go to test set (default: 0.85)')

trainGroup$add_argument('--nboot', metavar = 'b', type = 'integer', default = 100, 
												help = 'number of bootstrap iterations to perform (default: 100)')

trainGroup$add_argument('--nthread', metavar = 't', type = 'integer', default = 1, 
												help = 'number of parallel threads allowed (default: 1)')



# Options to control bagging:
bagGroup <- parser$add_argument_group('Bagging options')
bagGroup$add_argument('--Ntop', type = 'integer', default = 10,
											help = 'Number of top models to use during bagging (default: 15)')



inputArgs <- parser$parse_args()


set.seed(inputArgs$random_seed)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Constants ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Hyper-parameter search (performed for each bootstrap):
CV_K <- 5                  # Number of folds for k-fold cross-validation



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Data ---------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(parallel)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(readxl)
library(ModelMetrics)
library(caret)


TrainData <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))
HoldOutData <- readRDS(file.path('CalculatedData', 'SplitData_Holdout.rds'))

TaxonomyData <- read_excel(file.path('ExternalData', 'ICTV_MasterSpeciesList_2018b.xlsx'),
													 sheet = 'ICTV 2018b Master Species #34 v', col_types = "text")
UnclassifiedTaxonomy <- read_csv(file.path('InternalData', 'Taxonomy_UnclassifiedViruses.csv'),
																 col_types = cols(.default = 'c'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Functions ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
predict_from_taxonomy <- function(train, test, cutoff) {
	## Predict zoonotic status of a test set, given the proportion of zoonotic viruses
	## observed in a test set and return a summary of accuracy
	stopifnot(sum(train$LatestSppName %in% test$LatestSppName) == 0) # Check for data leaks
	
	# Get family proportions from training data:
	familyProps <- train %>% 
		group_by(.data$Family) %>% 
		summarise(PropHumanInFamily = sum(.data$InfectsHumans) / n())
	
	# Predict
	preds <- test %>% 
		left_join(familyProps, by = "Family") %>% 
		mutate(PropHumanInFamily = if_else(is.na(.data$PropHumanInFamily), 0, .data$PropHumanInFamily)) %>% 
		mutate(PredictedStatus = .data$PropHumanInFamily > cutoff)
	
	# Return
	preds %>% 
		mutate(Cutoff = cutoff)
}


train_tax_predictor <- function(data, trainInds) {
	## Find the optimal cutoff which balances sensitivity and specificity when predicting a test set
	## 'trainInds' defines the rows in data to use for training - all remaining data
	## are treated as the test set
	
	# Get this fold's test and train data
	train <- data[trainInds, , drop = FALSE]
	
	test <- data %>% 
		filter(! .data$LatestSppName %in% train$LatestSppName)
	
	# Trying all cutoffs, in 1% increments
	preds <- lapply(seq(0, 1, by = 0.01), predict_from_taxonomy, train = train, test = test) %>%
		bind_rows()
	
	# Score cutoffs:
	# - For our purposes, this cutoff should minimise the distance between sensitivity and specificity
	accuracy <- preds %>% 
		group_by(.data$Cutoff) %>% 
		summarise(Sensitivity = sum(.data$InfectsHumans & (.data$InfectsHumans == .data$PredictedStatus)) / sum(.data$InfectsHumans),
							Specificity = sum(!.data$InfectsHumans & (.data$InfectsHumans == .data$PredictedStatus)) / sum(!.data$InfectsHumans),
							Score = abs(.data$Sensitivity - .data$Specificity)) %>% 
		mutate(Score = -Score) # Make score more intuitive, i.e. higher means better
	
	# Return the best cutoff:
	# - In if multiple equivalent cutoffs are found, select randomly
	accuracy %>% 
		ungroup() %>% 
		top_n(n = 1, wt = .data$Score) %>%
		sample_n(size = 1)
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Add family info ----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
virusFams <- c(TaxonomyData$Family, UnclassifiedTaxonomy$Family)
names(virusFams) <- c(TaxonomyData$Species, UnclassifiedTaxonomy$UniversalName)

TrainData <- TrainData %>% 
	mutate(Family = virusFams[.data$LatestSppName]) %>% 
	distinct(.data$LatestSppName, .data$Family, .data$InfectsHumans)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Training: find optimal cutoff --------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# As in main models:
# 	- generate bootstrap training samples
# 	- then cross-validate to find the best cutoff for that sample
# 	- then predict a test set to get the final reported accuracy

train_bootstrap <- function(i, data) {
	## Create a bootstrap sample of data, tune model using cross-validation, and
	## return predictions for all viruses not in the bootstrap sample
	
	# Get a bootstrap sample (stratified to maintain class imbalances):
	# - Calibration data not used here, but removed to ensure subset sizes match main training script
	train <- data %>% 
		group_by(.data$InfectsHumans) %>% 
		sample_frac(size = inputArgs$trainProportion, replace = TRUE) %>% 
		ungroup()
	
	calibrate <- data %>% 
		filter(! .data$LatestSppName %in% train$LatestSppName) %>% 
		group_by(.data$InfectsHumans) %>% 
		sample_frac(size = 0.5, replace = FALSE) %>% 
		ungroup()
	
	test <- data %>% 
		filter(! .data$LatestSppName %in% c(train$LatestSppName, calibrate$LatestSppName))
	
	# Cross-validate to get optimal cutoff:
	# - Best cutoff is the one which occurs most often, i.e. works most consistently
	# - As above, a random choice is made if multiple equivalent cutoffs are found
	cvIndexes <- createFolds(train$InfectsHumans, k = CV_K)
	
	foldCutOffs <- lapply(cvIndexes, train_tax_predictor, data = train) %>% 
		bind_rows()
	
	bestCutoff <- foldCutOffs %>% 
		group_by(.data$Cutoff) %>% 
		summarise(occurences = n()) %>% 
		ungroup() %>% 
		top_n(n = 1, wt = .data$occurences) %>%  
		sample_n(size = 1)
	
	# Predict test set:
	preds <- predict_from_taxonomy(train = train, test = test, cutoff = bestCutoff$Cutoff)
	
	# Get AUC for this model
	auc <- auc(actual = preds$InfectsHumans,
						 predicted = preds$PropHumanInFamily)
	
	# Return
	preds %>% 
		mutate(Iteration = i,
					 Cutoff = bestCutoff$Cutoff,
					 AUC = auc)
}

bootstrapPreds <- mclapply(1:inputArgs$nboot, train_bootstrap,
													 data = TrainData,
													 mc.cores = inputArgs$nthread) %>% 
	bind_rows()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Bag predictions ----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
baggedPreds <- bootstrapPreds %>% 
	group_by(.data$LatestSppName) %>% 
	top_n(n = inputArgs$Ntop, wt = .data$AUC) %>% 
	summarise(BagScore = sum(.data$PredictedStatus) / n(),
						BaggedPrediction = .data$BagScore > 0.5)  # Zoonotic when more than half of models say so
	
# Add true status:
trueStatus <- TrainData %>% 
	distinct(.data$LatestSppName, .data$InfectsHumans)

baggedPreds <- baggedPreds %>% 
	left_join(trueStatus, by = 'LatestSppName')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output ----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.create(file.path('RunData', 'TaxonomyHeuristic'))

saveRDS(bootstrapPreds, file.path('RunData', inputArgs$runid, 'Test_BootstrapPredictions.rds'))
saveRDS(baggedPreds, file.path('RunData', inputArgs$runid, 'Test_BaggedPredictions.rds'))