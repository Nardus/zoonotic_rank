#! Rscript
# Part of Zoonosis prediction pipeline
#	 - Repeatedly train and validate models on subsets of the data

suppressPackageStartupMessages({
	library(rprojroot)
	library(seqinr)
	library(ape)
	library(plyr)  # Not used directly, but loaded by caret::train() - here ensuring it is loaded before dplyr to avoid issues
	library(dplyr)
	library(tidyr)
	library(readxl)
	library(stringr)
	library(caret)
	library(doParallel)
})

RootDir <- find_rstudio_root_file()
setwd(RootDir)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Constants ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# These options are not currently exposed:
BLAST_CACHE <- file.path(RootDir, 'cached_blast_searches')


# Hyper-parameter search (performed for each bootstrap):
CV_K <- 5                  # Number of folds for k-fold cross-validation
N_HYPER_PARAMS <- 500      # Number of hyper-parameter combinations to try in each iteration


# Tuning parameter values to try:
# - N = N_HYPER_PARAMS random combinations of these will be tested 
TUNING_PARAMETERS <- list(eta = c(0.001, 0.005, seq(0.01, 0.2, by = 0.02)),
													max_depth = seq(6, 15, by = 1),
													subsample = seq(0.6, 1.0, by = 0.1),
													colsample_bytree = seq(0.5, 1.0, by = 0.1),
													nrounds = seq(50, 250, by = 10), # Keeping this somewhat low to prevent over-fitting
													min_child_weight = c(5, 6, 8, 10),
													gamma = seq(0, 7, by = 0.5))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Utility functions --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
source(file.path('Utils', 'selection_utils.R'))
source(file.path('Utils', 'subset_utils.R'))
source(file.path('Utils', 'blast_utils.R'))
source(file.path('Utils', 'PNsummary_utils.R'))
source(file.path('Utils', 'taxonomy_utils.R'))
source(file.path('Utils', 'xgboost_utils.R'))
source(file.path('Utils', 'select_features.R'))


# Record a consistent set of variables for predictions of various datasets:
record_predictions <- function(model, raw_data, caret_data, iteration, dataset_name) {
	data.frame(Iteration = iteration,
						 Dataset = dataset_name,
						 UniversalName = raw_data$UniversalName,
						 Strain = raw_data$Strain,
						 InfectsHumans = raw_data[[PREDICT_COLUMN]], 
						 Predicted = predict(model, newdata = caret_data),
						 RawScore = predict(model, newdata = caret_data, type = 'prob')[[POSITIVE_NAME]])
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Process input options and load data as needed to match these options:
source(file.path('Utils', 'process_training_options.R'))


set.seed(INPUT$RandomSeed)

BLAST_FRAGMENTS <- INPUT$nthread * 4  # Number of parallel/sequential blast searches (see blastn_parallel())
																			# Currently, 4 fragments processed on each thread

cl <- makeCluster(INPUT$nthread)
registerDoParallel(cl)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Prepare data for caret/xgboost -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FinalData$InfectsHumans <- ifelse(FinalData$InfectsHumans, 'True', 'False')
FinalData$InfectsHumans <- factor(FinalData$InfectsHumans, levels = c('True', 'False'))  # Some caret functions take first factor level as positive

POSITIVE_NAME <- 'True'
PREDICT_COLUMN <- 'InfectsHumans'


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Phylogenetic neighbourhood calculations ----------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# These vary depdending on the test/train split: only training viruses are considered 'known' and
# thus available to be phylogenetic neigbour
#	However, we can perform a single blast search to define neighbours, since the only
#  things that will change for each training set are 
#			(a) which viruses are available to match
#			(b) the e-values of these matches (since database size has changed)
#	These factors are corrected for in each iteration below

if (INPUT$includePN) {
	# Blast all-against-all
	used_seqs <- subset_sequences(AllSeqs, unique(FinalData$Accessions))
	
	full_blast_results <- cached_blast(querySeqs = used_seqs, 
																		 dbSeqs = used_seqs,
																		 max_target_seqs = 10000,  # Very high, so this param does not become the limiting factor for matches
																		 nfragments = BLAST_FRAGMENTS, 
																		 nthreads = INPUT$nthread, 
																		 cache_dir = BLAST_CACHE)
	
	
	# Set up default data and options for PN summaries:
	pn_defaults <- list(pIdentityCutoffs = c(0, 60, 70, 80),
											maxNeighbours = 5,  
											positiveName = POSITIVE_NAME,
											nthreads = INPUT$nthread,
											minimal = TRUE)
	
	
	# Main logic for calculating PN features:
	get_pn_features <- function(queryAccessions, dbData, removeSelf = FALSE,
															blastResults = full_blast_results, allSequences = used_seqs, fullData = FinalData, 
															hostPhylogeny = HostPhylo, 
															reservoirPhylogeny = ReservoirPhylo,
															pnDefaults = pn_defaults) {
		
		# Extract relevant sequences
		train_seqs <- subset_sequences(allSequences, unique(dbData$Accessions))
		
		# Correct blast results for the new database
		queryAccessions <- strsplit(queryAccessions, split = '; ') %>% 
			unlist()
		
		corrected_blast_results <- blastResults %>% 
			filter(queryID %in% queryAccessions) %>% 
			subset_blast(full_db_sequences = allSequences, 
									 subset_db_sequences = train_seqs)
		
		
		# Post-process blast results
		corrected_blast_results <- corrected_blast_results %>%
			clean_blast_hits(removeSelf = removeSelf) %>%
			add_names_to_blast(fullData)
		
		
		# Calculate PN features
		pnParams <- c(pnDefaults,
									list(namedBlastRes = corrected_blast_results,
											 matchData = unique(dbData)))
		
		pnSummary <- do.call(summarise_pn, pnParams)
		
		pnSummary
	}
} else {
	get_pn_features <- function(...) {}
}


# Setting seed again to ensure all random splits are the same across different runs
# (different input options may have triggered differing numbers of sampling steps before this point)
set.seed(INPUT$RandomSeed)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Feature selection --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Find the best features predicting the response - this speeds up code below and also helps reduce
# overfitting

# Non-feature columns:
# - These may be needed in calculating features, but should not be used for training directly
RemoveCols <- c('UniversalName', 'Strain', 'LatestSppName', 'Accessions', 'OtherTaxonomicInfo',
								'Reservoir', 'Hosts', 'VectorBorne', 'Vector', 'PubmedResults', 'Reference_Reservoir',
								'Reference_Vector', 'Reference', 'DataSource', 'HumanOnly', 'IndirectDetectionOnly',
								'HumanVirus', 'Olival', 'SppName_ICTV_MSL2018v1', 'Remove', 'X')


selected_features <- select_best_features(data = FinalData, 
																					predict_column = PREDICT_COLUMN, 
																					positive_name = POSITIVE_NAME, 
																					removecols = RemoveCols, 
																					pn_feature_function = get_pn_features, 
																					random_seed = INPUT$RandomSeed, 
																					include_pn = INPUT$includePN, 
																					include_taxonomy = INPUT$includeTaxonomy, 
																					n_features = INPUT$topFeatures, 
																					train_proportion = INPUT$trainProportion)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Train and test model -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Repeated test/train splits


predictions <- data.frame()
performance <- data.frame()
fits <- list()
calculated_data <- list()

inner_replicates <- ceiling(INPUT$nboot / INPUT$nseeds)
iteration <- 0

for (s in 1:INPUT$nseeds) {
	xboost_seed <- sample(1e8, size = 1)  # This seed will be recorded automatically under fits[[iteration]]$finalModel$param$seed
	
	for (i in 1:inner_replicates) {
		iteration <- iteration + 1
		
		## Test/train splits:
		#		- All splits stratified to maintain class imbalance
		#		- After selecting training data, half of remaining data used to calibrate model, other half used to test model (score performance)
		trainData <- FinalData %>% 
			group_by(InfectsHumans) %>% 
			sample_frac(size = INPUT$trainProportion, replace = FALSE) %>%
			ungroup()
		
		calibrationData <- FinalData %>% 
			filter(! LatestSppName %in% trainData$LatestSppName) %>% 
			group_by(InfectsHumans) %>% 
			sample_frac(size = 0.5, replace = FALSE) %>% 
			ungroup()
		
		testData <- FinalData %>% 
			filter(! LatestSppName %in% c(trainData$LatestSppName, calibrationData$LatestSppName))
		
		
		## Add PN / taxonomy features (these are based on 'known' viruses, i.e. those present in training data)
		if (INPUT$includePN) {
			train_PN_summary <- get_pn_features(queryAccessions = trainData$Accessions, dbData = trainData, removeSelf = TRUE)
			
			calibration_PN_summary <- get_pn_features(queryAccessions = calibrationData$Accessions, dbData = trainData)  # Blast db is always training data only
			test_PN_summary <- get_pn_features(queryAccessions = testData$Accessions, dbData = trainData)								 # Blast db is always training data only
			
			# Joins below only valid if UniversalNames are unique:
			stopifnot(nrow(trainData) == length(unique(trainData$UniversalName)))
			stopifnot(nrow(calibrationData) == length(unique(calibrationData$UniversalName)))
			stopifnot(nrow(testData) == length(unique(testData$UniversalName)))
			
			trainData <- full_join(trainData, train_PN_summary, by = c(UniversalName = 'queryUniversalName'))
			calibrationData <- full_join(calibrationData, calibration_PN_summary, by = c(UniversalName = 'queryUniversalName'))
			testData <- full_join(testData, test_PN_summary, by = c(UniversalName = 'queryUniversalName'))
		}
		
		if (INPUT$includeTaxonomy) {
			train_proportions <- trainData %>%
				group_by(Taxonomy_Family) %>%
				summarise(Taxonomy_PropPositiveInFamily = sum(InfectsHumans == POSITIVE_NAME) / n()) %>%
				ungroup()
			
			trainData <- left_join(trainData, train_proportions, by = 'Taxonomy_Family')
			calibrationData <- left_join(calibrationData, train_proportions, by = 'Taxonomy_Family')
			testData <- left_join(testData, train_proportions, by = 'Taxonomy_Family')
		}
		
		
		## Final set of features to use:
		all_columns <- colnames(trainData)
		removed_cols <- all_columns[!all_columns %in% selected_features]
		
		
		## Set up training:
		if (any(testData$LatestSppName %in% trainData$LatestSppName)) stop('Data leak detected in test data')
		if (any(calibrationData$LatestSppName %in% trainData$LatestSppName)) stop('Data leak detected in calibration data')
		
		caretData <- trainData %>% 
			as_caret_data(labelcol = PREDICT_COLUMN, removecols = removed_cols)
		
		training_features <- colnames(caretData) %>% 
			subset(!. %in% c(PREDICT_COLUMN, removed_cols))
		
		
		## Train/optimize, using a random search for hyper-param values
		parameter_combos <- lapply(TUNING_PARAMETERS, sample, size = N_HYPER_PARAMS, replace = T) %>% 
			bind_rows() %>% 
			as.data.frame()
		
		searchCriteria <- trainControl(method = 'adaptive_cv',
																	 number = CV_K,
																	 search = 'random',
																	 summaryFunction = twoClassSummary,
																	 classProbs = TRUE,
																	 adaptive = list(min = 3, 
																	 								alpha = 0.05,
																	 								method = 'gls', 
																	 								complete = FALSE)) 
		
		trainedModels <- train(x = caretData[, training_features], 
													 y = caretData[[PREDICT_COLUMN]],
													 method = 'xgbTree',
													 tuneGrid = parameter_combos,
													 trControl = searchCriteria,
													 metric = 'ROC',
													 seed = xboost_seed,
													 nthread = 1)  # run training in parallel and xgboost sequentially (see benchmarks here: http://appliedpredictivemodeling.com/blog/2018/1/17/parallel-processing)
		
		
		## Record predictions for each data split:
		#		- calibration data:
		calibration_caretData <- as_caret_data(calibrationData, labelcol = PREDICT_COLUMN, removecols = removed_cols)
		
		
		calibration_preds <- record_predictions(trainedModels,
																						raw_data = calibrationData,
																						caret_data = calibration_caretData,
																						iteration = iteration,
																						dataset_name = 'calibration')
		
		#  - training set:
		train_preds <- record_predictions(trainedModels, 
																			raw_data = trainData, 
																			caret_data = caretData, 
																			iteration = iteration, 
																			dataset_name = 'training')
		
		#  - test set:
		test_caretData <- as_caret_data(testData, labelcol = PREDICT_COLUMN, removecols = removed_cols)
		test_preds <- record_predictions(trainedModels, 
																		 raw_data = testData, 
																		 caret_data = test_caretData, 
																		 iteration = iteration, 
																		 dataset_name = 'test') 
		
		
		predictions <- rbind(predictions, train_preds, calibration_preds, test_preds)
		
		
		## Store models and calculated data:
		# - The entire train object is rather large, but so far we only need the final model and the
		#   data used to train it (which may now contain calculated variables), so keep just these instead
		calculated_data[[as.character(iteration)]] <- list(train = caretData,
																											 calibration = calibration_caretData,
																											 test = test_caretData)
		
		fits[[as.character(iteration)]] <- list(Iteration = iteration, 
																						finalModel = trainedModels$finalModel,
																						trainingData = trainedModels$trainingData,
																						trainingNames = trainData$LatestSppName)
	}
}





# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Save results: -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
outPath <- paste0('RunData/', INPUT$RunID, '/')
if (!dir.exists('RunData')) dir.create('RunData')
dir.create(outPath)
setwd(outPath)


saveRDS(calculated_data, file = paste0(INPUT$RunID, '_CalculatedData.rds'))
saveRDS(fits, file = paste0(INPUT$RunID, '_ModelFits.rds'))
saveRDS(predictions, file = paste0(INPUT$RunID, '_Predictions.rds'))


## Clean up:
stopCluster(cl)
warnings()
