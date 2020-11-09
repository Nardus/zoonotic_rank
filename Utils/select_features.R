#
# Utility function for feature selection
# 

library(stringr)

# USAGE:
# select_features(FinalData, predict_column = PREDICT_COLUMN, positive_name = POSITIVE_NAME,
# 								removecols = removeCols, pn_feature_function = get_pn_features,
# 								random_seed = INPUT$randomSeed,
# 								include_pn = INPUT$includePN, inlcude_taxonomy = INPUT$includeTaxonomy,
# 								n_features = INPUT$nfeatures, train_proportion = INPUT$trainProportion)


# Select the top n features predicting a response variable
# - data: the final dataset, including features
# - predict_column: the column containing the response variable
# - positive_name: value in predict_column to be considered as the 'TRUE' or positive class
# - removecols: columns in data which are not to be considered features
# - pn_feature_function: a function capable of calculating PN features
# - random_seed: a random seed to use during training
# - include_pn: Should PN features be included?
# - inlcude_taxonomy: Should taxonomy features be included?
# - n_features: number of features to retain
# - n_repeats: number of random test-train splits to perform
# - train_proportion: proportion of data to select for training
select_best_features <- function(data, predict_column, positive_name, removecols, 
																 pn_feature_function, random_seed,
																 include_pn = FALSE, include_taxonomy = FALSE,
																 n_features = 100, n_repeats = 100, train_proportion = 0.7) {
	
	# These are the defaults for xgboost
	tuning_params <- data.frame(eta = 0.3,
															max_depth = 6,
															subsample = 1,
															colsample_bytree = 1,
															nrounds = 150, # Has no default in xgboost, taken from Babayan et al. (2018) instead
															min_child_weight = 1,
															gamma = 0)
	
	varimps <- data.frame()
	
	for (i in 1:n_repeats) {  
		# Select training data
		trainData <- data %>%
			group_by_at(predict_column) %>%
			sample_frac(size = train_proportion, replace = FALSE) %>%
			ungroup()
		
		# Calculate proportions / PN features (which vary by training set)
		if (include_pn) {
			stopifnot(nrow(trainData) == length(unique(trainData$UniversalName))) # Join below only valid if UniversalNames are unique
			
			train_PN_summary <- pn_feature_function(queryAccessions = trainData$Accessions, dbData = trainData, removeSelf = TRUE)
			trainData <- full_join(trainData, train_PN_summary, by = c(UniversalName = 'queryUniversalName'))
		}
		
		if (include_taxonomy) {
			train_proportions <- trainData %>%
				group_by(.data$Taxonomy_Family) %>%
				summarise(Taxonomy_PropPositiveInFamily = sum(.data[[predict_column]] == positive_name) / n()) %>%
				ungroup()
			
			trainData <- left_join(trainData, train_proportions, by = 'Taxonomy_Family')
		}
		
		# Training
		caretData <- trainData %>%
			as_caret_data(labelcol = predict_column, removecols = removecols)
		
		xgboost_seed <- sample(1e8, size = 1)
		
		trainedModel <- train(x = caretData[, !colnames(caretData) == predict_column],
													y = caretData[[predict_column]],
													method = 'xgbTree',
													tuneGrid = tuning_params,
													trControl = trainControl(method = 'none', number = 1),
													seed = xgboost_seed, nthread = 1)
		
		current_importance <- varImp(trainedModel)
		varimps <- bind_rows(varimps, data.frame(Iteration = i,
																						 Feature = rownames(current_importance$importance),
																						 Importance = current_importance$importance$Overall,
																						 stringsAsFactors = FALSE))
		
	}
	
	## Summarise importance and select feature:
	#   - Taxonomy features need special handling, since most are one-hot encoded
	#   - Importance in this case is the sum of importance values across all one-hot encoded columns
	keep_features <- varimps %>%
		mutate(Feature = str_replace(.data$Feature, "(Taxonomy_[[:alpha:]]+)[.]*[[:alpha:]_]*$", "\\1")) %>% 
		group_by(.data$Iteration, .data$Feature) %>% 
		summarise(Importance = sum(.data$Importance)) %>% 
		
		# Summarise across iterations
		group_by(.data$Feature) %>%
		summarise(Importance = mean(.data$Importance)) %>%
		ungroup() %>%
		
		# Select top features
		top_n(n = n_features, wt = .data$Importance) %>%
		pull(.data$Feature)
		
	keep_features
}
