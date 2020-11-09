require(dplyr)

## Generate a dummy-coded dataset compatible with caret
#  - 'labelcol': name of the column containing the labels (response) to be predicted
#  - 'removecols': names of columns to be removed
as_caret_data <- function(data, labelcol, removecols) {
	# Dummy encoding:
	# - All factor columns except the response should be dummy-coded
	labelData <- data[, labelcol]
	
	featureData <- data[, !colnames(data) %in% c(removecols, labelcol)]
	dummyEncoded <- dummyVars(~ ., data = featureData)
	dataMat <- predict(dummyEncoded, newdata = featureData)
	
	# Add label column and return:
	finalData <- data.frame(labelData, dataMat)
	colnames(finalData) <- c(labelcol, colnames(dataMat))
	
	finalData
}



# Feature (pre-)selection:
filter_correlated_features <- function(training_data, cor_cutoff = 0.95, keep_columns = c()) {
	# Reduce the number of features by selecting the most central representative for clusters of
	# highly correlated features. Note that for very small groups of correlated features 
	# (e.g. pairs), the choice of feature to keep may be arbitrary
	# 
	# 'cor_cutoff': Features showing more than or equal to this level of correlation will be
	#               reduced by selecting a single feature
	# 'keep_columns': non-feature columns which should be retained. These will not be considered
	#                 candidates for clustering
	# 
	# returns: a reduced dataframe containing only the selected features, along with any columns 
	#          named in 'keep_columns'. Note that column order is not preserved.
	require(apcluster)
	
	# Allow keep_columns to list columns not actually present in the data
	keep_columns <- keep_columns[keep_columns %in% colnames(training_data)]
	
	# Calculate correlations:
	correlations <- training_data %>% 
		select(-one_of(keep_columns)) %>% 
		cor(use = 'pairwise.complete') %>% 
		data.frame(From = rownames(.), .) %>% 
		gather(-From, key = 'To', value = 'Similarity') %>% 
		mutate(Similarity = abs(Similarity))
	
	# ID features which are definitely kept
	keep_features <- correlations %>% 
		filter(From != To) %>% 
		group_by(From) %>% 
		summarise(Keep = all(Similarity < cor_cutoff, na.rm = T)) %>% 
		filter(Keep) %>% 
		.$From %>% 
		as.character()
		
	# Cluster remaining features to find representatives
	if (any(!correlations$From %in% keep_features)) {  # Some features failed the filter, so need to cluster
		correlations <- correlations %>% 
			filter(!(From %in% keep_features) & !(To %in% keep_features)) %>% 
			spread(key = 'To', value = 'Similarity')
		
		similarity_matrix <- correlations %>% 
			select(-From) %>% 
			as.matrix()
		
		rownames(similarity_matrix) <- correlations$From
		
		clusters <- apcluster(s = similarity_matrix)
		
		# Simply keep the exemplar for each cluster:
		keep_features <- c(keep_features, names(clusters@exemplars))
	}
	
	# Return reduced training data
	training_data %>% 
		select(one_of(keep_columns, keep_features))
}
