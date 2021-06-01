####################################################################################################
## 
## Utility functions for result plots
## 
####################################################################################################
require(dplyr)
require(tidyr)
require(xgboost)
require(caret)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- General ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
save_pdf <- function (p, filename, width, height) {
	# Save a single plot as pdf
	filename <- file.path(OUT_DIR, filename)
	
	pdf(filename, width = width, height = height)
	print(p)
	dev.off()
}


list.append <- function(currlist, newobj) {
	# Append to an existing list
	# - Useful for growing a list of plots to save into a single file
	ind <- length(currlist) + 1
	currlist[[ind]] <- newobj
	currlist
}


se <- function(x, na.rm = FALSE) {
	# Calculate standard error of the mean
	if (na.rm) 
		x <- x[!is.na(x)]
	
	sd(x)/sqrt(length(x))
}


standardise <- function(x) {
	# Scale values to be within [0:1], where 0 represents the smallest value
	# - Note: simply dividing by the maximum does not work for negative values
	#         (that would scale values to be in [-1:1])
	((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}


# Other utils for processing feature names:
escape_special_chars <- function(x) {
	# Escape names containing special chars by adding backticks
	# Useful for recreating the column names used by readr
	
	allowedX <- gsub('[._]', '', x)  # Remove dots and underscores, which are allowed (should not be escaped)
	
	startsSpecial <- grepl('^[_[:digit:]]', x)  # Can't start with an underscore or a letter
	containsSpecial <- grepl('[[:punct:]]', allowedX)
	
	# Escape and return:
	ifelse(startsSpecial | containsSpecial,
				 paste0('`', x, '`'),
				 x)
}

# Emulate ggplot's colour selection
gg_colors <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_colors_shuffled <- function(n) {
	# Return n colours, shuffled such that similar colours do not occur next to each other
	# Based on code in apcluster:::heatmap.AggExResult.matrix
	# - Main idea is to cut the vector of colours in half, then alternate between colours
	#   from the first and second half
	nEven <- n + if (n %% 2) 1 else 0
	inds <- as.vector(t(matrix(1:nEven, nEven / 2)))[1:n]
	
	gg_colors(n)[inds]
}

round_up <- function(x, digits = 1) {
	# Round values upwards (useful when calculating scale limits)
	ceiling(x * 10^digits) / 10^digits
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Optimize cutoff ----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Find the optimal cutoff balancing sensitivity and specificity
# - Achieved by minimizing the distance between sensitivity and specificity
find_balanced_cutoff <- function(observed_labels, predicted_score, positive_value = "True", increment_size = 0.0001) {
	dist_stat <- function(cutoff, obs, prob) {
		stopifnot(length(obs) == length(prob))
		
		pred <- prob > cutoff
		sensitivity = sum(obs & (obs == pred)) / sum(obs)
		specificity = sum(!obs & (obs == pred)) / sum(!obs)
		
		abs(sensitivity - specificity)
	}
	
	obs <- observed_labels == positive_value
	try_cutoffs <- seq(0, 1, by = increment_size)
	
	distances <- vapply(try_cutoffs, dist_stat, 
											FUN.VALUE = numeric(1),
											obs = obs, prob = predicted_score)
	
	try_cutoffs[which.min(distances)]
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Parse variable types -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

## Parse variable names to ID genomic feature classes:
id_genomic_features <- function(featureNames) {
	matchExpressions <- c('Nucleotide content' = '^[ATGCN]$',
												'AT / GC content' = '^[ATGC]{2}$',
												'Dinucleotide bias' = '^[AUGC]p[ATUGC]',
												'Bridge dinucleotide bias' = '^br[AUGC]p[AUGC]',
												'Non-bridge dinucleotide bias' = 'NonBr[AUGC]p[AUGC]',
												'Amino acid bias' = '^[A-Z].Bias',
												'Codon bias' = '[ATGC]{3}.Bias',
												'Codon pair bias' = '[ATGC]{3}.[A-Z]..[ATGC]{3}.[A-Z].')
	
	featureNames <- unique(featureNames)
	
	## Create a dataframe of matches: 
	#  - Columns are names of matchExpressions
	#  - Cell are whether or not the expression matched, for each featureName (rows)
	matches <- lapply(matchExpressions, grepl, x = featureNames)
	matches <- do.call(tibble, matches)
	
	## Check and return the name of matched feature type
	measureTypes <- matches %>% 
		mutate(featureName = featureNames) %>% 
		gather(-featureName, key = 'MeasureType', value = 'Matches') %>% 
		group_by(featureName) %>% 
		mutate(Nmatches = sum(Matches)) %>% 
		ungroup()
	
	if (any(measureTypes$Nmatches > 1))   # Non-matches allowed
		stop('Some variables match more than one genomic summary type')
	
	measureTypes %>% 
		filter(Matches) %>% 
		select(-Matches, -Nmatches)
	
}


## Parse variable names to identify complete variable types (e.g GenomicDensity, MotifDistance, etc.):
# Returns the original data frame, with relevant columns added
id_variable_types <- function(data, varnameColumn) {
	# varnameColumn is the name of the column to be parsed to identify variable types
	derivedFeatures = c('GenomicDensity', 'GenomicDistance')
	genomicFeatures = c('VirusDirect', derivedFeatures)
	
	extract_part <- function(vec, pos) sapply(strsplit(vec, '_'), '[', pos)
	
	extract_multi <- function(vec, positions) {
		parts <- extract_part(vec, positions)  # Will return a matrix if pos is a vector
		joined <- apply(parts, 2, function(x) paste(na.omit(x), collapse = '_'))  # Only collapse if there's more than one part
		
		as.vector(joined)
	}
	
	
	# Split variables into component parts:
	varnames <- data[[varnameColumn]] %>% 
		unique() %>% 
		as.character()
	
	result <- data.frame(Variable = varnames, stringsAsFactors = FALSE) %>% 
		mutate(VariableType = extract_part(Variable, 1),
					 IsDerivedFeature = VariableType %in% derivedFeatures,
					 IsGenomicFeature = VariableType %in% genomicFeatures,
					 Part2 = extract_part(Variable, 2),
					 Part3 = extract_part(Variable, 3),
					 Part4 = extract_part(Variable, 4)) %>% 
		
		mutate(Gene = ifelse(IsDerivedFeature, Part2, NA),
					 Location = ifelse(IsDerivedFeature, Part4,
					 									 ifelse(IsGenomicFeature, Part3, NA)),
					 Measure = ifelse(IsDerivedFeature, Part3,
					 								  ifelse(!is.na(Part2), Part2, VariableType))) %>% 
		select(-Part2, -Part3, -Part4, -IsDerivedFeature)
	
	# Identify genomic feature classes
	genomicFeatureTypes <- result %>% 
		filter(IsGenomicFeature) %>% 
		.$Measure %>% 
		id_genomic_features()
	
	# Join results and return
	result <- result %>% 
		left_join(genomicFeatureTypes, by = c('Measure' = 'featureName')) %>% 
		mutate(MeasureType = ifelse(is.na(MeasureType), VariableType, MeasureType))
	
	resultNames <- colnames(result)
	resultNames[resultNames == 'Variable'] <- varnameColumn
	colnames(result) <- resultNames
	
	data %>% 
		left_join(result, by = varnameColumn)
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Parse Feature sets -------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# This assumes id_variable_types() has already been called on data
id_feature_set <- function(data) {
	data %>% 
		mutate(SetName = case_when(Gene == 'ISG' ~ 'Similarity to ISGs',
															 Gene == 'Housekeeping' ~ 'Similarity to housekeeping genes',
															 Gene == 'Remaining' ~ 'Similarity to remaining genes',
															 is.na(Gene) & VariableType == 'PN' ~ 'Phylogenetic neighbourhood',
															 is.na(Gene) & VariableType == 'VirusDirect' ~ 'Viral genomic features'))
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- SHAP values -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_shap_contributions <- function(trainedModel, newdata, collapseOn, keepBias = F) {
	# Extract SHAP scores for individual vars and observations
	# - Based on code in xgb.plot.shap()
	# - newdata is optional - when not given, shap contributions for the training data will be
	#   returned
	# - collapseOn is optional. If given, collapse should be a vector of column names 
	#   specifying one or more grouping variables to collapse / summarise contributions to
	
	if (missing(newdata)) {
		newdata <- as.matrix(trainedModel$trainingData[, -1])
	}
	
	colnames(newdata) <- escape_special_chars(colnames(newdata))  # predict.xgb checks that column names are identical to model$feature_names
	newdata <- newdata[ , trainedModel$finalModel$feature_names]
	
	if (!missing(newdata)) {
		newdata <- as.matrix(newdata)
	}
	
	model <- trainedModel$finalModel
	
	shap <- xgboost:::predict.xgb.Booster(model, newdata, predcontrib = TRUE, approxcontrib = FALSE, predinteraction = FALSE)
	
	shap <- data.frame(shap, check.names = FALSE)
	Bias <- shap$BIAS
	
	shap <- shap %>% 
		mutate(RowID = 1:n()) %>% 
		select(-BIAS) %>% 
		gather(-RowID, key = 'Feature', value = 'SHAP') %>% 
		mutate(Feature = gsub('`', '', Feature))
	
	# Collapse to a single value per gene for each feature class (if requested)
	#   - Note that while SHAP values are additive (by definition), this is 
	#   merely an assumption about the structure of the model
	if (!missing(collapseOn)) {
		shap <- shap %>% 
			id_variable_types('Feature') %>% 
			group_by_at(c('RowID', collapseOn)) %>% 
			summarise(SHAP = sum(SHAP))
	}
	
	# Return:
	if (keepBias)
		attr(shap, 'Bias') <- Bias
	
	return(shap)
}


get_shap_interaction_vals <- function(trainedModel, newdata, keepBias = F) {
	# Extract SHAP interaction values for individual vars and observations
	# - newdata is optional - when not given, shap contributions for the training data will be
	#   returned
	
	if (missing(newdata)) {
		newdata <- as.matrix(trainedModel$trainingData[, -1])
	}
	
	colnames(newdata) <- escape_special_chars(colnames(newdata))  # predict.xgb checks that column names are identical to model$feature_names
	newdata <- newdata[ , trainedModel$finalModel$feature_names]
	
	if (!missing(newdata)) {
		newdata <- as.matrix(newdata)
	}
	
	model <- trainedModel$finalModel
	
	shap <- xgboost:::predict.xgb.Booster(model, newdata, predcontrib = FALSE, approxcontrib = FALSE, predinteraction = TRUE)
	
	shap <- data.frame(shap, check.names = FALSE)
	Bias <- shap$BIAS
	
	shap <- shap %>% 
		mutate(RowID = 1:n()) %>% 
		select(-BIAS) %>% 
		gather(-RowID, key = 'Feature', value = 'SHAP') %>% 
		mutate(Feature = gsub('`', '', Feature))
	
	# Collapse to a single value per gene for each feature class (if requested)
	#   - Note that while SHAP values are additive (by definition), this is 
	#   merely an assumption about the structure of the model
	if (!missing(collapseOn)) {
		shap <- shap %>% 
			id_variable_types('Feature') %>% 
			group_by_at(c('RowID', collapseOn)) %>% 
			summarise(SHAP = sum(SHAP))
	}
	
	# Return:
	if (keepBias)
		attr(shap, 'Bias') <- Bias
	
	return(shap)
}




get_varimp_shap <- function(trainedModel, collapseOn) {
	# Get overall variable importance, based on the SHAP values for individual observations.
	# The overall (model- or dataset-wide) importance of a variable is defined as the mean
	# magnitude of SHAP values across all data points
	# - colapseOn is optional (see get_shap_contributions())
	
	shapContribs <- get_shap_contributions(trainedModel, collapseOn = collapseOn)
	
	if (missing(collapseOn)) {
		shapContribs <- shapContribs %>% 
			group_by(Feature)
		
	} else {
		# Feature no longer present, so calculate importance for the requested groups instead
		shapContribs <- shapContribs %>% 
			group_by_at(collapseOn)
	}
	
	varImpShap <- shapContribs %>% 
		summarise(MeanSHAP = mean(abs(SHAP))) %>% 
		mutate(Rank = rank(-MeanSHAP, ties.method = 'min')) %>% 
		mutate(Rank = ifelse(MeanSHAP == 0, NA, Rank)) %>%   # Unused features should not be included in top 50, etc.
		arrange(MeanSHAP)
	
	varImpShap$Feature <- factor(varImpShap$Feature, levels = varImpShap$Feature)
	
	varImpShap
}



# Get obs vs expected from variable importance
# `groupingVars` is character vector of variable names to use for grouping
# `top` indicates the number of top features to compare
get_obs_vs_expected <- function(variableImportanceSummary, groupingVars, top = 50) {
	res <- variableImportanceSummary %>% 
		ungroup() %>% 
		mutate(Nvars = length(unique(Feature))) %>% 
		
		group_by_at(vars(one_of(groupingVars))) %>% 
		summarise(Nvars = unique(Nvars),
							N = n(),
							Expected = N/Nvars * top,   # Frequency of var x topN
							Observed = sum(Rank <= top, na.rm = T),
							ObsExp = Observed/Expected)
	
	if(sum(res$Observed) > top) 
		warning(paste('Number of features with ranks <=', top, 'is larger than', top, '- ranks and output should be checked.'))
	
	return(res)
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Clusters of features -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_cluster_members <- function(apresult) {
	# Convert an APResult object to a dataframe describing cluster membership
	
	clusterMembers <- lapply(apresult@clusters, function(x) tibble(Member = names(x)))
	names(clusterMembers) <- names(apresult@exemplars)
	
	clusterMembers <- bind_rows(clusterMembers, .id = 'Exemplar')
	
	# Get arrangement in dendrogram
	dendro <- as.dendrogram(apresult)
	labelOrder <- labels(dendro)
	
	clusterMembers$Member <- factor(clusterMembers$Member, levels = labelOrder)
	
	# Add cluster numbers and return:
	clusterMembers %>% 
		arrange(.data$Member) %>%   # Takes ordering from factor levels
		mutate(Exemplar = factor(.data$Exemplar, unique(.data$Exemplar)),
					 Cluster = as.numeric(.data$Exemplar)) 
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Supervised clustering ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
calculate_shapley_stacking <- function(individual_shapley_vals, summarised_shapley_vals, 
																			 bias, non_feature_cols) {
	# Calculate stacking data to make a plot similar to the 'supervised clustering' plot in 
	# Lundberg et al. 2018 (arXiv:1802.03888)
	# Use 'non_feature_cols' to list columns in individual_shapley_vals which should 
	# not be considered model features
	
	# Reshape and arrange contribution data for plotting:
	# - Features arranged by global importance (most important plotted first / closest to y = 0),
	#   but sign (direction) takes precedence over this when plotting 
	# - Thus, most important negative feature will be plotted after least important positive one,
	#   but still first among the negative features
	shap_contributions <- individual_shapley_vals %>% 
		mutate(Bias = bias) %>% 
		gather(-one_of(non_feature_cols), -.data$Bias, key = 'Feature', value = 'SHAP') %>% 
		group_by(.data$LatestSppName, .data$Strain) %>% 
		mutate(ModelOutput = unique(.data$Bias) + sum(.data$SHAP)) %>% #,
		#RelativeSize = abs(.data$SHAP) / sum(abs(.data$SHAP)),
		#Contribution = ModelOutput * RelativeSize) %>% 
		ungroup() %>% 
		
		left_join(summarised_shapley_vals, by = 'Feature') %>% 
		arrange(Rank)
	
	# Calculate stacking sizes:
	# - First arrange so largest features are closest to 'ModelOutput' line:
	# 		this means positive and negative features need to be sorted separately
	shap_positive <- shap_contributions %>% 
		filter(.data$SHAP >= 0) %>% 
		arrange(-.data$SHAP)
	
	shap_negative <- shap_contributions %>% 
		filter(.data$SHAP < 0) %>% 
		arrange(.data$SHAP)
	
	
	shap_contributions <- bind_rows(shap_positive, shap_negative) %>%   # negative above line, positive below
		mutate(ShapPositive = .data$SHAP >= 0) %>% 
		group_by(.data$LatestSppName, .data$Strain, .data$ShapPositive) %>% 
		mutate(upper = cumsum(abs(.data$SHAP)),
					 lower = upper - abs(.data$SHAP),
					 ymax = if_else(.data$ShapPositive, 
					 							 .data$ModelOutput - .data$upper, 
					 							 .data$ModelOutput + .data$lower),   # Somewhat counter-intuitively, positive effects decrease ymax (and visa versa), since we are showing what log odds would have been without this effect
					 ymin = if_else(.data$ShapPositive, 
					 							 .data$ModelOutput - .data$lower, 
					 							 .data$ModelOutput + .data$upper))
	
	shap_contributions
}
