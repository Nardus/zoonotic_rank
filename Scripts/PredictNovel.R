#!/usr/bin/env Rscript
#
# Predict zoonotic status for novel viruses
#
suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Predict zoonotic status from genome sequences.")

parser$add_argument("format", type = "character", metavar = "format", choices = c("genbank", "fasta"),
                    help = "Format of supplied sequence data. Either 'genbank' or 'fasta'")

parser$add_argument("sequences", type = "character",
                    help = "Path to a sequence file in either fasta or genbank flatfile format.")

parser$add_argument("metadata", type = "character",
                    help = paste("CSV file with columns: Name, SequenceID, [CodingStart, CodingStop].",
                    						 
                    						 "The first two columns are always required.",

                                 "SequenceID should match names in the sequence file (up to the",
                                 "first space, e.g. 'KC660259.1' in the fasta sequence name",
                                 "'>KC660259.1 Rabies virus isolate 10/564, complete genome').",

                                 "If sequences are supplied in fasta format, the CodingStart and",
                    						 "CodingStop columns are required, and should specify the start and",
                                 "stop location of a coding sequence (using 1-based indexing).",
                                 "To specify multiple coding sequences in the same sequence, add",
                                 "additional lines, repeating Name and SequenceID.",

                                 "Similarly, segmented viruses can be specified by repeating the",
                                 "same Name for different SequenceIDs.",
                    						 
                    						 "To specify genes on the complementary strand, reverse the",
                    						 "coordinates (i.e. CodingStart > CodingStop)."))

parser$add_argument("out", type = "character",
                    help = paste("base name for output files. Output generated are",
                                 "{out}.predictions.csv (the actual predictions),",
                                 "{out}.explanations.csv (contributions of individual features to",
                                 "these estimates, for the individual models used in bagging),",
                                 "{out}.feature_dendrogram.pdf (dendrogram showing distances",
                                 "between input viruses and the training data, in terms of genome",
                                 "features), and",
                                 "{out}.similar_explanations.csv / {out}.similar_explanations.pdf",
                                 "(viruses from the training data with",
                                 "similar explanations to the submitted viruses, obtained using",
                                 "affinity propagation clustering based on the median of influence",
                                 "(TreeSHAP) values for each feature in across models. Clusters can",
                                 "be identified in this file as all viruses sharing the same exemplar."))

parser$add_argument("--exclude", type = "character",
										help = paste("Comma-separated list of viruses which should NOT be present in",
																 "the training data of models used for prediction. Names should",
																 "match version 2018b of ICTV taxonomy. Note that excluding even 1",
																 "virus can severely decrease the number of models available for", 
																 "averaged predictions."))

parser$add_argument("--n_models", type = "double", default = 100,
                    help = paste("Number of models to use when bagging predictions - the top",
                                 "'n_models' will be used (default: 100)."))

parser$add_argument("--cutoff", type = "double", default = 0.293,
										help = paste("Probability cutoff to use when making binary predictions: Viruses",
																 "with an average probability >= to this will be labelled as zoonotic",
																 "(default: 0.293)."))

parser$add_argument("--random_seed", type = "integer", default = trunc(runif(1, max = 1e5)),
                    help = "Random seed to use (default: a random integer between 0 and 1e5)")


inputArgs <- parser$parse_args()

# Check input:
if (!file.exists(inputArgs$sequences))
  stop(sprintf("%s not found", inputArgs$sequences))

if (!file.exists(inputArgs$metadata))
  stop(sprintf("%s not found", inputArgs$metadata))


message("Random seed set to ", inputArgs$random_seed)
set.seed(inputArgs$random_seed)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load utility functions and data ------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
suppressPackageStartupMessages({
  library(rprojroot)
  library(seqinr)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(readr)
	library(ggplot2)
  library(xgboost)
  library(caret)
	library(ModelMetrics)
  library(apcluster)
  library(ape)
	library(EnvStats)
})


RUN_ID <- "AllGenomeFeatures_LongRun"  	# RunID of trained models to use (only genome feature and mimicry-based runs implemented)
USE_DERIVED_FEATURES <- TRUE  					# Whether derived (mimicry) features should be calculated
MODEL_TYPE <- "xgbTree"       					# Type of model used in this run


ROOT_DIR <- find_rstudio_root_file()
genome_feature_script <- file.path(ROOT_DIR, "Utils", "GenomeFeatures.py")

source(file.path(ROOT_DIR, "Utils", "cds_parser.R"))
source(file.path(ROOT_DIR, "Utils", "xgboost_utils.R"))
source(file.path(ROOT_DIR, "Utils", "plot_utils.R"))
source(file.path(ROOT_DIR, "Utils", "prediction_utils.R"))
source(file.path(ROOT_DIR, 'Utils', 'calibration_utils.R'))
source(file.path(ROOT_DIR, 'Utils', 'rundata_utils.R'))

if (USE_DERIVED_FEATURES)
	source(file.path('Utils', 'derived_genome_feature_utils.R'))


# Trained models:
trainedModels <- paste0(RUN_ID, "_ModelFits.rds") %>%
  file.path(ROOT_DIR, "RunData", RUN_ID, .) %>%
  readRDS()

if (inputArgs$n_models > length(trainedModels))
  stop(sprintf("n_models is too high: only %d models have been trained in this model set", 
               length(trainedModels)))

# Predictions from trained models (used to assess accuracy and for calibration)
all_model_preds <- paste0(RUN_ID, "_Predictions.rds") %>% 
	file.path(ROOT_DIR, "RunData", RUN_ID, .) %>%
	readRDS()

test_preds <- all_model_preds %>% 
	filter(Dataset == 'test')

calibration_preds <- all_model_preds %>% 
	filter(Dataset == 'calibration')


# Training data and features
trainingData <- file.path(ROOT_DIR, "CalculatedData", "SplitData_Training.rds") %>%
  readRDS() %>%
  ungroup()

trainingFeatures <- file.path(ROOT_DIR, "CalculatedData", "GenomicFeatures-Virus.rds") %>%
	readRDS() %>%
  rename_at(vars(-UniversalName, -Strain), ~ paste0("VirusDirect_", .)) %>%
  filter(paste(UniversalName, Strain) %in% paste(trainingData$UniversalName, trainingData$Strain))


if (USE_DERIVED_FEATURES) {
	# Data needed to calculate derived features:
	humanFeatures <- readRDS(file.path('CalculatedData', 'GenomicFeatures-HumanCombined.rds'))
	
	# Data needed to identify training viruses by feature value:
	virusDistFeatures <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Distances.rds'))
	trainingFeatures <- trainingFeatures %>% 
		left_join(virusDistFeatures, by = c('UniversalName', 'Strain'))
}


# Metadata for novel viruses:
message("\n\nReading metadata")

expectedCols <- c("Name", "SequenceID", "Start", "Stop")
colTypes <- list(Name = "c", SequenceID = "c", Start = "d", Stop = "d")

if (inputArgs$format == "genbank") {
  expectedCols <- expectedCols[1:2]
  colTypes <- colTypes[1:2]
}

colTypes <- do.call(cols, colTypes)

# Ignore actual header so column names are predictable (but assume there is one):
metaData <- read_csv(inputArgs$metadata, col_names = expectedCols,
                     col_types = colTypes, skip = 1)


if (!is.null(inputArgs$exclude)) {
	excluded_spp <- str_split(inputArgs$exclude, pattern = ",")[[1]]
	
	if (!all(excluded_spp %in% trainingData$LatestSppName)) {
		not_found <- excluded_spp[!excluded_spp %in% trainingData$LatestSppName]
		
		stop("Some species marked for exclusion were either not in the training data or their names are invalid:\n\t", 
						paste(unique(not_found), collapse = "\n\t"))
	}
}
	

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Parse out coding sequences -----------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if (inputArgs$format == "fasta") {
    message("\nParsing coding sequences from fasta file\n")

  sequences <- read.fasta(inputArgs$sequences)

  if (!all(metaData$SequenceID %in% names(sequences)))
    stop(paste("Not all sequence IDs were found in the supplied fasta file. Ensure the second",
               "column of the supplied metadata contains valid IDs (see PredictNovel.R --help)"))

  sequences <- sequences[metaData$SequenceID]
  stopifnot(length(sequences) == nrow(metaData))

  coding <- mapply(extract_cds,
                   sequence = sequences,
                   start_coordinate = metaData$Start,
                   stop_coordinate = metaData$Stop,
  								 MoreArgs = list(allow_complementary = TRUE))

  # Output: Each CDS simply gets the same name as the parent sequence (but a different SeqID)
  new_ids <- paste(names(coding), seq_len(length(coding)), sep = "_")
  new_metadata <- metaData %>%
    mutate(SequenceID = new_ids)
  
  names(coding) <- new_ids
  
  temp_coding_file <- tempfile()

  write.fasta(coding, names = names(coding), file.out = temp_coding_file)
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Simplify metadata --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Always creating a tempfile here, since GenomeFeatures.py expects specific column names
temp_metadata_file <- tempfile()

metaData <- metaData %>% 
  distinct(Name, SequenceID)

write.csv(metaData, temp_metadata_file)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Calculate genome features ------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("\nCalculating genome features")
temp_out_coding <- tempfile()
temp_out_entire <- tempfile()

if (inputArgs$format == "fasta") {
  temp_extra_metadata_file <- tempfile()
  write.csv(new_metadata, temp_extra_metadata_file)
  
  cmd_string_coding <- paste(genome_feature_script, temp_extra_metadata_file, temp_coding_file, 
  													 "--out", temp_out_coding)
  
  cmd_string_entire <- paste(genome_feature_script, temp_metadata_file, inputArgs$sequences,
  													 "--noncoding", "--out", temp_out_entire)
} else {
  # Rely on GenomeFeatures.py to extract coding sequences from the genbank file
	cmd_string_coding <- paste(genome_feature_script, temp_metadata_file, inputArgs$sequences,
														 "--extract", "--out", temp_out_coding)
	
	cmd_string_entire <- paste(genome_feature_script, temp_metadata_file, inputArgs$sequences,
														 "--extract", "--noncoding", "--out", temp_out_entire)
}


# Run these calls:
system2("python3", cmd_string_coding)
system2("python3", cmd_string_entire)


# Read output
# Currently CPB_machine"s output has some issues:
# (all rows except the header end in a tab, causing data to be shifted relative to column names)
message("\nParsing genome features\n")

features_coding <- read.delim(temp_out_coding, na.strings = c("?", "NA"), row.names = NULL,
															stringsAsFactors = F)
colnames(features_coding) <- c(colnames(features_coding)[-1], "XX_REMOVE")
features_coding <- features_coding[, colnames(features_coding) != "XX_REMOVE"]


features_entire <- read.csv(temp_out_entire)


# Clean up tempfiles:
# - Capturing return code to preven errors when a file is missing
if (inputArgs$format == "fasta") {
  unlink(temp_coding_file)
  unlink(temp_extra_metadata_file)
}


unlink(temp_metadata_file)
unlink(temp_out_coding)
unlink(temp_out_entire)


# Merge feature sets
features_coding <- features_coding %>% 
	select(-TaxID, -Species) %>% 
	rename_at(vars(-SeqName), ~ paste(., 'Coding', sep = '_'))

features_entire <- features_entire %>% 
	rename_at(vars(-SeqName), ~ paste(., 'EntireSeq', sep = '_'))

all_direct_features <- features_coding %>% 
	full_join(features_entire, by = 'SeqName')

# Rename direct features to their final names: 
genomeFeatures <- all_direct_features %>% 
	rename_at(vars(-SeqName), ~ paste("VirusDirect", ., sep = '_')) %>% 
	rename(Name = SeqName)


# Calculate derived genome features ('distance' from human genes) if needed
if (USE_DERIVED_FEATURES) {
	feature_names <- colnames(humanFeatures)
	feature_names <- feature_names[! feature_names %in% c('GeneID', 'meanCPM', 'TranscriptID', 'GeneSet')]
	
	all_direct_features <- all_direct_features %>% 
		mutate(UniversalName = SeqName,
					 Strain = NA_character_)  # These columns expected by get_feature_dists()
	
	isg_features <- humanFeatures %>% 
		filter(GeneSet == 'ISG')
	
	isg_features <- get_feature_dists(virusFeatures = all_direct_features, 
																		geneFeatures = isg_features,
																		setprefix = 'ISG',
																		featureColNames = feature_names)
	
	housekeeping_features <- humanFeatures %>% 
		filter(GeneSet == 'Housekeeping')
	
	housekeeping_features <- get_feature_dists(virusFeatures = all_direct_features, 
																						 geneFeatures = housekeeping_features,
																						 setprefix = 'Housekeeping',
																						 featureColNames = feature_names)
	
	remaining_features <- humanFeatures %>% 
		filter(GeneSet == 'Remaining')
	
	remaining_features <- get_feature_dists(virusFeatures = all_direct_features, 
																					geneFeatures = remaining_features,
																					setprefix = 'Remaining',
																					featureColNames = feature_names)
	
	# Join:
	all_indirect_features <- isg_features %>% 
		full_join(housekeeping_features, by = c('UniversalName', 'Strain')) %>% 
		full_join(remaining_features, by = c('UniversalName', 'Strain')) %>% 
		select(-Strain)
	
	genomeFeatures <- genomeFeatures %>% 
		left_join(all_indirect_features, by = c('Name' = 'UniversalName'))
}




# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Predict ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("\n\nMaking predictions\n")

# Find top models:
# - The best set of models changes if exclusions are requested:
#    In that case, 
#      (1) only models with none of the excluded viruses in the training or optimisation data 
#          are considered
#      (2) accuracy is calculated without considering how well the model performed on any of 
#          the excluded viruses (which will be in the test set of these models)
if (!is.null(inputArgs$exclude)) {
	virus_names <- trainingData %>% 
		distinct(UniversalName, LatestSppName)
	
	test_preds <- test_preds %>% 
		left_join(virus_names, by = "UniversalName")
	
	valid_models <- test_preds %>% 
		group_by(Iteration) %>% 
		summarise(keep = all(excluded_spp %in% LatestSppName)) %>% 
		filter(keep) %>% 
		pull(Iteration)
	
	if (length(valid_models) < inputArgs$n_models)
		stop("Not enough models remain after excluding the requested viruses. ", 
				 "Adjust either 'n_models' or the exclusion list.")
	
	cat(length(valid_models), "models remain after excluding requested viruses\n")
	
	test_preds <- test_preds %>% 
		filter(Iteration %in% valid_models) %>% 
		filter(!LatestSppName %in% excluded_spp) %>% 
		select(-LatestSppName)
}


trainedAccuracy <- test_preds %>% 
	mutate(InfectsHumans = InfectsHumans == 'True') %>% 
	group_by(Iteration) %>% 
	summarise(AUC = auc(actual = InfectsHumans, predicted = RawScore)) %>% 
	ungroup() %>% 
  mutate(Rank = rank(-AUC, ties.method = "random")) 

topModels <- trainedAccuracy %>% 
  filter(Rank <= inputArgs$n_models) %>% 
  pull(Iteration)
  
bagModels <- trainedModels[topModels]

# Predict:
predictionData <- lapply(bagModels, prepare_prediction_data, newdata = genomeFeatures)

allPredictions <- mapply(get_prediction, fit = bagModels, newdata = predictionData,
                         MoreArgs = list(virusnames = genomeFeatures$Name,
                         								 model_type = MODEL_TYPE),
                         SIMPLIFY = FALSE, USE.NAMES = TRUE) %>%
  bind_rows()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Calibrate predictions ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
allPredictions <- allPredictions %>% 
	select(.data$Iteration, .data$Name, RawScore = .data$True) %>% 
	group_by(.data$Iteration) %>% 
	group_modify(~ calibrate_preds(.x, calibration_preds = calibration_preds))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Bagging & priority classification ----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Classify by priority
#  - Very high: entire CI above cutoff
#  - High: Median and upper limit above cutoff, but CI crosses cutoff
#  - Medium: Only upper limit is above cutoff, i.e. still potentially zoonotic, but less likely
#  - Low: entire CI below cutoff, i.e. very unlikely to be zoonotic
prioritize <- function(lower_bound, upper_bound, median, cutoff) {
	stopifnot(length(cutoff) == 1)
	stopifnot(cutoff > 0 & cutoff < 1)
	
	p <- if_else(lower_bound > cutoff, 'Very high',
							 if_else(median > cutoff, 'High',
							 				if_else(upper_bound > cutoff, 'Medium',
							 								'Low')))
	factor(p, levels = c('Low', 'Medium', 'High', 'Very high'))
}

baggedPredictions <- allPredictions %>% 
  group_by(Name) %>% 
  summarise(calibrated_score_mean = mean(CalibratedScore),
  					calibrated_score_lower = quantile(CalibratedScore, probs = 0.05/2),
  					calibrated_score_upper = quantile(CalibratedScore, probs = 1 - 0.05/2)) %>% 
  mutate(bagged_prediction = calibrated_score_mean >= inputArgs$cutoff,
  			 priority_category = prioritize(lower_bound = calibrated_score_lower,
  			 															  upper_bound = calibrated_score_upper,
  			 															  median = calibrated_score_mean,
  			 															  cutoff = inputArgs$cutoff)) %>%
  ungroup()



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Find viruses with similar explanations -----------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("Finding viruses with similar explanations\n")

# Explanations
bagModelData <- lapply(bagModels, function(x) x$trainingData)

trainExplanations <- mapply(get_shap_contributions, 
														trainedModel = bagModels, 
														newdata = bagModelData, 
														SIMPLIFY = FALSE, USE.NAMES = TRUE)

novelExplanations <- mapply(get_shap_contributions, 
                            trainedModel = bagModels, 
                            newdata = predictionData, 
                            SIMPLIFY = FALSE, USE.NAMES = TRUE)

## Add names:
#   - Each training set contains a different set of viruses
trainNames <- lapply(trainedModels, function(x) data.frame(RowID = 1:length(x$trainingNames),
																													 Name = x$trainingNames,
																													 stringsAsFactors = FALSE)) %>% 
	bind_rows(.id = 'Iteration')

#		- Names for new viruses: need to make sure there are no name clashes with training viruses
novelNames <- genomeFeatures %>%
  ungroup() %>%
  select(Name) %>%
  mutate(RowID = 1:n(),
         Name = paste(Name, "NOVEL", sep = "__"))  # Need to avoid name clashes between training and novel viruses


#		- Add to explanation data
trainExplanations <- trainExplanations %>%
  bind_rows(.id = "Iteration") %>%
  left_join(trainNames, by = c("Iteration", "RowID"))

novelExplanations <- novelExplanations %>%
  bind_rows(.id = "Iteration") %>%
  left_join(novelNames, by = "RowID")


## Cluster:
#		- Summarising explanations across iterations for each virus, matching what was done in figure 2
allExplanations <- trainExplanations %>%
  bind_rows(novelExplanations) %>%
	group_by(Name, Feature) %>% 
	summarise(SHAP_mean = mean(SHAP)) %>%  # Keeping directionality: if a feature changes sign for a single virus we do want it to impact effect size measured
	ungroup() %>% 
	spread(key = Feature, value = SHAP_mean)

shapMatrix <- allExplanations %>%
  select(-Name) %>%
  as.matrix()

rownames(shapMatrix) <- allExplanations$Name

clusters <- apcluster(negDistMat(r = 2), x = shapMatrix)


# Summarise clusters:
extract_clusters <- function(cluster, exemplar) {
  data.frame(Name = names(cluster),
             Exemplar = exemplar,
             stringsAsFactors = FALSE)
}

clusterMembers <- mapply(extract_clusters, cluster = clusters@clusters, exemplar = names(clusters@exemplars),
                         SIMPLIFY = FALSE) %>%
  bind_rows()

clusterNumbers <- clusterMembers %>%
  distinct(Exemplar) %>% 
  mutate(Cluster = 1:n())

clusterMembers <- clusterMembers %>%
  left_join(clusterNumbers, by = "Exemplar")

keepClusters <- clusterMembers %>%
    filter(str_ends(Name, "__NOVEL")) %>%
    .$Cluster

keepClusters <- clusterMembers %>%
  filter(Cluster %in% keepClusters)


# Finally, add back accessions so similar viruses are easier to identify, and clean up:
accessions <- trainingData %>% 
	select(Name = LatestSppName, Accessions)

similarViruses <- keepClusters %>%
  left_join(accessions, by = "Name") %>%
  mutate(Name = str_remove(Name, "__NOVEL$")) %>%   # Remove 'NOVEL' tags from user viruses (used to avoid name clashes above)
  arrange(Exemplar, Name) %>%
  select(Cluster, Exemplar, Name, Accessions)  # Arrange columns in a more logical order 


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Find distance to dataset in general --------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Above step may fail to find viruses with the exact *combination* of explanations triggered. This
# is not neccesarily bad, but we need to check that the virus is not too unusual compared to the
# training data to know how much to trust the predictions.
# What matters here is not the pure genetic distance, but distance in terms of genome features.

# TODO: unclear *which* genome features to use here - currently using direct only, but it may also be valid to use either those used in training or all possible features

message("Checking distance to training dataset\n")
latest_spp_names <- trainingData %>% 
	distinct(UniversalName, Strain, LatestSppName)

trainingFeatures <- trainingFeatures %>% 
	left_join(latest_spp_names, by = c("UniversalName", "Strain")) %>% 
	rename(Name = LatestSppName) %>% 
	select(-UniversalName, -Strain)

genomeFeatures <- genomeFeatures[, colnames(genomeFeatures) %in% colnames(trainingFeatures)]
allFeatures <- rbind(trainingFeatures, genomeFeatures)

allFeaturesMat <- allFeatures %>%
	select(starts_with('VirusDirect_'))

rownames(allFeaturesMat) <- allFeatures$Name

hc <- allFeaturesMat %>%
  scale() %>%       # Scaled to ensure all features have equal influence on the distance
  dist() %>%
  hclust(method = "ward.D2")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("Finalizing output\n")
outDir <- dirname(inputArgs$out)

if (!dir.exists(outDir))
  dir.create(outDir, recursive = TRUE)


## Genome features --------------------------------------------------------------------------------
write_excel_csv(genomeFeatures, path = sprintf("%s.genome_features.csv", inputArgs$out))


## Predictions ------------------------------------------------------------------------------------
write_excel_csv(baggedPredictions, path = sprintf("%s.predictions.csv", inputArgs$out))


## Explanations (at a per-model, per feature level): ----------------------------------------------
trainedAccuracy <- trainedAccuracy %>%
  mutate(Iteration = as.character(Iteration))

novelExplanations <- novelExplanations %>%
	mutate(Name = str_remove(Name, "__NOVEL$")) %>%
  
  left_join(trainedAccuracy, by = "Iteration") %>%
  id_variable_types("Feature") %>%
  select(Name,
         Model = Iteration, Model_AUC = AUC, Model_Rank = Rank,
         Feature, Variable_Type = VariableType, Gene, Measure, Measure_Type = MeasureType, Location,
         SHAP) %>% 
  arrange(Model_Rank, Name, Feature)


write_excel_csv(novelExplanations, path = sprintf("%s.explanations.csv", inputArgs$out))


## Similar viruses, based on median SHAP values across bagged models --------------------------
write_excel_csv(similarViruses, path = sprintf("%s.similar_explanations.csv", inputArgs$out))


## Dendrogram of these clusters: ------------------------------------------------------------------
save_tree <- function(tree, labelTips, filename, cex = 0.5, width = 14, height = 14) {
  ## Save a phylo object, highlighting specific tips
  # Remove tip labels for internal viruses:
  #tree$tip.label <- if_else(tree$tip.label %in% labelTips, tree$tip.label, "")

  # Colour tips of input viruses:
  branchCols <- rep("black", Nedge(tree))
  inputVirusTips <- which(tree$tip.label %in% labelTips)
  branchCols[which(tree$edge[, 2] %in% inputVirusTips)] <- "red"

  # Save
  pdf(filename, width = width, height = height)
    plot(tree, type = "fan", cex = cex, underscore = TRUE,
       edge.color = branchCols)
  invisible(dev.off())
}

aptree <- clusters %>%
  as.hclust() %>%
  as.phylo() %>%
  ladderize()

save_tree(aptree, cex = 0.8,
      labelTips = paste("Cluster", keepClusters$Cluster),
      filename = sprintf("%s.similar_explanations.pdf", inputArgs$out))


## Tree dendrogram showing relationship to training data: -----------------------------------------
hctree <- as.phylo(hc) %>%
  ladderize()

save_tree(hctree, labelTips = genomeFeatures$Name, cex = 0.2,
      filename = sprintf("%s.feature_dendrogram.pdf", inputArgs$out))


## TEMP: All cluster members ----------------------------------------------------------------------
write_excel_csv(clusterMembers, path = sprintf("%s.all_cluster_members.csv", inputArgs$out))


warnings()