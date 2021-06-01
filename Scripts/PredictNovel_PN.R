#!/usr/bin/env Rscript
#
# Predict zoonotic status based on the phylogenetic neighbourhood
#
suppressPackageStartupMessages({
	library(argparse)
})

parser <- ArgumentParser(description = paste("Predict zoonotic status from the phylogenetic neighbourhood.", 
																						 "NOTE THAT THIS MODEL HAS POOR PREDICTIVE ACCURACY", 
																						 "- this script is for illustrative purposes only"))

parser$add_argument("format", type = "character", metavar = "format", choices = c("genbank", "fasta"),
										help = "Format of supplied sequence data. Either 'genbank' or 'fasta'")

parser$add_argument("sequences", type = "character",
										help = "Path to a sequence file in either fasta or genbank flatfile format.")

parser$add_argument("metadata", type = "character",
										help = paste("CSV file with columns: Name, SequenceID, [CodingStart, CodingStop].",
																 
																 "SequenceID should match names in the sequence file (up to the",
																 "first space, e.g. 'KC660259.1' in the fasta sequence name",
																 "'>KC660259.1 Rabies virus isolate 10/564, complete genome').",
																 
																 "If sequences are supplied in fasta format, the final two columns",
																 "are required, and should specify the start and stop location of a",
																 "coding sequence (using 1-based indexing). To specify multiple",
																 "coding sequences in the same sequence, add additional lines,",
																 "repeating Name and SequenceID.",
																 
																 "Similarly, segmented viruses can be specified by repeating the",
																 "same Name for different SequenceIDs."))

parser$add_argument("out", type = "character",
										help = paste("base name for output files. Output will be in",
																 "{out}.predictions.csv"))

parser$add_argument("--n_models", type = "double", default = 100,
										help = paste("Number of models to use when bagging predictions - the top",
																 "'n_models' will be used (default: 100)."))

parser$add_argument("--cutoff", type = "double", default = 0.2991,
										help = paste("Probability cutoff to use when making binary predictions: Viruses",
																 "with an average probability >= to this will be labelled as zoonotic",
																 "(default: 0.2991)."))

parser$add_argument("--random_seed", type = "integer", default = trunc(runif(1, max = 1e5)),
										help = "Random seed to use (default: a random integer between 0 and 1e5)")

parser$add_argument("--nthread", type = "integer", default = 1,
										help = "Number of parallel threads allowed (default: 1)")


inputArgs <- parser$parse_args()


# Warn user that this might not be the right script
message("\n\nWARNING: the phylogenetic neighbourhood model is not recommended. ",
				"To get predictions from the best performing model, see ./PredictNovel.R\n\n")

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
	library(xgboost)
	library(caret)
	library(ModelMetrics)
	library(apcluster)
	library(ape)
	library(betacal)
})


RUN_ID <- "PN_LongRun"  										# RunID of trained models to use (only PN models implemented here)
MODEL_TYPE <- "xgbTree"       					    # Type of model used in this run
BLAST_FRAGMENTS <- inputArgs$nthread * 4    # Number of parallel/sequential blast searches (see blastn_parallel())
																				    # Currently, 4 fragments processed on each thread

ROOT_DIR <- find_rstudio_root_file()

source(file.path(ROOT_DIR, "Utils", "xgboost_utils.R"))
source(file.path(ROOT_DIR, "Utils", "subset_utils.R"))
source(file.path(ROOT_DIR, "Utils", "prediction_utils.R"))
source(file.path(ROOT_DIR, 'Utils', 'calibration_utils.R'))
source(file.path(ROOT_DIR, 'Utils', 'rundata_utils.R'))

source(file.path(ROOT_DIR, 'Utils', 'blast_utils.R'))
source(file.path(ROOT_DIR, 'Utils', 'PNsummary_utils.R'))


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


# Training data and sequences
trainingData <- file.path(ROOT_DIR, "CalculatedData", "SplitData_Training.rds") %>%
 	readRDS() %>%
 	ungroup()

training_seqs <- read.fasta(file.path('ExternalData', 'Sequences', 'CombinedSequences.fasta'), as.string = T)
training_seqs <- subset_sequences(training_seqs, unique(trainingData$Accessions))



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


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Parse out full genome sequences ------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("\n\nReading sequences")

if (inputArgs$format == "fasta") {
	novel_seqs <- read.fasta(inputArgs$sequences, as.string = T)
	
} else {
	# If input sequences are in genbank format, need to extract these to fasta
	# Not currently needed, so this is not implemented
	stop("Only fasta file input supported (genbank extraction not implemented)")
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Simplify metadata --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Need to gather all sequence IDs associated with each virus species

novel_virus_data <- metaData %>% 
	group_by(.data$Name) %>% 
	summarise(Accessions = paste(unique(.data$SequenceID), collapse = "; ")) %>% 
	ungroup() %>% 
	rename(UniversalName = .data$Name)

combined_data <- bind_rows(trainingData, novel_virus_data)  # used for name lookups


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Identify top models ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

trainedAccuracy <- test_preds %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True') %>% 
	group_by(.data$Iteration) %>% 
	summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$RawScore)) %>% 
	ungroup() %>% 
	mutate(Rank = rank(-.data$AUC, ties.method = "random")) 

topModels <- trainedAccuracy %>% 
	filter(.data$Rank <= inputArgs$n_models) %>% 
	pull(.data$Iteration)

bagModels <- trainedModels[topModels]


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

message("\n\nRunning blast search")

# Blast novel viruses against all training viruses
blast_cache <- paste(inputArgs$out, "cached_blast_searches", sep = ".")
dir.create(blast_cache, recursive = TRUE)
	
full_blast_results <- cached_blast(querySeqs = novel_seqs, 
																	 dbSeqs = training_seqs,
																	 max_target_seqs = 10000,  # Very high, so this param does not become the limiting factor for matches
																	 nfragments = BLAST_FRAGMENTS, 
																	 nthreads = inputArgs$nthread, 
																	 cache_dir = blast_cache)


message("\n\nCalculating PN features")

# Set up default data and options for PN summaries:
pn_defaults <- list(pIdentityCutoffs = c(0, 60, 70, 80),
										maxNeighbours = 5,  
										positiveName = TRUE,
										nthreads = inputArgs$nthread,
										minimal = TRUE)
	
	
# Main logic for calculating PN features:
get_pn_features <- function(queryAccessions, dbData, removeSelf = FALSE,
														blastResults = full_blast_results, allSequences = training_seqs, fullData = combined_data, 
														hostPhylogeny = HostPhylo, 
														reservoirPhylogeny = ReservoirPhylo,
														pnDefaults = pn_defaults) {
		
	# Extract relevant sequences
	train_seq_subset <- subset_sequences(allSequences, unique(dbData$Accessions))
	
	# Correct blast results for the new database
	queryAccessions <- strsplit(queryAccessions, split = '; ') %>% 
		unlist()
	
	corrected_blast_results <- blastResults %>% 
		filter(queryID %in% queryAccessions) %>% 
		subset_blast(full_db_sequences = allSequences, 
								 subset_db_sequences = train_seq_subset)
	
	
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



## Calculate PN features matching each training set
#   - Only need to do this for the best models
novel_features <- lapply(bagModels, function(x) {
	train_subset <- filter(trainingData, .data$LatestSppName %in% x$trainingNames)
	get_pn_features(queryAccessions = novel_virus_data$Accessions, dbData = train_subset)
})


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Predict ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("\n\nMaking predictions\n")

novel_names <- lapply(novel_features, function(x) x$queryUniversalName)  # Can't guarantee that order is constant

predictionData <- mapply(FUN = prepare_prediction_data, 
												 fit = bagModels,
												 newdata = novel_features, 
												 MoreArgs = list(labelcol = "queryUniversalName"),
												 SIMPLIFY = FALSE)

allPredictions <- mapply(FUN = get_prediction, 
												 fit = bagModels, 
												 newdata = predictionData,
												 virusnames = novel_names,
												 MoreArgs = list(model_type = MODEL_TYPE),
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
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
message("Finalizing output\n")

write_excel_csv(baggedPredictions, path = sprintf("%s.predictions.csv", inputArgs$out))

