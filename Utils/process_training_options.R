##
## Part of Zoonosis prediction pipeline:
## 	- Process input options passed to training script, loading optional feature data as needed
## 	
## 

library(argparse)

parser <- ArgumentParser(description = 'Train a model with specific feature sets. Validation is against subsets of the test set only - this script does not see the holdout data.')
parser$add_argument('RandomSeed', type = 'integer', 
										help = 'a random seed (used by both R and xgboost)')

parser$add_argument('RunID', type = 'character', 
										help = 'an identifier for this run, used as a base to name output files')


# Optional arguments:
# Use publication count to weight samples?
parser$add_argument('--useWeights', action = 'store_const', const = TRUE, default = FALSE,
										help = 'Use publication count to weight samples')

# Feature types to include in training:
featureGroup <- parser$add_argument_group('Feature sets')
featureGroup$add_argument('--includePN', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include features summarising the proportions of zoonotic and vector-borne viruses in the phylogenetic neigbourhood')
featureGroup$add_argument('--includeVirusFeatures', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include direct genomic features (calculated directly from virus genomes)')
featureGroup$add_argument('--includeISG', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include genomic features relative to human ISGs')
featureGroup$add_argument('--includeHousekeeping', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include genomic features relative to human housekeeping genes')
featureGroup$add_argument('--includeRemaining', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include genomic features relative to remaining human genes')
featureGroup$add_argument('--includeProteinMotifs', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include protein motif features')
featureGroup$add_argument('--includeTaxonomy', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include virus taxonomy (family and proportion zoonotic in family)')


# Options to control training
trainGroup <- parser$add_argument_group('Training options')
trainGroup$add_argument('--trainProportion', metavar = 'p', type = 'double', default = 0.7, 
												help = 'proportion of virus species to select for each training set - all remaining virus species go to test set (default: 0.85)')
trainGroup$add_argument('--topFeatures', metavar = 'f', type = 'double', default = 100, 
												help = 'Number of features to retain for final training')
trainGroup$add_argument('--nboot', metavar = 'b', type = 'integer', default = 100, 
												help = 'number of iterations to perform (default: 100)')
trainGroup$add_argument('--nseeds', metavar = 's', type = 'integer', default = 10, 
												help = 'number of xgboost random seeds to use (default: 10). If e.g. nboot = 100 and nseeds = 10, there will be 10 iterations using each seed.')

trainGroup$add_argument('--nthread', metavar = 't', type = 'integer', default = 1, 
												help = 'number of parallel threads allowed (default: 1)')


# Parse input:
# - These may be specified by the calling script (as a vector of strings named OVERRIDE_EXTERNAL_COMMANDS),
#   but will be taken from arguments used to invoke the calling script if this is left unspecified
if (exists("OVERRIDE_EXTERNAL_COMMANDS")) {
	INPUT <- parser$parse_args(OVERRIDE_EXTERNAL_COMMANDS)
} else {
	INPUT <- parser$parse_args()
}

if (INPUT$trainProportion <= 0 | INPUT$trainProportion >= 1) 
	stop('trainProportion should be between 0 and 1') # But exactly 0 or 1 makes no sense either


if (! any(INPUT$includePN, INPUT$includeVirusFeatures, INPUT$includeISG, 
					INPUT$includeHousekeeping, INPUT$includeRemaining, INPUT$includeProteinMotifs, 
					INPUT$includeTaxonomy))
	stop('At least one of the feature set flags must be set. See TrainAndValidate.R --help')


if (INPUT$useWeights)
	stop('Weighting training data by publication count is no longer implemented - will need changes to TrainAndValidate.R')

if (INPUT$nboot <= 0)
	stop('nboot must be positive')

if (INPUT$nseeds <= 0)
	stop('nseeds must be positive')

if (INPUT$nseeds > INPUT$nboot)
	stop('nseeds must be less than or equal to nboot')

if (INPUT$nboot %% INPUT$nseeds != 0)
	warning('nseeds is not a multiple of nboot. nboot will be increased')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Core data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
set.seed(INPUT$RandomSeed)

InputData <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load/Add optional data based on input options ----------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FinalData <- InputData  # Will be modified below

# Publication counts: Used for weighting cases (if requested):
if (INPUT$useWeights) {
	PublicationCounts <- read.csv(file.path('InternalData', 'PublicationCounts.csv'))
	
	FinalData <- InputData %>% 
		left_join(PublicationCounts, by = c('LatestSppName' = 'VirusName'))
}


if (INPUT$includePN) {
	# Load sequences
	AllSeqs <- read.fasta(file.path('ExternalData', 'Sequences', 'CombinedSequences.fasta'), as.string = T)
}


if (INPUT$includeVirusFeatures) {
	# Virus features (created by CalculateGenomicFeatures.R)
	VirusFeaturesDirect <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Virus.rds'))
	
	VirusFeaturesDirect <- VirusFeaturesDirect %>% 
		rename_at(vars(-UniversalName, -Strain), ~ paste('VirusDirect', ., sep = '_'))
	
	FinalData <- FinalData %>% 
		left_join(VirusFeaturesDirect, by = c('UniversalName', 'Strain'))
}


if (INPUT$includeISG | INPUT$includeHousekeeping | INPUT$includeRemaining) {
	# Genomic features (created by CalculateGenomicFeatures.R)
	GenomicDistances <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Distances.rds'))
	
	if (!INPUT$includeISG) GenomicDistances <- select(GenomicDistances, -contains('_ISG_'))
	if (!INPUT$includeHousekeeping) GenomicDistances <- select(GenomicDistances, -contains('_Housekeeping_'))
	if (!INPUT$includeRemaining) GenomicDistances <- select(GenomicDistances, -contains('_Remaining_'))
	
	FinalData <- FinalData %>% 
		left_join(GenomicDistances, by = c('UniversalName', 'Strain'))
}


if (INPUT$includeProteinMotifs) {
	# This option not fully implemented:
	warning('Inclusion of protein motifs requested - check that these have been subjected to screening in SelectFeatures.R. Otherwise, this option will have no effect.')
	
	MotifBias <- read.csv(file.path('CalculatedData', 'ProteinMotifs_Bias.csv'))
	MotifDistance <- read.csv(file.path('CalculatedData', 'ProteinMotifs_Distance.csv'))
	
	MotifBias <- MotifBias %>% 
		select(-Nmatches, -Npossible, -ExpectedFrequency) %>% 
		mutate(ELM_Accession = paste('MotifBias', ELM_Accession, sep = '_')) %>% 
		spread(key = 'ELM_Accession', value = 'MatchBias')
	
	MotifDistance <- MotifDistance %>% 
		mutate(ELM_Accession = paste('MotifDistance', ELM_Accession, sep = '_')) %>% 
		spread(key = 'ELM_Accession', value = 'MeanDistance')
	
	
	FinalData <- FinalData %>% 
		left_join(MotifBias, by = c('UniversalName', 'Strain')) %>% 
		left_join(MotifDistance, by = c('UniversalName', 'Strain'))
}


if (INPUT$includeTaxonomy) {
	# Include info on taxonomy (down to family level) as well as the proportion of the family 
	# which is zoonotic
	#		- Proportion zoonotic needs to be calculated on the fly during training, since it depends on the 
	#		  test data in each iteration
	
	# Merge taxonomy sources:
	# - This ensures higher level taxonomic info for unclassified viruses is up to date - only the
	#   family comes from Genbank / the original authors
	# - Missing levels are interpolated to ensure lineages are maintained
	TaxonomyData <- read_excel(file.path('ExternalData', 'ICTV_MasterSpeciesList_2018b.xlsx'),
														 sheet = 'ICTV 2018b Master Species #34 v', col_types = "text")
	UnclassifiedTaxonomy <- read.csv(file.path('InternalData', 'Taxonomy_UnclassifiedViruses.csv'),
																	 stringsAsFactors = FALSE)
	
	TaxonomyData <- TaxonomyData %>% 
		mutate(Species = str_replace(Species, '\u00A0$', ''))  # Remove trailing non-breaking spaces in names
	
	Lineages <- TaxonomyData %>% 
		distinct(Phylum, Subphylum, Class, Subclass, Order, Suborder, Family)
	
	MergedTaxonomy <- UnclassifiedTaxonomy %>% 
		select(Species = UniversalName, Family) %>% 
		filter(! Species %in% TaxonomyData$Species) %>% 
		left_join(Lineages, by = 'Family') %>% 
		bind_rows(TaxonomyData) %>% 
		filter(Species %in% FinalData$LatestSppName) %>% 
		add_artificial_levels() %>% 
		select(Phylum, Subphylum, Class, Order, Suborder, Family, Species) %>% 
		rename_at(vars(-Species), list(~paste0('Taxonomy_', .))) %>% 
		mutate_at(vars(-Species), as.factor)  # This ensures dummy columns are created for all levels,
																					# even when some are missing in a particular dataset
	
	# Join to other data:
	stopifnot(all(FinalData$LatestSppName %in% MergedTaxonomy$Species))
	
	FinalData <- FinalData %>% 
		left_join(MergedTaxonomy, by = c('LatestSppName' = 'Species'))
}
