## Zoonosis predictor: utility functions for calculating derived genome features
library(EnvStats)

## Calculate distances between virus and host features (i.e. distance from host optimum)
#  - weightBy is optional, but when given should be a vector with same length as geneFeatures, 
#  	 specifying a variable to use for weighting geneFeature observations. 
#  	 This should be in raw form to allow NA values in geneFeatures to be accommodated.
#  - featureColNames is a list of column names (in both virusFeatures and geneFeatures) for
#    which distances should be calculated
#  - setprefix is the abbreviation to be used to identify feature columns (using the format:
#     <..DistanceType..>_<..setprefix..>_<..feature..>)
get_feature_dists <- function(virusFeatures, geneFeatures, weightBy, featureColNames, setprefix) {
	featureDists <- virusFeatures %>% 
		select(UniversalName, Strain)
	
	# These columns do not vary (e.g. because there's only one codon):
	featureColNames <- featureColNames[! featureColNames %in% c('N', 'ATG.Bias')]
	
	# Get 'distances' for each remaining feature column:
	# (actually the density given the human distribution)
	for (featureCol in featureColNames) {
		virusValues <- virusFeatures[[featureCol]]
		geneValues <- geneFeatures[[featureCol]]
		
		if (!missing(weightBy)) {
			# Calculate weights
			# - Done after removing any values corresponding to missing observations
			featureWeights <- weightBy[!is.na(geneValues)]
			featureWeights <- featureWeights / sum(featureWeights)  # Weights must sum to 1
			geneValues <- geneValues[!is.na(geneValues)]
			
			# Get density, given the distribution of weights (e.g. observed values for the humans)
			geneDensity <- demp(virusValues, obs = geneValues,
													density.arg.list = list(weights = featureWeights))
		} else {
			# No weights: treat all host gene values as equally important
			geneValues <- geneValues[!is.na(geneValues)]
			geneDensity <- demp(virusValues, obs = geneValues)
		}
		
		densityName <- paste('GenomicDensity', setprefix, featureCol, sep = '_')

		featureDists[[densityName]] <- geneDensity
	}
	
	featureDists
}