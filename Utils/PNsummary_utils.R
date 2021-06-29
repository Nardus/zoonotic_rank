#
# Utility functions for summarising features in the phylogenetic neighbourhood of target viruses
# 
# 
library(rprojroot)
library(rjson)
library(ape)
library(dplyr)
library(tidyr)
library(parallel)

ROOT_DIR <- find_rstudio_root_file()


# Add UniversalNames for both the query and match accession numbers in a blast result table
# blastRes: a blast result table
#	data: a data frame containing columns 'UniversalName' and 'Accession', with entries for all
#				viruses in both the blast database and the query (multiple Accession entries may be 
#				colon-separated)
# 
# returns: the blast result table with columns 'queryUniversalName' and 'matchUniversalName' 
# 					added (note the case)
add_names_to_blast <- function(blastRes, data) {
	accessionDF <- data %>% 
		select(UniversalName, Accessions) %>% 
		mutate(Accessions = strsplit(Accessions, split = '; ')) %>% 
		unnest() %>% 
		rename(Accession = Accessions)
	
	AccessionNames <- accessionDF$UniversalName
	names(AccessionNames) <- accessionDF$Accession
	
	# Check that all accessions can be found
	if (! all(blastRes$queryID %in% accessionDF$Accession))
		stop('Not all query IDs were found in the Accession column of the supplied data')
	if (! all(blastRes$matchID %in% accessionDF$Accession))
		stop('Not all match IDs were found in the Accession column of the supplied data - ensure data contains all viruses in blast database')
	
	# Add names and return
	blastRes %>% 
		mutate(queryUniversalName = unname(AccessionNames[queryID]),
					 matchUniversalName = unname(AccessionNames[matchID]))
}


# Calculate the phylogenetic alpha diversity of observations in each PN
# observations: a vector of names (e.g. species) (in raw form - will be counted in the function)
# neigbourhoodIDs: a vector of identifiers specifying where the observations occur
# phylo: a phylogeny containing all observations as tip labels (can also contain other names,
#        these will be removed)
# 
# returns: a data frame with columns 'Neighbourhood' and 'Diversity'
calculate_diversity <- function(observations, neighbourhoodIDs, phylo) {
	require(rdiversity)

	# Check input
	if (length(observations) != length(neighbourhoodIDs))
		stop('observations and neighbourhoodIDs should have equal length')
	
	## Handle input of all NA's:
	if (all(is.na(observations))) {
		blankOut <- data.frame(Neighbourhood = neighbourhoodIDs,
													 Diversity = NA,
													 stringsAsFactors = FALSE)
		return(blankOut)
	}
	
	# Similarly, if there is only one non-NA species, set diversity to 1 (since there 
	# cannot by a phylogeny in this case)
	uniqueSpp <- unique(na.omit(observations))
	
	if (length(uniqueSpp) == 1) {
		blankOut <- data.frame(Neighbourhood = neighbourhoodIDs,
													 Diversity = 1,
													 stringsAsFactors = FALSE)
		return(blankOut)
	}
	
	
	## At least some observations are known:
	# Remove NA's so they don't get counted as a category
	naPositions <- is.na(observations)
	observations <- observations[!naPositions]
	neighbourhoodIDs <- neighbourhoodIDs[!naPositions]
	
	if (!all(observations %in% phylo$tip.label))
		stop('Not all observed species are present in the phylogeny')
	
	# Subset phylogeny to contain only the observed species
	discard <- phylo$tip.label[! phylo$tip.label %in% observations]
	subphylo <- drop.tip(phylo, discard)
	
	# Count species occurences
	abundance <- table(observations, neighbourhoodIDs) %>% 
		as.data.frame(stringsAsFactors = FALSE) %>% 
		spread(key = neighbourhoodIDs, value = Freq)
	
	rownames(abundance) <- abundance$observations
	abundance <- abundance[, -1, drop = FALSE]
	
	# Normalise to relative abundances
	abundance <- abundance / sum(abundance)
	
	# Re-arrange abundances to match phylo
	abundance <- abundance[subphylo$tip.label, , drop = FALSE]
	
	
	# Calculate phylogenetic alpha diversity
	metacom <- metacommunity(abundance, subphylo)
	normAlpha <- norm_alpha(metacom)
	result <- subdiv(normAlpha, qs = 1)
	
	# Return only names and diversity values
	if (ncol(abundance) == 1) {
		# Partition name not properly returned when there is only 1 partition
		partitions <- colnames(abundance)
	} else {
		partitions <- as.character(result$partition_name)
	}
	
	data.frame(Neighbourhood = partitions,
						 Diversity = result$diversity,
						 stringsAsFactors = FALSE)
}


# Summarise the support for each reservoir category based on those of neigbouring viruses
# blastResults: a table of blast results, annotated with 'queryUniversalName' and
# 							'matchReservoir' columns
# reservoirList: the reservoirs for which probabilities should be calculated - those not present
# 							 in the PN will be given a probability of 0
# maxNeigbours: the maximum number of neigbours to condsider; if more are available, the
# 							top hits will be chosen based on bit scores (in case of ties, the 
# 							percentIdentity will be considered, for any remaining ties, the first
# 							entry is taken)
# 
# returns: a data frame with columns 'queryUniversalName', 'Reservoir' and 'Support'
summarise_pn_reservoir <- function(blastResults, reservoirList, maxNeighbours = 5) {
	# Check input
	if (! all(na.omit(blastResults$matchReservoir) %in% reservoirList))
		stop('blastResults contains reservoirs not listed in reservoirList')
	
	if (length(reservoirList) != length(unique(reservoirList)))
		stop('reservoirList contains duplicates')  # this is not a problem (easy to deal with), but may indicate the wrong argument was passed
	
	
	# Restrict neighbourhood size
	blastResults <- blastResults %>% 
		group_by(queryUniversalName) %>% 
		mutate(MatchRanking = order(bitScore, percentIdentity, decreasing = T)) %>% 
		filter(MatchRanking <= maxNeighbours)
	
	
	# Add missing reservoirs, so we get support values for all of them:
	reservoirs <- expand.grid(queryUniversalName = unique(blastResults$queryUniversalName), 
														matchReservoir = reservoirList, stringsAsFactors = FALSE)
	
	blastResults <- blastResults %>% 
		right_join(reservoirs, by = c('queryUniversalName', 'matchReservoir'))
	
	
	# Calculate relative support
	blastResults$SupportAdded <- blastResults$percentIdentity / sum(blastResults$percentIdentity, na.rm = T)
	
	blastResults %>% 
		group_by(queryUniversalName, matchReservoir) %>% 
		summarise(Support = sum(SupportAdded, na.rm = T)) %>% 
		rename(Reservoir = matchReservoir) %>% 
		ungroup()
}


## Calculate distance-corrected proportions:
# 'percentIdentities': a vector of percentIdentity scores
# 'hasFeature': a logical vector (same length as blastResults), indicating whether
#               the match in the corresponding element of percentIdentities has
#               the feature to be counted
# returns: a single proportion value
calculate_corrected_proportion <- function(percentIdentities, hasFeature) {
	if (length(percentIdentities) != length(hasFeature))
		stop('percentIdentities and hasFeature must have equal length')
	
	if (min(percentIdentities) < 0 | max(percentIdentities) > 100)
		stop('percentIdentities not in [0, 100]')
	
	if (! is.logical(hasFeature))
		stop('hasFeature must a vector of logical (true/false) values')
	
	similarities <- percentIdentities[hasFeature] / 100
	N <- length(percentIdentities)  # The total number of neighbours
	
	sum(similarities) / N
}


# Calculate distance-corrected publication count:
# 'percentIdentities': a vector of percentIdentity scores
# 'publicationCounts': a vector of publication counts
# returns: a single publication count
summarise_publication_count <- function(percentIdentities, publicationCounts) {
	if (length(percentIdentities) != length(publicationCounts))
		stop('percentIdentities and publicationCounts must have equal length')
	
	similarities <- percentIdentities/100
	sum( similarities * publicationCounts, na.rm = TRUE) / length(percentIdentities)
}



## Calculate the average distance from RESERVOIRS in each PN to humans
# 'blastResults': named blast results, annotated with a 'matchReservoir' column
# 'reservoirDistances': a data frame of distances between all reservoir groups, with 
#                       columns FromReservoir, ToReservoir, and Distance
summarise_pn_reservoir_dist <- function(blastResults, reservoirDistances, targetReservoir = 'Primates') {
	reservoirDistances <- filter(reservoirDistances, ToReservoir == targetReservoir)
	
	# Check inputs:
	#  - Allow reservoirs to be unknown (NA), but all known reservoirs must have distance data
	#  - Otherwise, left join below will introduce new NA's which are unlikely to be detected downstream
	knownReservoirs <- na.omit(blastResults$matchReservoir)
	
	if (! all(knownReservoirs %in% reservoirDistances$FromReservoir))
		stop('Not all reservoirs in blastResults have associated distances to the target reservoir, this will introduce NAs')
	
	# Retrieve distances:
	blastResults %>% 
		left_join(reservoirDistances, by = c('matchReservoir' = 'FromReservoir')) %>% 
		group_by(queryUniversalName) %>% 
		summarise(PN_ReservoirDistToHuman = mean(Distance, na.rm = T)) %>% 
		ungroup()
}


## Calculate the average distance from HOSTS in each PN to humans
# 'reservoirDistances': a data frame of distances between all reservoir groups, with 
#                       columns FromReservoir, ToReservoir, and Distance
summarise_pn_host_dist <- function(blastResults, hostphylo, targetHost = 'Homo_sapiens') {
	if (! all(na.omit(blastResults$matchHost) %in% hostphylo$tip.label))
		stop('Not all hosts in blast results are present in the host phylogeny')
	if (! targetHost %in% hostphylo$tip.label)
		stop('Target host not present in supplied host phylogeny')
	
	distmat <- cophenetic.phylo(hostphylo)
	distmat <- distmat / max(distmat)  # Normalise
	distDF <- as_data_frame(distmat)
	distDF$From <- rownames(distmat)
	
	distDF <- distDF %>% 
		gather(-From, key = 'To', value = 'Distance') %>% 
		filter(To == targetHost)
	
	blastResults %>% 
		left_join(distDF, by = c('matchHost' = 'From')) %>% 
		group_by(queryUniversalName) %>% 
		summarise(PN_HostDistToHuman = mean(Distance, na.rm = T)) %>% 
		ungroup()
}


# Utility function, see summarise_pn() below
# pIdentCutoff: a single percent identity value at which to restrict the phylogenetic neigbourhood
# all other arguments: see see summarise_pn()
# 
get_distance_based_summaries <- function(pIdentCutoff, namedBlastRes, matchData, hostPhylogeny, 
																				 reservoirPhylogeny, reservoirDistances) {
	# Restrict PN
	restrictedBlastRes <- namedBlastRes %>% 
		filter(percentIdentity >= pIdentCutoff)
	
	# Join data
	reservoirBlast <- matchData %>% 
		select(matchUniversalName = UniversalName,
					 matchReservoir = Reservoir) %>% 
		right_join(restrictedBlastRes, by = 'matchUniversalName')
	
	hostBlast <- matchData %>% 
		mutate(Hosts = strsplit(Hosts, split = '; ')) %>% 
		unnest() %>% 
		select(matchUniversalName = UniversalName,
					 matchHost = Hosts) %>% 
		right_join(restrictedBlastRes, by = 'matchUniversalName')
	
	
	## Calculate summaries:
	# Reservoir diversity
	reservoirDiversity <- calculate_diversity(observations = reservoirBlast$matchReservoir,
																						neighbourhoodIDs = reservoirBlast$queryUniversalName,
																						phylo = reservoirPhylogeny) %>% 
		rename(queryUniversalName = Neighbourhood,
					 PN_ReservoirDiversity = Diversity)
	
	# Host diversity
	hostDiversity <- calculate_diversity(observations = hostBlast$matchHost, 
																			 neighbourhoodIDs = hostBlast$queryUniversalName,
																			 phylo = hostPhylogeny) %>% 
		rename(queryUniversalName = Neighbourhood,
					 PN_HostDiversity = Diversity)
	
	# Distance from reservoirs in PN to Primates
	reservoirHumanDist <- summarise_pn_reservoir_dist(reservoirBlast, reservoirDistances, targetReservoir = 'Primates')
	
	# Distance from hosts in PN to Humans
	hostHumanDist <- summarise_pn_host_dist(hostBlast, hostPhylogeny)
	
	# Other summaries (currently only size):
	OtherSummaries <- matchData %>% 
		select(matchUniversalName = UniversalName) %>% 
		right_join(restrictedBlastRes, by = 'matchUniversalName') %>% 
		group_by(queryUniversalName) %>% 
		summarise(PN_Size = n())
	
	
	# Join results
	results <- hostDiversity %>% 
		full_join(reservoirDiversity, by = 'queryUniversalName') %>% 
		full_join(reservoirHumanDist, by = 'queryUniversalName') %>% 
		full_join(hostHumanDist, by = 'queryUniversalName') %>% 
		full_join(OtherSummaries, by = 'queryUniversalName')
	
	# Add pIdentCutoff as a label to all summary columns and return
	results %>% 
		rename_at(vars(starts_with('PN_')), ~ paste(., pIdentCutoff, sep = '_'))
}


# Calculate all PN summaries and join results
# namedBlastRes: cleaned blast results, annotated with names
# matchData: a data frame containing 'UniversalName', 'Accession' and 'Reservoir' columns,
#								 with entries for ALL viruses in the blast database
# hostPhylogeny: a phylogeny containing all hosts named in the Hosts column of matchData
# reservoirPhylogeny: a phylogeny containing all names present in the Reservoir column of matchData
# reservoirDistances: a data frame with distances between all reservoir categories, containing 
#                     columns 'FromReservoir', 'ToReservoir', and 'Mean'
#	reservoirList: a list of all known reservoirs
#	vectorList: a list of all known vectors
#	pIdentityCutoffs: a vector of percent identity cutoff values at which to calculate PN summaries
#	nthreads: the number of threads available for parallel execution
#	maxNeigbours: the maximum number of hits to include in reservoir support calculations
#	minimal: calculate a basic summary of zoonotic viruses in the neigbourhood only; reservoir, etc
#	         is ignored in this case
# 
# returns: a data frame with columns 'queryUniversalName', and various columns starting PN_XX,
#          where XX is a summary type
summarise_pn <- function(namedBlastRes, matchData, hostPhylogeny, reservoirPhylogeny, 
												 reservoirDistances, reservoirList, vectorList, pIdentityCutoffs, nthreads, 
												 maxNeighbours = 5, positiveName = POSITIVE_NAME, minimal = FALSE, check = TRUE) {
	## Check data
	if (check) {
		if (! all(namedBlastRes$matchUniversalName %in% matchData$UniversalName))
			stop('Not all viruses in the matchUniversalName column of namedBlastRes are present in matchData')
		
		reservoirCounts <- rowSums(table(matchData$UniversalName, matchData$Reservoir))
		
		if (any(reservoirCounts > 1))
			stop('At least one virus species has more than one reservoir listed in matchData')
	}
	
	# Simplify matches to 1 per species
	#  - Queries and matches representing segmented viruses will have multiple entries, but should not be counted
	#    multiple times 
	#  - So far this is only an issue for reservoir support values - for proportions it doesn't matter, since the same feature would be repeated
	namedBlastRes <- namedBlastRes %>% 
		group_by(queryUniversalName, matchUniversalName) %>% 
		summarise(percentIdentity = mean(percentIdentity),
							bitScore = mean(bitScore))
	
	## Actual calculations
	if (minimal) {
		# Calculate basic summaries only
		if (length(unique(matchData$UniversalName)) != nrow(matchData))
			stop("matchData contains duplicates - this is not accounted for when making minimal summaries")  # Still valid for other summaries, since the different reservoirs would require such replication (see test above)
		
		result <- matchData %>% 
			select(matchUniversalName = UniversalName,
						 matchInfectsHumans = InfectsHumans,
						 matchVectorBorne = VectorBorne) %>% 
			right_join(namedBlastRes, by = 'matchUniversalName') %>% 
			mutate(matchInfectsHumans = matchInfectsHumans == positiveName,
						 matchIsVectorBorne = matchVectorBorne == 1) %>% 
			group_by(queryUniversalName) %>% 
			summarise(PN_HumanProportion_DC = calculate_corrected_proportion(percentIdentity, matchInfectsHumans),
								PN_HumanProportion_raw = sum(matchInfectsHumans, na.rm = T) / length(na.omit(matchInfectsHumans))) %>% 
			ungroup() %>% 
			select(queryUniversalName, starts_with('PN_'))
		
		
	} else {
		## These do not vary by distance, since we use maxNeighbours as a closer cut-off:
		# Reservoir support values
		reservoirSummary <- matchData %>% 
			select(matchUniversalName = UniversalName,
						 matchReservoir = Reservoir) %>% 
			right_join(namedBlastRes, by = 'matchUniversalName') %>% 
			summarise_pn_reservoir(reservoirList, maxNeighbours) %>% 
			mutate(Reservoir = paste0('PN_Reservoir', Reservoir, 'Support')) %>%  # Rename before they become columns
			spread(key = Reservoir, value = Support)
		
		# Vector support values
		vectorSummary <- matchData %>% 
			select(matchUniversalName = UniversalName,
						 matchReservoir = Vector) %>% 
			right_join(namedBlastRes, by = 'matchUniversalName') %>% 
			summarise_pn_reservoir(vectorList, maxNeighbours) %>% 
			mutate(Reservoir = paste0('PN_Vector', Reservoir, 'Support')) %>%
			spread(key = Reservoir, value = Support)
		
		
		# Distance-corrected proportions and publication count:
		# - Also including non-corrected versions, in case our correction is actually unhelpful
		proportionSummary <- matchData %>% 
			select(matchUniversalName = UniversalName,
						 matchInfectsHumans = InfectsHumans,
						 matchVectorBorne = VectorBorne,
						 matchPublications = PubmedResults) %>% 
			right_join(namedBlastRes, by = 'matchUniversalName') %>% 
			mutate(matchInfectsHumans = matchInfectsHumans == positiveName,
						 matchIsVectorBorne = matchVectorBorne == 1) %>% 
			group_by(queryUniversalName) %>% 
			summarise(PN_HumanProportion_DC = calculate_corrected_proportion(percentIdentity, matchInfectsHumans),
								PN_VectorBorneProportion_DC = calculate_corrected_proportion(percentIdentity, matchIsVectorBorne),
								PN_Publications_DC = summarise_publication_count(percentIdentity, matchPublications),
								
								PN_HumanProportion_raw = sum(matchInfectsHumans, na.rm = T) / length(na.omit(matchInfectsHumans)),
								PN_VectorBorneProportion_raw = sum(matchIsVectorBorne, na.rm = T) / length(na.omit(matchIsVectorBorne)),
								PN_Publications_raw = sum(matchPublications) / n()) %>% 
			ungroup() %>% 
			select(queryUniversalName, starts_with('PN_'))
		
		
		## All other summaries change depending on cut-off:
		distanceSummaries <- mclapply(pIdentityCutoffs, FUN = get_distance_based_summaries,
																	namedBlastRes = namedBlastRes, 
																	matchData = matchData, 
																	hostPhylogeny = hostPhylogeny, 
																	reservoirPhylogeny = reservoirPhylogeny, 
																	reservoirDistances = reservoirDistances,
																	mc.cores = nthreads)
		
		distanceSummaries <- Reduce(distanceSummaries, f = function(...) full_join(..., by = 'queryUniversalName'))
		
		## Join all summaries and return:
		result <- full_join(reservoirSummary, vectorSummary, by = 'queryUniversalName') %>% 
			full_join(proportionSummary, by = 'queryUniversalName') %>% 
			full_join(distanceSummaries, by = 'queryUniversalName')
	}
	
	result
}
