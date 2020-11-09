#! Rscript
#
# Tests the phylogenetic neighbourhood summary functions
# - Note: diversity calculations are tested separately, in test_diversity_calculations.R;
#         only the wrapper (input/output processing) function is tested here
#

library(testthat)
library(ape)
library(dplyr)
library(tidyr)
library(rdiversity)

ROOT_DIR <- rprojroot::find_rstudio_root_file()
setwd(ROOT_DIR)
source('./Utils/PNsummary_utils.R')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- add_names_to_blast() -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('add_names_to_blast() does not allow unmatched queryIDs', {
	# Want an error in this case, not simply silent NAs in result
	blastRes <- frame_data(~queryID, ~matchID, ~percentIdentity, ~bitScore,  # Only these columns are used by the function
												 'ID_1',    'ID_A',   77,               1e-4,
												 'ID_2',    'ID_B',   75,               1e-4,
												 'Missing', 'ID_B',   80,               1e-4 )
	
	nameData <- frame_data(~UniversalName, ~Accessions,
												 'Virus1',       'ID_1',
												 'Virus2',       'ID_2',
												 'VirusA',       'ID_A',
												 'VirusB',       'ID_B')
	
	expect_error(add_names_to_blast(blastRes, nameData))
})


test_that('add_names_to_blast() does not allow unmatched matchIDs', {
	# Want an error in this case, not simply silent NAs in result
	blastRes <- frame_data(~queryID, ~matchID,  ~percentIdentity, ~bitScore,
												 'ID_1',    'ID_A',    77,               1e-4,
												 'ID_2',    'ID_B',    75,               1e-4,
												 'ID_3',    'Missing', 80,               1e-4 )
	
	nameData <- frame_data(~UniversalName, ~Accessions,
												 'Virus1',       'ID_1',
												 'Virus2',       'ID_2',
												 'Virus3',       'ID_3',
												 'VirusA',       'ID_A',
												 'VirusB',       'ID_B')
	
	expect_error(add_names_to_blast(blastRes, nameData))
})


test_that('add_names_to_blast() matches accessions and names correctly', {
	blastRes <- frame_data(~queryID, ~matchID, ~percentIdentity, ~bitScore,
												 'ID_1',    'ID_A',   77,               1e-4,
												 'ID_3',    'ID_C',   75,               1e-4,
												 'ID_1',    'ID_A',   35,               1e-4,
												 'ID_4',    'ID_B',   100,              1e-4,
												 'ID_2',    'ID_C',   80,               1e-4 )
	
	nameData <- frame_data(~UniversalName, ~Accessions,
												 'Virus1',       'ID_1',
												 'Virus2',       'ID_2',
												 'Virus3',       'ID_3',
												 'Virus4',       'ID_4',
												 'VirusA',       'ID_A',
												 'VirusB',       'ID_B',
												 'VirusC',       'ID_C')
	
	result <- add_names_to_blast(blastRes, nameData)
	
	virusNames <- nameData$UniversalName
	names(virusNames) <- nameData$Accessions
	
	expect_true(all(result$queryUniversalName == virusNames[blastRes$queryID]))
	expect_true(all(result$matchUniversalName == virusNames[blastRes$matchID]))
})


test_that('add_names_to_blast() can parse multiple accessions associated with a single virus', {
	blastRes <- frame_data(~queryID, ~matchID,  ~percentIdentity, ~bitScore,
												 'ID_1',    'ID_A',    77,               1e-4,
												 'ID_2',    'ID_B',    75,               1e-4,
												 'ID_2',    'ID_C',    75,               1e-4,
												 'ID_3',    'ID_B',    80,               1e-4 )
	
	nameData <- frame_data(~UniversalName, ~Accessions,
												 'Virus1A',       'ID_1; ID_A',
												 'Virus2B',       'ID_2; ID_B',
												 'Virus3C',       'ID_3; ID_C')
	
	result <- add_names_to_blast(blastRes, nameData)
	
	expect_true(all(result$queryUniversalName == c('Virus1A', 'Virus2B', 'Virus2B', 'Virus3C')))
	expect_true(all(result$matchUniversalName == c('Virus1A', 'Virus2B', 'Virus3C', 'Virus2B')))
})



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- calculate_diversity() ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('calculate_diversity() rejects input vectors of different lengths', {
	tree <- rtree(10)
	species <- sample(tree$tip.label, size = 20, replace = TRUE)
	neighbourhoods <- sample(LETTERS[1:5], size = 15, replace = TRUE)
	
	expect_error(calculate_diversity(species, neighbourhoods, tree))
})


test_that("calculate_diversity() rejects phylogenies that don't contain all observations", {
	tree <- rtree(10)
	species <- c(sample(tree$tip.label, size = 20, replace = T),
							 'UnknowSpp', 'UnknownSpp')
	neighbourhoods <- sample(LETTERS[5:10], size = 22, replace = T)
	
	expect_error(calculate_diversity(species, neighbourhoods, tree))
})


test_that('calculate_diversity() returns expected results (general i/o test)', {
	# This tests that counts are correct, that phylogeny is handled correctly, etc.
	tree <- rtree(10)
	species <- sample(tree$tip.label, size = 20, replace = TRUE)  # May end up not using all species in tree, but function should handle this
	neighbourhoods <- sample(LETTERS[1:5], size = 20, replace = TRUE)
	
	result <- calculate_diversity(species, neighbourhoods, tree)
	
	# Calculate expected result (independent implementation)
	partition <- table(species, neighbourhoods) %>% 
		as.data.frame() %>% 
		spread(key = neighbourhoods, value = Freq)
	
	rownames(partition) <- partition$species
	partition <- select(partition, -species)
	partition <- partition / sum(partition)  # Normalise
	
	subtree <- drop.tip(tree, tree$tip.label[! tree$tip.label %in% species])
	
	metac <- metacommunity(partition, subtree)
	expected <- subdiv(norm_alpha(metac), qs = 1)
	
	expect_equal(result$Diversity, expected$diversity)
})


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- summarise_pn_reservoir() -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('summarise_pn_reservoir() fails when blast table contains extra reservoirs', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,  # Only these columns are used by the function
												 'Virus1',             'Match1',            77,               1e-4,      'SppA',
												 'Virus2',             'Match2a',           75,               1e-4,      'SppA',
												 'Virus2',             'Match2b',           35,               1e-4,      'SppB',
												 'Virus3',             'Match3a',           100,              1e-4,      'SppB',
												 'Virus3',             'Match3b',           80,               1e-4,      'SppWithNoMatch')

	reservoirs <- c('SppA', 'SppB')

	expect_error(summarise_pn_reservoir(blastRes, reservoirs))
})


test_that('summarise_pn_reservoir() respects maxNeighbours setting', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName,  ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',             80,               1e-4,      'SppA',
												 'Virus1',             'Match2',             75,               1e-4,      'SppB',
												 'Virus1',             'Match3',             60,               1e-4,      'SppC',
												 'Virus1',             'Match4',             100,              1e-5,      'HighIdent-LowScore',
												 'Virus1',             'Match5',             35,               1e-5,      'LowIdent-LowScore')

	reservoirs <- unique(blastRes$matchReservoir)
	result <- summarise_pn_reservoir(blastRes, reservoirs, maxNeighbours = 3)

	# Top matches are chosen by bit score, so expect these to be excluded from calculation:
	zeroScores <- result$Support[result$Reservoir %in% c('HighIdent-LowScore', 'LowIdent-LowScore')]

	expect_equal(zeroScores, c(0, 0))
})


test_that('summarise_pn_reservoir() returns values for all reservoirs', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'SppA',
												 'Virus2',             'Match2a',           75,               1e-4,      'SppA',
												 'Virus2',             'Match2b',           35,               1e-4,      'SppB',
												 'Virus3',             'Match3',            100,              1e-4,      'SppB')

	reservoirs <- c('SppA', 'SppB', 'SppWithNoMatch')

	result <- summarise_pn_reservoir(blastRes, reservoirs)

	expect_true('SppWithNoMatch' %in% result$Reservoir)

	# Expect 0's (not NA), and should have an entry for each UniversalName:
	expect_equal(result$Support[result$Reservoir == 'SppWithNoMatch'], c(0,0,0))
})


test_that('summarise_pn_reservoir() returns correct relative support when there are matches to all reservoirs', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'SppA',
												 'Virus1',             'Match2',            75,               1e-4,      'SppB',
												 'Virus1',             'Match3',            35,               1e-4,      'SppC')

	reservoirs <- c('SppA', 'SppB', 'SppC')
	expectedSupport <- blastRes$percentIdentity / sum(blastRes$percentIdentity)  # Each reservoir occurs once (and in order)

	result <- summarise_pn_reservoir(blastRes, reservoirs)

	expect_equal(result$Support, expectedSupport)
})


# summarise_pn_reservoir() should return correct relative support when only one reservoir matches
test_that('summarise_pn_reservoir() returns correct relative support when only one reservoir matches', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'SppA',
												 'Virus1',             'Match2',            75,               1e-4,      'SppA')

	reservoirs <- c('SppA', 'SppB', 'SppC')
	expectedSupport <- sum( blastRes$percentIdentity / sum(blastRes$percentIdentity) )

	result <- summarise_pn_reservoir(blastRes, reservoirs)

	expect_equal(result$Support, c(expectedSupport, 0, 0))
})



# TODO: Not sure about this requirement, but something similar may be needed?
# test_that('summarise_pn_reservoir() returns equal prob for all reservoirs if there are no good matches', {
# 	blastRes <- frame_data(~queryID,    ~matchID,    ~percentIdentity, ~eValue, 
# 												 'Virus1',    'GoodMatch1', 77,               1e-4,
# 												 'Virus2',    'GoodMatch2', 75,               1e-4,
# 												 'Virus2',    'GoodMatch1', 35,               1e-4,
# 												 'SelfMatch', 'SelfMatch',  100,              1e-4,
# 												 'qBad_eVal',  'mBad_eVal', 80,               0.1)
# 	
# 	featureData <- frame_data(~Accession,   ~Reservoir,
# 														'GoodMatch1', 'HostA',
# 														'GoodMatch2', 'HostB',
# 														'SelfMatch',  'HostA',
# 														'mBad_eVal',  'HostB')
# 	
# 	hostSummary <- summarise_pn_reservoir(blastRes, featureData, removeSelf = FALSE)
# 	
# 	testRows <- hostSummary$queryID %in% c('SelfMatch', 'qBad_eVal')
# 	expect_equal(object = hostSummary$HostA[testRows],
# 							 expected = c(0.5, 0.5))
# 	expect_equal(object = hostSummary$HostB[testRows],
# 							 expected = c(0.5, 0.5))
# })
# 


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- calculate_corrected_proportion() -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('calculate_corrected_proportion() rejects arguments with different length', {
	percentIdentities <- c(100, 70, 80)
	hasFeature <- c(T, F, F)
	
	expect_error(calculate_corrected_proportion(percentIdentities, hasFeature[1:2]))
	expect_error(calculate_corrected_proportion(percentIdentities[1:2], hasFeature))
})


test_that('calculate_corrected_proportion() allows only valid percentIdentity arguments', {
	percentIdentities <- c(100, 70, 200)
	hasFeature <- c(T, F, F)
	
	expect_error(calculate_corrected_proportion(percentIdentities, hasFeature))
})



test_that('calculate_corrected_proportion() allows only logical hasFeature argument', {
	percentIdentities <- c(100, 70, 80)
	hasFeature <- c(1, 0, 1)  # Binary or other values would cause indexing to fail
	
	expect_error(calculate_corrected_proportion(percentIdentities, hasFeature))
})


test_that('calculate_corrected_proportion() allows NAs in hasFeature argument', {
	percentIdentities <- c(100, 70, 80)
	hasFeature <- c(T, F, NA)
	
	# Don't want an error here
	expect_error(calculate_corrected_proportion(percentIdentities, hasFeature), regexp = NA)
})



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- summarise_pn_reservoir_dist() --------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('summarise_pn_reservoir_dist() fails when not all matched reservoirs have distances', {
	# Unmatched reservoirs will introduce NA's downstream
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'Res1',
												 'Virus1',             'Match2',            75,               1e-4,      'Missing' )
	
	reservoirDistances <- frame_data(~FromReservoir, ~ToReservoir, ~Distance,
																	 'Res1',         'Primates',    0.4)
	
	expect_error(summarise_pn_reservoir_dist(blastRes, reservoirDistances))
})


test_that('summarise_pn_reservoir_dist() fails when targetReservoir is not in distance table', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'Res1',
												 'Virus1',             'Match2',            75,               1e-4,      'Res2' )
	
	reservoirDistances <- frame_data(~FromReservoir, ~ToReservoir, ~Distance,
																	 'Res1',          'Primates',   0.4,  
																	 'Res2',          'Primates',   0.8 )
	
	expect_error(summarise_pn_reservoir_dist(blastRes, reservoirDistances, targetReservoir = 'Missing'))
})


test_that('summarise_pn_reservoir_dist() returns values for all reservoirs (no NAs)', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'Res1',
												 'Virus1',             'Match2',            75,               1e-4,      'Res2',
												 'Virus2',             'Match1',            100,              1e-4,      'Res1')
	
	reservoirDistances <- frame_data(~FromReservoir, ~ToReservoir, ~Distance,
																	 'Res1',          'Primates',   0.4,  
																	 'Res2',          'Primates',   0.8 )
	
	result <- summarise_pn_reservoir_dist(blastRes, reservoirDistances = reservoirDistances)
	
	expect_true(all(blastRes$queryUniversalName %in% result$queryUniversalName))
	expect_false(any(is.na(result)))
})


test_that('summarise_pn_reservoir_dist() returns expected values', {
	blastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchReservoir,
												 'Virus1',             'Match1',            77,               1e-4,      'Res1',
												 'Virus1',             'Match2',            75,               1e-4,      'Res2',
												 'Virus2',             'Match1',            100,              1e-4,      'Res1')
	
	reservoirDistances <- frame_data(~FromReservoir, ~ToReservoir, ~Distance,
																	 'Res1',          'Primates',   0.4,  
																	 'Res2',          'Primates',   0.8 )
	
	result <- summarise_pn_reservoir_dist(blastRes, reservoirDistances)
	virus1Result <- result$PN_ReservoirDistToHuman[result$queryUniversalName == 'Virus1']
	virus2Result <- result$PN_ReservoirDistToHuman[result$queryUniversalName == 'Virus2']
	
	expect_equal(virus1Result, mean(c(0.4, 0.8)))
	expect_equal(virus2Result, 0.4)
})


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- summarise_pn_host_dist() -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('summarise_pn_host_dist() fails when not all matched hosts present in phylogeny', {
	blastResults <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchHost,
														 'Virus1',             'Match1',            77,               1e-4,      'Host1',
														 'Virus1',             'Match2',            75,               1e-4,      'Host2',
														 'Virus2',             'Match1',            100,              1e-4,      'Missing')
	
	hostPhylo <- rtree(3, tip.label = c('Host1', 'Host2', 'TargetHost'))
	
	expect_error(summarise_pn_host_dist(blastResults, hostPhylo, targetHost = 'TargetHost'))
})


test_that('summarise_pn_host_dist() fails when targetHost is not present in phylogeny', {
	blastResults <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchHost,
														 'Virus1',             'Match1',            77,               1e-4,      'Host1',
														 'Virus1',             'Match2',            75,               1e-4,      'Host2',
														 'Virus2',             'Match1',            100,              1e-4,      'Host3')
	
	hostPhylo <- rtree(3, tip.label = c('Host1', 'Host2', 'Host3'))
	
	expect_error(summarise_pn_host_dist(blastResults, hostPhylo, targetHost = 'Missing'))
})


test_that('summarise_pn_host_dist() returns values for all hosts (no NAs)', {
	blastResults <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchHost,
														 'Virus1',             'Match1',            77,               1e-4,      'Host1',
														 'Virus1',             'Match2',            75,               1e-4,      'Host2',
														 'Virus2',             'Match1',            100,              1e-4,      'Host3')
	
	hostPhylo <- rtree(4, tip.label = c('Host1', 'Host2', 'Host3', 'TargetHost'))
	
	result <- summarise_pn_host_dist(blastResults, hostPhylo, targetHost = 'TargetHost')
	
	expect_true(all(blastResults$queryUniversalName %in% result$queryUniversalName))
	expect_false(any(is.na(result)))
})


test_that('summarise_pn_host_dist() returns expected values', {
	blastResults <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore, ~matchHost,
														 'Virus1',             'Match1',            77,               1e-4,      'Host1',
														 'Virus1',             'Match2',            75,               1e-4,      'Host2',
														 'Virus2',             'Match1',            100,              1e-4,      'Host3')
	
	hostPhylo <- rtree(4, tip.label = c('Host1', 'Host2', 'Host3', 'TargetHost'))
	
	result <- summarise_pn_host_dist(blastResults, hostPhylo, targetHost = 'TargetHost')
	resultDists <- result$PN_HostDistToHuman
	names(resultDists) <- result$queryUniversalName
	
	# Expected:
	hostDists <- cophenetic.phylo(hostPhylo)
	hostDists <- (hostDists / max(hostDists)) %>% 
		as.data.frame() %>% 
		mutate(From = rownames(.)) %>% 
		select(From, 
					 Dist = TargetHost) %>% 
		filter(From != 'TargetHost')
	
	expectDists <- hostDists$Dist
	names(expectDists) <- hostDists$From
	
	expect_equal(resultDists[['Virus1']], mean(expectDists[c('Host1', 'Host2')]))
	expect_equal(resultDists[['Virus2']], expectDists[['Host3']])
})



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- get_distance_based_summaries() -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Mainly calls the functions already tested above, so simply test it's handling of these results

## This function and the next one require a large number of separate objects - create them once 
## to reduce amount of code:
namedBlastRes <- frame_data(~queryUniversalName, ~matchUniversalName, ~percentIdentity, ~bitScore,
														'Virus1',             'Match1',            77,               1e-4,
														'Virus2',             'Match2a',           75,               1e-4,
														'Virus2',             'Match2b',           35,               1e-4,
														'Virus3',             'Match3',           100,              1e-4)

matchData <- frame_data(~UniversalName, ~Reservoir, ~Vector, ~Hosts,
												'Match1',       'Res1',    'Vec1',   'Host1; Host2',
												'Match2a',      'Res2',    'Vec2',   'Host2',
												'Match2b',      'Res3',    'Vec1',   'Host2',
												'Match3',       'Res3',    'Vec2',   'Host2; Host3')
matchData$PubmedResults <- 1
matchData$InfectsHumans <- 'True'
matchData$VectorBorne <- 0
positive_name <- 'True'


reservoirList <- unique(matchData$Reservoir)
vectorList <- unique(matchData$Vector)

hostPhylo <- rtree(n = 4, tip.label = c('Host1', 'Host2', 'Host3', 'Homo_sapiens'))
reservoirPhylo <- rtree(n = length(unique(matchData$Reservoir)), tip.label = unique(matchData$Reservoir))

reservoirDistances <- frame_data(~FromReservoir, ~ToReservoir, ~Distance,
																 'Res1',         'Primates',    0.4,
																 'Res2',         'Primates',    0.2,
																 'Res3',         'Primates',    0.6)


## Tests
test_that('get_distance_based_summaries() does not introduce NAs when joining results', {
	result <- get_distance_based_summaries(pIdentCutoff = 0,
																				 namedBlastRes = namedBlastRes,
																				 matchData = matchData,
																				 hostPhylogeny = hostPhylo,
																				 reservoirPhylogeny = reservoirPhylo,
																				 reservoirDistances = reservoirDistances)
	
	expect_false(any(is.na(result)))
})


test_that('get_distance_based_summaries() does not drop viruses when joining results', {
	result <- get_distance_based_summaries(pIdentCutoff = 0,
																				 namedBlastRes = namedBlastRes,
																				 matchData = matchData,
																				 hostPhylogeny = hostPhylo,
																				 reservoirPhylogeny = reservoirPhylo,
																				 reservoirDistances = reservoirDistances)
	
	expect_true(all(namedBlastRes$queryUniversalName %in% result$queryUniversalName))
})


test_that('get_distance_based_summaries() does not duplicate viruses when joining results', {
	result <- get_distance_based_summaries(pIdentCutoff = 0,
																				 namedBlastRes = namedBlastRes,
																				 matchData = matchData,
																				 hostPhylogeny = hostPhylo,
																				 reservoirPhylogeny = reservoirPhylo,
																				 reservoirDistances = reservoirDistances)
	
	expect_equal(length(result$queryUniversalName), length(unique(result$queryUniversalName)))
})



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- summarise_pn() -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Main function called by other scripts
# - Ensure input is checked for potential issues
# - Check handling of results from other functions

test_that('summarise_pn() fails when not all matches are present in matchData', {
	blastExtra <- rbind(namedBlastRes,
											c('Virus3', 'Missing', 100, 1e-4))
	
  
	expect_error(summarise_pn(blastExtra, matchData, 
														hostPhylogeny = hostPhylo, 
														reservoirPhylogeny = reservoirPhylo,
														reservoirDistances = reservoirDistances,
														reservoirList = reservoirList, 
														vectorList = vectorList,
														pIdentityCutoffs = 0,
														nthreads = 1,
														positiveName = positive_name))
})



test_that('summarise_pn() fails when data contains species with multiple reservoirs', {
	matchDataExtra <- rbind(matchData,
													c('Match3', 'ExtraReservoir', 'Vec2', 'Host2', 1, 'InfectsHumans', 'False'))
	
	reservoirDistanceExtra <- rbind(reservoirDistances,
																	c('ExtraReservoir', 'Primates', 0.6))

	reservoirPhyloExtra <- rtree(n = length(unique(matchDataExtra$Reservoir)), tip.label = unique(matchDataExtra$Reservoir))
	

	expect_error(summarise_pn(namedBlastRes, matchDataExtra, 
														hostPhylogeny = hostPhylo, 
														reservoirPhylogeny = reservoirPhyloExtra,
														reservoirDistances = reservoirDistanceExtra,
														reservoirList = reservoirList, 
														vectorList = vectorList,
														pIdentityCutoffs = 0,
														nthreads = 1,
														positiveName = positive_name))
})


test_that('summarise_pn() correctly parses viruses with multiple hosts', {
	result <- summarise_pn(namedBlastRes, matchData, 
												 hostPhylogeny = hostPhylo, 
												 reservoirPhylogeny = reservoirPhylo,
												 reservoirDistances = reservoirDistances,
												 reservoirList = reservoirList, 
												 vectorList = vectorList,
												 pIdentityCutoffs = 0,
												 nthreads = 1,
												 positiveName = positive_name)
	
	# Check if host dists could be calculated:
	resultHostDists <- result$PN_HostDistToHuman_0
	names(resultHostDists) <- result$queryUniversalName
	
	# Expected:
	hostDists <- cophenetic.phylo(hostPhylo)
	hostDists <- (hostDists / max(hostDists)) %>% 
		as.data.frame() %>% 
		mutate(From = rownames(.)) %>% 
		select(From, 
					 Dist = Homo_sapiens) %>% 
		filter(From != 'Homo_sapiens')
	
	expectDists <- hostDists$Dist
	names(expectDists) <- hostDists$From
	
	
	expect_equal(resultHostDists[['Virus1']], mean(expectDists[c('Host1', 'Host2')])) # One match, which has two hosts
	expect_equal(resultHostDists[['Virus2']], expectDists[['Host2']]) # Two matches, but they have the same host
})


test_that('summarise_pn() returns values for all viruses (no NAs)', {
	result <- summarise_pn(namedBlastRes, matchData, 
												 hostPhylogeny = hostPhylo, 
												 reservoirPhylogeny = reservoirPhylo,
												 reservoirDistances = reservoirDistances,
												 reservoirList = reservoirList, 
												 vectorList = vectorList,
												 pIdentityCutoffs = 0,
												 nthreads = 1,
												 positiveName = positive_name)

	expect_false(any(is.na(result)))
})


test_that('summarise_pn() does not drop viruses', {
	result <- summarise_pn(namedBlastRes, matchData, 
												 hostPhylogeny = hostPhylo, 
												 reservoirPhylogeny = reservoirPhylo,
												 reservoirDistances = reservoirDistances,
												 reservoirList = reservoirList, 
												 vectorList = vectorList,
												 pIdentityCutoffs = 0,
												 nthreads = 1,
												 positiveName = positive_name)
	
	expect_true(all(namedBlastRes$queryUniversalName %in% result$queryUniversalName))
})


test_that('summarise_pn() does not duplicate viruses', {
	# This could happen during joining
	result <- summarise_pn(namedBlastRes, matchData, 
												 hostPhylogeny = hostPhylo, 
												 reservoirPhylogeny = reservoirPhylo,
												 reservoirDistances = reservoirDistances,
												 reservoirList = reservoirList, 
												 vectorList = vectorList,
												 pIdentityCutoffs = 0,
												 nthreads = 1,
												 positiveName = positive_name)
	
	expect_equal(length(result$queryUniversalName), length(unique(result$queryUniversalName)))
})