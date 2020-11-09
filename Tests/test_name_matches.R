#! Rscript
#
# Tests virus matching across all datasets ##
#
#

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Dependencies and data ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(testthat)
library(readxl)
library(dplyr)

options(stringsAsFactors = FALSE)
Sys.setenv(TZ="Europe/London")  # Needed by readxl (to avoid a warning)

WoolhouseData <- read_xlsx('../ExternalData/WoolhouseBrierley_2018.xlsx', trim_ws = TRUE)

BabayanData <- read.csv('../ExternalData/BabayanEtAl_Viruses.csv')

#OlivalViruses <- read.csv('./ExternalData/Olival2017viruses.csv')
OlivalAssociations <- read.csv('../ExternalData/Olival2017associations.csv')
OlivalHosts <- read.csv('../ExternalData/Olival2017hosts.csv')


NameMatches <- read.csv('../InternalData/NameMatches_All.csv', sep = ',', strip.white = TRUE)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Tests --------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Virus name matching
# All names in Woolhouse data should have matching UniversalNames
test_that('All Woolhouse data should have UniversalName matches', {
	testData <- WoolhouseData %>% 
		left_join(NameMatches, by = c(Species = 'Woolhouse'))
	
	expect_true(!any(is.na(testData$UniversalName)))
})


# All names in Babayan data should have matching UniversalNames
test_that('All Babayan data should have UniversalName matches', {
	testData <- BabayanData %>% 
		left_join(NameMatches, by = c(Virus.name = 'Babayan'))
	
	expect_true(!any(is.na(testData$UniversalName)))
})


# All names in Olival should have have matching UniversalNames
test_that('All Olival data should have matching UniversalNames', {
	testData <- OlivalHosts %>% 
		left_join(OlivalAssociations, ., by = 'hHostNameFinal') %>% 
		mutate_at(c('vVirusNameCorrected', 'hHostNameFinal'), ~ gsub('_', ' ', .)) %>% 
		left_join(NameMatches, by = c(vVirusNameCorrected = 'Olival')) %>% 
		filter(! OtherTaxonomicInfo %in% c('Unknown', 'Split'))
	
	expect_true(!any(is.na(testData$UniversalName)))
})



## Row duplication:
# 	(Olival data gets summarised to spp level, so duplication shouldn't matter)
# Joining should not cause duplicate rows in Woolhouse data
test_that('Joining should not duplicate Woolhouse data', {
	testData <- WoolhouseData %>% 
		left_join(NameMatches, by = c(Species = 'Woolhouse'))
	
	expect_equal(nrow(WoolhouseData), nrow(testData))
})

# Joining should not cause duplicate rows in Babayan data
test_that('Joining should not duplicate Babayan data', {
	testData <- BabayanData %>% 
		left_join(NameMatches, by = c(Virus.name = 'Babayan'))
	
	expect_equal(nrow(BabayanData), nrow(testData))
})