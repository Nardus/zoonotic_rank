#! Rscript
#
# Tests the subset generation functions in selection_utils.R
#

library(testthat)

ROOT_DIR <- rprojroot::find_rstudio_root_file()
setwd(ROOT_DIR)

source('./Utils/selection_utils.R')

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- sample_strains -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('sample_strains() returns all columns in original data', {
	testData <- data.frame(UniversalName = LETTERS[1:10],
												 Strain = LETTERS[1:10],
												 ColumnA = runif(10),
												 ColumnB = runif(10))
	
	processedData <- sample_strains(testData)
	
	expect_equal(colnames(processedData), colnames(testData))
})


test_that('sample_strains() does not return grouped data frame', {
	# Unclear how downstream code would be affected by this, so best to remove grouping
	testData <- data.frame(UniversalName = LETTERS[1:10],
												 Strain = LETTERS[1:10],
												 ColumnA = runif(10))
	
	processedData <- sample_strains(testData)
	
	expect_false(is.grouped_df(processedData))
})


test_that('sample_strains() returns just one strain per species', {
	testData <- data.frame(UniversalName = sample(LETTERS[1:5], 26, replace = T),
												 Strain = LETTERS,
												 ColumnA = runif(26))
	
	processedData <- sample_strains(testData)
	
	expect_equal(length(processedData$UniversalName), 5)
})


test_that('sample_strains() allows species with no (NA) strains', {
	testData <- data.frame(UniversalName = c('A', 'A', 'B'),
												 Strain = c('Strain A', 'Strain B', NA),
												 ColumnA = runif(3), stringsAsFactors = FALSE)
	
	processedData <- sample_strains(testData)
	
	expect_equal(processedData$UniversalName, c('A', 'B'))
	expect_true(is.na(processedData[processedData$UniversalName == 'B', 'Strain']))
})


test_that('sample_strains() does not allow repeated spp-strain combinations', {
	# This would make it impossible to trace which row of the data was actually selected
	testData <- data.frame(UniversalName = c('A', 'A', 'B'),
												 Strain = c('Strain A', 'Strain A', NA),
												 ColumnA = runif(3))
	
	expect_error(sample_strains(testData), 'Invalid input')
})



test_that('sample_strains() returns only valid species-strain combinations', {
	species <- replicate(paste(sample(LETTERS, 5), collapse = ''), n = 200)  # generate unique spp ids
	
	testData <- data.frame(UniversalName = sample(species, 10000, replace = TRUE),
												 Strain = replicate(paste(sample(LETTERS, 10), collapse = ''), n = 10000))
	
	inputIDs <- paste(testData$UniversalName, testData$Strain)
	
	processedData <- sample_strains(testData)
	selectedIDs <- paste(processedData$UniversalName, processedData$Strain)
	
	expect_true(all(selectedIDs %in% inputIDs))
})



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- downsample ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
test_that('downsample() returns all columns in original data', {
	testData <- data.frame(ClassColumn = rep('A', 10),
												 ColumnB = LETTERS[1:10],
												 ColumnC = runif(10))
	
	processedData <- downsample(testData, ClassColumn)
	
	expect_equal(colnames(processedData), colnames(testData))
})


test_that('downsample() does not return a grouped dataframe', {
	testData <- data.frame(ClassColumn = rep('A', 10),
												 ColumnB = LETTERS[1:10],
												 ColumnC = runif(10))
	
	processedData <- downsample(testData, ClassColumn)
	
	expect_false(is.grouped_df(processedData))
})


test_that('downsample() returns equal class frequencies', {
	testData <- data.frame(ClassColumn = c(rep('A', 10), rep('B', 4)),
												 ColumnB = LETTERS[1:14],
												 ColumnC = runif(14))
	
	processedData <- downsample(testData, ClassColumn)
	freqA = sum(processedData$ClassColumn == 'A')
	freqB = sum(processedData$ClassColumn == 'B')
	
	expect_equal(freqA, freqB)
})


test_that('downsample() reduces frequency by the expected amount', {
	testData <- data.frame(ClassColumn = c(rep('A', 10), rep('B', 4)),
												 ColumnB = LETTERS[1:14],
												 ColumnC = runif(14))
	
	processedData <- downsample(testData, ClassColumn)
	
	expect_equal(sum(processedData$ClassColumn == 'A'), 4)
})

