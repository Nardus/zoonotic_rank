library(testthat)

ROOT_DIR <- rprojroot::find_rstudio_root_file()
setwd(ROOT_DIR)

source('./Utils/cds_parser.R')

SEQ_CHARS <- strsplit("AAATGCTACACTAGCGATGGTGA", split = "")[[1]]
SEQUENCE <- seqinr::as.SeqFastadna(SEQ_CHARS, name = "TestSeq")
VALID_START <- 3
VALID_STOP <- 11
TRAILING_STOP <- 14
INTERNAL_STOP <- 20


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Input checks -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

test_that("extract_cds should accept only single seqinr::SeqFastadna objects", {
	expect_error(extract_cds(list(SEQUENCE), VALID_START, VALID_STOP), 'Invalid sequence input')
})

test_that("extract_cds should not allow string-based sequence objects", {
	# seqinr's read.fasta has an 'as.string' option, but this means subsetting does not work
	testSeq <- seqinr::as.SeqFastadna("AAATGCTACACTAGCGATGGTGA", name = "Test")
	expect_error(extract_cds(testSeq, VALID_START, VALID_STOP), 'Sequence should be a vector')
})


test_that("extract_cds should not allow reversed coordinates", {
	expect_error(extract_cds(SEQUENCE, VALID_STOP, VALID_START), 'Start coordinate must be lower than stop coordinate')
})


test_that("extract_cds should not allow equal coordinates", {
	expect_error(extract_cds(SEQUENCE, VALID_START, VALID_START), 'Start coordinate must be lower than stop coordinate')
})


test_that("extract_cds should warn about possible 0-based coordinates", {
	# In reality, this would almost certainly cause an error anyway, since the extracted
	# sequence would not be divisible by 3. However, the warning may help pinpoint the
	# cause of this.
	testSeq <- seqinr::as.SeqFastadna(c("A", "T", "G", "C", "T", "A"), name = "Test")
	
	expect_warning(extract_cds(SEQUENCE, 0, 6), '0-based')
})


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Logic checks -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

test_that("extract_cds should return the expected region", {
	expectedSeq <- seqinr::as.SeqFastadna(SEQ_CHARS[VALID_START:VALID_STOP], name = 'TestSeq')
	returnedSeq <- extract_cds(SEQUENCE, VALID_START, VALID_STOP)
	
	expect_equal(returnedSeq, expectedSeq)
})


test_that("extract_cds should retain sequence name", {
	returnedSeq <- extract_cds(SEQUENCE, VALID_START, VALID_STOP)
	expectedName <- attr(SEQUENCE, 'name')
	returnedName <- attr(returnedSeq, 'name')
	
	expect_equal(returnedName, expectedName)
})


test_that("extract_cds should not return sequences which are not divisible by 3", {
	expect_error(extract_cds(SEQUENCE, VALID_START, VALID_STOP-2), "not a multiple of 3")
})


test_that("extract_cds should not return sequences with internal stop codons", {
	expect_error(extract_cds(SEQUENCE, VALID_START, INTERNAL_STOP), "stop codon")
})


test_that("extract_cds should allow sequences ending in a stop codon", {
	expect_error(extract_cds(SEQUENCE, VALID_START, TRAILING_STOP), NA) # Expect no error
})


test_that("extract_cds should strip trailing stop codon", {
	expectedSeq <- SEQUENCE[VALID_START:VALID_STOP]
	returnedSeq <- extract_cds(SEQUENCE, VALID_START, TRAILING_STOP)
	returnedSeq <- as.character(returnedSeq)
	
	expect_equal(returnedSeq, expectedSeq)
})
