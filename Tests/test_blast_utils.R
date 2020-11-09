#! Rscript
#
# Tests the wrapper functions for interacting with blast in Utils/blast_utils.R
#

RootDir <- rprojroot::find_rstudio_root_file()
setwd(RootDir)

source('./Utils/blast_utils.R')

# A sequence set for testing 'real' blast searches:
BIG_SEQUENCE_SET <- read.fasta('./ExternalData/Sequences/CombinedSequences.fasta', as.string = T)
BIG_TEST_DB <- make_blast_db(BIG_SEQUENCE_SET)

#
## GENERAL
#
test_that('blast commandline tools are available and reachable', {
	exitCode <- system2('blastn', '-version', stdout = F)
	expect_equal(exitCode, 0)
})


#
## MAKE_BLAST_DB()
#
test_that('make_blast_db() returns a valid path', {
	fastaLines <- c(">Seq1", "ATTCAA",
									">Seq2", "ATGCCC")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	dbpath <- make_blast_db(sequences)
	
	dbBaseName <- basename(dbpath)
	dbdir <- sub(dbBaseName, '', dbpath)
	
	# Can't make assumptions about the file extensions used by different versions of blast, so
	# just check that ot created a directory of db files:
	expect_true(dir.exists(dbdir))
	unlink_blast_db(dbpath)
})


test_that('make_blast_db() rejects input when sequence name contains spaces', {
	# When this happens, SeqInR names the sequence incorrectly, so we won't be able to identify 
	# it in matches
	fastaLines <- c(">Name_with_underscores", "ATGC",
									">Name with spaces", "ATGC")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	expect_error(make_blast_db(sequences))
})


test_that('make_blast_db() rejects input with duplicated names', {
	fastaLines <- c(">Duplicate", "ATTCAA",
									">Duplicate", "ATGCCC")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	expect_error(make_blast_db(sequences))
})


test_that('make_blast_db() correctly names sequences in created database', {
	fastaLines <- c(">Name_with_underscores", "ATGC",
									">NameWithStr@ngeCharacters!", "ATGC")  # Blast won't take unicode, but these should work at least
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	dbpath <- make_blast_db(sequences)
	
	exitCode1 <- system2('blastdbcmd', stdout = FALSE,
											 args = paste0('-db "', dbpath, '" -entry "Name_with_underscores" '))
	exitCode2 <- system2('blastdbcmd', stdout = FALSE,
											 args = paste0('-db "', dbpath, '" -entry "NameWithStr@ngeCharacters!" '))
	
	expect_equal(exitCode1, 0)
	expect_equal(exitCode2, 0)
	unlink_blast_db(dbpath)
})


#
## BLASTN()
#

test_that('blastn() rejects query with duplicate names', {
	# The results of such a query would be difficult to handle, and may 
	# invalidate summary stats calculated later
	fastaLines <- c(">Seq1", "ATTCAA",
									">Seq2", "ATGCCC")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	testDB <- make_blast_db(sequences)
	
	queryLines <- c(">Duplicate", "ATACC",
									">Duplicate", "ATGTCC")
	fastaFile <- textConnection(queryLines)
	querySet <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	expect_error(blastn(querySet, testDB, 10))
	unlink_blast_db(testDB)
})


test_that('blastn() returns a dataframe of results', {
	fastaLines <- c(">Seq1", "ATTCAAATTCAAGCGC",
									">Seq2", "ATGCCCGTTAC",
									">Seq3", "ATGATTCACGTTCGACCA")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	testDB <- make_blast_db(sequences)
	
	querySet <- sequences["Seq3"]
	
	result <- blastn(querySet, testDB, 10)
	expect_equal(class(result), 'data.frame')  # This is pretty strict - what if it's a tibble?
	unlink_blast_db(testDB)
})


test_that('blastn() respects the max_target_seqs option passed', {
	fastaLines <- c(">SimilarSeq_1", "ATTCAAATTCAAGCGCCCCGTTACATGATTCACGTTCGACCA",
									">SimilarSeq_2", "ATTCAAGTTCAAGCGACCCGTGACATGATTCAAGTTCGACCA",
									">SimilarSeq_3", "AGTCAAATTCAAGGGCCCCGTTACTTGATTCAAGTTCGACCA",
									">SimilarSeq_4", "ATTCCAATTTAAGCGCCCCGTAACATGATTAACGTTCGACCA")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	testDB <- make_blast_db(sequences)
	
	querySet <- sequences["SimilarSeq_3"]
	
	result <- blastn(querySet, testDB, 2)
	expect_equal(nrow(result), 2)
	unlink_blast_db(testDB)
})


#
## CLEAN_BLAST_HITS()
#
test_that('clean_blast_hits() removes low e-value hits', {
	# A query set with a range of e values
	fastaLines <- c(">Rabies_lyssavirus_N", "cacctctacaatggatgccgacaagattgtattcaaagtcaataatcaggtggtctctttgaagcctgagattatcgtggatcaatatgagtacaagtaccctgccatcaaagatttgaaaaagccctgtataactctaggaaaggctcccgatttaaataaagcatacaagtcagttttatcatgcatgagcgccgccaaacttgatcctgacgatgtatgttcctatttggcggcggcaatgcagttttttgaggggacatgtccggaagactggaccagctatggaatcgtgattgcacgaaaaggagataagatcaccccaggttctctggtggagataaaacgtactgatgtagaagggaattgggctctgacaggaggcatggaactgacaagagaccccactgtccctgagcatgcgtccttagtcggtcttctcttgagtctgtataggttgagcaaaatatccgggcaaagcactggtaactataagacaaacattgcagacaggatagagcagatttttg",
									">Australian_bat_lyssavirus_N", "atggagtctgataagattgcctttaagatcaacaatcaattggtgtctgttaagccggaggtgatagtagatcagtatgaatataagtaccctgcaatcaaagatcagaggaagcctagcattactcttggaaaggccccagatttaaataaagcgtacaagtctatattatctggcatgaacgccgcgaagttggacccggacgatgtttgctcctacctagctgcagctatggagttatttgaggggatctgcccagaggactggacaagttacgggattttgattgccagaaaaggagacaaaatcacaccggctactttagttgacataaggagaacagatattcagggcagctgggctctggcaggggggcaggactttaccagagaccctacaatcgcagagcatgcatctctggtgggtcttcttctgagcctctacagattgagcaaaatttcaggtcaaaacacgggaaactacaaaaccaacatcgcagacaggattgagcagatttttgag",
									">Bad_hit_distant", "aacaggggatgacatctgtcaaaattttcacagtatcctggagggaaggagcaaactacgagtcagttcaagttggtgattttttaaagagagtagatgatgatgagaaaatgagtagtacatctgtcaaaattttaaatgcaatcaagccattgaacagtcacccatcatacatttcaaatgcattcagaattcccaagacttgttccacaaatgcattcagaattcccaagacttgttccacaaatgcattcagaattcccaagacttgttccactttggaaggagcaaaggtttggaaggagcaaaggtttggaaggatcaccccaagcaaaggtttctttagttgacataggaaggagcaaagg")
	fastaFile <- textConnection(fastaLines)
	querySet <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	result <- blastn(querySet, BIG_TEST_DB, 10)
	
	# Check that the test is working:
	if (! all(result$eValue[result$queryID == 'Bad_hit_distant'] >= 1E-3))
		stop('This unit test is broken - check sequence db and input sequences')
	
	cleanResult <- clean_blast_hits(result)
	
	expect_false("Bad_hit_distant" %in% cleanResult$queryID)
})


test_that('clean_blast_hits() removes hits to self when asked', {
	fastaLines <- c(">SimilarSeq_1", "ATTCAAATTCAAGCGCCCCGTTACATGATTCACGTTCGACCA",
									">SimilarSeq_2", "ATTCAAGTTCAAGCGACCCGTGACATGATTCAAGTTCGACCA",
									">SimilarSeq_3", "AGTCAAATTCAAGGGCCCCGTTACTTGATTCAAGTTCGACCA",
									">SimilarSeq_4", "ATTCCAATTTAAGCGCCCCGTAACATGATTAACGTTCGACCA")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	testDB <- make_blast_db(sequences)
	
	querySet <- sequences["SimilarSeq_3"]
	
	result <- blastn(querySet, testDB, 10)
	cleanresult <- clean_blast_hits(result, removeSelf = TRUE)
	
	expect_false(any(cleanresult$queryID == cleanresult$matchID))
	unlink_blast_db(testDB)
})


test_that('clean_blast_hits() keeps hits to self when asked', {
	fastaLines <- c(">SimilarSeq_1", "ATTCAAATTCAAGCGCCCCGTTACATGATTCACGTTCGACCA",
									">SimilarSeq_2", "ATTCAAGTTCAAGCGACCCGTGACATGATTCAAGTTCGACCA",
									">SimilarSeq_3", "AGTCAAATTCAAGGGCCCCGTTACTTGATTCAAGTTCGACCA",
									">SimilarSeq_4", "ATTCCAATTTAAGCGCCCCGTAACATGATTAACGTTCGACCA")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	testDB <- make_blast_db(sequences)
	
	querySet <- sequences["SimilarSeq_3"]
	
	result <- blastn(querySet, testDB, 10)
	cleanresult <- clean_blast_hits(result, removeSelf = FALSE)
	
	expect_true(any(cleanresult$queryID == cleanresult$matchID))
	unlink_blast_db(testDB)
})


test_that('clean_blast_hits() removes hits to self by default', {
	fastaLines <- c(">SimilarSeq_1", "ATTCAAATTCAAGCGCCCCGTTACATGATTCACGTTCGACCA",
									">SimilarSeq_2", "ATTCAAGTTCAAGCGACCCGTGACATGATTCAAGTTCGACCA",
									">SimilarSeq_3", "AGTCAAATTCAAGGGCCCCGTTACTTGATTCAAGTTCGACCA",
									">SimilarSeq_4", "ATTCCAATTTAAGCGCCCCGTAACATGATTAACGTTCGACCA")
	fastaFile <- textConnection(fastaLines)
	sequences <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	testDB <- make_blast_db(sequences)
	
	querySet <- sequences["SimilarSeq_3"]
	
	result <- blastn(querySet, testDB, 10)
	cleanresult <- clean_blast_hits(result)
	
	expect_false(any(cleanresult$queryID == cleanresult$matchID))
	unlink_blast_db(testDB)
})


#
## SUBSET_BLAST()
#
test_that('subset_blast() returns corrected result', {
	# subsetting a blast result should produce the same effect as directly blasting against the 
	# smaller database
	
	small_db_seqs <- BIG_SEQUENCE_SET[1:500]
	small_test_db <- make_blast_db(small_db_seqs)
	
	query_seqs <- BIG_SEQUENCE_SET[c(1, 10, 20)]
	
	blast_result_full <- blastn(query_seqs, BIG_TEST_DB, max_target_seqs = 1000)
	blast_result_subset <- blastn(query_seqs, small_test_db, max_target_seqs = 1000)
	unlink_blast_db(small_test_db)
	
	# We expect the match sets (blast_result_full vs blast_result_subset) to differ in the low-quality
	# matches found, but that doesn't matter as these would be filtered out anyway:
	blast_result_full <- filter(blast_result_full, eValue < 1e-2)
	blast_result_subset <- filter(blast_result_subset, eValue < 1e-2)
	
	# Does our subset_blast() function replicate this?
	subset_result <- subset_blast(full_result = blast_result_full,
																full_db_sequences = BIG_SEQUENCE_SET,
																subset_db_sequences = small_db_seqs)
	
	expect_equal(subset_result, blast_result_subset, tolerance = 1e-6)
})


test_that('subset_blast() raises error when sequences not in original db', {
	fastaLines <- c(">Rabies_lyssavirus_N", "cacctctacaatggatgccgacaagattgtattcaaagtcaataatcaggtggtctctttgaagcctgagattatcgtggatcaatatgagtacaagtaccctgccatcaaagatttgaaaaagccctgtataactctaggaaaggctcccgatttaaataaagcatacaagtcagttttatcatgcatgagcgccgccaaacttgatcctgacgatgtatgttcctatttggcggcggcaatgcagttttttgaggggacatgtccggaagactggaccagctatggaatcgtgattgcacgaaaaggagataagatcaccccaggttctctggtggagataaaacgtactgatgtagaagggaattgggctctgacaggaggcatggaactgacaagagaccccactgtccctgagcatgcgtccttagtcggtcttctcttgagtctgtataggttgagcaaaatatccgggcaaagcactggtaactataagacaaacattgcagacaggatagagcagatttttg",
									">Australian_bat_lyssavirus_N", "atggagtctgataagattgcctttaagatcaacaatcaattggtgtctgttaagccggaggtgatagtagatcagtatgaatataagtaccctgcaatcaaagatcagaggaagcctagcattactcttggaaaggccccagatttaaataaagcgtacaagtctatattatctggcatgaacgccgcgaagttggacccggacgatgtttgctcctacctagctgcagctatggagttatttgaggggatctgcccagaggactggacaagttacgggattttgattgccagaaaaggagacaaaatcacaccggctactttagttgacataaggagaacagatattcagggcagctgggctctggcaggggggcaggactttaccagagaccctacaatcgcagagcatgcatctctggtgggtcttcttctgagcctctacagattgagcaaaatttcaggtcaaaacacgggaaactacaaaaccaacatcgcagacaggattgagcagatttttgag",
									">Bad_hit_distant", "aacaggggatgacatctgtcaaaattttcacagtatcctggagggaaggagcaaactacgagtcagttcaagttggtgattttttaaagagagtagatgatgatgagaaaatgagtagtacatctgtcaaaattttaaatgcaatcaagccattgaacagtcacccatcatacatttcaaatgcattcagaattcccaagacttgttccacaaatgcattcagaattcccaagacttgttccacaaatgcattcagaattcccaagacttgttccactttggaaggagcaaaggtttggaaggagcaaaggtttggaaggatcaccccaagcaaaggtttctttagttgacataggaaggagcaaagg")
	fastaFile <- textConnection(fastaLines)
	querySet <- read.fasta(fastaFile, as.string = TRUE)
	close(fastaFile)
	
	result <- blastn(querySet, BIG_TEST_DB, 10)
	
	expect_error(subset_blast(result, BIG_SEQUENCE_SET, querySet), "Not a subset")
})



unlink_blast_db(BIG_TEST_DB)

