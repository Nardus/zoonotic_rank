# 
# Utility functions for calculating the phylogenetic neighbourhood
# Requires BLAST commandline tools (for installation instructions, see 
#   https://www.ncbi.nlm.nih.gov/books/NBK279671/)
#
require(seqinr)
require(parallel)
require(digest)
require(dplyr)


## Check for issues with sequences passed to the functions below
#  sequences: A list of sequences as produced by seqinr::read.fasta('path', as.string = T)
#  
#  returns: NULL, but raises an error if issues are found
check_sequences <- function(sequences) {
	if (! all(lapply(sequences, class) == "SeqFastadna"))
		stop('Supplied sequence object invalid (expected a list of SeqFastadna objects)')
	
	seqnames <- names(sequences)
	annotations <- unlist(lapply(sequences, attr, 'Annot'))
	annotations <- sub(pattern = '>', replacement = '', annotations)
	
	if (any(seqnames != annotations))
		stop('Sequence names may have been parsed incorrectly by seqinr. This is known to happen when sequence names contain spaces in the fasta file.')
	
	if (length(unique(seqnames)) != length(seqnames))
		stop('The supplied sequence object contains duplicate names. This will cause problems downstream.')
}



## Wrapper function to make a temporary blast database
# sequences: A list of sequences as produced by seqinr::read.fasta('path', as.string = T)
#
# returns: a path to the database created
make_blast_db <- function(sequences) {
	check_sequences(sequences)
	
	seqpath <- tempfile(pattern = 'make_blast_db-sequences-', fileext = '.fasta')
	dbpath <- tempfile(pattern = 'make_blast_db-blastdb-')
	
	write.fasta(sequences, names(sequences), file.out=seqpath, open = 'w', nbchar = 100, as.string = T)
	
	# Make db
	argString <- paste('-dbtype nucl -parse_seqids',
										 '-in', seqpath,
										 '-out', dbpath)
	returnCode <- system2('makeblastdb', argString, stdout = FALSE)
	
	if (returnCode != 0)
		stop("Failed to create blast database")
	
	# Return the path to this db
	unlink(seqpath)
	return(dbpath)
}


## Clean up temporary files created by make_blast_db()
# dbpath: the full path to the blast db to remove
unlink_blast_db <- function(dbpath) {
	dbBaseName <- basename(dbpath)
	dbdir <- sub(dbBaseName, '', dbpath)
	dbfiles <- list.files(path = dbdir, pattern = dbBaseName)
	dbfiles <- paste0(dbdir, dbfiles)
	
	file.remove(dbfiles)
}


## Run a blastn search and load results
# queryset: A list of sequences as produced by seqinr::read.fasta('path', as.string = T)
# dbpath: Path to the blast database to query
# Other options (max_target_seqs, num_threads, max_hsps, reward, evalue, 
# word_size, gapopen, gapextend): Blast options (see https://www.ncbi.nlm.nih.gov/books/NBK279684/)
#
# returns: a dataframe
blastn <- function(queryset, dbpath, max_target_seqs, 
									 num_threads = 2, max_hsps = 1, reward = 2, evalue = 10, word_size = 8,
									 gapopen = 2, gapextend = 2) {
	
	check_sequences(queryset)
	
	if (num_threads < 1 | num_threads %% 1 != 0)
		stop("Invalid num_threads argument")
	
	# The blast options:
	blastopts <- paste('-max_target_seqs', max_target_seqs, '-num_threads', num_threads, 
										 '-max_hsps', max_hsps, '-reward', reward, '-evalue', evalue, 
										 '-word_size', word_size, '-gapopen', gapopen, '-gapextend', gapextend,
										 '-task blastn',
										 '-outfmt "10 qaccver saccver pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"')
	
	# Temporarily save the query set:
	tempQ <- tempfile(pattern = 'blastn-queryset-', fileext = '.fasta')
	write.fasta(queryset, names(queryset), file.out=tempQ, open = 'w', nbchar = 100, as.string = T)
	
	# Run the blast query:
	tempOut <- tempfile(pattern = 'blastn-outfile-', fileext = '.out')
	
	argString <- paste('-db', dbpath, 
										 '-query', tempQ,
										 '-out', tempOut,
										 blastopts)
	returnCode <- system2('blastn', argString, stdout = FALSE)
	
	if (returnCode != 0)
		stop("Running blastn command failed")
	
	
	# Load results:
	result <- read.csv(file = tempOut,
										 col.names = c('queryID', 'matchID', 'percentIdentity', 'qlen', 'alignedLen', 
										 							'nMismatch', 'gapOpen', 'qStart', 'qEnd', 'sStart', 'sEnd',
										 							'eValue', 'bitScore'), header=F, stringsAsFactors = FALSE)
	
	# Clean temp files and return:
	unlink(c(tempQ, tempOut))
	return(result)
}


## Run parallel blastn searches and load results
#  Splits the query set into equally sized fragments, and blasts these in parallel, before 
#  merging results
#  
# queryset: A list of sequences as produced by seqinr::read.fasta('path', as.string = T)
# dbpath: Path to the blast database to query
# nfragments: number of fragments into which queryset should be split.
# nthreads: maximum number of threads that can be used
#             - If there are more threads available than nfragments,
# 					    some will remain unused, with the free threads divided equally among 
#               blast instances.
# 					  - If there are more fragments than threads, each blast instance will
# 					    run on a single thread, with other fragments processed as threads
# 					    become available
# ...: further options passed to blastn()
blastn_parallel <- function(queryset, dbpath, nfragments, nthreads, ...) {
	# Figure out how to divide threads
	if (nfragments >= nthreads) {
		fragmentThreads <- nthreads
		blastThreads <- 1
	} else {
		fragmentThreads <- nfragments
		blastThreads <- nthreads %/% nfragments
	}
	
	# Split the query sequences
	inds <- 1:length(queryset)
	inds <- sample(inds)   # Randomly distribute sequences, so long sequences (from the same 
												 # family, for example), don't end up on the same node. This will work
												 # best if blast jobs are small (i.e. there are more fragments than 
												 # threads, which will also alow cores to be re-assigned if some jobs 
												 # take longer)
	fragmentInds <- split(inds, cut(seq_along(inds), nfragments, labels = FALSE)) 
	
	# Set up parallel blast searches
	queryFragments <- lapply(fragmentInds, function(x) queryset[x])
	blastResults <- mclapply(queryFragments, FUN = blastn, 
													 dbpath = dbpath, num_threads = blastThreads, ..., 
													 mc.cores = fragmentThreads)
	
	# Combine results and return
	blastResults <- do.call("rbind", blastResults)
	return(blastResults)
}



## Clean blastn results by removing bad hits and (optionally) hits to self
# blastresult: a dataframe produced by blastn()
# eCutOff: the e-value cut-off for included (counted) hits
# removeSelf: whether mathces of to the query sequence to itself should be removed (occur when 
#							the blast database contained the query sequence)
clean_blast_hits <- function(blastresult, eCutOff = 1e-3, removeSelf = T) {
	if (removeSelf) {
		selfMatches <- (blastresult$queryID == blastresult$matchID) & (blastresult$percentIdentity == 100)
		blastresult <- blastresult[!selfMatches, ]
	}
	
	# Remove bad hits
	blastresult <- blastresult[blastresult$eValue <= eCutOff, ]
	
	return(blastresult)
}



## Cached blast: Check if this exact combination of query and test data has been seen
## before, and load the existing result instead if this is the case
#		querySeqs / dbSeqs: A list of sequences as produced by seqinr::read.fasta('path', as.string = T),
#		used for the query and the search database, respecively
#		max_target_seqs: blast option (see 'man blastn' on commandline)
#		nfragments: number of fragments into which queryset should be split.
# 	nthreads: maximum number of threads that can be used (see blastn_parallel above)
# 	cache_dir: a directory used for storing cached results
cached_blast <- function(querySeqs, dbSeqs, max_target_seqs, nfragments, nthreads, cache_dir) {
	# Check for existing blast results:
	# - Results are stored in RDS format under cache_dir/dbMD5_queryMD5
	# - Both sets of sequences are sorted by name to ensure order does not matter.
	# 	However, this may cause unexpected results during blasting when max_target_seqs is low
	# 	and there are many equivalent matches - arguably we want a random selection of the best
	# 	results, not an alphabetically biased one. Thus, we need keep the original input order
	# 	for blasting to allow such issues to be dealt with in higher-level code (if needed)
	dbSeqsOrdered <- dbSeqs[order(names(dbSeqs))]
	querySeqsOrdered <- querySeqs[order(names(querySeqs))]
	
	dbMD5 <- digest(dbSeqsOrdered)
	queryMD5 <- digest(querySeqsOrdered)
	
	cacheFile <- file.path(cache_dir, paste(dbMD5, queryMD5, sep = '_'))
	
	if (file.exists(cacheFile)) {
		cachedResult <- readRDS(cacheFile)
		return(cachedResult)
	}
	
	# No existing result, so do the search:
	blastDB <- make_blast_db(dbSeqs)
	
	blastResults <- blastn_parallel(querySeqs, blastDB, 
																	max_target_seqs = max_target_seqs,
																	nfragments = nfragments, 
																	nthreads = nthreads)
	unlink_blast_db(blastDB)
	
	# Save results in cache and return:
	if (!file.exists(cache_dir))
		dir.create(cache_dir)
	
	saveRDS(blastResults, file = cacheFile)
	
	blastResults
}



## Correct an existing blast result tables to match the effect that would have been obtained when 
## blasting against only a subset of the original blast database
#		full_result: the original blast result, obtained for the full database
#		full_db_sequences: sequences used to create the original, full database (in seqinr::read.fasta's format)
#		subset_db_sequences: sequences matching the new database, which must be a subset of full_db_sequences
subset_blast <- function(full_result, full_db_sequences, subset_db_sequences) {
	# Check input
	# - New db sequences should be a subset of original:
	if (!all(names(subset_db_sequences) %in% names(full_db_sequences)))
		stop("Not a subset: not all sequence names in 'subset_db_sequences' are present in the original database")
	
	# - DB length calculations below assume a list of character sequences:
	if (!is.list(subset_db_sequences) | !class(subset_db_sequences[[1]]) %in% c('SeqFastadna', 'character'))
		stop("subset_db_sequences is not a list of SeqFastadna objects or character strings. Check input")
	
	# Remove matches no longer present in search database:
	subset_result <- full_result %>% 
		filter(.data$matchID %in% names(subset_db_sequences))
	
	# Correct e-value
	# E = (qlen x db_len) / 2^bitscore
	new_db_length =  sapply(subset_db_sequences, nchar) %>% 
		sum() %>% 
		as.double()  # Using double/numeric to prevent integer overflow in next step
	
	subset_result %>% 
		mutate(eValue = (.data$qlen * new_db_length) / 2^.data$bitScore)
}

