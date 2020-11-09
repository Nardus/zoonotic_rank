# 
# Utility functions for subsetting different types of data
# - These are needed because several columns in the dataset contain multiple entries per 
#   individual (and thus per row), e.g. accession number and host
#
require(dplyr)


# Subset sequences using a vector of strings, which may contain multiple entries
# sequences: a list of sequences as produced by seqinr::read.fasta()
# names: a vector of sequence names to keep; each element may contain multiple entries separated by the 
# 			 character(s) listed in 'separator'
# separator: a character or string to use when splitting elements in 'names'
# 
# returns: a seqinr-compatible list of sequences
subset_sequences <- function(sequences, names, separator = '; ') {
	splitnames <- names %>% 
		strsplit(split = separator) %>% 
		unlist() %>% 
		unique()
	
	# Check that all names are present
	if (! all(splitnames %in% names(sequences)))
		stop('Not all names are present in the sequence list')
	
	# Return subset
	sequences[splitnames]
}


# Subset a phylogeny to contain only the specified tips
# phylo: an ape::phylo object
# tips: a vector of tip labels; each element may contain multiple entries separated by the 
# 			character(s) listed in 'separator' (tips do not need to be unique - they will be 
# 			checked after splitting the strings)
# separator: a character or string to use when splitting elements in 'tips'
#
# returns: a phylo object
subset_phylo <- function(phylo, tips, separator = '; ') {
	splitnames <- tips %>% 
		strsplit(split = separator) %>% 
		unlist() %>% 
		unique() %>% 
		na.omit()
	
	# Check that all names are present
	if (! all(splitnames %in% phylo$tip.label)) {
		missingTips <- unique(splitnames[!splitnames %in% phylo$tip.label])
		missingString <- paste(missingTips, collapse = ', ')
		stop(paste('The following tips were not found in the phylogeny:', missingString))
	}
	
	# Return subset
	discard <- phylo$tip.label[!(phylo$tip.label %in% splitnames)]
	drop.tip(phylo, discard)
}
