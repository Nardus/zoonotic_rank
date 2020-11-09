## Calculate genomic features
# Calculates the distance between the genome feature values of each virus and the distribution 
# of the same feature observed in the host genome
# - 'Distances' are calculated as the area under the curve (or density) when the host values are
#   treated as an emperical distribution

library(tidyr)
library(dplyr)
library(seqinr)
library(EnvStats)
library(matrixStats)

library(rprojroot)
ROOT_DIR <- find_rstudio_root_file()
setwd(ROOT_DIR)

source(file.path('Utils', 'derived_genome_feature_utils.R'))


CPM_CUTOFF <- 1 # mean CPM values above this are considered as evidence that the gene is expressed


VIRUS_DATA <- './CalculatedData/FinalData_Cleaned.rds'
VIRUS_SEQS_GENBANK <- './ExternalData/Sequences/'
VIRUS_SEQS_FASTA <- './ExternalData/Sequences/CombinedSequences.fasta'

ISG_IDENTITY <- './InternalData/Shaw2017_raw/ISG_PublishedData_Web.csv'
EXPRESSION_DATA <- './InternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv'  # TODO: should be using TPM instead of these CPM values

HOUSEKEEPING_IDENTITY <- './CalculatedData/HumanGeneSets/HousekeepingGeneIDs.csv'

TRANSCRIPT_DATA <- './CalculatedData/HumanGeneSets/TranscriptData.csv'
TRANSCRIPT_SEQS <- './CalculatedData/HumanGeneSets/TranscriptSequences.fasta'


## Regular expressions for column classes:
DINUCLEOTIDE <- '^[AUGC]p[ATUGC]'
BRIDGE_DINUCLEOTIDE <- '^br[AUGC]p[AUGC]'
NON_BRIDGE_DINUCLEOTIDE <- 'NonBr[AUGC]p[AUGC]'
AA_BIAS <- '^[A-Z].Bias'
CODON_BIAS <- '[ATGC]{3}.Bias'

# Not currently used:
NUCLEOTIDE_CONTENT <- '^[ATGCN]_'
AT_GC_CONTENT <- '^[ATGC]{2}_'
CODON_PAIR_BIAS <- '[ATGC]{3}.[A-Z]..[ATGC]{3}.[A-Z].'
CPS <- 'Cps'  # TODO: Not sure what these two columns are

ALL_FEATURES <- paste(NUCLEOTIDE_CONTENT, AT_GC_CONTENT, DINUCLEOTIDE, BRIDGE_DINUCLEOTIDE, 
											NON_BRIDGE_DINUCLEOTIDE, AA_BIAS, CODON_BIAS, CODON_PAIR_BIAS, CPS, 
											sep = '|')

USED_FEATURES <- paste(DINUCLEOTIDE, BRIDGE_DINUCLEOTIDE, NON_BRIDGE_DINUCLEOTIDE, 
											 AA_BIAS, CODON_BIAS, sep = '|')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Functions ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Calculate genomic features for a single set of sequences
# Run python3 ./GenomeFeatures.py --help in the terminal for arguments / details
calculate_genomic <- function(metaData, sequenceLocation, extract = FALSE, noncoding = FALSE) {
	## Write data and call python script
	write.csv(metaData, './_temp.csv')
	cmdString <- paste('./Utils/GenomeFeatures.py _temp.csv', sequenceLocation)
	if (extract) cmdString <- paste(cmdString, '--extract')
	if (noncoding) cmdString <- paste(cmdString, '--noncoding')
	
	system2('python3', cmdString)
	
	## Read in result
	if (!noncoding) {
		# Currently CPB_machine's output has some issues:
		# - all rows except the header end in a tab, causing data to be shifted relative to column names:
		result <- read.delim('./_temp_GenomeFeatures.csv', na.strings = c('?', 'NA'), row.names = NULL,
												 stringsAsFactors = F)
		colnames(result) <- c(colnames(result)[-1], 'XX_REMOVE')
		result <- result[, colnames(result) != 'XX_REMOVE']
	} else {
		result <- read.csv('./_temp_GenomeFeatures.csv')
	}
	
	
	## Clean up temp files and return
	unlink('./_temp.csv')
	unlink('./_temp_GenomeFeatures.csv')
	return(result)
}


# Check for transcript IDs that are missing sequences
check_missing_seqs <- function(transcriptData, transcriptSeqs) {
	seqData <- names(transcriptSeqs) %>% 
		data.frame(SeqName = ., stringsAsFactors = F) %>% 
		separate(SeqName, into = c('GeneID', 'TranscriptID'), sep = '_')
	
	if (any(is.na(transcriptData$TranscriptID)))
		warning('Some genes are missing matching transcript data')
	
	if (! all(transcriptData$TranscriptID %in% seqData$TranscriptID) )
		warning('Some transcript sequences are missing')
	
	transcriptData %>% 
		filter(! TranscriptID %in% seqData$TranscriptID)
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Data ---------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Virus data
virusMetaData <- readRDS(VIRUS_DATA)


## ISG data
isgIdentityData <- read.csv(ISG_IDENTITY, stringsAsFactors = FALSE)

## Expression data from ISG experiment
ExpressionData <- read.csv(EXPRESSION_DATA, stringsAsFactors=FALSE)
colnames(ExpressionData) <- c('GeneID', colnames(ExpressionData)[-1]) # First column unnamed in file


## Housekeeping genes:
housekeepingData <- read.csv(HOUSEKEEPING_IDENTITY)


## Transcript data
transcriptData <- read.csv(TRANSCRIPT_DATA)
transcriptSequences <- read.fasta(TRANSCRIPT_SEQS)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Viruses ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Calculating features for both coding sequences and across the entire genome, since it remains 
# unclear which features are important at which stage of virus replication
# - This also helps to capture data from coding sequences which have to be skipped because they
#   are not divisible by 3 (often because the sequence is incomplete)

# Each virus strain gets a unique summary:
virusMetaData <- virusMetaData %>% 
	select(UniversalName, Strain, Accessions) %>% 
	mutate(virusID = paste(UniversalName, Strain, sep = '_'))


# Each strain may be associated with multiple sequences (if it's segmented)
# By giving these sequences the same name, summaries are combined:
featureDat <- virusMetaData %>% 
	select(Name = virusID,
				 SequenceID = Accessions) %>% 
	mutate(SequenceID = strsplit(SequenceID, split = '; ')) %>% 
	unnest()

virusCoding <- calculate_genomic(featureDat, VIRUS_SEQS_GENBANK, extract = TRUE) %>% 	
	select(SeqName, matches(USED_FEATURES)) %>% 
	select(-X.Bias, -ATG.Bias, -TGG.Bias) %>%   # These have 0 variance (while X bias is not a real feature)
	rename_at(vars(-SeqName), ~ paste(., 'Coding', sep = '_'))

file.rename('./_temp_Problems.csv', './CalculatedData/GenomicFeatures_CDSExtractionFailed-Virus.csv')

virusEntireSeq <- calculate_genomic(featureDat, VIRUS_SEQS_FASTA, noncoding = TRUE) %>% 	
	select(SeqName, matches(USED_FEATURES)) %>% 
	rename_at(vars(-SeqName), ~ paste(., 'EntireSeq', sep = '_'))


virusFeatures <- virusCoding %>% 
	full_join(virusEntireSeq, by = 'SeqName') %>% 
	left_join(virusMetaData, by = c('SeqName' = 'virusID'))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Human genes --------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Calculate features for all expressed genes:
## - Using canonical transcript for each gene (or the longest transcript
##   if none are marked as canonical)
## - This has already been determined during download by Misc/DownloadGeneSets.py, 
##   meaning only the relevant sequences are included in transcriptSequences
## - Not all of the detected transcripts represent protein-coding genes:
##  		- For non-coding genes, calculate a subset of valid features
transcriptData <- transcriptData %>% 
	mutate(SeqName = paste(GeneID, TranscriptID, sep = '_')) %>% 
	filter(SeqName %in% names(transcriptSequences)) %>% 
	distinct()


# TODO: Below currently fails for:
# - 7 mitochondrially encoded transcripts
# - 1 endogenous retrovirus transcript
# - 1 human gene marked as 3' incomplete on Ensembl
codingFeatures <- transcriptData %>%
	filter(Biotype == 'protein_coding') %>% 
	select(SequenceID = SeqName) %>% 
	mutate(Name = SequenceID) %>% 
	calculate_genomic(TRANSCRIPT_SEQS) %>% 
	filter(!is.na(A)) %>%  # See TODO above
	select(SeqName, matches(USED_FEATURES)) %>% 
	rename_at(vars(-SeqName), ~ paste(., 'Coding', sep = '_'))

file.rename('./_temp_Problems.csv', './CalculatedData/GenomicFeatures_CDSExtractionFailed-Human.csv')

noncodingFeatures <- transcriptData %>% 
	select(SequenceID = SeqName) %>% 
	mutate(Name = SequenceID) %>% 
	calculate_genomic(TRANSCRIPT_SEQS, noncoding = TRUE) %>% 
	select(SeqName, matches(USED_FEATURES)) %>% 
	rename_at(vars(-SeqName), ~ paste(., 'EntireSeq', sep = '_'))


# Combining these creates NA's for the features that don't apply to non-protein-coding RNAs
# - This will cause these genes to be ignored when calculating distances for CDS-related features,
#   which is what we want
combinedFeatures <- codingFeatures %>% 
	full_join(noncodingFeatures, by = 'SeqName') %>% 
	separate(SeqName, into = c('GeneID', 'TranscriptID'), sep = '_')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Gene sets ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Split genes into: ISG, non-ISG housekeeping genes, Remaining genes
# In all cases, only genes with expression data from Shaw et al. are considered (i.e. genes
# detected as being expressed in the primary skin fibroblasts used by them)

## Summarise expression data:
ExpressionData <- ExpressionData %>% 
	gather(-GeneID, key = 'key', value = 'CPM') %>%   
	separate(key, into = c('Condition', 'Replicate')) %>% 
	
	group_by(GeneID, Condition) %>% 
	summarise(meanCPM = mean(CPM)) %>%  
	ungroup() %>% 
	
	filter(meanCPM >= CPM_CUTOFF) %>% 
	filter(! startsWith(GeneID, 'CVR'))


MockExpression <- ExpressionData %>% 
	filter(Condition == 'C1')

StimulatedExpression <- ExpressionData %>% 
	filter(Condition == 'C2')


## IDs for each set:
isgIDs <- isgIdentityData %>% 
	filter(Species == 'Homo sapiens') %>% 
	filter(Expression == 'up_regulated') %>% 
	filter(! startsWith(ENSEMBL.ID, 'CVR')) %>%  
	.$ENSEMBL.ID

housekeepingIDs <- housekeepingData %>% 
	filter(! GeneID %in% isgIDs) %>%   									# Removes 321 housekeeping genes which are interferon-stimulated
	filter(GeneID %in% MockExpression$GeneID) %>%  	# Removes 97 genes not detected in Shaw et al's data (at least not in mock-stimulated data)
	.$GeneID

remainingIDs <- MockExpression %>% 
	filter(! GeneID %in% isgIDs) %>% 
	filter(! GeneID %in% housekeepingIDs) %>% 
	filter(! startsWith(GeneID, 'CVR')) %>% 
	.$GeneID


## Expression and feature data for each set
## - Housekeeping and remaining expression taken from the mock-stimulated cells, so results are
## 	 comparable to those when using ISGs
isgFeatures <- StimulatedExpression %>% 
	filter(GeneID %in% isgIDs) %>% 
	left_join(combinedFeatures, by = 'GeneID') %>% 
	ungroup()

housekeepingFeatures <- MockExpression %>% 
	filter(GeneID %in% housekeepingIDs) %>% 
	left_join(combinedFeatures, by = 'GeneID') %>% 
	ungroup()
	
remainingFeatures <- MockExpression %>% 
	filter(GeneID %in% remainingIDs) %>% 
	left_join(combinedFeatures, by = 'GeneID') %>% 
	ungroup()


# Check joins (expect one row of features for each gene - multiple transcripts need to have been
# removed/summarised by this point):
stopifnot(nrow(isgFeatures) == length(unique(isgFeatures$GeneID)))
stopifnot(nrow(housekeepingFeatures) == length(unique(housekeepingFeatures$GeneID)))
stopifnot(nrow(remainingFeatures) == length(unique(remainingFeatures$GeneID)))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Distances / densities --------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
featureColNames <- colnames(virusFeatures)[grepl(ALL_FEATURES, colnames(virusFeatures))]

isgDists <- get_feature_dists(virusFeatures, 
															geneFeatures = isgFeatures, 
															setprefix = 'ISG',
															featureColNames = featureColNames)

housekeepingDists <- get_feature_dists(virusFeatures, 
																			 geneFeatures = housekeepingFeatures, 
																			 setprefix = 'Housekeeping',
																			 featureColNames = featureColNames)

remainingDists <- get_feature_dists(virusFeatures, 
																		geneFeatures = remainingFeatures, 
																		setprefix = 'Remaining',
																		featureColNames = featureColNames)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Save results -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Raw virus features:
featureColNames <- colnames(virusFeatures)[grepl(ALL_FEATURES, colnames(virusFeatures))]
featureColNames <- featureColNames[! featureColNames %in% c('N', 'ATG.Bias')] # These don't vary

virusFeatures <- virusFeatures %>% 
	select(UniversalName, Strain, featureColNames)

saveRDS(virusFeatures, './CalculatedData/GenomicFeatures-Virus.rds')


# Raw human features, needed to calculate features for novel viruses:
isgFeatures$GeneSet <- 'ISG'
housekeepingFeatures$GeneSet <- 'Housekeeping'
remainingFeatures$GeneSet <- 'Remaining'

humanFeatures <- bind_rows(isgFeatures, housekeepingFeatures, remainingFeatures) %>% 
	select(-.data$Condition)

saveRDS(humanFeatures, './CalculatedData/GenomicFeatures-HumanCombined.rds')


# For machine learning, using the distances:
featureDists <- isgDists %>% 
	full_join(housekeepingDists, by = c('UniversalName', 'Strain')) %>% 
	full_join(remainingDists, by = c('UniversalName', 'Strain'))

saveRDS(featureDists, './CalculatedData/GenomicFeatures-Distances.rds')
