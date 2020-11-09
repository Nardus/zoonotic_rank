#! Rscript
# Part 1b of Zoonosis prediction pipeline
#	 - Perform final merges and clean data
#	 - Zoonotic status data is merged separately, by 'MergeZoonoticStatusData.R'
#	 - This script combines these data with final accessions, etc. and performs final
#	   data cleaning steps


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Constants ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ALLOW_INDIRECT_DETECTION <- FALSE  # Whether serological and enzyme-based detections should be
																	 # included when determining zoonotic status and known hosts

KEEP_HUMAN_VIRUSES <- TRUE        # Should exclusively human viruses be retained? If set to TRUE,
																	 # 'IsZoonotic' will be replaced with 'InfectsHumans'

KEEP_UNCLASSIFIED <- FALSE        # Should unclassified viruses (i.e. those not present in 
																	# MSL 2018b) be kept?


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Set up working directory -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
RootDir <- rprojroot::find_rstudio_root_file()
setwd(RootDir)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Dependencies and data ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(rjson)

options(stringsAsFactors = FALSE)


InternalData <- read.csv(file.path('InternalData', 'AllInternalData_Checked.csv'), na.strings = c('', NA))
FinalAccessions <- read.csv(file.path('InternalData', 'Final_Accessions_Unique_Spp.csv'))
ZoonoticStatusData <- readRDS(file.path('CalculatedData', 'ZoonoticStatus_Merged.rds'))

NameMatches <- read.csv(file.path('InternalData', 'NameMatches_All.csv'), sep = ',', strip.white = TRUE)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Initial cleanup / prepare for joins --------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Name Matches: 
#   - All calculations should be at the level of species, but some previously recognised species
#     may have been merged in the mean time. Thus, use the latest name where available, but keep
#     unclassified viruses as separate species
NameMatches <- NameMatches %>% 
	mutate(LatestSppName = if_else(is.na(SppName_ICTV_MSL2018b), UniversalName, SppName_ICTV_MSL2018b))

NameMatchesUniversal <- NameMatches %>% 
	distinct(UniversalName, LatestSppName)

NameMatchesOlival <- NameMatches %>% 
	distinct(Olival, LatestSppName)



## Internal reservoir data:
#		- Remove unnecessary columns
#		- Replace ? in Vector.borne column and 'Unknown'/'None' in Vector column with NA
#		- Fix capatilization of vector names
InternalData <- InternalData %>% 
	mutate(Reservoir = ifelse(Reservoir == 'Orphan', NA, Reservoir),
				 Vector.borne = ifelse(Vector.borne == '?', NA, Vector.borne),
				 Vector.borne = as.numeric(Vector.borne),
				 Vector = str_to_title(Vector),
				 Vector = ifelse(Vector %in% c('Unknown', 'None'), NA, Vector))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Merge all interal data sources ---------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# NOTE:
#   - The InternalData contains strains associated with different reservoirs, but reservoirs
#     turned out to not be predictive of zoonotic status (see Mollentze & Streicker, 2020:
#     Viral zoonotic risk is homogenous among taxonomic orders of mammalian and avian reservoir hosts)
#   - Here reducing data to a single representative sequence per virus species (with species defined
#     by LatestSppName)

# Name matches / zoonotic status
# - Only viruses for which zoonotic status is known or which have a human reservoir
FinalData <- InternalData %>% 
	left_join(NameMatchesUniversal, by = 'UniversalName') %>% 
	full_join(ZoonoticStatusData, by = 'LatestSppName') %>% 
	filter(!is.na(IsZoonotic) | Reservoir == 'Human')


# Single representative per species:
FinalData <- FinalData %>% 
	select(-Genbank.accession, -WholeGenome, -Reservoir) %>% 
	inner_join(FinalAccessions, by = c('UniversalName', 'Strain', 'LatestSppName')) %>% 
	rename(Reservoir = Reservoir_of_this_strain)

stopifnot(nrow(FinalData) <= nrow(FinalAccessions))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Main cleanup -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Modify zoonotic status of viruses thought to be zoonotic based only on indirect 
# detection methods (if requested):
if (!ALLOW_INDIRECT_DETECTION) {
	FinalData <- FinalData %>% 
		mutate(Indirect = if_else(is.na(IndirectDetectionOnly), FALSE, IndirectDetectionOnly),
					 IsZoonotic = if_else(Indirect, FALSE, IsZoonotic)) %>% 
		select(-Indirect)
}

# Remove human viruses (if requested):
if (!KEEP_HUMAN_VIRUSES) {
	FinalData <- FinalData %>% 
		filter(!Reservoir == 'Human' | is.na(Reservoir)) %>% 
		filter(!(is.na(Reservoir) & HumanOnly) | is.na(HumanOnly)) %>%  # Serves as second check in case reservoir is unknown
		filter(! is.na(IsZoonotic))
	
} else {
	# If kept, this means 'IsZoonotic' actualy measures something like 'can infect humans'
	#		- The test "(FALSE | NA)" yields NA, so we have deal with NA's first to avoid losing data
	#		- We can do this by assuming viruses are not able to infect humans unless proven otherwise,
	#		  but have to also exclude viruses for which none of the relevant variables are known as 
	#		  'not assessed'
	infection_status <- FinalData %>% 
		filter(!(is.na(Reservoir) & is.na(HumanOnly) & is.na(IsZoonotic))) %>% 
		mutate(HumanVirus = if_else(is.na(Reservoir), FALSE, Reservoir == 'Human'),
					 HumanOnly = if_else(is.na(HumanOnly), FALSE, HumanOnly),
					 IsZoonotic = if_else(is.na(IsZoonotic), FALSE, IsZoonotic)) %>% 
		group_by(LatestSppName) %>% 
		summarise(InfectsHumans = any(HumanVirus) | any(HumanOnly) | any(IsZoonotic)) %>% 
		ungroup()
	
	FinalData <- FinalData %>% 
		left_join(infection_status, by = 'LatestSppName') %>% 
		select(-IsZoonotic) %>%  # Remove column to avoid confusion and since it should not be a feature in models
		filter(! is.na(InfectsHumans))
}


## Other cleanup steps:
#		- Remove viruses marked as having multiple reservoirs 
#	  - Remove viruses marked as having no reservoir (used for accute transforming 
#		  retroviruses)
#		- Remove Hepatitis delta virus
#		- Remove viruses which have no available sequence information
#		- Remove unclassified viruses (if requested)
FinalData <- FinalData %>% 
	filter(! Reservoir %in% c('MULTIPLE', 'NONE')) %>% 
	filter(LatestSppName != 'Hepatitis delta virus') %>% 
	filter(! is.na(Genbank.accession))


if (!KEEP_UNCLASSIFIED) {
	unclassified <- NameMatches %>% 
		filter(is.na(SppName_ICTV_MSL2018b)) %>% 
		distinct(SppName_ICTV_MSL2018b, LatestSppName)
	
	FinalData <- FinalData %>% 
		filter(! LatestSppName %in% unclassified$LatestSppName)
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Remove NA's from Strain column -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# - Strain is simply an identifier used for joins, etc. 
# - Here replacing NA's by repeating UniversalName

FinalData <- FinalData %>% 
	mutate(Strain = if_else(is.na(Strain), UniversalName, Strain))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Rename columns to more convenient names ----------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FinalData <- FinalData %>% 
	rename(Accessions = Genbank.accession,
				 VectorBorne = Vector.borne)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output final data --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
saveRDS(FinalData, file.path('CalculatedData', 'FinalData_Cleaned.rds'))
write.csv(FinalData, file.path('CalculatedData', 'FinalData_Cleaned.csv'), row.names = FALSE)
