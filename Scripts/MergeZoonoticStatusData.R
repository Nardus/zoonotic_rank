##
## - Merge zoonotic status information from Olival2017, Woolhouse2018 and internally collected data
## - Apply corrections identified while searching for reservoirs of the viruses in these datasets
##


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Constants ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Animal virus families below produced by checking all families in the 2018 MSL against
# Fields virology, NCBI taxonomy (Host field), and ViralZone (as well as literature searches
# when these sources had no match)
ANIMAL_VIRUS_FAMILIES <- c('Adenoviridae', 'Anelloviridae', 'Arenaviridae', 'Arteriviridae',
													 'Asfarviridae', 'Astroviridae', 'Birnaviridae', 'Bornaviridae',
													 'Caliciviridae', 'Circoviridae', 'Coronaviridae', 'Filoviridae',
													 'Flaviviridae', 'Genomoviridae', 'Hantaviridae', 'Hepadnaviridae',
													 'Hepeviridae', 'Herpesviridae', 'Matonaviridae', 'Nairoviridae', 
													 'Orthomyxoviridae', 'Papillomaviridae', 'Paramyxoviridae', 
													 'Parvoviridae', 'Peribunyaviridae', 'Phasmaviridae', 'Phenuiviridae', 
													 'Picobirnaviridae', 'Picornaviridae', 'Pneumoviridae', 
													 'Polyomaviridae', 'Poxviridae', 'Reoviridae', 'Retroviridae', 
													 'Rhabdoviridae', 'Sunviridae', 'Tobaniviridae', 'Togaviridae')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Set up working directory -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(rprojroot)
RootDir <- find_rstudio_root_file()
setwd(RootDir)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Dependencies and data ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

options(stringsAsFactors = FALSE)
Sys.setenv(TZ="Europe/London")  # Needed by readxl (to avoid a warning)


UnclassifiedTaxonomy <- read.csv(file.path('InternalData', 'Taxonomy_UnclassifiedViruses.csv'),
																 stringsAsFactors = FALSE)

CurrentTaxonomyData <- read_excel('./ExternalData/ICTV_MasterSpeciesList_2018b.xlsx',
																	sheet = 'ICTV 2018b Master Species #34 v', col_types = 'text')

WoolhouseData <- read_xlsx('./ExternalData/WoolhouseBrierley_2018.xlsx', trim_ws = TRUE)
WoolhouseTaxonomyData <- read_excel('./ExternalData/ICTV_MasterSpeciesList_2016v1.3.xlsx',
																		sheet = 'ICTV 2016 Master Species #31')

OlivalAssociations <- read.csv('./ExternalData/Olival2017associations.csv')
OlivalViruses <- read.csv('./ExternalData/Olival2017viruses.csv')

BabayanZoonoticStatus <- read.csv('./InternalData/SourcesOfZoonoses_BabayanZoonotic.csv')

NameMatches <- read.csv('./InternalData/NameMatches_All.csv', sep = ',', strip.white = TRUE)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Initial cleanup / prepare for joins --------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Viruses marked for removal when manually matching names (duplicates):
#		- Keep virus if nothing is listed in Remove column
NameMatches <- NameMatches %>% 
	mutate(Remove = if_else(is.na(Remove), FALSE, Remove))


# Taxonomy:
UnclassifiedTaxonomy <- UnclassifiedTaxonomy %>% 
	select(LatestSppName = UniversalName, Family, Genus)

CurrentTaxonomyData <- CurrentTaxonomyData %>% 
	select(LatestSppName = Species, Family, Genus) %>% 
	mutate(LatestSppName = str_replace(LatestSppName, '\u00A0$', ''))  # Remove trailing non-breaking spaces in names

TaxonomyAll <- bind_rows(UnclassifiedTaxonomy, CurrentTaxonomyData)


## Woolhouse data:
#   - Collapse reference columns (currently ignoring host range references...)
#		- Add name matching data
# 			- Remove duplicates, etc. found during name matching
#   - Select relevant columns
WoolhouseData <- WoolhouseData %>% 
	rename(WoolhouseName = Species) %>% 
	unite('Reference', `Reference (discovery)`, `Reference (Transmission route)`, `Ref2 (TR)`, 
				sep = '; ') %>% 
	mutate(Reference = gsub('; NA', '', Reference),  # Removing NA's introduced by unite
				 IndirectDetectionOnly = `Serological detection only` == 'Y') %>% 
	left_join(NameMatches, by = c(WoolhouseName = 'Woolhouse')) %>% 
	filter(! Remove) %>%
	distinct(UniversalName, WoolhouseName, IndirectDetectionOnly, Reference,
					 TransmissionLevelWoolhouse = `Transmission level`)


## Olival data:
#		- Fix virus names
#		- Mark indirect detections
#		- Add name matching data
#		   - Remove duplicates, etc. found during name matching
#		- Select relevant columns
OlivalAssociations <- OlivalAssociations %>% 
	mutate(OlivalName = gsub('_', ' ', vVirusNameCorrected),
				 IndirectDetectionOnly = DetectionQuality %in% c(0, 1, NA)) %>%
	left_join(NameMatches, by = c(OlivalName = 'Olival')) %>% 
	filter(! Remove) %>% 
	distinct(UniversalName, IndirectDetectionOnly, Reference, hHostNameFinal)


OlivalViruses <- OlivalViruses %>% 
	mutate(OlivalName = gsub('_', ' ', vVirusNameCorrected)) %>%
	left_join(NameMatches, by = c(OlivalName = 'Olival')) %>% 
	filter(! Remove) %>% 
	distinct(UniversalName, ReverseZoonoses)



## Babayan et al's zoonotic status data:
#		- Add universal names
#		- Add name matching data
#		    - Remove duplicates, etc. found during name matching
#		- Select relevant columns
BabayanZoonoticStatus <- BabayanZoonoticStatus %>% 
	left_join(NameMatches, by = c('GenbankID' = 'Babayan')) %>% 
	filter(! Remove) %>% 
	distinct(UniversalName, 
					 TransmissionLevelBabayan = zoonotic, 
					 Reference = Ref) %>% 
	mutate(IndirectDetectionOnly = FALSE)  # Always false for these transmission levels



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Summarise zoonotic status: Olival et al ----------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Re-derived from Olival et al.'s association data so we can capture references

## Zoonotic viruses
#   - Viruses infecting humans considered zoonotic unless we know them to be a human virus (in which case 
#     they would have been marked as a reverse zooonosis)
reverseZoonoses <- OlivalViruses %>% 
	filter(ReverseZoonoses == 1) %>% 
	.$UniversalName

OlivalZoonotic <- OlivalAssociations %>% 
	group_by(UniversalName) %>% 
	mutate(Nhosts = n()) %>% 
	ungroup() %>% 
	
	filter(hHostNameFinal == 'Homo_sapiens') %>% 
	mutate(IsZoonotic = !(UniversalName %in% reverseZoonoses),
				 HumanOnly =  !IsZoonotic) %>% 
	select(-hHostNameFinal)

# Collapse duplicates:
OlivalZoonotic <- OlivalZoonotic %>% 
	group_by(UniversalName, IsZoonotic, HumanOnly) %>% 
	summarise(IndirectDetectionOnly = all(IndirectDetectionOnly),
						Reference = paste(na.omit(Reference), collapse = '; ')) %>% 
	ungroup()


# Add non-zoonotic/non-human infecting viruses assessed by Olival et al.
OlivalNonZoonotic <- OlivalAssociations %>% 
	filter(! UniversalName %in% OlivalZoonotic$UniversalName) %>% 
	mutate(IsZoonotic = FALSE,
				 HumanOnly = FALSE,
				 IndirectDetectionOnly = NA) %>% 
	select(-hHostNameFinal, -Reference) %>% 
	distinct()  # Removes duplicates caused by multiple host associations

OlivalData <- bind_rows(OlivalZoonotic, OlivalNonZoonotic) %>% 
	mutate(DataSource = 'Olival2017')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Summarise zoonotic status: Woolhouse et al -------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# In theory, Woolhouse2018's data contain all RNA viruses known to be able to infect humans
#   - "Viruses were eligible for inclusion in our catalogue only if they were classified as 
#       species by the ICTV as of July 21st 2017"
#   - The closest matching master species list is MSL 2016 V1.3 from 21 March 2017
# Thus: RNA virus species in this MSL which are not present in the Woolhouse database are considered
# non-zoonotic (animal only, or 'A'), while presence may indicate zoonotic ('Z') or 
# human-only ('H'), depending on the `TransmissionLevelWoolhouse` column
# 

## Need to match to the 2016 MSL:
#  - UniversalNames are newer, but if they don't match, the name used the original data should
WoolhouseData <- WoolhouseData %>% 
	mutate(MatchingName = if_else(UniversalName %in% WoolhouseTaxonomyData$Species, UniversalName, WoolhouseName))

stopifnot(all(WoolhouseData$MatchingName %in% WoolhouseTaxonomyData$Species))

## Extract all RNA viruses recognized in 2016 MSL:
RNAgenomes <- c("ssRNA(-)", "ssRNA(+/-)", "ssRNA(+)", "dsRNA", "ssRNA", "ssRNA-RT")

WoolhouseTaxonomyData <- WoolhouseTaxonomyData %>% 
	filter(`Genome Composition` %in% RNAgenomes) %>% 
	select(Species, Family)


## Derive zoonotic status:
#   - The new (non-zoonotic) viruses don't have a UniversalName (yet) - using the name they had in 
#     the 2016 MSL
WoolhouseData <- WoolhouseData %>% 
	full_join(WoolhouseTaxonomyData, by = c('MatchingName' = 'Species')) %>% 
	mutate(UniversalName = ifelse(is.na(UniversalName), MatchingName, UniversalName)) %>% 
	select(-WoolhouseName, -MatchingName) %>% 
	
	mutate(TransmissionLevelWoolhouse = if_else(is.na(TransmissionLevelWoolhouse), 'A',
																							if_else(TransmissionLevelWoolhouse %in% c('2', '3', '4a'), 'Z',
																											if_else(TransmissionLevelWoolhouse == '4b', 'H', 'UNHANDLED_CODE')))) %>% 
	mutate(IsZoonotic = TransmissionLevelWoolhouse == 'Z',
				 HumanOnly = TransmissionLevelWoolhouse == 'H',
				 DataSource = 'Woolhouse2018')

stopifnot(!any(WoolhouseData$TransmissionLevelWoolhouse == 'UNHANDLED_CODE'))


## Remove families known not be in ANIMAL_VIRUS_FAMILIES
#   - ANIMAL_VIRUS_FAMILIES refers to families as of 2018:
#		- Since families may have been different in 2016, we can only really check for absence
#		- This also means letting viruses with family listed as 'Unassigned' through so they can be 
#		  checked futher down
KnownNonAnimal <- CurrentTaxonomyData %>% 
	filter(!Family %in% ANIMAL_VIRUS_FAMILIES) %>% 
	.$Family %>% 
	unique()

WoolhouseData <- WoolhouseData %>% 
	filter(!Family %in% KnownNonAnimal) %>% 
	select(-Family, -TransmissionLevelWoolhouse)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Summarise zoonotic status: Babayan et al ---------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These data contains the same transmission levels (A, Z and H)
#		- Missing statuses removed, since we cannot assume they mean 'not zoonotic' (in that
#		 case they should have been marked 'A')

BabayanData <- BabayanZoonoticStatus %>% 
	filter(!is.na(TransmissionLevelBabayan) & !TransmissionLevelBabayan %in% c('?', '')) %>% 
	group_by(UniversalName) %>% 
	summarise(IsZoonotic = any(TransmissionLevelBabayan == 'Z'),
						HumanOnly = all(TransmissionLevelBabayan == 'H'),
						IndirectDetectionOnly = all(IndirectDetectionOnly), 
						Reference = if_else(IsZoonotic, paste(na.omit(Reference[TransmissionLevelBabayan == 'Z']), collapse = ';'),
																if_else(HumanOnly, paste(na.omit(Reference[TransmissionLevelBabayan == 'H']), collapse = ';'), 'NA')),
						DataSource = 'Internal')





# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Join and clean up -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FinalData <- OlivalData %>% 
	full_join(WoolhouseData) %>% 
	full_join(BabayanData)
	

## Remove duplicate rows where only the references or DataSource differ (by collapsing these):
FinalData <- FinalData %>% 
	group_by_at(vars(-IndirectDetectionOnly, -Reference, -DataSource)) %>% 
	summarise(IndirectDetectionOnly = all(IndirectDetectionOnly),
						Reference = paste(na.omit(Reference), collapse = '; '),
						DataSource = paste(DataSource, collapse = '; ')) %>% 
	ungroup()


stopifnot(! any(is.na(FinalData$IndirectDetectionOnly) & FinalData$IsZoonotic))  # Not clear how to handle this case, but it doesn't occur anyway


## Remove discrepancies where the only conflicting evidence is absence-only
# - Being listed as non-zoonotic simply means no evidence was found, but that can be 
#   overridden by another dataset finding evidence for being zoonotic
# - Similarly, 'HumanOnly' can be overriden by another dataset, under the assumption
#   that they have evidence to consider it a human virus
# - However, these categories are mutually exclusive, so evidence that the virus has a 
#   zoonotic reservoir overrides any data saying it is exclusive to humans
FinalData <- FinalData %>% 
	group_by_at(vars(-IsZoonotic, -HumanOnly, -IndirectDetectionOnly, -Reference, -DataSource)) %>% 
	summarise(IsZoonoticCollapsed = any(IsZoonotic),
						HumanOnlyCollapsed = any(HumanOnly) & !IsZoonoticCollapsed,
						
						
						IndirectDetectionOnly = if_else(any(IsZoonotic), all(IndirectDetectionOnly[IsZoonotic]), 
																						if_else(any(HumanOnly), all(IndirectDetectionOnly[HumanOnly]), 
																										all(IndirectDetectionOnly))),
						
						Reference = paste(na.omit(Reference), collapse = '; '),  # Keeping all references here
						
						DataSource = if_else(any(IsZoonotic), paste(DataSource[IsZoonotic], collapse = '; '),
																 if_else(any(HumanOnly), paste(DataSource[HumanOnly], collapse = '; '),
																 				 paste(DataSource, collapse = '; ')))) %>% 
	ungroup() %>% 
	rename(IsZoonotic = IsZoonoticCollapsed,
				 HumanOnly = HumanOnlyCollapsed)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Check for remaining disagreements ----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Disagreements between datasets will cause duplicated records

Disagreements <- FinalData %>% 
	group_by(UniversalName) %>% 
	summarise(N = n()) %>% 
	filter(N > 1)

stopifnot(nrow(Disagreements) == 0)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Add latest species names and collapse further ----------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Remove any remaining viruses from non-animal virus families:
# - Now have few enough viruses to make listing name changes (e.g. between Woolhouse's MSL and the 
#   latest one) feasible - these are included in NameMatches_All.csv
# - Using latest family assignments from either MSL 2018b or NCBI taxonomy (in the case of unclassified viruses)
stopifnot(all(ANIMAL_VIRUS_FAMILIES %in% TaxonomyAll$Family))
stopifnot(all(FinalData$UniversalName %in% NameMatches$UniversalName))

# Add latest spp name:
FinalData <- NameMatches %>% 
	distinct(UniversalName, SppName_ICTV_MSL2018b) %>% 
	right_join(FinalData, by = 'UniversalName') %>% 
	mutate(LatestSppName = if_else(is.na(SppName_ICTV_MSL2018b), UniversalName, SppName_ICTV_MSL2018b))

stopifnot(all(FinalData$LatestSppName %in% TaxonomyAll$LatestSppName))

# Add family and filter:
FinalData <- FinalData %>% 
	left_join(TaxonomyAll, by = 'LatestSppName') %>% 
	filter(!is.na(Family)) %>%
	filter(Family %in% ANIMAL_VIRUS_FAMILIES | Genus == 'Deltavirus') %>% 
	select(-Family, -SppName_ICTV_MSL2018b)



## Collapse to latest species names
#  - Needed in case the latest taxonomy update merged species further
FinalData <- FinalData %>% 
	group_by_at(vars(-UniversalName, -Genus,
									 -IsZoonotic, -HumanOnly, -IndirectDetectionOnly, -Reference, -DataSource)) %>% 
	summarise(IsZoonoticCollapsed = any(IsZoonotic),
						HumanOnlyCollapsed = any(HumanOnly) & !IsZoonoticCollapsed,
						
						
						IndirectDetectionOnly = if_else(any(IsZoonotic), all(IndirectDetectionOnly[IsZoonotic]), 
																						if_else(any(HumanOnly), all(IndirectDetectionOnly[HumanOnly]), 
																										all(IndirectDetectionOnly))),
						
						Reference = paste(na.omit(Reference), collapse = '; '),  # Keeping all references here
						
						DataSource = if_else(any(IsZoonotic), paste(DataSource[IsZoonotic], collapse = '; '),
																 if_else(any(HumanOnly), paste(DataSource[HumanOnly], collapse = '; '),
																 				paste(DataSource, collapse = '; ')))) %>% 
	ungroup() %>% 
	rename(IsZoonotic = IsZoonoticCollapsed,
				 HumanOnly = HumanOnlyCollapsed)


## Checks:
stopifnot(! any(is.na(FinalData$IsZoonotic)))
stopifnot(! any(FinalData$IsZoonotic & FinalData$HumanOnly))  # Mutually exclusive


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
saveRDS(FinalData, file.path('CalculatedData', 'ZoonoticStatus_Merged.rds'))