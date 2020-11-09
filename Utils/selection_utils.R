#
# Utility functions for selecting/sampling test & train data
# 

library(dplyr)


# Randomly select one strain for each species
# 'sppCol': The column representing species (as an unquoted name, as in dplyr)
# 'strainCol': The column containing strain id's (as an unquoted name)
# 'check': Should virus identifiers be checked for validity?
# returns: the original dataframe, with just one strain associated with each species
sample_strains <- function(data, sppCol = UniversalName, strainCol = Strain, check = TRUE) {
	# Capture column names
	sppCol <- enquo(sppCol)
	strainCol <- enquo(strainCol)
	sppColName <- quo_name(sppCol)
	strainColName <- quo_name(strainCol)
	
	# Check that all spp-strain combinations are unique (so the selected rows are identifiable later)
	if (check) {
		ids <- paste(data[[sppColName]], data[[strainColName]])
		if (length(ids) != length(unique(ids))) stop('Invalid input: Not all species-strain combinations are unique')
	}
	
	
	# Sample strains
	selected <- data %>% 
		group_by(!! sppCol) %>% 
		summarise(!!strainColName := sample(!! strainCol, 1)) %>% 
		ungroup()
	
	selected %>% 
		left_join(data, by = c(sppColName, strainColName))
}



# Downsample to balance class frequencies
# 'classCol': The column containing classes, as an unquoted name
downsample <- function(data, classCol) {
	classCol <- enquo(classCol)
	
	data %>% 
		group_by(!! classCol) %>% 
		add_count() %>% 
		sample_n(size = min(.$n)) %>% 
		select(-n) %>% 
		ungroup()
}