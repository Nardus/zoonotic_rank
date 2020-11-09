## Utilities for handling taxonomy data:

SEARCH_ORDER <- c('Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 
									'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 
									'Species')

## Ensure missing taxonomic levels are handled correctly:
#  - For missing levels, simply use the next available level
#  			- e.g. if Family and Genus is available, but Subfamily is not, all individuals in the same 
#  			       genus should also be in the same (artificial) subfamily (but checking Family is needed
#  			       in case genus is not unique to one family)
#  - All added levels start with 'artificial_', allowing them to be detected if needed
add_artificial_levels <- function(data, search_order = SEARCH_ORDER) {
	
	# Extract taxonomic info - each row forms a lineage:
	lineages <- data[, search_order]
	lowest_level <- search_order[length(search_order)]
	
	# Find positions of gaps in each row (lineage):
	for (i in 1:nrow(data)) {
		if (! is.na(lineages[i, lowest_level])) {  # If lowest level is NA, taxonomy is unknown, so ignore
			gaps <- which(is.na(lineages[i, ]))
			nongaps <- which(! is.na(lineages[i, ]))
			
			for (gapInd in gaps) {
				
				# Find tax level closest to, but above, the gap:
				if (any(nongaps < gapInd)) {
					# We have a level above to use
					upperInd <- max(nongaps[nongaps < gapInd])
					levelAbove <- lineages[i, upperInd]
				} else {
					# This is the highest known level
					levelAbove <- 'toplevel'
				}
				
				# Find nearest known tax level below gap:
				lowerInd <- min(nongaps[nongaps > gapInd])
				levelBelow <- lineages[i, lowerInd]
				
				# Create new level for gap:
				levelName <- SEARCH_ORDER[gapInd]
				newLevel <- paste(paste0('artificial_', levelName),
													levelAbove, levelBelow,
													sep = '_')
				data[i, levelName] <- newLevel
			}
		}
	}
	
	data
}
