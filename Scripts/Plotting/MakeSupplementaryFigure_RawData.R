##
## Plot raw data remaining after cleanining steps (i.e. the data actually used for training)
## 

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tidytext)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))


## Read data
AllData <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

# Infection type/class (calculated by MakeFigure1.R):
zoo_status <- readRDS(file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))

# Taxonomy (produced by MakeFigure1.R)
merged_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))


## Summarise data and plot:
genome_type_order <- unique(merged_taxonomy$GenomeType)
dna <- genome_type_order[grepl('DNA', genome_type_order)]
rna <- genome_type_order[!grepl('DNA', genome_type_order)]
genome_type_order <- c(sort(dna), sort(rna))


usedSpecies <- AllData %>% 
	distinct(.data$LatestSppName, .data$InfectsHumans) %>% 
	left_join(merged_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
	left_join(zoo_status, by = 'LatestSppName') %>% 
	group_by(.data$Family, .data$GenomeType, .data$Class) %>% 
	summarise(Count = n()) %>% 
	group_by(.data$Family, .data$GenomeType) %>% 
	mutate(TotalCount = sum(.data$Count)) %>% 
	ungroup() %>% 
	mutate(GenomeType = factor(.data$GenomeType, levels = genome_type_order))


raw_plot <- ggplot(usedSpecies, aes(x = reorder_within(Family, -TotalCount, GenomeType), y = Count, fill = Class)) +
	geom_col(width = 0.85, size = BAR_LINE_SIZE) +
	labs(x = 'Family', y = 'Number of species', fill = NULL) +
	facet_grid(. ~ GenomeType, scales = 'free_x', space = 'free') +
	
	scale_fill_manual(values = ZOONOTIC_STATUS_COLOURS, guide = FALSE) +
	scale_y_continuous(limits = c(0, 150), expand = expand_scale(add = c(0, 2))) +
	scale_x_reordered(expand = expand_scale(add = c(1.7, 1.7))) +
	PLOT_THEME +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
				axis.ticks.x = element_blank(),
				strip.text = element_text(angle = 45))


ggsave(file.path('Plots', 'Supplement_RawData.pdf'), raw_plot, width = 7, height = 4)

## Save an unsummarised version of this figure's data:
AllData %>% 
	distinct(.data$LatestSppName, .data$InfectsHumans) %>% 
	left_join(merged_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
	left_join(zoo_status, by = 'LatestSppName') %>% 
	write_excel_csv(file.path('FigureData', 's1_fig.csv'))
