#
# Plot an overview of ranks for training set viruses
# 

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(readxl)
library(ggplot2)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))

all_ranks <- read_csv(file.path('Plots', 'Intermediates', 'combined_virus_ranks.csv'))

merged_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))
zoonotic_status_raw <- readRDS(file.path('CalculatedData', 'ZoonoticStatus_Merged.rds'))


## Get cutoff from training data
training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

CUTOFF <- sum(training_data$InfectsHumans)/nrow(training_data)


## Prepare data
genome_type_order <- unique(merged_taxonomy$GenomeType)
dna <- genome_type_order[grepl('DNA', genome_type_order)]
rna <- genome_type_order[!grepl('DNA', genome_type_order)]
genome_type_order <- rev(c(sort(dna), sort(rna)))  # RNA viruses at top of plot


training_ranks <- all_ranks %>% 
	filter(.data$dataset == 'Training data') %>% 
	left_join(merged_taxonomy, by = c('species' = 'Species')) %>% 
	arrange(.data$probability_mean) %>% 
	mutate(species = factor(.data$species, levels = .data$species),
				 Family = factor(.data$Family, levels = sort(unique(.data$Family), decreasing = TRUE)),
				 priority = factor(.data$priority, levels = c('Low', 'Medium', 'High', 'Very high')),
				 infection_class = if_else(.data$infection_class == 'No known human infections', 
				 													'No human infections', .data$infection_class),
				 infection_class = factor(.data$infection_class, 
				 												 levels = c('Human virus', 'Zoonotic', 'No human infections')),
				 GenomeType = factor(.data$GenomeType, levels = genome_type_order))


## Main panel
main_panel <- ggplot(training_ranks, aes(x = species, y = probability_mean, colour = infection_class)) +
	geom_errorbar(aes(ymin = probability_lower, ymax = probability_upper), width = 0.5, alpha = 0.5) +
	geom_point(size = 0.5) +
	geom_hline(yintercept = CUTOFF, linetype = 2, colour = 'grey10') +
	scale_colour_manual(values = ZOONOTIC_STATUS_COLOURS) +
	scale_y_continuous(limits = c(0, 1), expand = expand_scale(add = c(0.02, 0.02))) +
	labs(x = NULL, y = SCORE_LABEL, colour = 'Current status') +
	PLOT_THEME +
	theme(panel.grid.major.x = element_blank(),
				axis.ticks.x = element_blank(),
				axis.text.x = element_blank(),
				plot.margin = margin(t = 5.5, r = 5.5, b = 1, l = 5.5))


## Indicate families
family_indicator_plot <- ggplot(training_ranks, aes(x = species, y = Family, fill = priority)) +
	geom_tile() + 
	facet_grid(rows = vars(GenomeType), scales = 'free', space = 'free') +
	scale_fill_manual(values = PRIORITY_COLOURS) +
	PLOT_THEME +
	labs(x = 'Virus species (ranked)', fill = 'Priority') +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				axis.text.y = element_text(lineheight = 0.65, face = 'italic'),
				panel.spacing = unit(0.05, 'lines'),
				panel.background = element_rect(fill = 'grey25'),
				strip.background = element_rect(fill = 'white'),
				strip.text.y = element_text(angle = 0),
				plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))


# Combine
final_plot <- plot_grid(main_panel, family_indicator_plot,
												nrow = 2, rel_heights = c(1, 5),
												align = 'v', axis = 'lr')



ggsave2(file.path('Plots', 'Supplement_TrainingSetRanks.pdf'), final_plot, width = 7, height = 6)

