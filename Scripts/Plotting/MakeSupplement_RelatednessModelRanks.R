#
# Plot an overview of ranks for training set viruses, using ranks from relatedness-based models
# - Very similar to "MakeSupplement_TrainingSetRanks.R", which produces the equivalent plot for
#   ranks from the main model
# 

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(readxl)
library(ggplot2)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'plot_utils.R'))

training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

merged_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))
zoo_status <- readRDS(file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))

taxonomy_predictions <- readRDS(file.path("RunData", "Taxonomy_LongRun", "Taxonomy_LongRun_Bagged_predictions.rds"))
pn_predictions <- readRDS(file.path("RunData", "PN_LongRun", "PN_LongRun_Bagged_predictions.rds"))


## Prepare data
genome_type_order <- unique(merged_taxonomy$GenomeType)
dna <- genome_type_order[grepl('DNA', genome_type_order)]
rna <- genome_type_order[!grepl('DNA', genome_type_order)]
genome_type_order <- rev(c(sort(dna), sort(rna)))  # RNA viruses at top of plot

zoo_status <- zoo_status %>% 
	rename(infection_class = .data$Class)

name_matches <- training_data %>% 
	select(.data$UniversalName, .data$LatestSppName)

stopifnot(n_distinct(name_matches$UniversalName) == nrow(name_matches))


## Add priority categories
prioritize <- function(lower_bound, upper_bound, median, cutoff) {
	stopifnot(length(cutoff) == 1)
	stopifnot(cutoff > 0 & cutoff < 1)
	
	p <- if_else(lower_bound > cutoff, 'Very high',
							 if_else(median > cutoff, 'High',
							 				if_else(upper_bound > cutoff, 'Medium',
							 								'Low')))
	factor(p, levels = c('Low', 'Medium', 'High', 'Very high'))
}


## Plots
plot_ranks <- function(prediction_data) {
	prediction_data <- prediction_data %>% 
		left_join(name_matches, by = "UniversalName") %>% 
		left_join(merged_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
		left_join(zoo_status, by = "LatestSppName")
	
	cutoff <- find_balanced_cutoff(observed_labels = prediction_data$InfectsHumans,
																 predicted_score = prediction_data$BagScore)
	
	cat("Using cutoff:", cutoff, "\n")
	
	ranks <- prediction_data %>% 
		mutate(priority = prioritize(lower_bound = .data$BagScore_Lower,
																 upper_bound = .data$BagScore_Upper,
																 median = .data$BagScore,
																 cutoff = cutoff)) %>% 
		
		arrange(.data$BagScore) %>% 
		mutate(species = factor(.data$LatestSppName, levels = .data$LatestSppName),
					 Family = factor(.data$Family, levels = sort(unique(.data$Family), decreasing = TRUE)),
					 priority = factor(.data$priority, levels = c('Low', 'Medium', 'High', 'Very high')),
					 infection_class = factor(.data$infection_class, 
					 												 levels = c('Human virus', 'Zoonotic', 'No human infections')),
					 GenomeType = factor(.data$GenomeType, levels = genome_type_order))
	
	
	## Main panel
	main_panel <- ggplot(ranks, aes(x = species, y = BagScore, colour = infection_class)) +
		geom_errorbar(aes(ymin = BagScore_Lower, ymax = BagScore_Upper), width = 0.5) +
		geom_step(group = 1, colour = 'grey10', size = 0.6) +
		geom_hline(yintercept = cutoff, linetype = 2, colour = 'grey10') +
		scale_colour_manual(values = ZOONOTIC_STATUS_COLOURS) +
		scale_y_continuous(expand = expand_scale(add = c(0.02, 0.02))) +
		labs(x = NULL, y = SCORE_LABEL_2LINE, colour = 'Current status') +
		PLOT_THEME +
		theme(panel.grid.major.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.text.x = element_blank(),
					plot.margin = margin(t = 5.5, r = 5.5, b = 1, l = 5.5))
	
	
	## Indicate families
	family_indicator_plot <- ggplot(ranks, aes(x = species, y = Family, fill = priority)) +
		geom_tile() + 
		facet_grid(rows = vars(GenomeType), scales = 'free', space = 'free') +
		scale_fill_manual(values = PRIORITY_COLOURS, drop = FALSE) +
		PLOT_THEME +
		labs(x = 'Virus species (ranked)', fill = 'Zoonotic potential') +
		theme(axis.text.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.text.y = element_text(lineheight = 0.65, face = 'italic'),
					panel.spacing = unit(0.05, 'lines'),
					panel.background = element_rect(fill = 'grey25'),
					strip.background = element_rect(fill = 'white'),
					strip.text.y = element_text(angle = 0),
					plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))
	
	
	plot_grid(main_panel, family_indicator_plot,
						nrow = 2, rel_heights = c(1, 5),
						align = 'v', axis = 'lr')
}


taxonomy_plot <- plot_ranks(taxonomy_predictions)
pn_plot <- plot_ranks(pn_predictions)


## Save
combined_plot <- plot_grid(taxonomy_plot, pn_plot,
													 nrow = 2, labels = c("A", "B"))

ggsave2(file.path('Plots', 'Supplement_RelatednessModelRanks.pdf'), combined_plot, width = 7, height = 8.5)