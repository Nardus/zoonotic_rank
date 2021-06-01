#
# Plot ranks for training set viruses 'misclassified' as high priority
# - These might be mislabelled in the literature instead
# 

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))

all_ranks <- read_csv(file.path('Plots', 'Intermediates', 'combined_virus_ranks.csv'))
merged_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))

zoonotic_status_raw <- readRDS(file.path('CalculatedData', 'ZoonoticStatus_Merged.rds'))

training_ranks <- all_ranks %>% 
	filter(dataset == 'Training data') %>% 
	left_join(merged_taxonomy, by = c('species' = 'Species')) %>% 
	arrange(probability_mean) %>% 
	mutate(species = factor(species, levels = species))


CUTOFF <- readRDS(file.path('Plots', 'Intermediates', 'figure1_prob_cutoff.rds'))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot high priority viruses not currently known to infect humans ----------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
missing_zoonoses <- training_ranks %>% 
	filter(infection_class == 'No known human infections' & priority == 'Very high')

# Do we have serological evidence of human infection for these viruses?
serological <- zoonotic_status_raw %>% 
	mutate(serological_detection = if_else(is.na(IndirectDetectionOnly), FALSE, IndirectDetectionOnly),
				 serological_detection = str_to_sentence(serological_detection),
				 serological_detection = factor(serological_detection, levels = c('True', 'False'))) %>% 
	select(species = LatestSppName, serological_detection)

missing_zoonoses <- missing_zoonoses %>% 
	left_join(serological, by = 'species') %>% 
	arrange(probability_mean) %>% 
	mutate(species = factor(species, levels = species))


# Main panel
sero_colours <- c('True' = '#fdae61', 'False' = '#66c2a5')

missing_zoonoses_plot <- ggplot(missing_zoonoses, aes(x = species, y = probability_mean, fill = serological_detection)) +
	geom_errorbar(aes(ymin = probability_lower, ymax = probability_upper), width = 0.5, colour = LINE_COLOUR) +
	geom_point(colour = LINE_COLOUR, shape = 21, size = 2.5) +
	geom_hline(yintercept = CUTOFF, linetype = 2, colour = 'grey10') +
	scale_fill_manual(values = sero_colours, guide = FALSE) +
	scale_y_continuous(limits = c(0, 1), expand = expand_scale(add = c(0.02, 0.02))) +
	labs(x = NULL, y = SCORE_LABEL) +
	coord_flip() +
	PLOT_THEME +
	theme(panel.grid.major.y = element_line(colour = 'grey92'),
				axis.text.y = element_text(face = 'italic'),
				plot.margin = margin(t = 5.5, r = 1, b = 5.5, l = 5.5))


# Indicate families
family_indicator_plot <- ggplot(missing_zoonoses, aes(x = Family, y = species, fill = serological_detection)) +
	geom_tile() + 
	scale_fill_manual(values = sero_colours) +
	PLOT_THEME +
	labs(fill = 'Serological\nevidence') +
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.65, face = 'italic'),
				plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))


# Combine
final_plot <- plot_grid(missing_zoonoses_plot, family_indicator_plot,
												ncol = 2, rel_widths = c(3, 1.3),
												align = 'h', axis = 'tb')



ggsave2(file.path('Plots', 'Supplement_HighPriority_MissingZoonoses.pdf'), final_plot, width = 7, height = 3)
