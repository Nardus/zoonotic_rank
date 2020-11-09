##
## Supplementary figures: Show effects of individual features
## 

library(scales)
library(tidyverse)
library(apcluster)
library(cowplot)
library(ggbeeswarm)


source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'plot_utils.R'))

# Data
# Created by MakeFigure2.R:
clusters <- readRDS(file.path('Plots', 'Intermediates', 'figure2_feature_clusters.rds'))
clusterdata_features <- clusters$cluster_data

shapley_vals_individual <- readRDS(file.path('Plots', 'Intermediates', 'figure2_virus_shapley_vals.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Effect direction ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Shorten variable names
clean_levels <- levels(clusterdata_features$Member) %>% 
	if_else(startsWith(., 'GenomicDensity_ISG'), paste(., 'similarity (I)'), .) %>% 
	if_else(startsWith(., 'GenomicDensity_Housekeeping'), paste(., 'similarity (H)'), .) %>% 
	if_else(startsWith(., 'GenomicDensity_Remaining'), paste(., 'similarity (R)'), .) %>% 
	str_remove('GenomicDensity_ISG_') %>% 
	str_remove('GenomicDensity_Housekeeping_') %>% 
	str_remove('GenomicDensity_Remaining_') %>% 
	str_remove('VirusDirect_') %>% 
	str_replace('.Bias', '') %>% 
	#str_replace('\\.', ' ') %>% 
	str_replace('pU', 'pT') %>% 
	str_replace('_Coding', '') %>% 
	str_replace('_EntireSeq', ' (e)') %>% 
	str_replace('NonBr', 'nbr') %>% 
	str_replace('U', 'T')  # Not all genomes RNA

names(clean_levels) <- levels(clusterdata_features$Member)


# Ensure virus and feature plotting order matches clustering:
shapley_vals_individual <- shapley_vals_individual %>% 
	mutate(FeatureLabel = clean_levels[as.character(.data$Feature)],
				 FeatureLabel = factor(.data$FeatureLabel, levels = clean_levels)) %>% 
	left_join(clusterdata_features, by = c(Feature = 'Member'))


# Scale features values to a common range
shapley_vals_individual <- shapley_vals_individual %>% 
	group_by(.data$Feature) %>% 
	mutate(FeatureValue_Scaled = rescale(.data$FeatureValue, to = c(0, 1))) %>% 
	ungroup()


## For clarity, only show clusters in which at least one feature does something:
#shapley_vals_individual <- shapley_vals_individual %>% 
#	group_by(.data$Cluster) %>% 
#	filter(!all(abs(.data$SHAP) < 0.2))


## Wrap cluster panels into multiple columns:
##  - since clusters vary in size, try to balance the number of features shown in each column
feats_per_column <- length(unique(shapley_vals_individual$Feature))/4

wrapping <- clusterdata_features %>% 
	group_by(.data$Label) %>% 
	summarise(Nfeatures = n()) %>% 
	arrange(.data$Label) %>% 
	mutate(Column = NA_character_)

col_num <- 1
feats_added <- 0

for (i in 1:nrow(wrapping)) {
	if (feats_added > feats_per_column) {
		col_num <- col_num + 1
		feats_added <- 0
	}
	
	wrapping$Column[i] <- col_num
	feats_added <- feats_added + wrapping$Nfeatures[i]
}


shapley_vals_individual <- shapley_vals_individual %>% 
	left_join(wrapping, by = 'Label')


## Plot
plot_direction_column <- function(col_num, all_data) {
	col_data <- all_data %>% 
		filter(.data$Column == col_num)
	
	ggplot(col_data, aes(x = FeatureLabel, y = FeatureValue_Scaled, colour = SHAP_mean)) +
		geom_point() +
		geom_quasirandom(alpha = 0.5, size = 0.5) +
		
		scale_colour_gradient2(labels = function(x) sprintf("%+1.1f", x),
													 low = "#2166ac", mid = "#E8E8E8", high = "#b2182b",
													 limits = round_up(range(all_data$SHAP_mean)),
													 guide = guide_colourbar(direction = "horizontal", 
													 												title.position = "top",
													 												barwidth = unit(3, 'cm'))) +
		
		labs(y = 'Feature value\n(scaled)', colour = 'Effect on log odds') +
		coord_flip() +
		facet_grid(Label ~ Column, scales = 'free_y', space = 'free') +
		PLOT_THEME +
		theme(strip.background = element_rect(fill = 'white', colour = LINE_COLOUR),
					strip.text.x = element_blank(),
					axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
					axis.title.y = element_blank(),
					panel.spacing = unit(0, 'lines'),
					plot.margin = margin(t = 5.5, r = 0.5, b = 5.5, l = 0.5))
}

direction_columns <- lapply(1:4, plot_direction_column, all_data = shapley_vals_individual)


## Combine columns:
# - Extract legend and add to bottom of last column
common_legend <- get_legend(direction_columns[[1]])
direction_columns <- lapply(direction_columns, function(p) {p + theme(legend.position = 'none')})

last_column <- plot_grid(direction_columns[[4]], common_legend,
												 ncol = 1, rel_heights = c(10, 1))

# - Add margin at edges of plot, while keeping margins between panels small:
first_column <- direction_columns[[c(1)]] +
	theme(plot.margin = margin(t = 5.5, r = 0.5, b = 5.5, l = 5.5))

last_column <- last_column +
	theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.5))

# - Combine
direction_plot <- plot_grid(first_column, direction_columns[[2]], 
														direction_columns[[3]], last_column, 
														rel_widths = c(1.1, 1, 1, 1),
														nrow = 1, align = 'v', axis = 't')


ggsave2(file.path('Plots', 'SupplementaryFigure_EffectDirection.png'), direction_plot, 
				width = 7, height = 7, units = 'in')
