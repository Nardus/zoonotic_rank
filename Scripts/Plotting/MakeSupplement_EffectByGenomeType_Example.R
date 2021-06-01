library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

source(file.path('Utils', 'plot_utils.R'))
source(file.path('Scripts/Plotting/PlottingConstants.R'))


shapley_vals_individual <- readRDS(file.path('Plots', 'Intermediates', 'figure2_virus_shapley_vals.rds'))
featureset_importance <- readRDS(file.path('Plots', 'Intermediates', 'figure2_feature_set_importance.rds'))
genome_types <- readRDS('Plots/Intermediates/figure1_merged_taxonomy.rds')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Example of an effect which works across genome types --------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
top_feature <- featureset_importance$Feature[featureset_importance$Rank == 1]

stopifnot(top_feature == 'GenomicDensity_Housekeeping_CTG.Bias_Coding')
top_feature_name <- 'CTG bias similarity (housekeeping)'

example_vals <- shapley_vals_individual %>% 
	filter(Feature == top_feature) %>% 
	left_join(genome_types, by = c('LatestSppName' = 'Species')) %>% 
	mutate(Class = if_else(Class == 'No human infections', as.character(Class), 'Infects humans'))

#ks.test(example_vals$FeatureValue[example_vals$Class == 'Infects humans'], example_vals$FeatureValue[example_vals$Class != 'InfectsHumans'])


axis_lims <- range(example_vals$FeatureValue)  # Make sure plots are aligned
cut_points <- quantile(example_vals$FeatureValue, probs = c(0.25, 0.75)) #c(0.55, 1.95)
colours <- ZOONOTIC_STATUS_COLOURS[1:2]
names(colours) <- c('No human infections', 'Infects humans')

hist_plot <- ggplot(example_vals, aes(x = FeatureValue, fill = Class)) +
	geom_histogram(bins = 150) +
	geom_vline(xintercept = cut_points, linetype = 2, colour = LINE_COLOUR) +
	facet_grid(rows = vars(Class)) +
	scale_fill_manual(values = colours, guide = FALSE) +
	xlim(axis_lims) +
	labs(y = 'Count') +
	PLOT_THEME +
	theme(axis.text.x = element_blank(),
				axis.title.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.spacing = unit(3, 'pt'),
				plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))


shap_plot <- ggplot(example_vals, aes(x = FeatureValue, y = SHAP_mean, colour = Class)) +
	geom_point(shape = 1) +
	geom_hline(yintercept = 0, linetype = 2, colour = LINE_COLOUR) +
	geom_vline(xintercept = cut_points, linetype = 2, colour = LINE_COLOUR) +
	scale_colour_manual(values = colours, guide = FALSE) +
	xlim(axis_lims) +
	labs(y = 'Contribution to log odds\n(SHAP value)') +
	PLOT_THEME +
	theme(axis.text.x = element_blank(),
				axis.title.x = element_blank(),
				axis.ticks.x = element_blank(),
				plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))




dist_plot <- ggplot(example_vals, aes(x = GenomeType, y = FeatureValue, fill = Class)) +
	geom_boxplot(colour = LINE_COLOUR, position = position_dodge(width = 0.9)) +
	geom_hline(yintercept = cut_points, linetype = 2, colour = LINE_COLOUR) +
	coord_flip() +
	#facet_grid(rows = vars(GenomeType), scales = 'free_y') +
	scale_fill_manual(values = colours, guide = FALSE) +
	ylim(axis_lims) +
	labs(x = 'Genome type', y = top_feature_name) +
	PLOT_THEME +
	theme(panel.spacing = unit(0, 'pt'),
				strip.text = element_blank())


p <- plot_grid(hist_plot, shap_plot, dist_plot,
							 nrow = 3, rel_heights = c(1.5, 1, 2),
							 align = 'v', axis = 'lr')

ggsave2(file.path('Plots', 'Supplement_EffectByGenomeType_Example.pdf'), p, width = 7, height = 8, units = 'in')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values for legend --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat("Cutpoints: c1 =", cut_points[1], "| c2 =", cut_points[2])

example_vals %>% 
	group_by(Class) %>% 
	summarise(below_c1 = sum(FeatureValue < cut_points[1]),
						above_c2 = sum(FeatureValue > cut_points[2])) %>% 
	ungroup() %>% 
	mutate(below_c1_prop = below_c1/sum(below_c1),
				 above_c2_prop = above_c2/sum(above_c2)) %>% 
	print()
