#
# Plot a labelled version of figure 1D, showing AUC for each family
# 

library(dplyr)
library(readxl)
library(readr)
library(ggplot2)
library(cowplot)
library(pROC)  # Always load this before ModelMetrics (using ModelMetrics::auc() below, not pROC::auc())
library(ModelMetrics)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'rundata_utils.R'))

BEST_RUN_ID <- 'AllGenomeFeatures_LongRun'


training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))
merged_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))

predictions_final_bagged <- read_run_rds(BEST_RUN_ID, '_Bagged_predictions.rds')


# Arrange genome types alphabetically, but group DNA / RNA genome types
get_genome_order <- function(all_genome_types) {
	genome_type_order <- unique(as.character(all_genome_types))
	dna <- genome_type_order[grepl('DNA', genome_type_order)]
	rna <- genome_type_order[!grepl('DNA', genome_type_order)]
	
	rev(c(sort(dna), sort(rna)))
}

merged_taxonomy <- merged_taxonomy %>% 
	mutate(GenomeType = factor(.data$GenomeType,
														 levels = get_genome_order(.data$GenomeType)),
				 BroadGenomeType = factor(.data$BroadGenomeType,
				 												 levels = get_genome_order(.data$BroadGenomeType)))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- AUC by family: overview --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_auc_ci <- function(actual, predicted) {
	if (length(unique(actual)) == 1) {
		# AUC can't be calculated
		ci_result <- c(NA_real_, NA_real_, NA_real_)
	} else {
		ci_result <- ci.auc(response = actual, predictor = predicted,
												levels = c(FALSE, TRUE), direction = "<")
	}
	
	# ModelMetrics::auc() will work even in very small families, where pROC::ci.auc() returns NA
	# - Given the very small sample size, these may not be reliable though
	direct_auc <- auc(actual = actual, predicted = predicted)
	
	data.frame(AUC_median = ci_result[2],
						 AUC_lower = ci_result[1],
						 AUC_upper = ci_result[3],
						 AUC_direct = direct_auc)
}


families <- training_data %>% 
	distinct(.data$LatestSppName, .data$UniversalName) %>% 
	left_join(merged_taxonomy, by = c('LatestSppName' = 'Species'))

family_counts <- training_data %>% 
	left_join(merged_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
	group_by(.data$Family) %>% 
	summarise(N = n(),
						Npositive = sum(.data$InfectsHumans),
						Nnegative = sum(!.data$InfectsHumans))

family_auc <- predictions_final_bagged %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True') %>% 
	left_join(families, by = 'UniversalName') %>% 
	group_by(.data$Family, .data$BroadGenomeType) %>% 
	group_modify(~ get_auc_ci(actual = .$InfectsHumans, predicted = .$BagScore)) %>% 
	ungroup() %>% 
	left_join(family_counts, by = c('Family')) %>% 
	mutate(AUC_median = if_else(is.na(.data$AUC_median), -Inf, .data$AUC_median))


# Same legend used across all panels:
COLOUR_RANGE <- range(family_auc$Npositive)


## Correlation:
make_overview_plot <- function(data, infs = FALSE, colour_range = COLOUR_RANGE, 
															 show_shape_legend = FALSE, show_colour_legend = FALSE) {
	# colour_range must be specified, since we need both panels below to have the same legend
	
	shape_guide <- if (show_shape_legend) guide_legend() else FALSE
	colour_guide <- if (show_colour_legend) guide_colourbar() else FALSE
	
	
	p <- ggplot(data, aes(x = log10(N), y = AUC_median))
	
	if (!infs) {
		p <- p + 
			geom_hline(yintercept = 0.5, linetype = 2, colour = LINE_COLOUR) +
			geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper), colour = 'grey80', width = 0.025) +
			coord_cartesian(ylim = c(0, 1))
	}
	
	p + 
		geom_point(aes(colour = Npositive), size = 1.6, shape = 19) +  # Add colour behind plain (DNA) symbols
		geom_point(aes(shape = BroadGenomeType, fill = Npositive), size = 1.5) +
		scale_shape_manual(values = c('ssDNA' = 4, 'dsDNA' = 3, 'dsDNA-RT' = 8,
																	'ssRNA' = 21, 'dsRNA' = 23, 'ssRNA-RT' = 22),
											 guide = shape_guide) +
		scale_colour_viridis_c(direction = -1, limits = colour_range, guide = colour_guide) +
		scale_fill_viridis_c(direction = -1, limits = colour_range, guide = colour_guide) +
		xlim(0, 2.2) +
		labs(x = expression(log[10]*group('(', textstyle(Number~of~species), ')')), 
				 y = 'AUC', shape = 'Genome',
				 colour = 'Known human-\ninfecting species', fill = 'Known human-\ninfecting species') +
		PLOT_THEME +
		theme(legend.key.size = unit(0.6, 'lines'),
					legend.margin = margin(t = 5.5, r = 0, b = 2.5, l = -2.5))
}

## Main plot
noninfs_plot <- family_auc %>% 
	filter(!is.infinite(.data$AUC_median)) %>% 
	make_overview_plot() +
	theme(axis.text.x = element_blank(),
				axis.title.x = element_blank(),
				plot.margin = margin(t = 5.5, r = 0, b = 1, l = 5.5))

## Show families for which AUC can't be calulated
infs_plot <- family_auc %>% 
	filter(is.infinite(.data$AUC_median)) %>% 
	group_by(N) %>% 
	mutate(AUC_median = rank(.data$AUC_median, ties.method = 'random')) %>%   # Simply stack points: y-axis will be blank  
	make_overview_plot(infs = TRUE) +
	scale_y_continuous(expand = expand_scale(add = 0.6)) +
	theme(axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				plot.margin = margin(t = 1, r = 0, b = 5.5, l = 5.5))

## Combine
legend_colour <- make_overview_plot(family_auc, show_colour_legend = TRUE) %>% 
	get_legend()

legend_shape <- make_overview_plot(family_auc, show_shape_legend = TRUE) %>% 
	get_legend()

overview_plot <- plot_grid(noninfs_plot, infs_plot, 
													 nrow = 2, rel_heights = c(2, 1.2),
													 align = 'v', axis = 'lr')
overview_plot <- plot_grid(overview_plot, legend_colour, legend_shape,
													 ncol = 3, rel_widths = c(7, 1, 1),
													 align = 'h', axis = 't')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Details (families labelled) ----------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Mark families for which AUC cannot be calculated
# - Otherwise, it's difficult to distinguish families which actually just have low AUC
family_auc <- family_auc %>% 
	mutate(AUC_median = if_else(is.infinite(.data$AUC_median), +Inf, .data$AUC_median),
				 Npositive = if_else(is.infinite(.data$AUC_median), NA_integer_, .data$Npositive))
	

# Arrange genome types alphabetically, but group DNA / RNA genome types
genome_types <- merged_taxonomy %>% 
	distinct(.data$Family, .data$GenomeType)

family_auc <- family_auc %>% 
	left_join(genome_types, by = 'Family')


# Plot
detail_plot <- ggplot(family_auc, aes(x = Family, y = AUC_median, fill = Npositive, colour = Npositive)) +
	geom_col() +
	geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper), colour = 'grey50', width = 0.4) +
	geom_hline(yintercept = 0.5, linetype = 2, colour = LINE_COLOUR) +
	facet_grid(cols = vars(GenomeType), scales = 'free_x', space = 'free') +
	scale_fill_viridis_c(direction = -1, na.value = 'grey95', limits = COLOUR_RANGE, guide = FALSE) +
	scale_colour_viridis_c(direction = -1, limits = COLOUR_RANGE, guide = FALSE) +
	scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.03))) +
	scale_x_discrete(expand = expand_scale(add = c(1.7, 1.7))) +
	labs(y = 'Within-family ranking performance (AUC)') +
	PLOT_THEME +
	theme(legend.position = 'top',
				axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = 'italic'),
				strip.text = element_text(angle = 45),
				panel.spacing = unit(0.1, 'lines'))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine and save ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
p <- plot_grid(overview_plot, detail_plot, 
							 nrow = 2, rel_heights = c(1.2, 2),
							 labels = c('A', 'B'))

ggsave2(file.path('Plots', 'Supplement_family_auc.pdf'), p, width = 7, height = 5)

write_excel_csv(family_auc, file.path('FigureData', 's3_fig.csv'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values mentioned in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat('\nDoes AUC increase with size of family?\n')
fit <- family_auc %>% 
	filter(!is.infinite(.data$AUC_median)) %>% 
	lm(AUC_median ~ log10(N), data = .)

print(summary(fit))

cat('\nDoes AUC increase with number of positives in family?\n')
fit2 <- family_auc %>% 
	filter(!is.infinite(.data$AUC_median)) %>% 
	lm(AUC_median ~ log10(Npositive), data = .)

print(summary(fit2))


cat('\nNumber of human infecting viruses strongly correlated with family size: Spearman rank cor =', 
		cor(family_counts$N, family_counts$Npositive, method = 'spearman'), '| Pearson cor = ',
		cor(family_counts$N, family_counts$Npositive, method = 'pearson'),
		'\n')


# Predictable families:
cat('\nFamilies with AUC > 0.9 (in order of decreasing AUC):\n')
family_auc %>% 
	filter(.data$AUC_median > 0.9) %>% 
	arrange(desc(.data$AUC_median)) %>% 
	pull(.data$Family) %>% 
	print()


stopifnot(length(unique(family_auc$Family)) == nrow(family_auc))

cat('\n\n', 
		sum(family_auc$AUC_median > 0.5), 'out of', 
		sum(!is.infinite(family_auc$AUC_median)),
		'families for which AUC can be calculated have AUC > 0.5\n')


cat('\nFamilies with AUC <= 0.5 (in order of decreasing AUC):\n')
family_auc %>% 
	filter(!is.infinite(.data$AUC_median) & .data$AUC_median <= 0.5) %>% 
	arrange(desc(.data$AUC_median)) %>% 
	print()
