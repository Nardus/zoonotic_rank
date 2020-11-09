## =================================================================================================
## Plot / compare virus clusters with taxonomy: is model simply re-constructing taxonomy?
## =================================================================================================
# Note: this script requires output from MakeFigure2.R, which should be run first 
# 			(this ensures all plots are always referring to the same feature importances & clusters)
set.seed(1521312)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(readxl)
library(cluster)
library(dendextend)
library(parallel)
library(EnvStats)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'taxonomy_utils.R'))
source(file.path('Utils', 'plot_utils.R'))

N_CORES <- 8

# Number of iterations to use when calculating empirical null distributions:
NULL_ITERATIONS <- 1000

# Distance based feature sets to include in feature-based clustering:
USE_ISG_FEATURES <- TRUE
USE_HOUSEKEEPING_FEATURES <- TRUE
USE_REMAINING_FEATURES <- TRUE


## Load data:
# Using only the official taxonomy here:
taxonomy_data <- read_excel(file.path('ExternalData', 'ICTV_MasterSpeciesList_2018b.xlsx'),
														sheet = 'ICTV 2018b Master Species #34 v', col_types = "text") %>% 
	mutate(Species = str_replace(.data$Species, '\u00A0$', ''))   # Remove trailing non-breaking spaces in names


data_training <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))
genome_features <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Virus.rds'))
dist_features <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Distances.rds'))


fig2_data <- readRDS(file.path('Plots', 'Intermediates', 'figure2_virus_clusters.rds'))
shapley_vals_individual <- readRDS(file.path('Plots', 'Intermediates', 'figure2_virus_shapley_vals.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Make taxonomy dendrograms ------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Extract included viruses:
included_viruses <- fig2_data$metadata %>% 
	filter(.data$LatestSppName %in% data_training$LatestSppName) %>% 
	rename(Species = .data$LatestSppName,
				 InfectionClass = .data$Class)

stopifnot(all(included_viruses$Species %in% taxonomy_data$Species))


# Interpolate missing levels
taxonomy_data <- taxonomy_data %>% 
	filter(.data$Species %in% included_viruses$Species) %>% 
	select(one_of(SEARCH_ORDER)) %>% 
	add_artificial_levels() %>% 
	left_join(included_viruses, by = 'Species')



# Do clustering
taxonomy_cluster_obj <- list()

for (inf_type in names(fig2_data$cluster_obj)) {
	if (inf_type == 'Human-infecting') {
		classes <- c('Human virus', 'Zoonotic')
	} else {
		classes <- 'No human infections'
	}
	
	sub_taxonomy <- taxonomy_data %>% 
		filter(.data$InfectionClass %in% classes)
	
	# Calculate distance:
	#   - Since all columns are categorical, using gower distance equates to dividing the number
	#     of shared levels by the total number of levels, i.e. we are getting the proportion
	#     of taxonomic levels which are shared between all species
	taxonomy_matrix <- sub_taxonomy %>% 
		select(one_of(SEARCH_ORDER)) %>% 
		mutate_all(as.factor) %>% 
		as.data.frame()
	
	rownames(taxonomy_matrix) <- sub_taxonomy$Species
	
	taxonomy_distances <- daisy(taxonomy_matrix, metric = 'gower')
	
	# Clustering
	sub_cluster <- agnes(taxonomy_distances, method = 'average')
	taxonomy_cluster_obj[[inf_type]] <- sub_cluster
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Alternative feature dendrogram (actual feature values) -------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Merge feature sets (if needed)
genome_features <- genome_features %>% 
	rename_at(vars(-.data$UniversalName, -.data$Strain), ~ paste0('VirusDirect_', .))

if (any(USE_ISG_FEATURES, USE_HOUSEKEEPING_FEATURES, USE_REMAINING_FEATURES)) {
	keep_sets <- c()
	
	if (USE_ISG_FEATURES) keep_sets <- c(keep_sets, 'ISG')
	if (USE_HOUSEKEEPING_FEATURES) keep_sets <- c(keep_sets, 'Housekeeping')
	if (USE_REMAINING_FEATURES) keep_sets <- c(keep_sets, 'Remaining')
	
	dist_names <- colnames(dist_features)
	dist_sets <- str_match(dist_names, '[[:alpha:]]+_([[:alpha:]]+)_')[, 2]
	
	keep_names <- dist_names[dist_sets %in% keep_sets]
	
	dist_features <- dist_features %>% 
		select(.data$UniversalName, .data$Strain, one_of(keep_names))
	
	genome_features <- genome_features %>% 
		left_join(dist_features, by = c('UniversalName', 'Strain'))
}


# Fix virus names to match other dendrograms:
name_matches <- data_training %>% 
	distinct(.data$UniversalName, .data$Strain, .data$LatestSppName) %>% 
	right_join(included_viruses, by = c('LatestSppName' = 'Species')) %>% 
	select(-.data$GenomeType)

genome_features <- genome_features %>% 
	left_join(name_matches, by = c('UniversalName', 'Strain')) %>% 
	filter(.data$LatestSppName %in% included_viruses$Species) %>% 
	select(-.data$UniversalName, -.data$Strain)


# Used features: any feature having an absolute shapley value > 0 for at least 1 virus:
# - Almost certainly this will be all 150 features making it past the initial screen, but
#   checking in case this changes
used_features <- shapley_vals_individual %>% 
	filter(abs(.data$SHAP_mean) > 0) %>% 
	pull(.data$Feature) %>% 
	unique()


# Do the clustering:
feature_cluster_obj <- list()

for (inf_type in names(fig2_data$cluster_obj)) {
	if (inf_type == 'Human-infecting') {
		classes <- c('Human virus', 'Zoonotic')
	} else {
		classes <- 'No human infections'
	}
	
	feature_subset <- genome_features %>% 
		filter(.data$InfectionClass %in% classes)
	
	genome_features_mat <- feature_subset %>% 
		select(one_of(used_features)) %>% 
		mutate_all(scale) %>% 
		as.matrix()
	
	rownames(genome_features_mat) <- feature_subset$LatestSppName
	
	# Distance and clustering:
	feature_distances <- genome_features_mat %>% 
		scale() %>% 
		daisy(metric = 'euclidean')
	
	temp_cluster <- agnes(genome_features_mat, method = 'average')
	feature_cluster_obj[[inf_type]] <- temp_cluster
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Compare clustering: depth ------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_bakers_gamma <- function(dend1, dend2, niter = NULL_ITERATIONS, ncores = N_CORES, 
															plot_theme = PLOT_THEME, calculated_values = NULL, arrow_ymax = 25,
															percentiles = c(0, 0.01, 0.1, 0.5, 0.9, 0.95, 0.99, 1)) {
	# Baker's Gamma: For each pair of points in each tree, calculate the highest possible k (number of 
	# clusters) for which the pair are still in the same cluster. Then calculate the Spearman 
	# correlation between these values across both trees. 
	# - Since this is a correlation, it ranges between [-1, 1]; values near 0 mean the dendrograms 
	# 	are not similar, while a value of 1 would mean they are perfectly correlated
	#
	# If this function has been called before, supply calculated_values to save time. If not 
	#  available, the null distribution and observed values will be calculated.
	# Returns: A list containing calculated values and the plot
	
	get_null <- function(i, dend1, dend2) {
		dend1_random <- sample.dendrogram(dend1, replace = FALSE)
		dend2_random <- sample.dendrogram(dend2, replace = FALSE)
		cor_bakers_gamma(dend1_random, dend2_random)
	}
	
	if (is.null(calculated_values)) {
		# Null distribution:
		gamma_null <- mclapply(1:niter, 
													 FUN = get_null, 
													 dend1 = dend1,
													 dend2 = dend2,
													 mc.cores = ncores)
		
		gamma_null <- data.frame(Type = 'Null distribution',
														 Values = as.numeric(gamma_null))
		
		# Observed
		gamma_observed <- cor_bakers_gamma(dend1, dend2)
		
		gamma_observed <- data.frame(Type = 'Observed',
																 Values = gamma_observed) %>% 
			mutate(ymin = demp(.data$Values, obs = gamma_null$Values))
		
	} else {
		gamma_null <- calculated_values$gamma_null
		gamma_observed <- calculated_values$gamma_observed
	}
	
	# Plot:
	percentile_labels <- percentiles
	percentiles <- quantile(gamma_null$Values, probs = percentile_labels)
	
	gamma_plot <- ggplot(gamma_null, aes(x = Values, fill = Type, colour = Type)) +
		geom_density() +
		geom_segment(aes(xend = Values, y = arrow_ymax, yend = ymin), data = gamma_observed,
								 size = 1, arrow = arrow(length = unit(0.5, 'lines'), type = 'closed')) +
		scale_x_continuous(name = 'Quantiles', breaks = percentiles, labels = percentile_labels, 
											 position = "top", minor_breaks = NULL,
											 sec.axis = sec_axis(trans = ~ ., name = "Baker's Gamma correlation coefficient")) +
		scale_fill_manual(values = c('Null distribution' = '#BBBBBB', 'Observed' = '#AA3377')) +
		scale_colour_manual(values = c('Null distribution' = '#BBBBBB', 'Observed' = '#AA3377')) +
		labs(y = 'Density', colour = NULL, fill = NULL) +
		plot_theme +
		theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))
	
	
	pval <- mean(gamma_null$Values >= gamma_observed$Values)
	cat(sprintf("Baker's gamma: %1.3f; Emperical p-value: %1.3f\n", gamma_observed$Values, pval))
	
	# Return:
	list(plot = gamma_plot,
			 calculated_values = list(gamma_null = gamma_null,
			 												  gamma_observed = gamma_observed))
}

## Plots
gamma_shapley_pos <- plot_bakers_gamma(dend1 = as.dendrogram(taxonomy_cluster_obj[["Human-infecting"]]),
																			 dend2 = as.dendrogram(fig2_data$cluster_obj[["Human-infecting"]]))

gamma_shapley_neg <- plot_bakers_gamma(dend1 = as.dendrogram(taxonomy_cluster_obj[["Non human-infecting"]]),
																			 dend2 = as.dendrogram(fig2_data$cluster_obj[["Non human-infecting"]]))


gamma_feature_pos <- plot_bakers_gamma(dend1 = as.dendrogram(taxonomy_cluster_obj[["Human-infecting"]]),
																			 dend2 = as.dendrogram(feature_cluster_obj[["Human-infecting"]]),
																			 arrow_ymax = 12.5, percentiles = c(0.01, 0.5, 0.99))

gamma_feature_neg <- plot_bakers_gamma(dend1 = as.dendrogram(taxonomy_cluster_obj[["Non human-infecting"]]),
																			 dend2 = as.dendrogram(feature_cluster_obj[["Non human-infecting"]]),
																			 arrow_ymax = 12.5, percentiles = c(0.01, 0.5, 0.99))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Compare clustering: order ------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_bk <- function(dend1, dend2, tax_data, n_iter = NULL_ITERATIONS, 
									 plot_theme = PLOT_THEME, calculated_values = NULL) {
	# Similarity of sub-clusters: B_k (the Fowlkes-Mallows index calculated for increasing numbers of 
	# clusters, k) measures the similarity the dendrograms as we move from the root outwards.
	#  - NOTE: not all families connect to a single root in the virus taxonomy, which means the
	#    minimum number of clusters is >> 1. Small values of k therefore make no sense, and the 
	#    clusters become arbitrary (since the dendrogram must be bifurcating, but in reality all
	#    these lineages split immediately at the root)
	# 
	# Returns: a list containing calculated values
	min_possible_k <- length(unique(tax_data$Realm))
	
	used_k_vals <- seq(min_possible_k, (nleaves(dend1) - 1), by = 1)
	
	Bk_null_values <- Bk_permutations(dend1, dend2, k = used_k_vals, R = n_iter, warn = TRUE)
	
	Bk_observed_values <- Bk(dend1, dend2, k = used_k_vals, include_EV = FALSE, warn = TRUE) %>% 
		unlist() %>% 
		data.frame(k = names(.),
							 Bk = .)
	
	# Return
	list(used_k_vals = used_k_vals,
			 Bk_null_values = Bk_null_values,
			 Bk_observed_values = Bk_observed_values)
}

# Calculations
tax_data_pos <- filter(taxonomy_data, InfectionClass %in% c('Human virus', 'Zoonotic'))
tax_data_neg <- filter(taxonomy_data, InfectionClass == 'No human infections')


bk_shapley_pos <- get_bk(dend1 = as.dendrogram(taxonomy_cluster_obj[["Human-infecting"]]),
												 dend2 = as.dendrogram(fig2_data$cluster_obj[["Human-infecting"]]),
												 tax_data = tax_data_pos)

bk_shapley_neg <- get_bk(dend1 = as.dendrogram(taxonomy_cluster_obj[["Non human-infecting"]]),
												 dend2 = as.dendrogram(fig2_data$cluster_obj[["Non human-infecting"]]),
												 tax_data = tax_data_neg)


bk_feature_pos <- get_bk(dend1 = as.dendrogram(taxonomy_cluster_obj[["Human-infecting"]]),
												 dend2 = as.dendrogram(feature_cluster_obj[["Human-infecting"]]),
												 tax_data = tax_data_pos)

bk_feature_neg <- get_bk(dend1 = as.dendrogram(taxonomy_cluster_obj[["Non human-infecting"]]),
												 dend2 = as.dendrogram(feature_cluster_obj[["Non human-infecting"]]),
												 tax_data = tax_data_neg)


## Make a combined plot
# Summarise Null distribution:
#  - k is discrete (integer only), so we do not really have confidence bands to show, merely
#    pointwise confidence intervals
#  - Calculating a one-sided CI, since Bk values can only be more than the null - being less 
#    clustered than random does not make sense
#    
#    
clean_null_distribution <- function(null_data, null_class, dendrogram_type) {
	bind_rows(null_data) %>% 
		gather(key = 'k', value = 'Bk') %>%
		mutate(k = as.numeric(.data$k)) %>% 
		mutate(Class = null_class,
					 Dendrogram = dendrogram_type)
}

null_shapley_neg <- clean_null_distribution(bk_shapley_neg$Bk_null_values, 'No human infections', 'SHAP-based')
null_shapley_pos <- clean_null_distribution(bk_shapley_pos$Bk_null_values, 'Infects humans', 'SHAP-based')
null_feature_neg <- clean_null_distribution(bk_feature_neg$Bk_null_values, 'No human infections', 'Feature-based')
null_feature_pos <- clean_null_distribution(bk_feature_pos$Bk_null_values, 'Infects humans', 'Feature-based')

Bk_null <- bind_rows(null_shapley_neg, null_shapley_pos, null_feature_neg, null_feature_pos) %>%
	group_by(.data$Class, .data$Dendrogram) %>% 
	mutate(n_tests = length(unique(.data$k)),
				 alpha_lvl = 0.05 / .data$n_tests) %>%   # Bonferroni correction
	group_by(.data$Class, .data$Dendrogram, .data$k) %>%
	summarise(Median = quantile(.data$Bk, probs = 0.5), 
						LowerQ = 0,  # Min possible value for Bk (not -Inf)
						UpperQ = quantile(.data$Bk, probs = 1 - unique(.data$alpha_lvl))) %>% 
	mutate(Type = 'Null',
				 Type = factor(.data$Type, levels = c('Null', 'Significant', 'NonSignificant')))


## Summarise observed values
obs_shapley_neg <- data.frame(Class = 'No human infections', Dendrogram = 'SHAP-based',
															bk_shapley_neg$Bk_observed_values)
obs_shapley_pos <- data.frame(Class = 'Infects humans', Dendrogram = 'SHAP-based',
															bk_shapley_pos$Bk_observed_values)
obs_feature_neg <- data.frame(Class = 'No human infections', Dendrogram = 'Feature-based',
															bk_feature_neg$Bk_observed_values)
obs_feature_pos <- data.frame(Class = 'Infects humans', Dendrogram = 'Feature-based',
															bk_feature_pos$Bk_observed_values)

Bk_observed <- bind_rows(obs_shapley_neg, obs_shapley_pos, obs_feature_neg, obs_feature_pos) %>% 
	mutate(k = as.numeric(as.character(.data$k))) %>%     
	left_join(Bk_null, by = c('Class', 'Dendrogram', 'k')) %>% 
	mutate(Type = if_else(.data$Bk > .data$UpperQ, 'Significant', 'NonSignificant'),
				 Type = factor(.data$Type, levels = c('Null', 'Significant', 'NonSignificant')))


## Cut-offs corresponding to main taxonomic levels:
tax_level_markers <- taxonomy_data %>% 
	gather(key = 'Level', value = 'Value', one_of(SEARCH_ORDER)) %>% 
	mutate(Class = if_else(InfectionClass %in% c('Zoonotic', 'Human virus'), 'Infects humans', 'No human infections')) %>% 
	group_by(Class, Level) %>% 
	summarise(k = length(unique(Value))) %>% 
	filter(Level %in% c('Realm', 'Phylum', 'Family', 'Subfamily', 'Genus', 'Species')) %>% 
	mutate(Level = if_else(Level == 'Subfamily', 'Sub-\nfamily', Level))


# Plot:
bk_plot <- ggplot(Bk_observed, aes(x = k, colour = Type)) +
	geom_segment(aes(xend = k), y = -Inf, yend = 0.87, linetype = 2, colour = LINE_COLOUR, data = tax_level_markers) +
	geom_text(aes(x = k, label = Level), y = 0.885, colour = LINE_COLOUR, data = tax_level_markers,
						 angle = 90, hjust = 0, size = 2.5, lineheight = 0.7) +
	geom_linerange(aes(ymin = LowerQ, ymax = UpperQ), size = 1, data = Bk_null) +
	geom_step(aes(y = Bk), group = 1, size = 0.8) +
	facet_grid(vars(Dendrogram), vars(Class), scales = 'free_x') +
	scale_colour_manual(values = c('Null' = '#BBBBBB', 
																 'NonSignificant' = '#4477AA',
																 'Significant' = '#AA3377'),
											labels = c('Null' = 'Null distribution',
																 'NonSignificant' = expression('Observed'~group('(', q >= 0.05, ')')),
																 'Significant' = expression('Observed'~group('(', q < 0.05, ')'))),
											drop = FALSE) +
	labs(x = 'Number of clusters (k)', y = 'Clustering similarity (Fowlkes-Mallows index)',
			 colour = NULL) +
	ylim(c(0, 1)) +
	scale_x_discrete(expand = expand_scale(mult = 0.025)) +
	PLOT_THEME +
	theme(legend.position = 'top')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Make tanglegrams ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Used for display only, as a quick and accessible way to show dissimilarity while also showing
# the two dendrograms used for the more quantitative analyses above

calculate_tanglegram <- function(dend1, dend2, metadata, colour_zoonoses = FALSE, aligned_dendrograms = NULL) {
	# Prepare a tanglegram object for plotting, colouring branches by family and highlighting 
	# zoonotic viruses
	# If colour_zoonoses = TRUE, the strands connecting dendrograms will be coloured red (zoonotic) /
	# blue (non-zoonotic). Otherwise, strands will be coloured by family (the default).
	# If a pre-aligned dendrogram object is given (aligned_dendrograms is not NULL), dend1 and dend2 
	# will be ignored, and colours added to the existing dendrom object
	
	# Rotate sub-clusters:
	if (is.null(aligned_dendrograms)) 
		aligned_dendrograms <- untangle(intersect_trees(dend1, dend2), method = 'step2side')
	
	# Colour branches by family:
	families <- metadata$Family
	names(families) <- metadata$Species
	
	#  - Order families by occurence in taxonomy dendrogram (so colours are neater):
	families <- families[labels(aligned_dendrograms[[1]])] %>%
		unique()
	
	family_colours <- gg_colors_shuffled(length(families))
	names(family_colours) <- families
	
	for (fam in families) {
		lbls <- metadata$Species[metadata$Family == fam]
		lbls <- lbls[lbls %in% labels(aligned_dendrograms[[2]])]
		
		aligned_dendrograms <- set(aligned_dendrograms, 'by_labels_branches_col',
															 value = lbls, TF_values = family_colours[fam])
	}
	
	# Also colour connectors:
	if (! colour_zoonoses) {
		# Match branch colours
		#  - Getting these from tree so the order is correct
		connector_colours <- get_leaves_branches_col(aligned_dendrograms[[1]])
	} else {
		zoo_lbls <- metadata$InfectionClass
		names(zoo_lbls) <- metadata$Species
		connector_class <- zoo_lbls[labels(aligned_dendrograms[[1]])]
		connector_colours <- ZOONOTIC_STATUS_COLOURS[connector_class]
	}
	
	# Return:
	list(aligned_dendrograms = aligned_dendrograms,
			 connector_colours = connector_colours)
}

# Prepare metadata:
all_metadata <- taxonomy_data %>% 
	select(.data$Species, .data$Family, .data$InfectionClass)

# Plot
shapley_tangle_pos <- calculate_tanglegram(dend1 = as.dendrogram(taxonomy_cluster_obj[["Human-infecting"]]),
																					 dend2 = as.dendrogram(fig2_data$cluster_obj[["Human-infecting"]]),
																					 metadata = all_metadata, colour_zoonoses = TRUE)

shapley_tangle_neg <- calculate_tanglegram(dend1 = as.dendrogram(taxonomy_cluster_obj[["Non human-infecting"]]),
																					 dend2 = as.dendrogram(fig2_data$cluster_obj[["Non human-infecting"]]),
																					 metadata = all_metadata, colour_zoonoses = TRUE)


feature_tangle_pos <- calculate_tanglegram(dend1 = as.dendrogram(taxonomy_cluster_obj[["Human-infecting"]]),
																					 dend2 = as.dendrogram(feature_cluster_obj[["Human-infecting"]]),
																					 metadata = all_metadata, colour_zoonoses = TRUE)

feature_tangle_neg <- calculate_tanglegram(dend1 = as.dendrogram(taxonomy_cluster_obj[["Non human-infecting"]]),
																					 dend2 = as.dendrogram(feature_cluster_obj[["Non human-infecting"]]),
																					 metadata = all_metadata, colour_zoonoses = TRUE)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
out_dir <- file.path('Plots', 'Intermediates')

if (!dir.exists(out_dir))
	dir.create(out_dir, recursive = TRUE)


## Plots:
# - Not saving gamma plots (report stats in text instead)

# Bk plots: 
ggsave2(file.path('Plots', 'Supplement_bk_plots.pdf'), bk_plot, width = 7, height = 6, units = 'in')



# Tanglegrams
#  - These are too complex for cowplot - need to combine them manually in latex
pdf(file.path('Plots', 'Supplement_tanglegram_features_positive.pdf'), width = 3.5, height = 4)
tanglegram(feature_tangle_pos$aligned_dendrograms, 
					 color_lines = feature_tangle_pos$connector_colours, 
					 faster = TRUE, match_order_by_labels = TRUE,
					 margin_inner = 0.4, lwd = 1)
dev.off()


pdf(file.path('Plots', 'Supplement_tanglegram_features_negative.pdf'), width = 3.5, height = 4)
tanglegram(feature_tangle_neg$aligned_dendrograms, 
					 color_lines = feature_tangle_neg$connector_colours, 
					 faster = TRUE, match_order_by_labels = TRUE,
					 margin_inner = 0.4, lwd = 1)
dev.off()


pdf(file.path('Plots', 'Supplement_tanglegram_shap_positive.pdf'), width = 3.5, height = 4)
tanglegram(shapley_tangle_pos$aligned_dendrograms, 
					 color_lines = shapley_tangle_pos$connector_colours, 
					 faster = TRUE, match_order_by_labels = TRUE,
					 margin_inner = 0.4, lwd = 1)
dev.off()

pdf(file.path('Plots', 'Supplement_tanglegram_shap_negative.pdf'), width = 3.5, height = 4)
tanglegram(shapley_tangle_neg$aligned_dendrograms, 
					 color_lines = shapley_tangle_neg$connector_colours, 
					 faster = TRUE, match_order_by_labels = TRUE,
					 margin_inner = 0.4, lwd = 1)
dev.off()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Save calculated objects --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
save_objs <- function(pos, neg, file_name, path = out_dir) {
	temp <- list(pos = pos, neg = neg)
	saveRDS(temp, file.path(path, file_name))
}

save_objs(pos = gamma_feature_pos, neg = gamma_feature_neg,
					'supplementary_figS2_gamma_feature.rds')
save_objs(pos = gamma_shapley_pos, neg = gamma_shapley_neg,
					'supplementary_figS2_gamma_shapley.rds')

save_objs(pos = bk_feature_pos, neg = bk_feature_neg,
					'supplementary_figS2_bk_feature.rds')
save_objs(pos = bk_shapley_pos, neg = bk_shapley_neg,
					'supplementary_figS2_bk_shapley.rds')

save_objs(pos = feature_tangle_pos, neg = feature_tangle_neg,
					'supplementary_figS2_feature_tangle.rds')
save_objs(pos = shapley_tangle_pos, neg = shapley_tangle_neg,
					'supplementary_figS2_shapley_tangle.rds')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values reported in text --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat('Highest observed correlation across cutpoints (b_k):\n')


cat('Corresponding taxonomy levels:\n')
taxonomy_k <- taxonomy_data %>% 
	gather(key = 'Level', value = 'Value', one_of(SEARCH_ORDER)) %>% 
	mutate(Class = if_else(InfectionClass %in% c('Zoonotic', 'Human virus'), 'Infects humans', 'No human infections')) %>% 
	group_by(Class, Level) %>% 
	summarise(k = length(unique(Value))) %>% 
	ungroup() %>% 
	arrange(Class, k) %>% 
	mutate(k = as.factor(k))

print(taxonomy_k)


cat('\n\nFeature-based, positive\n')
bk_feature_pos$Bk_observed_values %>% 
	filter(.data$Bk == max(.data$Bk)) %>% 
	left_join(taxonomy_k, by = 'k') %>% 
	print()

cat('\nFeature-based, negative\n')
bk_feature_neg$Bk_observed_values %>% 
	filter(.data$Bk == max(.data$Bk)) %>% 
	left_join(taxonomy_k, by = 'k') %>%
	print()

cat('\n\nSHAP-based, positive\n')
bk_shapley_pos$Bk_observed_values %>% 
	filter(.data$Bk == max(.data$Bk)) %>% 
	left_join(taxonomy_k, by = 'k') %>%
	print()

cat('\nSHAP-based, negative\n')
bk_shapley_neg$Bk_observed_values %>% 
	filter(.data$Bk == max(.data$Bk)) %>% 
	left_join(taxonomy_k, by = 'k') %>%
	print()
