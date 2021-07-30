## =================================================================================================
## Plot variable importance / model explanations
## =================================================================================================
set.seed(1521312)

library(scales)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(readxl)
library(apcluster)
library(cluster)
library(ggplot2)
library(cowplot)
library(egg)
library(ape)
library(ggtree)
library(phylogram)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'rundata_utils.R'))
source(file.path('Utils', 'plot_utils.R'))

BEST_RUN_ID <- 'AllGenomeFeatures_LongRun'

## Load data
# Metadata
zoo_status <- readRDS(file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))

taxonomy <- read_excel(file.path('ExternalData', 'ICTV_MasterSpeciesList_2018b.xlsx'),
											 sheet = 'ICTV 2018b Master Species #34 v', col_types = "text") %>% 
	mutate(Species = str_replace(.data$Species, '\u00A0$', '')) %>%   # Remove trailing non-breaking spaces in names
	select(.data$Species, .data$Family, GenomeType = .data$`Genome Composition`)

# Training data / features
data_training_viruses <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))
data_genome_features <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Virus.rds'))
data_dist_features <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Distances.rds'))

data_genome_features <- data_genome_features %>% 
 	rename_at(vars(-.data$UniversalName, -.data$Strain), ~ paste0('VirusDirect_', .))


# Cutoff:
CUTOFF <- readRDS(file.path('Plots', 'Intermediates', 'figure1_prob_cutoff.rds'))


# Features available for selection in the best run:
direct_features <- data_genome_features %>% 
	select(-.data$UniversalName, -.data$Strain) %>% 
	colnames()

mimicry_features <- data_dist_features %>% 
	select(-.data$UniversalName, -.data$Strain) %>% 
	colnames()

available_features <- c(direct_features, mimicry_features)


# Model fits:
model_fits <- read_run_rds(BEST_RUN_ID, '_ModelFits.rds')
iteration_data <- read_run_rds(BEST_RUN_ID, '_CalculatedData.rds')


# Bagged predictions:
bagged_preds <- read_run_rds(BEST_RUN_ID, '_Bagged_predictions.rds')

# Get latest spp names (training predictions, etc. currently contain only UniversalNames):
spp_names <- data_training_viruses %>% 
	distinct(.data$UniversalName, .data$Strain, .data$LatestSppName)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Top features, by feature set ---------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Get variable importance
#  - Variable importance / effect on train set predictions in each run (using training set, so these are 
#    conceptually similar to coefficients in a glm)
#  - Reporting the variation across iterations, not across viruses
#  			- Thus, for each level of summary reported (e.g. variable type, variable cluster, etc.),
#  			  calculate average importance across all train viruses in a given iteration
#  			- Then report mean + sd of variation of this total importance across runs...
train_data_splits <- lapply(iteration_data, function(x) x$train)
used_features <- colnames(train_data_splits$`1`)
used_features <- used_features[used_features != 'InfectsHumans']

# Get per-virus variable importance:
varimp_raw <- mapply(FUN = get_shap_contributions, 
										 trainedModel = model_fits, newdata = train_data_splits,
										 SIMPLIFY = FALSE)


## Calculate importance across all models
#  - Since we are primarily interested in the bagged model, calculate the mean effect size
#    across all iterations
stopifnot(length(available_features) == length(unique(available_features)))
n_total_features <- length(available_features)


featureset_importance <- bind_rows(varimp_raw) %>% 
	id_variable_types('Feature') %>% 
	
	group_by(.data$Feature, .data$VariableType, .data$Gene) %>% 
	summarise(mean_importance = mean(abs(.data$SHAP))) %>% 
	
	ungroup() %>% 
	mutate(Rank = rank(-.data$mean_importance, ties.method = 'random', na.last = 'keep'),
				 Rank = ifelse(is.na(.data$Rank), Inf, .data$Rank))
	

## Clean up names for plot:
featureset_importance <- featureset_importance %>% 
	ungroup() %>% 
	id_feature_set() %>% 
	mutate(SetName = factor(.data$SetName, levels = FEATURE_SET_ORDER, labels = FEATURE_SET_LABELS),
				 SetName_numeric = as.numeric(.data$SetName),
				 SetName_numeric = .data$SetName_numeric - min(.data$SetName_numeric) + 1)



featureset_imp_plot <- ggplot(featureset_importance, aes(x = SetName, y = Rank, colour = SetName)) +
	geom_blank() +  # stops next layer from converting x-axis to numeric
	geom_segment(aes(x = SetName_numeric - 0.5, xend = SetName_numeric + 0.5, yend = Rank), colour = 'grey80') +
	
	geom_boxplot(outlier.size = 0.5, size = 1, fill = 'white', alpha = 0.3) +
	
	scale_colour_manual(values = FEATURE_SET_COLOURS_LB,
											guide = FALSE) +
	scale_y_reverse() +
	labs(x = 'Feature set', y = 'Feature rank') +
	PLOT_THEME




# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Change in rank -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Show change in rank for different representations of the same feature
rank_change <- featureset_importance %>% 
	select(-VariableType, -Gene) %>% 
	id_variable_types('Feature') %>% 
	mutate(Measure = paste(Measure, Location)) %>% 
	
	add_count(Measure, name = 'N_sets') %>% 
	filter(N_sets >= 2)

# Keep only features shared with the unreferenced genome feature set 
#  - All comparisons will be relative to this set
in_unreferenced <- rank_change %>% 
	filter(SetName == 'Viral genomic\nfeatures') %>% 
	pull(Measure)

rank_change <- rank_change %>% 
	filter(Measure %in% in_unreferenced)

# Split into separate pairwise comparisons
combos <- rank_change %>% 
	filter(SetName != 'Viral genomic\nfeatures') %>% 
	arrange(SetName) %>% 
	pull(SetName) %>% 
	as.character() %>% 
	unique() %>% 
	mapply(c, ., 'Viral genomic\nfeatures', SIMPLIFY = FALSE, USE.NAMES = FALSE)

names(combos) <- LETTERS[1:length(combos)]

rank_change <- lapply(combos, function(x) filter(rank_change, SetName %in% x)) %>% 
	bind_rows(.id = 'Combo') %>% 
	group_by(Combo) %>% 
	add_count(Measure, name = 'N_sets') %>% 
	filter(N_sets == 2)


# Identify which set is dominant for each feature (in each pairwise comparison)
rank_change <- rank_change %>% 
	group_by(Combo, Measure) %>% 
	mutate(Dominant = SetName[which.min(Rank)],
				 SetLabel = if_else(SetName == 'Viral genomic\nfeatures', 'Unref.', 'Sim.'))


rank_change_plot <- ggplot(rank_change, aes(x = SetLabel, y = Rank, colour = SetName, group = Measure)) +
	geom_line(aes(colour = Dominant)) +
	geom_point(size = 1.1) +
	
	facet_grid(cols = vars(Combo)) +
	
	scale_colour_manual(values = FEATURE_SET_COLOURS_LB,
											guide = FALSE) +
	scale_x_discrete(expand = expand_scale(add = 0.22)) +
	scale_y_reverse() +
	labs(x = 'Feature type', y = 'Feature rank  ') +  # Extra spaces to leave room for panel label
	PLOT_THEME +
	theme(strip.text = element_blank(),
				panel.spacing = unit(2.5, 'pt'),
				plot.margin = margin(t = 3, r = 5.5, b = 5.5, l = 5.5))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Feature importance - ranked clusters -------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Since individual features may be highly correlated, show importance for clusters of correlated 
#  features

## Cluster features by correlation:
# - Clustering only features which made it past the initial screening (and only viruses in the final dataset):
stopifnot(length(data_training_viruses$UniversalName) == length(unique(data_training_viruses$UniversalName)))

final_features <- left_join(data_genome_features, data_dist_features, by = c('UniversalName', 'Strain')) %>% 
	select(.data$UniversalName, .data$Strain, one_of(used_features)) %>% 
	filter(.data$UniversalName %in% data_training_viruses$UniversalName)


# - Do clustering:
feature_matrix <- final_features %>% 
	select(-.data$UniversalName, -.data$Strain) %>% 
	t()

clusterobj_features = apcluster(corSimMat(method = 'spearman'), feature_matrix, 
																details = T, includeSim = T)

clusterdata_features <- get_cluster_members(clusterobj_features) 


# - Make labels for each cluster:
clusterlabs_features <- clusterdata_features %>% 
	group_by(.data$Cluster, .data$Exemplar) %>% 
	summarise(Members = paste(.data$Member, collapse = ', ')) %>% 
	id_variable_types('Exemplar') %>% 
	rename(Exemplar_Measure = .data$Measure,
				 Exemplar_MeasureType = .data$MeasureType,
				 Exemplar_VariableType = .data$VariableType,
				 Exemplar_Gene = .data$Gene,
				 Exemplar_Location = .data$Location) %>% 
	select(-.data$IsGenomicFeature)


# Feature importance: cluster-level summary (across all viruses):
make_clustered_iteration_summary <- function(varimps, clusterdata = clusterdata_features) {
	# Pre-summarise individual iterations before we convert to a single dataframe (to reduce memory usage):
	# Calculating importance of each cluster in each iteration:
	clusterdata <- clusterdata %>% 
		mutate(Member = as.character(.data$Member))
	
	varimps %>% 
		left_join(clusterdata, by = c('Feature' = 'Member')) %>% 
		id_variable_types('Feature') %>% 
		group_by(.data$Exemplar, .data$Cluster, .data$RowID) %>%  
		summarise(TotalAbsSHAP = sum(abs(.data$SHAP))) %>%       # Total effect size of all features in each cluster, for each virus
		group_by(.data$Exemplar, .data$Cluster) %>% 
		summarise(MeanAbsoluteSHAP = mean(.data$TotalAbsSHAP))   # Mean effect size of each cluster, across all viruses
}

varimp_summary_clusters <- lapply(varimp_raw, make_clustered_iteration_summary) %>% 
	bind_rows(.id = 'Iteration') 


# Distribution of importance values across iterations:
varimp_summary_clusters <- varimp_summary_clusters %>% 
	group_by(.data$Exemplar, .data$Cluster) %>% 
	summarise(AbsSHAP_mean = mean(.data$MeanAbsoluteSHAP),
						AbsSHAP_se = sd(.data$MeanAbsoluteSHAP)/sqrt(n())) %>%   # std error
	ungroup() %>% 
	
	# Tidy labels for plotting:
	left_join(clusterlabs_features, by = c('Exemplar', 'Cluster')) %>% 
	mutate(Exemplar_Gene = if_else(.data$Exemplar_Gene == 'ISG', 'ISG', str_to_lower(.data$Exemplar_Gene)),
				 ExemplarLabel = str_replace(.data$Exemplar_Measure, '\\.', ' '),
				 ExemplarLabel = str_replace(.data$ExemplarLabel, 'Bias', 'bias'),
				 ExemplarLabel = str_replace(.data$ExemplarLabel, 'U', 'T'),  # Not all genomes RNA
				 ExemplarLabel = str_replace(.data$ExemplarLabel, '^br', 'Bridge '),
				 ExemplarLabel = str_replace(.data$ExemplarLabel, '^NonBr', 'Non-bridge '),
				 ExemplarLabel = if_else(.data$Exemplar_VariableType == 'GenomicDensity', 
				 												paste0(.data$ExemplarLabel, ' similarity (', .data$Exemplar_Gene, ')'),
				 												.data$ExemplarLabel))
	


## Rename clusters so the most important cluster is cluster 1, second-most important is cluster 2, etc:
cluster_ranks <- varimp_summary_clusters %>% 
	distinct(.data$Cluster, .data$Exemplar, .data$AbsSHAP_mean) %>% 
	mutate(Label = rank(-.data$AbsSHAP_mean),
				 Label = factor(.data$Label)) %>% 
	select(-.data$AbsSHAP_mean)

clusterdata_features <- clusterdata_features %>% 
	left_join(cluster_ranks, by = c('Cluster', 'Exemplar'))

varimp_summary_clusters <- varimp_summary_clusters %>% 
	left_join(cluster_ranks, by = c('Cluster', 'Exemplar')) %>% 
	arrange(.data$Label) %>%   # Exemplar labels need to be in same order as these labels
	mutate(ExemplarLabel = factor(.data$ExemplarLabel, levels = unique(.data$ExemplarLabel)))


## Add feature class (if cluster is purely one class)
cluster_types <- clusterdata_features %>% 
	id_variable_types('Member') %>% 
	id_feature_set() %>% 
	group_by(.data$Cluster) %>% 
	mutate(n_features_in_cluster = n()) %>% 
	group_by(.data$Cluster, .data$SetName) %>% 
	summarise(set_proportion = n()/unique(.data$n_features_in_cluster))


varimp_summary_proportions <- varimp_summary_clusters %>% 
	left_join(cluster_types, by = 'Cluster') %>% 
	mutate(AbsSHAP_mean = .data$AbsSHAP_mean * .data$set_proportion,
				 SetName = factor(.data$SetName, levels = FEATURE_SET_ORDER))


## Restrict to top 25
varimp_summary_clusters_top25 <- varimp_summary_clusters %>% 
	filter(.data$AbsSHAP_mean > 0) %>% 
	filter(as.numeric(.data$Label) <= 25)

varimp_summary_proportions <- varimp_summary_proportions %>% 
	filter(.data$AbsSHAP_mean > 0) %>% 
	filter(as.numeric(.data$Label) <= 25)

## Main plot:
cluster_labs <- levels(varimp_summary_clusters$Label)
exemplar_labs <- levels(varimp_summary_clusters$ExemplarLabel)

ranks_plot <- ggplot(varimp_summary_clusters_top25, aes(x = as.numeric(Label), y = AbsSHAP_mean)) +
	geom_col(aes(fill = SetName), position = 'stack', width = 1, data = varimp_summary_proportions) +
	geom_col(fill = NA, colour = LINE_COLOUR, width = 1) +
	geom_errorbar(aes(ymin = AbsSHAP_mean - AbsSHAP_se, ymax = AbsSHAP_mean + AbsSHAP_se), width = 0.4, size = 0.35) +
	
	labs(x = 'Feature cluster', y = 'Combined effect magnitude', fill = NULL) +
	
	scale_x_reverse(breaks = 1:length(cluster_labs), labels = cluster_labs,
									sec.axis = sec_axis(trans = ~ ., breaks = 1:length(cluster_labs), labels = exemplar_labs,
																			name = 'Exemplar'),
									expand = expand_scale(mult = c(0, 0))) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
	scale_fill_manual(values = FEATURE_SET_COLOURS, guide = FALSE) +
	coord_flip() +
	PLOT_THEME



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Viruses with similar explanations ----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Cluster viruses with similar explanations
#   (termed 'supervised clustering' in Lundberg et al. [arXiv:1802.03888])
#  - Clustering viruses based on average SHAP values across all iterations
#  - Clustering human infecting / non-infecting separately: we expect SHAP values to separate
#    the classes; what is of interest though is whether there are common patterns/explanations
#    learned across completely unrelated human viruses

## Summarise varimps to virus level
# - Accumulating values, since we cannot process the entire df of varimps from all iterations at 
#   once when number of iterations is high
shap_totals <- expand.grid(LatestSppName = unique(spp_names$LatestSppName),
													 Feature = used_features,
													 Total_SHAP = 0,
													 N = 0)

for (iteration in names(varimp_raw)) {
	train_names <- data.frame(LatestSppName = model_fits[[iteration]]$trainingNames,
														stringsAsFactors = FALSE) %>%
		mutate(RowID = 1:n())

	varimp_named <- varimp_raw[[iteration]] %>%
		left_join(train_names, by = 'RowID') %>%
		select(.data$LatestSppName, .data$Feature, .data$SHAP)

	shap_totals <- shap_totals %>%
		left_join(varimp_named, by = c('LatestSppName', 'Feature')) %>%
		mutate(Total_SHAP = if_else(is.na(.data$SHAP), .data$Total_SHAP,  # This virus not seen in this iteration
																.data$Total_SHAP + .data$SHAP),
					 N = if_else(is.na(.data$SHAP), .data$N, .data$N + 1)) %>%
		select(-.data$SHAP)
}

# Calculate mean
shapley_vals_individual <- shap_totals %>%
	group_by(.data$LatestSppName, .data$Feature) %>%
	filter(.data$N != 0) %>%  # Remove unused features
	summarise(SHAP_mean = .data$Total_SHAP / .data$N)  # Keeping directionality: if a feature changes sign for a single virus we do want it to impact mean effect size measured


## Clustering:
# - For viruses, using hierarchical clustering, since we don't need discrete clusters, but do
#   need a very good ordering
shapley_vals_wide <- shapley_vals_individual %>% 
	spread(key = .data$Feature, value = .data$SHAP_mean) %>% 
	ungroup()


clust_viruses <- list('Human-infecting' = NA, 'Non human-infecting' = NA)
clusterdata_viruses <- data.frame()

for (inf_type in names(clust_viruses)) {
	if (inf_type == 'Human-infecting') {
		classes <- c('Human virus', 'Zoonotic')
	} else {
		classes <- 'No human infections'
	}
	
	keep <- zoo_status %>% 
		filter(.data$Class %in% classes) %>% 
		pull(.data$LatestSppName)
	
	shapley_df <- shapley_vals_wide %>% 
		filter(.data$LatestSppName %in% keep)
	
	shapley_matrix <- shapley_df %>% 
		select(-.data$LatestSppName) %>% 
		as.matrix()
	
	rownames(shapley_matrix) <- shapley_df$LatestSppName
	
	temp_clust <- shapley_matrix %>% 
		scale() %>% 
		daisy(metric = 'euclidean') %>% 
		agnes(method = 'average')
	
	clust_viruses[[inf_type]] <- temp_clust
	
	clusterdata_viruses <- rbind(clusterdata_viruses,
															 data.frame(InfectsHumans = inf_type,
															 					  Member = factor(temp_clust$order.lab, 
															 					 								 levels = temp_clust$order.lab)))
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot genome type for similar viruses -------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
genome_types <- taxonomy %>% 
	distinct(.data$GenomeType, .data$Species)

make_cluster_plot <- function(inf_type, label_y = TRUE,
															clust_data = clusterdata_viruses, clust_obj = clust_viruses, 
															bag_data = bagged_preds, name_data = spp_names,
															genome_type = genome_types, zoo_status_df = zoo_status) {
	## Main panel:
	genome_type_arranged <- clust_data %>% 
		filter(.data$InfectsHumans == inf_type) %>% 
		left_join(genome_type, by = c('Member' = 'Species')) %>% 
		left_join(zoo_status_df, by = c('Member' = 'LatestSppName')) %>% 
		mutate(Member = factor(.data$Member, levels = levels(clust_data$Member)))
	
	genome_plot <- ggplot(genome_type_arranged, aes(x = Member, y = GenomeType, colour = Class, fill = Class)) +
		geom_tile() +
		scale_colour_manual(values = ZOONOTIC_STATUS_COLOURS, guide = FALSE) +
		scale_fill_manual(values = ZOONOTIC_STATUS_COLOURS, guide = FALSE) +
		labs(y = 'Genome') +
		PLOT_THEME +
		theme(axis.text.x = element_blank(),
					axis.title.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.text.y = element_text(size = 6.5))
	
	
	## Matching genealogy:
	virus_dendro <- as.dendrogram(clust_obj[[inf_type]]) %>% 
		as.phylo()
	
	dendro_panel <- ggtree(virus_dendro, aes(colour = Class), ladderize = FALSE) %<+% zoo_status +
		scale_x_reverse(expand = c(0.005, 0.005)) +
		scale_y_continuous(expand = c(0, 0)) +
		coord_flip() +
		scale_colour_manual(values = c(ZOONOTIC_STATUS_COLOURS[1:3]), na.value = 'grey45', guide = FALSE) +
		theme_void() +
		theme(plot.margin = margin(t = 5.5, b = 0.3))
	
	
	## Predicted scores:
	inf_type2 <- if_else(inf_type == 'Human-infecting', 'True', 'False')
	
	bag_data <- bag_data %>%
		filter(.data$InfectsHumans == inf_type2) %>%
		left_join(name_data, by = "UniversalName") %>% 
		left_join(genome_type, by = c('LatestSppName' = 'Species')) %>% 
		left_join(zoo_status_df, by = "LatestSppName") %>% 
		mutate(LatestSppName = factor(.data$LatestSppName, levels = levels(clusterdata_viruses$Member)))

	score_plot <- ggplot(bag_data, aes(x = LatestSppName, y = BagScore, colour = Class, fill = Class)) +
		geom_col(width = 1, size = 0) +
		geom_hline(yintercept = CUTOFF, linetype = '22', colour = LINE_COLOUR) +
		scale_colour_manual(values = ZOONOTIC_STATUS_COLOURS, guide = FALSE) +
		scale_fill_manual(values = ZOONOTIC_STATUS_COLOURS, guide = FALSE) +
		labs(x = NULL, y = SCORE_LABEL_2LINE) +
		ylim(0, 1) +
		PLOT_THEME +
		theme(axis.text.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.text.y = element_text(size = 6.5))
	
	## Adjust axis labels and margins for displaying side-by-side
	if (label_y) {
		genome_plot <- genome_plot +
			theme(plot.margin = margin(t = 0, r = -5.5, b = 0, l = 5.5))
		
		score_plot <- score_plot +
			theme(plot.margin = margin(t = 1, r = -5.5, b = -1, l = 5.5))
		
	} else {
		genome_plot <- genome_plot +
			theme(axis.text.y = element_blank(),
						axis.title.y = element_blank(),
						axis.ticks.y = element_blank(),
						plot.margin = margin(t = 0, r = 5.5, b = 0, l = 1))
		
		score_plot <- score_plot +
			theme(axis.text.y = element_blank(),
						axis.title.y = element_blank(),
						axis.ticks.y = element_blank(),
						plot.margin = margin(t = 1, r = 5.5, b = -1, l = 1))
	}
	
	## Combine panels and return
	plot_grid(dendro_panel, genome_plot, score_plot,
						nrow = 3, rel_heights = c(0.4, 1, 0.6),
						align = 'v', axis = 'lr')
}


clust_plot_positive <- make_cluster_plot('Human-infecting')

clust_plot_negative <- make_cluster_plot('Non human-infecting', label_y = FALSE)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine plots ------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Top row (cluster plots)
#  - Adding x-axis label here so it can be centered below both plots
axis_plot <- ggplot() + 
	xlab('Virus species (sorted by explanation similarity)') +
	PLOT_THEME +
	theme(plot.margin = margin())

axis_lab <- get_plot_component(axis_plot, "xlab-b")

top_row <- plot_grid(clust_plot_positive, clust_plot_negative,
										 ncol = 2, rel_widths = c(0.4, 0.7),  # widths roughly match ratio of pos:neg, but need extra space for y-labels in clust_plot_positive
										 labels = c('A', ''))

top_row <- plot_grid(top_row, axis_lab, nrow = 2, 
										 rel_heights = c(1, 0.15))


# Bottom row (feature importance)
bottom_left <- plot_grid(featureset_imp_plot, rank_change_plot,
												 nrow = 2, rel_heights = c(2, 1.2),
												 align = 'v', axis = 'lr',
												 labels = c('B', 'C'), vjust = c(1.5, 0.1))

bottom_row <- plot_grid(bottom_left, ranks_plot,
												ncol = 2, rel_widths = c(1, 1.2),
												labels = c('', 'D'))


## Final figure
combined_plot <- plot_grid(top_row, bottom_row,
													 nrow = 2, rel_heights = c(2, 3))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
out_dir <- file.path('Plots', 'Intermediates')

if (!dir.exists(out_dir))
	dir.create(out_dir, recursive = TRUE)

# Plot
ggsave2(file.path('Plots', 'Figure2.pdf'), combined_plot, width = 7.5, height = 4.6, units = 'in')



# Calculated data - needed for supplementary figures:
list(cluster_obj = clusterobj_features,
		 cluster_data = clusterdata_features) %>% 
	saveRDS(file.path(out_dir, 'figure2_feature_clusters.rds'))

saveRDS(varimp_raw, file.path(out_dir, 'figure2_varimp_raw.rds'))


feature_vals <- final_features %>%
	gather(-.data$UniversalName, -.data$Strain, key = 'Feature', value = 'FeatureValue') %>%
	left_join(spp_names, by = c('UniversalName', 'Strain')) %>%
	left_join(zoo_status, by = 'LatestSppName')

shapley_vals_individual <- shapley_vals_individual %>% 
	left_join(feature_vals, by = c('LatestSppName', 'Feature'))

saveRDS(shapley_vals_individual, file.path(out_dir, 'figure2_virus_shapley_vals.rds'))


saveRDS(featureset_importance, file.path(out_dir, 'figure2_feature_set_importance.rds'))


list(cluster_obj = clust_viruses,
		 cluster_data = clusterdata_viruses,
		 metadata = left_join(zoo_status, genome_types, by = c('LatestSppName' = 'Species'))) %>% 
	saveRDS(file.path(out_dir, 'figure2_virus_clusters.rds'))


## Save values underlying all panels in a human-readable format:
write_excel_csv(shapley_vals_individual, file.path('FigureData', 'fig2_a.csv'))

featureset_importance %>% 
	select(-SetName_numeric) %>% 
	mutate(SetName = str_replace_all(SetName, '\n', ' ')) %>% 
	write_excel_csv(file.path('FigureData', 'fig2_b.csv'))

rank_change %>% 
	arrange(Combo, Rank) %>% 
	mutate(SetName = str_replace_all(SetName, '\n', ' ')) %>% 
	select(Combo, Feature, Measure, SetName, SetLabel, mean_importance, Rank) %>% 
	write_excel_csv(file.path('FigureData', 'fig2_c.csv'))

varimp_summary_clusters %>% 
	select(Cluster, Exemplar, ExemplarLabel, AbsSHAP_mean, AbsSHAP_se) %>% 
	write_excel_csv(file.path('FigureData', 'fig2_d_clusters.csv'))

varimp_summary_proportions %>% 
	select(Cluster, SetName, set_proportion) %>% 
	write_excel_csv(file.path('FigureData', 'fig2_d_featuresets.csv'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values discussed in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat('Median feature ranks:\n')
featureset_importance %>% 
	group_by(VariableType, Gene) %>% 
	summarise(Median_rank = median(Rank)) %>% 
	print()


cat('\n\n', length(unique(rank_change$Measure)), 'measures occur in >1 set', 
		'(i.e., multiple representations of the same feature), making up',
		length(unique(rank_change$Feature)), 'features in total\n')


# See "Descriptive_Stats_FeatureClusters.R" for more summary statistics