## 
## Supplementary figure: Show feature clusters
## 
set.seed(1521312)

library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(apcluster)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'plot_utils.R'))

# Feature cluster data from figure 2:
fig2_clusters <- readRDS('Plots/Intermediates/figure2_feature_clusters.rds')
clusterobj_features <- fig2_clusters$cluster_obj
clusterdata_features <- fig2_clusters$cluster_data

varimp_raw <- readRDS('Plots/Intermediates/figure2_varimp_raw.rds')


## Summarise variable importance to feature level
# - Get importance across all viruses in each iteration
summarise_importance <- function(varimps) {
	varimps %>% 
		group_by(.data$Feature) %>% 
		summarise(MeanAbsSHAP = mean(abs(.data$SHAP)))
}

varimp_summary <- lapply(varimp_raw, summarise_importance) %>% 
	bind_rows(.id = 'Iteration')

# - Get average across iterations:
varimp_summary <- varimp_summary %>% 
	group_by(.data$Feature) %>% 
	summarise(MeanAbsSHAP = mean(MeanAbsSHAP))


## Get 2-d coordinates for each feature:
feature_coords <- (1 - clusterobj_features@sim) %>% 
	as.dist() %>% 
	cmdscale() %>% 
	data.frame(Member = rownames(.), .) %>% 
	rename(X = .data$X1, Y = .data$X2)

exemplar_coords <- feature_coords %>% 
	rename(Exemplar = .data$Member,
				 Exemplar_X = .data$X,
				 Exemplar_Y = .data$Y)


## Join with cluster info:
cluster_coords <- clusterdata_features %>% 
	rename(ClusterLabel = .data$Label) %>% 
	left_join(feature_coords, by = 'Member') %>% 
	left_join(exemplar_coords, by = 'Exemplar') %>% 
	left_join(varimp_summary, by = c('Member' = 'Feature')) 

message(sum(is.na(cluster_coords$ClusterLabel)), ' features belonging to clusters in ',
				'which all features were removed during pre-screening will not be shown')

## Clean labels
cluster_coords <- cluster_coords %>% 
	filter(!is.na(.data$ClusterLabel)) %>% 
	id_variable_types("Member") %>% 
	id_feature_set() %>% 
	rename(FeatureSet = .data$SetName) %>% 
	mutate(FeatureSet = factor(.data$FeatureSet, levels = FEATURE_SET_ORDER),
				 Label = str_remove(.data$Member, 'GenomicDensity_'),
				 Label = str_remove(.data$Label, 'VirusDirect_'),
				 Label = str_remove(.data$Label, 'ISG_'),
				 Label = str_remove(.data$Label, 'Housekeeping_'),
				 Label = str_remove(.data$Label, 'Remaining_'),
				 Label = str_remove(.data$Label, '_Coding'),
				 Label = str_remove(.data$Label, '.Bias'),
				 Label = str_replace(.data$Label, '_EntireSeq', '.e'),
				 Label = str_replace(.data$Label, 'U', 'T'),  # Not all genomes RNA
				 Label = str_replace(.data$Label, 'NonBr', 'n'),
				 Label = str_replace(.data$Label, 'br', 'b')) %>% 
	arrange(.data$MeanAbsSHAP)  # Ensure used features are plotted on top


## Exemplar points plotted again to ensure they are on top:
exemplar_coords <- cluster_coords %>% 
	filter(.data$Member == .data$Exemplar)


## Plot:
apply_plot_styling <- function(p) {
	p +
		scale_shape_manual(values = c('Viral genomic features' = 21, 
																	'Similarity to ISGs' = 23,
																	'Similarity to housekeeping genes' = 24,
																	'Similarity to remaining genes' = 22),
											 na.value = 13) +
		scale_fill_viridis_c(breaks = seq(0, 1, by = 0.025), direction = -1) +
		
		labs(shape = 'Feature class',
				 fill = "Mean effect magnitude") +
		PLOT_THEME +
		theme(axis.text = element_blank(),
					axis.title = element_blank(),
					axis.ticks = element_blank(),
					plot.margin = margin(t = 10.5, r = 5.5, b = 5.5, l = 10.5))
}

unscaled_plot <- ggplot(cluster_coords, aes(x = X, y = Y, fill = MeanAbsSHAP, group = ClusterLabel)) +
	geom_segment(aes(xend = Exemplar_X, yend = Exemplar_Y), colour = 'grey70') +
	geom_point(aes(shape = FeatureSet), size = 6, colour = 'white') +
	geom_label(aes(label = Label), size = 2, fill = 'white', colour = 'grey10', alpha = 0.8,
						 label.padding = unit(0.05, 'lines'), label.size = 0) +
	
	geom_point(aes(shape = FeatureSet), size = 7, colour = 'grey40', data = exemplar_coords) +
	geom_label(aes(label = Label), size = 2, fill = 'white', alpha = 0.8, fontface = 'bold',
						 label.padding = unit(0.05, 'lines'), label.size = 0, data = exemplar_coords) +
	
	facet_wrap(~ ClusterLabel, scales = 'free') +
	scale_x_continuous(expand = expand_scale(mult = c(0.25, 0.25))) +
	scale_y_continuous(expand = expand_scale(mult = c(0.25, 0.25))) +
	coord_cartesian(clip = 'off')


unscaled_plot <- apply_plot_styling(unscaled_plot) +
	guides(shape = guide_legend(override.aes = list(size = 2), ncol = 2),
				 fill = guide_colorbar(barwidth = unit(4.7, 'cm'))) + 
	theme(panel.border = element_blank(),
				panel.spacing.x = unit(0, 'lines'),
				strip.text = element_text(size = 8, hjust = 0, face = 'bold'),
				legend.position = c(0.85, 0.14), 
				legend.direction = 'horizontal',
				legend.justification = c("right", "top"),
				legend.box.just = "right")



## Also plot a version in which all clusters are on the same scale:
equal_scale_plot <- ggplot(cluster_coords, aes(x = X, y = Y, fill = MeanAbsSHAP, group = ClusterLabel)) +
	geom_segment(aes(xend = Exemplar_X, yend = Exemplar_Y), colour = 'grey70') +
	geom_point(aes(shape = FeatureSet), size = 2, colour = 'white') +
	geom_point(aes(shape = FeatureSet), size = 2.3, colour = 'grey40', data = exemplar_coords) +
	facet_wrap(~ ClusterLabel)

equal_scale_plot <- apply_plot_styling(equal_scale_plot) +
	theme(panel.spacing.x = unit(-0.18, 'lines'),
				panel.spacing.y = unit(0, 'lines'),
				strip.text = element_text(size = 8, face = 'bold'),
				legend.position = 'none')


## Combine and save
combined_plot <- plot_grid(unscaled_plot, equal_scale_plot,
													 nrow = 2, rel_heights = c(2, 1),
													 align = 'v', axis = 'lr',
													 labels = c('A', 'B'), hjust = -0.15, vjust = 1.1)

ggsave2(file.path('Plots', 'SupplementaryFigure_FeatureClusters.pdf'), combined_plot, 
				width = 7, height = 9, units = 'in')


