## 
## Descriptive statistics for feature clusters
## 
library(dplyr)
library(tidyr)
library(apcluster)
library(stringr)

N_TOP <- 25  # Number of clusters displayed in figure 2 (not currently used, but stats for this given at end)

source(file.path('Utils', 'plot_utils.R'))

## Load clusters
fig2_clusters <- readRDS('Plots/Intermediates/figure2_feature_clusters.rds')
clusterdata <- fig2_clusters$cluster_data

top_clusters <- clusterdata %>% 
	filter(as.numeric(as.character(.data$Label)) <= N_TOP)

final_features <- clusterdata %>% 
	id_variable_types('Member')

n_features <- length(final_features$Member)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Number of clusters
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
length(unique(clusterdata$Label))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Cluster composition
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

composition <- clusterdata %>% 
	id_variable_types('Member') %>% 
	mutate(MeasureType = str_remove(.data$MeasureType, 'Non-bridge '),
				 MeasureType = str_remove(.data$MeasureType, 'Bridge '),
				 MeasureType = str_to_sentence(.data$MeasureType)) %>% 
	group_by(.data$Label) %>% 
	summarise(MeasureTypes = paste(sort(unique(.data$MeasureType)), collapse = ', '),
						Ntypes = length(unique(.data$MeasureType)))


# Overview
counts <- composition %>% 
	group_by(.data$MeasureTypes) %>% 
	summarise(N = n()) %>% 
	arrange(desc(.data$N))

counts

# Proportion of clusters containing more than 1 type:
nclusters <- sum(counts$N)

composition %>% 
	filter(.data$Ntypes > 1) %>% 
	nrow()/nclusters


# Clusters containing all 3 types:
composition %>% 
	filter(.data$Ntypes == 3)


# Most important clusters: ('Label' indicates cluster rank)
View(composition)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Exemplars
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Calculate obs/exp for a given measure:
#  - Note that clusters with < 3 members are excluded here, since for them the choice of exemplar 
#  would be arbitrary (i.e. not representative of a true cluster center)
get_exemplar_ratio <- function(cluster_data, grouping_cols, all_feature_data = final_features) {
	dinuc_names <- c("Dinucleotide bias", "Bridge dinucleotide bias", "Non-bridge dinucleotide bias")
	cluster_data <- cluster_data %>% 
		add_count(.data$Cluster, name = 'N_members') %>% 
		filter(.data$N_members >= 3)
	
	n_clusters <- length(unique(cluster_data$Cluster))
	
	retained_freq <- all_feature_data %>% 
		filter(.data$Member %in% cluster_data$Member) %>%   # Remove features belonging to the small clusters removed above
		mutate(MeasureType_collapsed = if_else(.data$MeasureType %in% dinuc_names,
																					 "Dinucleotide bias", .data$MeasureType)) %>% 
		group_by_at(vars(grouping_cols)) %>% 
		summarise(N_retained = n())
	
	cluster_data %>% 
		distinct(.data$Exemplar) %>% 
		id_variable_types('Exemplar') %>% 
		mutate(MeasureType_collapsed = if_else(.data$MeasureType %in% dinuc_names,
																					 "Dinucleotide bias", .data$MeasureType)) %>% 
		group_by_at(vars(grouping_cols)) %>% 
		summarise(N = n()) %>% 
		left_join(retained_freq, by = grouping_cols) %>% 
		mutate(Prop = .data$N / n_clusters,
					 Expected_Prop = .data$N_retained / n_features,
					 Obs_exp = .data$Prop / .data$Expected_Prop) %>% 
		arrange(-.data$Prop)
}


## Feature origin (gene):
get_exemplar_ratio(clusterdata, 'Gene')  # All clusters

## Are the clusters we show in the figure representative of this pattern?
get_exemplar_ratio(top_clusters, 'Gene')


## Measure types:
# (all OERs ~ 1?)
get_exemplar_ratio(clusterdata, 'MeasureType_collapsed')
