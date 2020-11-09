## =================================================================================================
## Plot virus ranks and 'missed zoonoses'
## =================================================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ape)
library(cluster)
library(ggplot2)
library(ggtree)
library(scales)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'rundata_utils.R'))


BEST_RUN_ID <- 'AllGenomeFeatures_LongRun'


## Load predictions (novel)
pred_col_types <- cols(Name = 'c', 
											 bagged_prediction = 'c', 
											 priority_category = 'c', 
											 .default = 'd')

predictions_novel <- read_csv(file.path('Predictions', 'NovelViruses.predictions.csv'),
															col_types = pred_col_types)

predictions_sarbeco <- read_csv(file.path('Predictions', 'sarbecovirus.predictions.csv'),
																col_types = pred_col_types)


## Get cutoff from training data
training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

CUTOFF <- sum(training_data$InfectsHumans)/nrow(training_data)


## Intermediates from other plots:
zoo_status <- readRDS(file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))
training_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))
known_priorities <- readRDS(file.path('Plots', 'Intermediates', 'figure1_trainingset_priorities.rds'))


## Taxonomy for novel viruses:
novel_taxonomy <- read_excel(file.path('ExternalData', 'NovelViruses', 'ICTV_MasterSpeciesList_2019.v1.xlsx'),
														 sheet = 'ICTV2019 Master Species List#35', col_types = 'text')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Ranks --------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Rank
novel_ranks <- predictions_novel %>% 
	mutate(Rank = rank(-.data$calibrated_score_mean, ties.method = 'min')) %>% 
	arrange(-.data$Rank) %>% 
	mutate(Name = factor(.data$Name, levels = .data$Name),
				 Priority = factor(.data$priority_category, 
				 									levels = rev(c('Low', 'Medium', 'High', 'Very high'))))

zoom_cutoff <- 25
zoom_start <- novel_ranks %>% 
	filter(Rank == 1) %>% 
	pull(Name)

zoom_stop <- novel_ranks %>% 
	filter(Rank == zoom_cutoff) %>% 
	pull(Name)


## Plots
rank_plot <- ggplot(novel_ranks, aes(x = Name, y = calibrated_score_mean, colour = Priority)) +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf, colour = NA), fill = 'grey90') +
	
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.6) +
	geom_point(shape = 20) +  # colour = 'grey10'
	geom_hline(yintercept = CUTOFF, colour = 'grey10', linetype = 2) +
	scale_y_continuous(limits = c(0, 1), expand = expand_scale(add = 0.005)) +
	scale_colour_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	labs(x = NULL, y = SCORE_LABEL) +
	PLOT_THEME +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				plot.margin = margin(t = 5.5, r = 5.5, b = 1, l = 5.5))


indicator_plot <- ggplot(novel_ranks, aes(x = Name, y = Priority, fill = Priority)) +
	geom_blank() +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
	geom_tile() +
	scale_x_discrete(expand = expand_scale(add = 0)) +
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	labs(x = 'Virus species', y = 'Priority') +
	PLOT_THEME +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.grid = element_blank(),
				plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))


rank_plot_combined <- plot_grid(rank_plot, indicator_plot,
																nrow = 2, rel_heights = c(3, 1.5),
																align = 'v', axis = 'lr')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Zoom in and label to viruses ---------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
novel_families <- novel_taxonomy %>% 
	select(.data$Species, .data$Family)

top_viruses <- novel_ranks %>% 
	filter(Rank <= zoom_cutoff) %>% 
	mutate(Species = as.character(.data$Name)) %>%  # Needed to avoid losing factor ordering during join to taxonomy
	left_join(novel_families, by = 'Species') %>% 
	mutate(Family = str_replace(.data$Family, 'viridae', '-\nviridae'))


top_virus_plot <-	ggplot(top_viruses, aes(x = Name, y = calibrated_score_mean, fill = Priority)) +
	geom_hline(yintercept = CUTOFF, colour = 'grey60') +
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.5, colour = LINE_COLOUR) +
	geom_point(colour = LINE_COLOUR, shape = 21, size = 2.5) +
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_y_continuous(limits = c(0, 1), expand = expand_scale(add = c(0.02, 0.02))) +
	labs(x = NULL, y = SCORE_LABEL) +
	coord_flip() +
	PLOT_THEME +
	theme(panel.grid.major.y = element_line(colour = 'grey92'),
				axis.text.y = element_text(face = 'italic'),
				plot.margin = margin(t = 5.5, r = 1, b = 5.5, l = 5.5))


family_indicator_plot <- ggplot(top_viruses, aes(x = Family, y = Name, fill = Priority)) +
		geom_tile(colour = NA) + 
		scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
		PLOT_THEME +
		theme(axis.text.y = element_blank(),
					axis.title.y = element_blank(),
					axis.ticks.y = element_blank(),
					axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.65, face = 'italic'),
					plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))


top_virus_plot_combined <- plot_grid(top_virus_plot, family_indicator_plot,
																		 ncol = 2, rel_widths = c(3, 1.3),
																		 align = 'h', axis = 'tb')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Coronavirus example ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Extract coronavirus predictions
training_corona <- known_priorities %>% 
	left_join(training_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
	filter(.data$Family == 'Coronaviridae') %>% 
	mutate(Name = case_when(.data$LatestSppName == 'Severe acute respiratory syndrome-related coronavirus' ~ 'SARS-CoV-1', 
													.data$LatestSppName == 'Middle East respiratory syndrome-related coronavirus' ~ 'MERS-CoV',
													TRUE ~ .data$LatestSppName),
				 InfectsHumans = .data$Class %in% c('Human\nvirus', 'Zoonotic')) %>% 
	select(Species = .data$LatestSppName, 
				 .data$Name, 
				 .data$InfectsHumans,
				 .data$BagScore, 
				 .data$BagScore_Lower, 
				 .data$BagScore_Upper, 
				 .data$Priority)

novel_corona <- novel_ranks %>% 
	left_join(novel_taxonomy, by = c('Name' = 'Species')) %>% 
	filter(.data$Family == 'Coronaviridae') %>% 
	mutate(Species = .data$Name,
				 Name = paste0(.data$Name, '*')) %>% 
	select(.data$Species, 
				 .data$Name,
				 BagScore = .data$calibrated_score_mean,
				 BagScore_Lower = .data$calibrated_score_lower, 
				 BagScore_Upper = .data$calibrated_score_upper, 
				 .data$Priority)

sars2 <- predictions_sarbeco %>% 
	filter(str_detect(.data$Name, 'SARS-CoV-2')) %>% 
	mutate(Species = 'Severe acute respiratory syndrome-related coronavirus',
				 Name = 'SARS-CoV-2',
				 InfectsHumans = TRUE) %>% 
	select(.data$Species, 
				 .data$Name,
				 .data$InfectsHumans,
				 BagScore = .data$calibrated_score_mean,
				 BagScore_Lower = .data$calibrated_score_lower, 
				 BagScore_Upper = .data$calibrated_score_upper, 
				 Priority = .data$priority_category)


# Merge and add latest taxonomy
stopifnot(all(training_corona$Species %in% novel_taxonomy$Species))

all_corona <- bind_rows(training_corona, novel_corona, sars2) %>% 
	left_join(novel_taxonomy, by = 'Species')


## Generate dendrogram from taxonomy
used_levels <- c('Family', 'Subfamily', 'Genus', 'Subgenus', 'Species')

corona_taxonomy <- all_corona %>% 
	select(one_of(used_levels)) %>% 
	mutate_all(as.factor) %>% 
	as.data.frame()

rownames(corona_taxonomy) <- all_corona$Name

corona_distances <- daisy(corona_taxonomy, metric = 'gower')
corona_clust <- agnes(corona_distances, method = 'average')

corona_phylo <- corona_clust %>% 
	as.hclust() %>% 
	as.phylo() %>% 
	di2multi()

taxonomy_plot <- ggtree(corona_phylo) +
	#geom_tiplab(size = 2, align = TRUE) +
	scale_x_continuous(expand = c(0, 0.005)) +
	scale_y_discrete(expand = c(0.011, 0.011)) +
	theme(plot.background = element_blank(),
				plot.margin = margin(t = 5.5, r = -5.45, b = 5.5, l = 10))

label_order <- taxonomy_plot$data %>% 
	arrange(.data$y) %>% 
	pull(.data$label) %>% 
	na.omit()


## Main panel
# Order should match dendrogram:
all_corona <- all_corona %>% 
	mutate(Name = factor(.data$Name, levels = label_order))

# Make species names italic (but not abbreviations)
format_spp <- function(name) {
	name <- as.character(name)
	if_else(name %in% c('MERS-CoV', 'SARS-CoV-1', 'SARS-CoV-2'), 
					sprintf('"%s"', name),
					sprintf('italic("%s")', name))
}

spp_labels <- format_spp(all_corona$Name)
names(spp_labels) <- as.character(all_corona$Name)
spp_labels <- sapply(spp_labels, function(x) parse(text = x))

# Mark human-infecting
human_infecting_names <- all_corona$Name[all_corona$InfectsHumans & !is.na(all_corona$InfectsHumans)]
human_infecting_marker_df <- data.frame(Name = human_infecting_names, ymin = Inf, ymax = 0.9,
																				Priority = NA)
novel_marker_df <- data.frame(Name = novel_corona$Name, BagScore = 1, Priority = NA)

# Plot
corona_rank_plot <- ggplot(all_corona, aes(x = Name, y = BagScore, fill = Priority)) +
	geom_blank() +
	geom_tile(aes(y = 0.5, height = 1), fill = 'grey80',
						 data = human_infecting_marker_df) +
	
	geom_hline(yintercept = CUTOFF, colour = 'grey60') +
	geom_errorbar(aes(ymin = BagScore_Lower, ymax = BagScore_Upper), width = 0.6, colour = LINE_COLOUR) +
	geom_point(colour = LINE_COLOUR, shape = 21, size = 1.6) +
	
	geom_segment(aes(xend = Name, y = ymin, yend = ymax),
							 arrow = arrow(length = unit(0.2, 'lines')),
							 lineend = 'butt', linejoin = 'mitre', size = 0.5,
							 colour = LINE_COLOUR,
							 data = human_infecting_marker_df) +
	
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_x_discrete(labels = spp_labels, position = 'top') +
	scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
	scale_size_continuous(range = c(2, 6)) +
	labs(x = NULL, y = SCORE_LABEL) +
	coord_flip() +
	PLOT_THEME +
	theme(panel.grid.major.y = element_line(colour = 'grey92'),
				axis.text.y = element_text(size = 5))



# Genus panel
corona_genera <- all_corona %>% 
	mutate(Genus = str_remove(.data$Genus, 'coronavirus'),
				 name_numeric = as.numeric(.data$Name)) %>%  # Factor already in correct order to match dendrogram plot
	group_by(.data$Genus) %>% 
	summarise(y_mid = median(.data$name_numeric),
						y_min = min(.data$name_numeric),
						y_max = max(.data$name_numeric))

genus_panel <- ggplot(corona_genera, aes(x = 1, xend = 1, y = y_min - 0.6, yend = y_max + 0.6, colour = Genus)) +
	geom_segment(size = 5) +
	geom_text(aes(y = y_mid, label = Genus), colour = 'grey20', fontface = 'italic', angle = 270, size = 2.5) +
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_continuous(position = 'right',
										 breaks = as.numeric(all_corona$Name),
										 labels = all_corona$Name,
										 expand = c(0, 0)) +
	scale_colour_manual(values = c('#d9d9d9', '#f7f7f7', '#d9d9d9', '#f7f7f7'), guide = FALSE) +
	theme_nothing() +
	theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = -5),
				panel.border = element_rect(colour = 'grey20'))


## Combine
corona_plot_combined <- plot_grid(taxonomy_plot, corona_rank_plot, genus_panel,
																	ncol = 3, rel_widths = c(0.4, 3, 0.15),
																	align = 'h', axis = 'tb')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
if (!dir.exists('Plots'))
	dir.create('Plots')

bottom_row <- plot_grid(top_virus_plot_combined, corona_plot_combined, 
												ncol = 2, rel_widths = c(1.5, 1),
												labels = c('B', 'C'))

final_plot <- plot_grid(rank_plot_combined, bottom_row,
												nrow = 2, rel_heights = c(1, 2),
												labels = c('A', ''))

ggsave2(file.path('Plots', 'Figure3.pdf'), final_plot, width = 7, height = 6)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Save a combined priority table -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
known <- known_priorities %>%
	mutate(dataset = 'Training data',
				 Class = str_replace_all(Class, '\\n', ' ')) %>% 
	left_join(training_taxonomy, by = c('LatestSppName' = 'Species')) %>%   # TODO: These family names may be out of date
	select(species = LatestSppName, 
				 family = Family,
				 dataset,
				 infection_class = Class, 
				 probability_mean = BagScore,
				 probability_lower = BagScore_Lower,
				 probability_upper = BagScore_Upper,
				 priority = Priority)

novel <- novel_ranks %>%
	mutate(dataset = 'Novel') %>% 
	left_join(novel_families, by = c("Name" = "Species")) %>% 
	select(species = Name, 
				 family = Family,
				 dataset,
				 probability_mean = calibrated_score_mean,
				 probability_lower = calibrated_score_lower,
				 probability_upper = calibrated_score_upper,
				 priority = Priority)

bind_rows(known, novel) %>% 
	arrange(desc(probability_mean)) %>% 
	write_excel_csv(file.path('Plots', 'Intermediates', 'combined_virus_ranks.csv'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values mentioned in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

cat('\nNovel viruses in very high priority category:\n')
high_priority <- novel_ranks %>% 
	filter(.data$Priority == 'Very high') %>% 
	left_join(novel_families, by = c('Name' = 'Species'))
	
print(table(high_priority$Family))

#print(high_priority)

cat('\n',  
		sum(novel_ranks$Priority == 'Low'),
		'novel viruses are classified as low priority\n\n')
