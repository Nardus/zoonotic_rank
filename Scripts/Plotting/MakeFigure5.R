## =================================================================================================
## Plot case study - coronaviruses
## =================================================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)
library(ape)
library(cluster)
library(ggplot2)
library(ggtree)
library(scales)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))


CUTOFF <- readRDS(file.path('Plots', 'Intermediates', 'figure1_prob_cutoff.rds'))


## Intermediates from other plots:
training_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))
known_priorities <- readRDS(file.path('Plots', 'Intermediates', 'figure1_trainingset_priorities.rds'))
novel_ranks <- readRDS(file.path('Plots', 'Intermediates', 'figure3_novel_ranks.rds'))

## Taxonomy for novel viruses:
novel_taxonomy <- read_excel(file.path('ExternalData', 'NovelViruses', 'ICTV_MasterSpeciesList_2019.v1.xlsx'),
														 sheet = 'ICTV2019 Master Species List#35', col_types = 'text')

## Sarbecovirus predictions:
pred_col_types <- cols(Name = 'c', 
											 bagged_prediction = 'c', 
											 priority_category = 'c', 
											 .default = 'd')

predictions_sarbeco <- read_csv(file.path('Predictions', 'sarbecovirus.predictions.csv'),
																col_types = pred_col_types)

sarbeco_tree <- read.newick(file.path('ExternalData', 'sarbecovirus', 'sarbeco_ml_phylogeny.treefile'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Main coronavirus example ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Extract coronavirus predictions
training_corona <- known_priorities %>% 
	left_join(training_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
	filter(.data$Family == 'Coronaviridae') %>% 
	filter(.data$LatestSppName != 'Severe acute respiratory syndrome-related coronavirus') %>%  # Replaced by an earlier isolate below
	mutate(Name = case_when(.data$LatestSppName == 'Middle East respiratory syndrome-related coronavirus' ~ 'MERS-CoV',
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
				 Name = case_when(.data$Name == 'NL63-related bat coronavirus strain BtKYNL63-9b' ~ 'NL63-related bat CoV',
				 								 .data$Name == 'Rhinolophus ferrumequinum alphacoronavirus HuB-2013' ~ 'R. ferrumequinum CoV HuB-2013',
				 								 .data$Name == 'Nyctalus velutinus alphacoronavirus SC-2013' ~ 'N. velutinus CoV SC-2013',
				 								 .data$Name == 'Myotis ricketti alphacoronavirus Sax-2011' ~ 'M. ricketti CoV Sax-2011',
				 								 TRUE ~ .data$Name)) %>% 
	mutate(Name = paste0(.data$Name, '*')) %>% 
	select(.data$Species, 
				 .data$Name,
				 BagScore = .data$calibrated_score_mean,
				 BagScore_Lower = .data$calibrated_score_lower, 
				 BagScore_Upper = .data$calibrated_score_upper, 
				 .data$Priority)

sars_variants <- predictions_sarbeco %>% 
	filter(str_detect(.data$Name, 'SARS-CoV-') | str_detect(.data$Name, 'RaTG13')) %>% 
	mutate(Species = 'Severe acute respiratory syndrome-related coronavirus',
				 Name = case_when(startsWith(.data$Name, 'HSZ-Cc') ~ 'SARS-CoV (HSZ-Cc)',
				 								 startsWith(.data$Name, 'Wuhan-Hu-1') ~ 'SARS-CoV-2 (Wuhan-Hu-1)',
				 								 startsWith(.data$Name, 'RaTG13') ~ 'SARS-related CoV (RaTG13)'),
				 InfectsHumans = .data$Name != 'SARS-related CoV (RaTG13)') %>% 
	select(.data$Species, 
				 .data$Name,
				 .data$InfectsHumans,
				 BagScore = .data$calibrated_score_mean,
				 BagScore_Lower = .data$calibrated_score_lower, 
				 BagScore_Upper = .data$calibrated_score_upper, 
				 Priority = .data$priority_category)


# Merge and add latest taxonomy
stopifnot(all(training_corona$Species %in% novel_taxonomy$Species))

all_corona <- bind_rows(training_corona, novel_corona, sars_variants) %>% 
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
	scale_y_discrete(expand = c(0.013, 0.013)) +
	theme(plot.background = element_blank(),
				plot.margin = margin(t = 5.5, r = -5.45, b = 5.5, l = 10))

label_order <- taxonomy_plot$data %>% 
	arrange(.data$y) %>% 
	pull(.data$label) %>% 
	na.omit()


## Main panel
# Order should match dendrogram:
all_corona <- all_corona %>% 
	mutate(Name = factor(.data$Name, levels = label_order),
				 Priority = factor(.data$Priority, levels = names(PRIORITY_COLOURS)))

# Make species names italic (but not abbreviations)
format_spp <- function(name) {
	shortened_names <- c('NL63-related bat CoV*', 'R. ferrumequinum CoV HuB-2013*', 
											 'N. velutinus CoV SC-2013*', 'M. ricketti CoV Sax-2011*', 
											 'MERS-CoV', 'SARS-CoV (HSZ-Cc)', 
											 'SARS-CoV-2 (Wuhan-Hu-1)', 'SARS-related CoV (RaTG13)')
	
	partial_italic <- c('R. ferrumequinum CoV HuB-2013*' = 'italic("R. ferrumequinum")~"CoV HuB-2013*"',
											'N. velutinus CoV SC-2013*' = 'italic("N. velutinus")~"CoV SC-2013*"',
											'M. ricketti CoV Sax-2011*' = 'italic("M. ricketti")~"CoV Sax-2011*"')
	
	name <- as.character(name)
	case_when(name %in% names(partial_italic) ~ partial_italic[name],
						name %in% shortened_names ~ sprintf('"%s"', name),
						TRUE ~ sprintf('italic("%s")', name))
}

spp_labels <- format_spp(all_corona$Name)
names(spp_labels) <- as.character(all_corona$Name)
spp_labels <- sapply(spp_labels, function(x) parse(text = x))

# Mark human-infecting
human_infecting_names <- all_corona$Name[all_corona$InfectsHumans & !is.na(all_corona$InfectsHumans)]
human_infecting_marker_df <- data.frame(Name = human_infecting_names, ymin = Inf, ymax = 0.9)
novel_marker_df <- data.frame(Name = novel_corona$Name, BagScore = 1)

# Plot
corona_rank_plot <- ggplot(all_corona, aes(x = Name, y = BagScore)) +
	geom_blank() +
	geom_tile(aes(y = 0.5, height = 1), fill = 'grey80',
						data = human_infecting_marker_df) +
	
	geom_hline(yintercept = CUTOFF, colour = 'grey60') +
	geom_errorbar(aes(ymin = BagScore_Lower, ymax = BagScore_Upper), width = 0.6, colour = LINE_COLOUR) +
	geom_point(aes(fill = Priority), colour = LINE_COLOUR, shape = 21, size = 1.6) +
	
	geom_segment(aes(xend = Name, y = ymin, yend = ymax),
							 arrow = arrow(length = unit(0.2, 'lines')),
							 lineend = 'butt', linejoin = 'mitre', size = 0.5,
							 colour = LINE_COLOUR,
							 data = human_infecting_marker_df) +
	
	scale_fill_manual(values = PRIORITY_COLOURS, drop = FALSE) +
	scale_x_discrete(labels = spp_labels, position = 'top') +
	scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
	scale_size_continuous(range = c(2, 6)) +
	labs(x = NULL, y = SCORE_LABEL, fill = sprintf("%s:", PRIORITY_LABEL)) +
	coord_flip() +
	PLOT_THEME +
	theme(panel.grid.major.y = element_line(colour = 'grey92'),
				axis.text.y = element_text(size = 5), 
				legend.position = 'top',
				legend.box.margin = margin(t = 2.5, r = 5.5, b = 0, l = 90))



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
# ---- Sarbecoviruses -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Shorten names
clean_sars_names <- function(name, include_host = TRUE) {
	host_name <- str_replace(name, '^[[:alnum:]-_]+\\|([[:alnum:]-_]+)[\\S]+', '\\1') %>% 
		str_replace("_", " ") %>% 
		str_replace_all("Bat-", "") %>% 
		str_replace("^R ", "Rhinolophus ") %>% 
		str_replace("Ferrumequinum", "ferrumequinum") %>% 
		str_replace("SARS-CoV-1", "SARS-CoV")
	
	name <- str_extract(name, '^[[:alnum:]-_]+')
	
	if (include_host) {
		name <- if_else(host_name %in% c("pangolin", "SARS-CoV-2", "SARS-CoV"),
										sprintf("'%s'~('%s')", name, host_name),
										sprintf("'%s'~(italic('%s'))", name, host_name))
	} else {
		sars_names <- host_name[grep('SARS', host_name)]
		
		name[grep('SARS', host_name)] <- sars_names
	}
	
	name
}

# Clean names, and remove duplicated sequences
sarbeco_clean <- predictions_sarbeco %>% 
	mutate(ShortName = clean_sars_names(.data$Name))

parsed_labels <- sarbeco_clean$ShortName
names(parsed_labels) <- sarbeco_clean$ShortName
parsed_labels <- sapply(parsed_labels, function(x) parse(text = x))

## Add ranks and priority:
sarbeco_priorities <- sarbeco_clean %>% 
	mutate(Rank = rank(-.data$calibrated_score_mean, ties.method = 'min')) %>% 
	arrange(-.data$Rank) %>% 
	mutate(Priority = factor(.data$priority_category,
													 levels = c('Low', 'Medium', 'High', 'Very high')))


## Phylogeny
# Clean up phylogeny
sarbeco_tree_clean <- sarbeco_tree
sarbeco_tree_clean$tip.label <- clean_sars_names(sarbeco_tree_clean$tip.label)

unused_tips <- sarbeco_tree_clean$tip.label[!sarbeco_tree_clean$tip.label %in% sarbeco_priorities$ShortName]
sarbeco_tree_clean <- drop.tip(sarbeco_tree_clean, unused_tips)

sarbeco_tree_clean <- drop.tip(sarbeco_tree_clean, clean_sars_names('BtKY72|Bat-R_spp|Kenya|KY352407|2007-10')) # removing outgroup


# Plot
phylo_plot <- ggtree(sarbeco_tree_clean) +
	geom_rootedge(rootedge = 0.01) +
	geom_tiplab(size = 2, align = TRUE) +
	geom_treescale(x = 0.08, y = 2, width = 0.05, fontsize = 2, offset = -1.15) +
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_discrete(expand = expand_scale(add = c(0.4, 0.6))) +
	theme(plot.background = element_blank(),
				plot.margin = margin(t = 5.5, r = -5.5, b = 5.5, l = 5.5))

label_order <- phylo_plot$data %>% 
	arrange(.data$y) %>% 
	pull(.data$label) %>% 
	na.omit() %>% 
	as.vector()


## Main panel
sarbeco_priorities <- sarbeco_priorities %>% 
	filter(Name != 'BtKY72|Bat-R_spp|Kenya|KY352407|2007-10') %>% 
	mutate(ShortName = factor(.data$ShortName, levels = label_order))

sars_names <- sarbeco_priorities$ShortName[grep('SARS', sarbeco_priorities$Name)]
sars_marker_df <- data.frame(ShortName = sars_names, ymin = Inf, ymax = 0.9,
														 Priority = 'Medium')

sarbeco_rank_plot <- ggplot(sarbeco_priorities, aes(x = ShortName, y = calibrated_score_mean, fill = Priority)) +
	geom_blank() +
	geom_tile(aes(y = 0.5, height = 1), fill = 'grey80',
						data = sars_marker_df) +
	
	geom_hline(yintercept = CUTOFF, colour = 'grey60') +
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.6, colour = LINE_COLOUR) +
	geom_point(colour = LINE_COLOUR, shape = 21, size = 1.6) +
	
	geom_segment(aes(xend = ShortName, y = ymin, yend = ymax), 
							 arrow = arrow(length = unit(0.2, 'lines')), 
							 lineend = 'butt', linejoin = 'mitre', size = 0.5,
							 colour = LINE_COLOUR,
							 data = sars_marker_df) +
	
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_x_discrete(labels = parsed_labels, position = 'top') +
	scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
	scale_size_continuous(range = c(2, 6)) +
	labs(x = NULL, y = SCORE_LABEL) +
	coord_flip() +
	PLOT_THEME +
	theme(panel.grid.major.y = element_line(colour = 'grey92'),
				axis.text.y = element_text(size = 5))


## Combine
sarbeco_plot_named <- plot_grid(phylo_plot, sarbeco_rank_plot, 
																ncol = 2, rel_widths = c(1, 3),
																align = 'h', axis = 'tb')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
out_dir <- file.path('Plots', 'Intermediates')

if (!dir.exists(out_dir))
	dir.create(out_dir, recursive = TRUE)


final_plot <- plot_grid(corona_plot_combined, sarbeco_plot_named,
												ncol = 2, rel_widths = c(1, 1),
												labels = c('A', 'B'))

ggsave2(file.path('Plots', 'Figure5.pdf'), final_plot, width = 7.5, height = 5)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values mentioned in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat("\n\nSarbecovirus ranks:")
print(table(predictions_sarbeco$priority_category))

cat('\nSarbecovirus priorities:\n')
print(table(sarbeco_priorities$Priority))

cat("\nHigh-ranking sarbecoviruses:")
predictions_sarbeco %>% 
	filter(!.data$priority_category %in% c("Low", "Medium")) %>% 
	print()

cat("\nCoronaviruses ranking at least as high as SARS-CoV-2:")
sars2_score <- all_corona$BagScore[all_corona$Name == "SARS-CoV-2 (Wuhan-Hu-1)"]

all_corona %>% 
	filter(.data$BagScore >= sars2_score) %>% 
	select(.data$Species, .data$Name, .data$InfectsHumans, .data$BagScore) %>% 
	arrange(.data$BagScore) %>% 
	print()
