## =================================================================================================
## Plot case study - all non-training viruses
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
source(file.path('Utils', 'rundata_utils.R'))


BEST_RUN_ID <- 'AllGenomeFeatures_LongRun'

CUTOFF <- readRDS(file.path('Plots', 'Intermediates', 'figure1_prob_cutoff.rds'))


## Load predictions (novel)
pred_col_types <- cols(Name = 'c', 
											 bagged_prediction = 'c', 
											 priority_category = 'c', 
											 .default = 'd')

predictions_novel <- read_csv(file.path('Predictions', 'NovelViruses.predictions.csv'),
															col_types = pred_col_types)

training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))
novel_data <- read_csv(file.path('ExternalData', 'NovelViruses', 'NovelViruses.csv'), col_types = cols(.default = 'c'))

## Intermediates from other plots:
zoo_status <- readRDS(file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))
training_taxonomy <- readRDS(file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))
known_priorities <- readRDS(file.path('Plots', 'Intermediates', 'figure1_trainingset_priorities.rds'))


## Taxonomy for novel viruses:
novel_taxonomy <- read_excel(file.path('ExternalData', 'NovelViruses', 'ICTV_MasterSpeciesList_2019.v1.xlsx'),
														 sheet = 'ICTV2019 Master Species List#35', col_types = 'text')

## Hosts of novel viruses:
novel_hosts <- read_csv(file.path("InternalData", "NovelVirus_Hosts_Curated.csv"), 
												col_types = cols(.default = "c"))

novel_virus_metadata <- read_csv(file.path("ExternalData", "NovelViruses", "NovelViruses.csv"),
																 col_types = cols(.default = "c"))

novel_hosts <- novel_hosts %>% 
	right_join(novel_virus_metadata, by = c("accession" = "SequenceID")) %>% 
	select(-.data$accession, -.data$notes) %>% 
	distinct()   # Remove duplicates from seqmented viruses


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Clean-up -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# For holdout/novel virus case study, show only viruses sampled in mammals, birds, taxonomic 
# orders containing known vectors of arboviruses, or with sampled host unknown. These are the only
# realistic candidates which might be zoonotic.
novel_hosts <- novel_hosts %>% 
	filter(is.na(.data$class) | .data$class %in% c("Mammalia", "Aves") | .data$order %in% c("Diptera", "Ixodida")) %>% 
	filter(!.data$Name == "Vaccinia virus")

predictions_novel <- predictions_novel %>%
	filter(.data$Name %in% novel_hosts$Name)

stopifnot(nrow(novel_hosts) == nrow(predictions_novel))

# Genome types
get_genome_order <- function(all_genome_types) {
	genome_type_order <- unique(as.character(all_genome_types))
	dna <- genome_type_order[grepl('DNA', genome_type_order)]
	rna <- genome_type_order[!grepl('DNA', genome_type_order)]
	
	rev(c(sort(dna), sort(rna)))
}

novel_taxonomy <- novel_taxonomy %>% 
	mutate(GenomeType = factor(.data$`Genome Composition`,
														 levels = get_genome_order(.data$`Genome Composition`)))



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

stopifnot(all(novel_ranks$Name %in% novel_hosts$Name))

host_data <- novel_ranks %>% 
	left_join(novel_hosts, by = "Name") %>% 
	mutate(Name = factor(.data$Name, levels = levels(novel_ranks$Name)))


# Determine zoom region:
# - Zoomed plot shows only non-human-associated viruses, so need to find the area covered by 
#   the top 25 viruses sampled from other species
human <- host_data %>% 
	filter(.data$host == "Homo sapiens")

zoom_cutoff <- 25
zoom_start <- novel_ranks %>% 
	filter(!.data$Name %in% human$Name) %>% 
	top_n(n = 1, wt = -.data$Rank) %>%  # The top-ranked non-human-associated virus
	pull(.data$Name)

zoom_stop <- novel_ranks %>% 
	filter(!.data$Name %in% human$Name) %>% 
	top_n(n = zoom_cutoff, wt = -.data$Rank) %>% 
	top_n(n = 1, wt = .data$Rank) %>%  # The lowest-ranked virus among the top 25
	pull(.data$Name)


## Plots
rank_plot <- ggplot(novel_ranks, aes(x = Name, y = calibrated_score_mean, colour = Priority)) +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf, colour = NA), fill = 'grey90') +
	
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.6) +
	geom_step(group = 1, colour = LINE_COLOUR, size = 0.5) +
	
	geom_rug(sides = 't', colour = LINE_COLOUR, data = human) +
	
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
	labs(x = NULL, y = PRIORITY_LABEL_2LINE) +
	PLOT_THEME +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				axis.title.y = element_text(margin = margin(r = 8)),
				panel.grid = element_blank(),
				plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5))


# Show hosts as another indicator plot:
source_data <- host_data %>% 
	group_by(.data$order) %>% 
	mutate(N = n()) %>% 
	ungroup() %>% 
	mutate(facet = case_when(is.na(.data$class) ~ "D",
				 									 .data$class == "Mammalia" ~ "A",
				 									 .data$class == "Aves" ~ "B", 
				 									 TRUE ~ "C"),
				 source = case_when(is.na(.data$class) ~ "Unknown",
				 									  .data$class == "Aves" ~ "Bird",
				 									  .data$class == "Mammalia" & N < 5 ~ "Other mammal",
				 									  TRUE ~ .data$order)) %>% 
	select(.data$Name, .data$source, .data$facet, .data$Rank, .data$Priority)

# - Order: place "Other mammal" last in the mammal facet:
source_name_order <- sort(unique(source_data$source), decreasing = TRUE)
source_name_order <- c("Other mammal", 
											 source_name_order[source_name_order != "Other mammal"])

source_data <- source_data %>% 
	mutate(source = factor(.data$source, levels = source_name_order))

# - Plot
source_plot <- ggplot(source_data, aes(x = Name, y = source, fill = Priority)) +
	geom_blank() +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
	geom_tile() +
	
	facet_grid(rows = vars(facet), scales = "free", space = "free") + 
	
	scale_x_discrete(drop = FALSE, expand = expand_scale(add = 0)) +
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	labs(x = NULL, y = "Sampled host or\nvector group") +
	PLOT_THEME +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.grid = element_blank(),
				strip.background = element_blank(),
				strip.text = element_blank(),
				panel.spacing = unit(-0.5, 'pt'),
				plot.margin = margin(t = 0, r = 5.5, b = -2, l = 5.5))


# Combine
rank_plot_combined <- plot_grid(rank_plot, indicator_plot, source_plot, 
																ncol = 1, rel_heights = c(3, 1.15, 3),
																align = 'v', axis = 'lr')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Zoom in and label top non-human-associated viruses ---------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
novel_families <- novel_taxonomy %>% 
	select(.data$Species, .data$Family, .data$GenomeType)

top_viruses <- novel_ranks %>% 
	filter(!.data$Name %in% human$Name) %>% 
	top_n(n = zoom_cutoff, wt = -.data$Rank) %>% 
	mutate(Species = as.character(.data$Name)) %>%  # Needed to avoid losing factor ordering during join to taxonomy
	left_join(novel_families, by = 'Species')

top_host_orders <- host_data %>% 
	filter(.data$Name %in% top_viruses$Name) %>% 
	mutate(facet = case_when(.data$class == "Mammalia" ~ "Mammal",
													 .data$class == "Aves" ~ "Bird", 
													 TRUE ~ "Vector"),
				 source = if_else(is.na(.data$order), "Unknown", .data$order),
				 facet = factor(.data$facet, levels = c("Mammal", "Bird", "Vector")),
				 source = factor(.data$source, levels = sort(unique(.data$source))))


top_virus_plot <-	ggplot(top_viruses, aes(x = Name, y = calibrated_score_mean, fill = Priority)) +
	geom_hline(yintercept = CUTOFF, colour = 'grey60', linetype = 2) +
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.5, colour = LINE_COLOUR) +
	geom_point(colour = LINE_COLOUR, shape = 21, size = 1.5) +
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
	xlab("Family") +
	PLOT_THEME +
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.text.x = element_text(face = 'italic', angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.65),
				panel.grid.major.y = element_line(colour = 'grey92'),
				plot.margin = margin(t = 5.5, r = 1, b = 5.5, l = 0))


host_indicator_plot <- ggplot(top_host_orders, aes(x = source, y = Name, fill = Priority)) +
	geom_tile(colour = NA) + 
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	facet_grid(cols = vars(facet), scales = "free", space = "free") + 
	xlab("Sampled host or vector order") +
	PLOT_THEME +
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.65),
				panel.grid.major.y = element_line(colour = 'grey92'),
				strip.background = element_blank(),
				panel.spacing = unit(-0.5, 'pt'),
				plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
out_dir <- file.path('Plots', 'Intermediates')

if (!dir.exists(out_dir))
	dir.create(out_dir, recursive = TRUE)


bottom_row <- plot_grid(top_virus_plot, family_indicator_plot, host_indicator_plot,
												ncol = 3, rel_widths = c(2, 1.3, 1.3),
												align = 'h', axis = 'tb',
												labels = c('B', '', ''))


final_plot <- plot_grid(rank_plot_combined, bottom_row,
												nrow = 2, rel_heights = c(1, 1.25),
												labels = c('A', ''))

ggsave2(file.path('Plots', 'Figure3.pdf'), final_plot, width = 7.5, height = 6)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Save a combined priority table -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
known_accessions <- training_data %>% 
	select(.data$LatestSppName, .data$Accessions)

novel_accessions <- novel_data %>% 
	group_by(.data$Name) %>% 
	summarise(Accessions = paste(sort(.data$SequenceID), collapse = "; "))

known <- known_priorities %>%
	mutate(dataset = 'Training data',
				 Class = str_replace_all(Class, '\\n', ' ')) %>% 
	left_join(training_taxonomy, by = c('LatestSppName' = 'Species')) %>% 
	left_join(known_accessions, by = "LatestSppName") %>% 
	select(species = LatestSppName, 
				 family = Family,
				 dataset,
				 infection_class = Class, 
				 accession = Accessions,
				 probability_mean = BagScore,
				 probability_lower = BagScore_Lower,
				 probability_upper = BagScore_Upper,
				 priority = Priority)

novel <- novel_ranks %>%
	mutate(dataset = 'Novel') %>% 
	left_join(novel_families, by = c("Name" = "Species")) %>% 
	left_join(novel_accessions, by = "Name") %>% 
	select(species = Name, 
				 family = Family,
				 dataset,
				 accession = Accessions,
				 probability_mean = calibrated_score_mean,
				 probability_lower = calibrated_score_lower,
				 probability_upper = calibrated_score_upper,
				 priority = Priority)

bind_rows(known, novel) %>% 
	arrange(desc(probability_mean)) %>% 
	mutate(dataset = if_else(dataset == 'Novel', 'Holdout', dataset)) %>% 
	write_excel_csv(file.path('Plots', 'Intermediates', 'combined_virus_ranks.csv'))

write_rds(novel_ranks, file.path('Plots', 'Intermediates', 'figure3_novel_ranks.rds'))

write_rds(host_data, file.path('Plots', 'Intermediates', 'figure3_host_data.rds'))
write_rds(novel_taxonomy, file.path('Plots', 'Intermediates', 'figure3_novel_taxonomy.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values mentioned in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat('\nNovel viruses found in humans:\n')
human_tab <- table(human$Priority)

print(human_tab)
print(human_tab/sum(human_tab))


cat("\nNovel viruses found in humans which are predicted low priority:")
human %>%
	filter(.data$Priority == "Low") %>%
	pull(.data$Name) %>%
	as.character() %>%
	print()


cat("\nFamilies not in training data but associated with humans:")
train_fams <- training_data %>%
	left_join(training_taxonomy, by = c("LatestSppName" = "Species")) %>%
	pull(.data$Family) %>%
	unique()

new_fams <- human %>%
	left_join(novel_taxonomy, by = c("Name" = "Species")) %>% 
		filter(!.data$Family %in% train_fams)

table(new_fams$Family, new_fams$Priority) %>%
	print()


cat('\nNovel viruses with no link to humans in the very high priority category:\n')
novel_ranks %>%
	filter(!.data$Name %in% human$Name) %>%
	pull(.data$Priority) %>%
	table() %>%
	print()

high_priority <- novel_ranks %>%
	filter(!.data$Name %in% human$Name) %>%
	filter(.data$Priority %in% c('Very high')) %>%
	left_join(novel_families, by = c('Name' = 'Species'))

table(high_priority$Family) %>%
	sort() %>%
	print()
