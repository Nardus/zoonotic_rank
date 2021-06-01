# Plot distribution of hosts by priority category
#  - Extends figure 3A/B

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scales)
library(mgcv)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'gam_plot_utils.R'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Data ---------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Load predictions (novel)
pred_col_types <- cols(Name = 'c', 
											 bagged_prediction = 'c', 
											 priority_category = 'c', 
											 .default = 'd')

predictions_novel <- read_csv(file.path('Predictions', 'NovelViruses.predictions.csv'),
															col_types = pred_col_types)

## Taxonomy for novel viruses:
novel_taxonomy <- read_excel(file.path('ExternalData', 'NovelViruses', 'ICTV_MasterSpeciesList_2019.v1.xlsx'),
														 sheet = 'ICTV2019 Master Species List#35', col_types = 'text')


## Hosts of novel viruses:
novel_hosts <- read_csv(file.path("InternalData", "NovelVirus_Hosts_Curated.csv"), 
												col_types = cols(.default = "c"))

novel_virus_metadata <- read_csv(file.path("ExternalData", "NovelViruses", "NovelViruses.csv"),
																 col_types = cols(.default = "c"))

novel_hosts <- novel_hosts %>% 
	left_join(novel_virus_metadata, by = c("accession" = "SequenceID")) %>% 
	select(-.data$accession, -.data$notes) %>% 
	distinct()   # Remove duplicates from seqmented viruses


## Filters applied in figure 3:
novel_hosts <- novel_hosts %>% 
	filter(is.na(.data$class) | .data$class %in% c("Mammalia", "Aves") | .data$order %in% c("Diptera", "Ixodida"))

predictions_novel <- predictions_novel %>%
	filter(.data$Name %in% novel_hosts$Name) %>% 
	filter(!.data$Name == "Vaccinia virus")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Clean-up -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Arrange genome types alphabetically, but group DNA / RNA genome types
get_genome_order <- function(all_genome_types) {
	genome_type_order <- unique(as.character(all_genome_types))
	dna <- genome_type_order[grepl('DNA', genome_type_order)]
	rna <- genome_type_order[!grepl('DNA', genome_type_order)]
	
	rev(c(sort(dna), sort(rna)))
}

novel_taxonomy <- novel_taxonomy %>% 
	mutate(GenomeType = factor(.data$`Genome Composition`,
														 levels = get_genome_order(.data$`Genome Composition`)))


## Prepare data
rank_data <- predictions_novel %>% 
	left_join(novel_taxonomy, by = c("Name" = "Species")) %>% 
	left_join(novel_hosts, by = "Name") %>% 
	rename(host_order = .data$order,
				 host_class = .data$class)

stopifnot(nrow(rank_data) == nrow(predictions_novel))

rank_data <- rank_data %>% 
	mutate(priority_category = factor(.data$priority_category, levels = names(PRIORITY_COLOURS)),
				 facet = case_when(is.na(.data$host_order) ~ "",
													 .data$host_order %in% c("Diptera", "Ixodida") ~ "Vector", 
													 TRUE ~ .data$host_class),
				 facet = factor(.data$facet, levels = c("Mammalia", "Aves", "Vector", "")),
				 
				 host_label = case_when(is.na(.data$host_order) ~ "Unknown",
																.data$host == "Homo sapiens" ~ "Human",
																.data$host_order == "Primates" ~ "Non-human primate",
																TRUE ~ .data$host_order)) %>% 
	add_count(.data$host_label, name = "n_by_host") %>% 
	add_count(.data$Family, name = "n_by_family") %>% 
	arrange(.data$n_by_host) %>% 
	mutate(host_label = factor(.data$host_label, levels = unique(.data$host_label))) %>% 
	arrange(.data$n_by_family) %>% 
	mutate(family_label = factor(.data$Family, levels = unique(.data$Family)))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot host overview -------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
p_host_count <- ggplot(rank_data, aes(x = host_label, fill = priority_category)) + 
	geom_bar() +
	coord_flip() +
	facet_grid(rows = vars(facet), scales = "free", space = "free") +
	scale_fill_manual(values = PRIORITY_COLOURS) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.01))) + 
	labs(y = "Count", x = "Sampled host order", fill = sprintf("%s:", PRIORITY_LABEL)) +
	PLOT_THEME +
	theme(strip.text = element_blank(),
				legend.position = "top",
				panel.spacing = unit(0.1, "lines"))


p_host_percent <- ggplot(rank_data, aes(x = host_label, fill = priority_category)) + 
	geom_bar(position = "fill") + 
	coord_flip() +
	facet_grid(rows = vars(facet), scales = "free", space = "free") +
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.01))) +
	labs(y = "Proportion", x = "Sampled host order", fill = "Priority") +
	PLOT_THEME +
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				strip.text.y = element_text(angle = 0),
				panel.spacing = unit(0.1, "lines"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot family overview --------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Showing non-human-associated only
rank_data2 <- rank_data %>% 
	filter(.data$host_label != "Human")

p_family_count <- ggplot(rank_data2, aes(x = family_label, fill = priority_category)) + 
	geom_bar() +
	coord_flip() +
	facet_grid(rows = vars(GenomeType), scales = "free", space = "free") +
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.01))) + 
	labs(y = "Count", x = "Family", fill = "Priority") +
	PLOT_THEME +
	theme(axis.text.y = element_text(face = "italic"),
				strip.text = element_blank(),
				panel.spacing = unit(0.1, 'lines'))


p_family_percent <- ggplot(rank_data2, aes(x = family_label, fill = priority_category)) + 
	geom_bar(position = "fill") + 
	coord_flip() +
	facet_grid(rows = vars(GenomeType), scales = "free", space = "free") +
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.01))) +
	labs(y = "Proportion", x = "Family", fill = "Priority") +
	PLOT_THEME +
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				strip.text.y = element_text(angle = 0),
				panel.spacing = unit(0.1, 'lines'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Save plots ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
extracted_legend <- get_legend(p_host_count)

p_host_count <- p_host_count +
	guides(fill = FALSE)

left_column <- plot_grid(p_host_count, 
												 p_family_count, 
												 ncol = 1, align = "v", axis = "lr",
												 labels = c("A", "B"))

right_column <- plot_grid(p_host_percent,
													p_family_percent,
													ncol = 1, align = "v", axis = "lr")

final_p <- plot_grid(left_column, right_column)

final_p <- plot_grid(extracted_legend, final_p, ncol = 1,
										 rel_heights = c(1, 20))


ggsave2("Plots/Supplement_NovelVirus_Hosts.pdf", final_p, width = 7, height = 8)
