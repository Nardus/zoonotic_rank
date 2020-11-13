#
# Plot ranks for different sarbecovirus isolates
#  - Similar to figure 3D, but zoomed-in to a single species
#

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ape)
library(ggplot2)
library(treeio)
library(ggtree)
library(scales)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Data ---------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

pred_col_types <- cols(Name = 'c', 
											 bagged_prediction = 'c', 
											 priority_category = 'c', 
											 .default = 'd')

predictions_sarbeco <- read_csv(file.path('Predictions', 'sarbecovirus.predictions.csv'),
																col_types = pred_col_types)

sarbeco_tree <- read.newick(file.path('ExternalData', 'sarbecovirus', 'sarbeco_ml_phylogeny.treefile'))


## Get cutoff from training data
training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

CUTOFF <- sum(training_data$InfectsHumans)/nrow(training_data)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Shorten names ------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
clean_sars_names <- function(name, include_host = TRUE) {
	host_name <- str_replace(name, '^[[:alnum:]-_]+\\|([[:alnum:]-_]+)[\\S]+', '\\1') %>% 
		str_replace("_", " ") %>% 
		str_replace_all("Bat-", "") %>% 
		str_replace("^R ", "Rhinolophus ") %>% 
		#str_replace("^R ", "R. ") %>% 
		#str_replace("^Aselliscus", "A.") %>% 
		#str_replace("^Chaerephon", "C.") %>% 
		str_replace("Ferrumequinum", "ferrumequinum")
	
	name <- str_extract(name, '^[[:alnum:]-_]+')
	
	if (include_host) {
		return(paste0(name, ' (', host_name, ')'))
	} else {
		sars_names <- host_name[grep('SARS', host_name)]
		
		name[grep('SARS', host_name)] <- sars_names
	}
	
	name
}

# Clean names, and remove duplicated sequences
sarbeco_clean <- predictions_sarbeco %>% 
	mutate(ShortName = clean_sars_names(.data$Name))

## Add ranks and priority:
sarbeco_priorities <- sarbeco_clean %>% 
	mutate(Rank = rank(-.data$calibrated_score_mean, ties.method = 'min')) %>% 
	arrange(-.data$Rank) %>% 
	mutate(Priority = factor(.data$priority_category,
													 levels = c('Low', 'Medium', 'High', 'Very high')))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Phylogeny ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Clean up phylogeny
sarbeco_tree_clean <- sarbeco_tree
sarbeco_tree_clean$tip.label <- clean_sars_names(sarbeco_tree_clean$tip.label)

unused_tips <- sarbeco_tree_clean$tip.label[!sarbeco_tree_clean$tip.label %in% sarbeco_priorities$ShortName]
sarbeco_tree_clean <- drop.tip(sarbeco_tree_clean, unused_tips)

sarbeco_tree_clean <- drop.tip(sarbeco_tree_clean, clean_sars_names('BtKY72|Bat-R_spp|Kenya|KY352407|2007-10')) # removing outgroup


# Plot
phylo_plot <- ggtree(sarbeco_tree_clean) +
	geom_rootedge(rootedge = 0.01) +
	geom_tiplab(size = 2, align = TRUE) +
	geom_treescale(x = 0.02, y = 3, width = 0.05) +
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_discrete(expand = c(0.011, 0.011)) +
	theme(plot.background = element_blank(),
				plot.margin = margin(t = 5.5, r = -5.5, b = 5.5, l = 5.5))

label_order <- phylo_plot$data %>% 
	arrange(.data$y) %>% 
	pull(.data$label) %>% 
	na.omit()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Main panel ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sarbeco_priorities <- sarbeco_priorities %>% 
	filter(Name != 'BtKY72|Bat-R_spp|Kenya|KY352407|2007-10') %>% 
	mutate(ShortName = factor(.data$ShortName, levels = label_order))

sars_names <- sarbeco_priorities$ShortName[grep('SARS', sarbeco_priorities$Name)]
sars_marker_df <- data.frame(ShortName = sars_names, ymin = Inf, ymax = 0.9,
														 Priority = 'Medium')

sarbeco_rank_plot <- ggplot(sarbeco_priorities, aes(x = ShortName, y = calibrated_score_mean, fill = Priority)) +
	geom_hline(yintercept = CUTOFF, linetype = 2, colour = 'grey10') +
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.6, colour = LINE_COLOUR) +
	geom_point(colour = LINE_COLOUR, shape = 21, size = 1.6) +
	
	geom_segment(aes(xend = ShortName, y = ymin, yend = ymax), 
							 arrow = arrow(length = unit(0.2, 'lines')), 
							 lineend = 'butt', linejoin = 'mitre', size = 0.5,
							 colour = LINE_COLOUR,
							 data = sars_marker_df) +
	
	scale_fill_manual(values = PRIORITY_COLOURS, guide = FALSE) +
	scale_x_discrete(position = 'top') +
	scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
	scale_size_continuous(range = c(2, 6)) +
	labs(x = NULL, y = SCORE_LABEL) +
	coord_flip() +
	PLOT_THEME +
	theme(panel.grid.major.y = element_line(colour = 'grey92'),
				axis.text.y = element_text(size = 5))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine and save ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sarbeco_plot_named <- plot_grid(phylo_plot, sarbeco_rank_plot, 
																ncol = 2, rel_widths = c(2, 3),
																align = 'h', axis = 'tb')


ggsave2(file.path('Plots', 'Supplement_Sarbecovirus_ranks.pdf'), sarbeco_plot_named, width = 7, height = 5)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values in text -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

cat('\nSarbecovirus priorities:\n')
print(table(sarbeco_priorities$Priority))
