## =================================================================================================
## Plot general data overview and model performance for main text
## =================================================================================================
set.seed(1521312)

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)
library(ggplot2)
library(cowplot)
library(ModelMetrics)
library(ggsignif)
library(plotROC)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'rundata_utils.R'))


## Constants
# Order here determines plotting order:
RUN_NAMES <- c(Taxonomy = 'Taxonomic',
							 PN = 'Phylogenetic\nneighbourhood',
							 VirusDirect = 'Viral genomic\nfeatures',
							 ISG = 'Similarity to\nISGs',
							 Housekeeping = 'Similarity to\nhousekeeping\ngenes',
							 Remaining = 'Similarity to\nremaining\ngenes',
							 AllGenomeFeatures = 'All genome\ncomposition\nfeature sets',
							 FeatureSelection_Top125 = 'All feature\nsets')
RUN_IDS <- names(RUN_NAMES)

BEST_ID <- 'AllGenomeFeatures_LongRun'            # RunID for final model
BEST_NAME <- 'All genome\ncomposition\nfeature sets'	  			  # Name for this model (should match one the run_names above so colours match)

GUESS_NAME <- 'Taxonomy-based\nheuristic'  # Name to use for 'informed guessing model'


# Pairs of run names for which p-values should be plotted:
# - Each level contains combinations which can be plotted at the same height without overlapping
P_VAL_COMPARISONS <- list(
	Level1 = list(c('Taxonomy-based\nheuristic', 'Taxonomic'),
								c('Viral genomic\nfeatures', 'Similarity to\nISGs')),
	Level2 = list(c('Taxonomy-based\nheuristic', 'Phylogenetic\nneighbourhood'),
								c('Viral genomic\nfeatures', 'Similarity to\nhousekeeping\ngenes')),
	Level3 = list(c('All feature\nsets', 'All genome\ncomposition\nfeature sets'),
								c('Viral genomic\nfeatures', 'Similarity to\nremaining\ngenes')),
	Level4 = list(c('Viral genomic\nfeatures', 'All genome\ncomposition\nfeature sets'))
)

P_VAL_HEIGHTS = c(Level1 = 0.835, Level2 = 0.870, Level3 = 0.905, Level4 = 0.940)



## Other utility functions:
merge_and_name <- function(preds, guesses) {
	# Merge prediction and guess data, adding a runid and fixing names / plotting order:
	combined_run_ids <- c('FeatureSelection_Top150', 'AllGenomeFeatures', 'VirusDirect_ISG', 'VirusDirect_ISG_Housekeeping')
	
	guesses %>% 
		mutate(RunID = 'Guess') %>% 
		bind_rows(preds) %>% 
		mutate(RunName = if_else(.data$RunID == 'Guess', GUESS_NAME, RUN_NAMES[.data$RunID]),
					 RunName = factor(.data$RunName, levels = c(GUESS_NAME, RUN_NAMES)),
					 Facet = if_else(.data$RunID == 'Guess', 'A', 
					 								if_else(.data$RunID %in% combined_run_ids, 'C', 'B')))  # Separate Guess, individual runs and combined run
}


## Read raw data
AllData <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

TaxonomyData <- read_excel(file.path('ExternalData', 'ICTV_MasterSpeciesList_2018b.xlsx'),
													 sheet = 'ICTV 2018b Master Species #34 v', col_types = "text")
UnclassifiedTaxonomy <- read_csv(file.path('InternalData', 'Taxonomy_UnclassifiedViruses.csv'),
																 col_types = cols(.default = 'c'))


## Read predictions:
guesses_test_boot <- readRDS(file.path('RunData', 'TaxonomyHeuristic', 'Test_BootstrapPredictions.rds'))
guesses_test_bag <- readRDS(file.path('RunData', 'TaxonomyHeuristic', 'Test_BaggedPredictions.rds'))

predictions_test_boot <- read_multiple_runs(RUN_IDS, '_Predictions.rds') %>% 
	filter(.data$Dataset == 'test')


## Read final model predictions:
predictions_final_boot <- read_run_rds(BEST_ID, '_Predictions.rds') %>% 
	mutate(RunName = BEST_NAME)

predictions_final_bagged <- read_run_rds(BEST_ID, '_Bagged_predictions.rds') %>% 
	mutate(RunName = BEST_NAME)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Classify viruses by infection type / class -------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
reservoir_human <- AllData %>% 
	group_by(.data$LatestSppName) %>% 
	summarise(HumanReservoir = 'Human' %in% .data$Reservoir)

zoo_status <- AllData %>% 
	left_join(reservoir_human, by = 'LatestSppName') %>% 
	mutate(Class = if_else(.data$HumanReservoir, 'Human virus',
												 if_else(.data$InfectsHumans, 'Zoonotic', 
												 				'No human infections'))) %>% 
	mutate(Class = factor(.data$Class, levels = rev(c('Human virus', 'Zoonotic', 'No human infections')))) %>% 
	distinct(.data$LatestSppName, .data$Class)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Merge taxonomies (in case unclassified viruses are included) -------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
TaxonomyData <- TaxonomyData %>% 
	mutate(Species = str_replace(.data$Species, '\u00A0$', '')) %>%   # Remove trailing non-breaking spaces in names
	select(.data$Species, .data$Family, GenomeType = .data$`Genome Composition`)

# Summarise genome types:
# - Showing the majority genome type (Phenuiviridae have both ssRNA(+/-) and ssRNA(-) isolates)
genome_types <- TaxonomyData %>% 
	group_by(.data$Family, .data$GenomeType) %>% 
	summarise(N = n()) %>% 
	top_n(n = 1, wt = N) %>% 
	ungroup() %>% 
	select(-.data$N) %>% 
	mutate(BroadGenomeType = str_replace(.data$GenomeType, '\\([-+/]+\\)', ''),   # Removing direction to simplify plots
				 BroadGenomeType = factor(.data$BroadGenomeType, levels = c('ssDNA', 'dsDNA', 'dsDNA-RT',
				 																													 'ssRNA', 'dsRNA', 'ssRNA-RT')))

stopifnot(nrow(genome_types) == length(unique(genome_types$Family)))  # Each family should now have just one genome type


## Merge
TaxonomyData <- TaxonomyData %>% 
	select(-.data$GenomeType)

merged_taxonomy <- UnclassifiedTaxonomy %>% 
	select(Species = .data$UniversalName, .data$Family) %>% 
	filter(! .data$Species %in% TaxonomyData$Species) %>% 
	bind_rows(TaxonomyData) %>% 
	left_join(genome_types, by = 'Family')

stopifnot(nrow(merged_taxonomy) == length(unique(merged_taxonomy$Species)))  # Each species should have a single entry


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Distribution of AUC values -----------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Clean up and merge bootstrap predictions
predictions_test_boot <- predictions_test_boot %>% 
	filter(.data$Dataset == 'test') %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True')

guesses_test_boot <- guesses_test_boot %>% 
	rename(RawScore = .data$PropHumanInFamily)


combined_test_boot <- merge_and_name(preds = predictions_test_boot,
																		 guesses = guesses_test_boot)

# Calculate AUC
auc_test_boot <- combined_test_boot %>% 
	group_by(.data$Facet, .data$RunName, .data$Iteration) %>% 
	summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$RawScore))


# Plot
run_colours <- rep('grey40', length(unique(auc_test_boot$RunName)))
names(run_colours) <- unique(auc_test_boot$RunName)
run_colours[BEST_NAME] <- '#44AA99' # Teal from Paul Tol's 'muted' palette

auc_plot <- ggplot(auc_test_boot, aes(x = RunName, y = AUC, colour = RunName)) +
	geom_blank() + # Ensures next line does not make x-scale continuous
	geom_vline(xintercept = c(1.56, 3.47, 7.5), colour = 'grey95', size = 1) +
	
	geom_hline(yintercept = 0.5, linetype = 2, colour = LINE_COLOUR) +
	geom_violin(colour = NA, fill = 'grey90') +
	geom_boxplot(width = 0.4, fill = NA) +
	
	# Annotate groups of features:
	annotate(geom = 'text', y = 1, size = 2.5, colour = 'grey30',
					 x = c(1, 2.5, 5.5, 8.57), 
					 label = c('Heuristic', 'Related viruses', 'Genome composition features', 'Combinations')) +
	coord_cartesian(ylim = c(0.38, 0.95), clip = 'off') +
	
	scale_colour_manual(values = run_colours, guide = FALSE) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
	xlab('Feature set') +
	PLOT_THEME +
	theme(plot.margin = margin(t = 14, r = 5.5, b = 5.5, l = 5.5),
				panel.border = element_rect(fill = NA, size = 0.3)) # Turning off clipping makes this appear thicker than in other plots


# Add p-value for specific comparisons:
for (lvl in names(P_VAL_COMPARISONS)) {
	p_vals <- c()
	
	for (comp in P_VAL_COMPARISONS[[lvl]]) {
		dta <- auc_test_boot %>% 
			filter(.data$RunName %in% comp)
		
		test_res <- kruskal.test(x = dta$AUC, g = dta$RunName)
		p_val <- test_res$p.value
		p_val <- if_else(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
		
		auc_plot <- auc_plot +
			geom_signif(comparisons = list(comp), 
									y_position = P_VAL_HEIGHTS[[lvl]],
									annotations = p_val,
									tip_length = 0.01, textsize = 2.5)
	}
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Accuracy by class --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
latest_names <- AllData %>% 
	distinct(.data$UniversalName, .data$LatestSppName)

prob_cutoff <- sum(AllData$InfectsHumans) / nrow(AllData)
message('Using cutoff ', prob_cutoff, ' for bagged accuracy plot')

get_confusion_df <- function(data, grouping, cutoff = prob_cutoff) {
	cm <- table(Observed = data$InfectsHumans == 'True', 
							Predicted = data$BagScore > cutoff,
							Class = data$Class)
	
	data.frame(cm)
}

bagged_confusion <- predictions_final_bagged %>% 
	left_join(latest_names, by = c('UniversalName')) %>% 
	left_join(zoo_status, by = 'LatestSppName') %>% 
	
	group_by(.data$RunName) %>% 
	group_modify(get_confusion_df) %>% 
	
	# Summarise to actual input classes, but label by human/zoonotic/other counts:
	mutate(Show = .data$Class != 'No human infections') %>% 
	group_by(.data$RunName, .data$Observed) %>% 
	mutate(Class_size = sum(.data$Freq)) %>% 
	group_by(.data$RunName, .data$Observed, .data$Predicted) %>% 
	summarise(Cell_count = sum(.data$Freq),
						Cell_prop = .data$Cell_count / unique(.data$Class_size),
						Main_label = paste(.data$Cell_count, unique(.data$Class_size), sep = "/"),
						Sub_label = if_else(unique(.data$Observed == 'FALSE'), NA_character_,
																paste0(.data$Class[.data$Show], ': ', .data$Freq[.data$Show], collapse = '\n'))) %>% 
	ungroup() %>% 
	mutate(Observed = factor(.data$Observed, levels = c('TRUE', 'FALSE'),
													 labels = c('Yes', 'No')),
				 Predicted = factor(.data$Predicted, levels = c('TRUE', 'FALSE'),
				 									  labels = c('Yes', 'No')),
				 RunName = str_replace(.data$RunName, '\n', ' '))



cm_plot <- ggplot(bagged_confusion, aes(x = Observed, y = Predicted, fill = Cell_prop, colour = Observed == Predicted)) +
	geom_tile() +
	geom_text(aes(label = Main_label), nudge_y = 0.05, size = 3) +
	geom_text(aes(label = Sub_label), nudge_y = -0.20, size = 2, lineheight = 0.8) +
	scale_fill_distiller(palette = 'Greens', direction = 1) +
	scale_colour_manual(values = c('TRUE' = 'grey95', 'FALSE' = 'grey30'), guide = FALSE) +
	labs(x = 'Prediction', y = 'Infects humans') +
	labs(x = 'Known to infect humans', y = 'Predicted to\ninfect humans', fill = 'Proportion') +
	PLOT_THEME +
	theme(legend.position = 'off',
				panel.spacing = unit(1.5, 'pt'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Priority categories ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Classify by priority
#  - Very high: entire CI above cutoff
#  - High: Median and upper limit above cutoff, but CI crosses cutoff
#  - Medium: Only upper limit is above cutoff, i.e. still potentially zoonotic, but less likely
#  - Low: entire CI below cutoff, i.e. very unlikely to be zoonotic
prioritize <- function(lower_bound, upper_bound, median, cutoff = prob_cutoff) {
	stopifnot(length(cutoff) == 1)
	stopifnot(cutoff > 0 & cutoff < 1)
	
	p <- if_else(lower_bound > cutoff, 'Very high',
							 if_else(median > cutoff, 'High',
							 				if_else(upper_bound > cutoff, 'Medium',
							 								'Low')))
	factor(p, levels = c('Low', 'Medium', 'High', 'Very high'))
}


known_priorities <- predictions_final_bagged %>% 
	left_join(latest_names, by = 'UniversalName') %>% 
	left_join(zoo_status, by = 'LatestSppName') %>% 
	
	mutate(Priority = prioritize(lower_bound = .data$BagScore_Lower,
															 upper_bound = .data$BagScore_Upper,
															 median = .data$BagScore),
				 Class = factor(.data$Class, levels = rev(levels(.data$Class)),
				 							  labels = c('Human\nvirus', 'Zoonotic',
				 							 					   'No known\nhuman\ninfections')))


# Plot
priority_plot <- ggplot(known_priorities, aes(x = Class, fill = Priority)) +
	geom_bar(position = 'fill') +
	scale_fill_manual(values = PRIORITY_COLOURS) +
	labs(x = NULL, y = 'Proportion') +
	PLOT_THEME +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- ROC curves for best classifier -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
predictions_final_boot <- predictions_final_boot %>% 
	filter(.data$Dataset == 'test') %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True')

predictions_final_bagged <- predictions_final_bagged %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True') %>% 
	rename(RawScore = .data$BagScore)


# Plot
roc_plot <- ggplot(predictions_final_bagged, aes(d = InfectsHumans, m = RawScore, colour = RunName)) +
	# Bootstraps:
	geom_roc(aes(fill = factor(Iteration)), n.cuts = 0, colour = 'grey90', size = 0.2,
					 data = predictions_final_boot) +
	
	# Bagged:
	geom_abline(colour = LINE_COLOUR, linetype = 2) +
	geom_roc(cutoffs.at = c(0.15, seq(0.2, 0.7, by = 0.1)), labelround = 2, labelsize = 2.5,
					 hjust = -0.35, vjust = 1.25, fontface = 'bold') +
	
	# Settings:
	scale_colour_manual(values = run_colours, guide = FALSE) +
	guides(fill = FALSE) +
	coord_equal() +
	scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
	labs(x = 'Proportion false positives', 
			 y = 'Proportion true positives',
			 colour = 'Model') +
	PLOT_THEME +
	theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 10.5))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine plots and save ---------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

bottom_row <- plot_grid(roc_plot, cm_plot, priority_plot,
												ncol = 3, rel_widths = c(1.2, 1.2, 1.7),
												labels = c('B', 'C', 'D'), hjust = -0.25)

combined_plot <- plot_grid(auc_plot, bottom_row,
													 nrow = 2, rel_heights = c(2, 1),
													 labels = c('A', ''), hjust = -0.25)


## Output:
out_dir <- file.path('Plots', 'Intermediates')

if (!dir.exists(out_dir))
	dir.create(out_dir, recursive = TRUE)


ggsave2(file.path('Plots', 'Figure1.pdf'), combined_plot, width = 7, height = 5)


## Save data transformations/calculations for use in other plots:
saveRDS(zoo_status, file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))
saveRDS(merged_taxonomy, file.path('Plots', 'Intermediates', 'figure1_merged_taxonomy.rds'))
saveRDS(known_priorities, file.path('Plots', 'Intermediates', 'figure1_trainingset_priorities.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values reported in text --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

cat('\n\nAUC: bootstrap stats\n')
boot_stats_auc <- auc_test_boot %>% 
	group_by(.data$RunName) %>% 
	summarise(Median_AUC = median(.data$AUC),
						SD_AUC = sd(.data$AUC),
						N_iterations = n(),
						Min_AUC = min(.data$AUC),
						Max_AUC = max(.data$AUC),
						Lower95_AUC = quantile(.data$AUC, probs = 0.025),
						Upper95_AUC = quantile(.data$AUC, probs = 0.975))

print(boot_stats_auc)


cat('\nAUC: bagging\n')
predictions_final_bagged  %>% 
	group_by(.data$RunName) %>% 
	summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$RawScore)) %>% 
	print()


cat('\nBinary classification: bagged models\n')
predictions_final_bagged %>% 
	left_join(latest_names, by = c('UniversalName')) %>% 
	left_join(zoo_status, by = 'LatestSppName') %>% 
	
	mutate(Prediction = .data$RawScore > prob_cutoff) %>% 
	group_by(.data$RunName, .data$Class) %>% 
	summarise(Accuracy = sum(.data$Prediction == .data$InfectsHumans) / n())


cat('\nPriorities for training viruses:\n')
tbl <- table(InfectsHumans = known_priorities$InfectsHumans, Priority = known_priorities$Priority)

cat('\n--Counts\n')
print(tbl)

cat('\n--Proportions\n')
print(tbl/rowSums(tbl))
