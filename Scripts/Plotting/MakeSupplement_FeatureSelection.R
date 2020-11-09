##
## Supplementary figure: number of top features to include
## 

library(dplyr)
library(tidyr)
library(readr)
library(ModelMetrics)
library(plotly)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'rundata_utils.R'))


RUN_NAMES <- c(FeatureSelection_Top10 = '10',
							 FeatureSelection_Top50 = '50',
							 FeatureSelection_Top100 = '100',
							 FeatureSelection_Top125 = '125',
							 FeatureSelection_Top150 = '150',
							 FeatureSelection_Top175 = '175',
							 FeatureSelection_Top200 = '200')
RUN_IDS <- names(RUN_NAMES)


## Read run data
predictions_test_boot <- read_multiple_runs(RUN_IDS, '_Predictions.rds') %>% 
	filter(.data$Dataset == 'test') %>% 
	mutate(RunName = RUN_NAMES[.data$RunID],
				 RunName = factor(.data$RunName, levels = RUN_NAMES))


## Calculate AUC
auc_test_boot <- predictions_test_boot %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True') %>% 
	group_by(RunName, Iteration) %>% 
	summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$RawScore))


## Plot
p_all <- ggplot(auc_test_boot, aes(x = RunName, y = AUC)) +
	geom_hline(yintercept = 0.5, linetype = 2, colour = LINE_COLOUR) +
	geom_violin(colour = NA, fill = 'grey90') +
	geom_boxplot(width = 0.4, fill = NA) +
	labs(x = 'Number of features retained', y = 'Ranking performance (AUC)') +
	scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
	PLOT_THEME +
	theme(strip.text = element_blank())


ggsave(file.path('Plots', 'Supplement_FeatureSelection.pdf'), p_all, width = 7, height = 5)


# Mean values:
auc_test_boot %>% 
	group_by(RunName) %>% 
	summarise(MeanAUC = mean(AUC),
						MedianAUC = median(AUC)) %>% 
	print()
