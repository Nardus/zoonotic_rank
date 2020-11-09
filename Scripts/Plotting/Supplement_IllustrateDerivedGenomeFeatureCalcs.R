#
# Illustrate the calculation of derived genome features (for inclusion in methods section)
# 

library(dplyr)
library(ggplot2)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))


# Feature to use as illustration
GENE_SET <- "Housekeeping"
FEATURE_PLAIN <- "CTG.Bias_Coding"
FEATURE_DENS <- paste("GenomicDensity", GENE_SET, FEATURE_PLAIN, sep = '_')
FEATURE_LABEL <- "CTG bias"

# Example viruses
# - Chosen to show a range of values in the feature above
EXAMPLE_VIRUSES <- c('Woolly monkey sarcoma virus', 'Bovine nidovirus 1', 
										 'Zika virus', 'Louping ill virus', 
										 'Human alphaherpesvirus 2', 'Foot-and-mouth disease virus',
										 'Severe acute respiratory syndrome-related coronavirus',
										 'Sindbis virus', 'Cervid alphaherpesvirus 1', 
										 'Canid alphaherpesvirus 1')


## Load data
# Latest spp names
training_data <- readRDS(file.path('CalculatedData', 'SplitData_Training.rds'))

latest_names <- training_data %>% 
	distinct(.data$UniversalName, .data$Strain, .data$LatestSppName)

stopifnot(all(EXAMPLE_VIRUSES %in% training_data$LatestSppName))


# Zoonotic status
zoo_status <- readRDS(file.path('Plots', 'Intermediates', 'figure1_zoo_status.rds'))


# Calculated features
virus_features_raw <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Virus.rds'))
virus_features_dens <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Distances.rds'))
human_features <- readRDS(file.path('CalculatedData', 'GenomicFeatures-HumanCombined.rds'))



## Select example feature / viruses
human_vals <- human_features %>% 
	filter(.data$GeneSet == GENE_SET) %>% 
	pull(FEATURE_PLAIN) %>% 
	data.frame(raw = ., dens = NA_real_)

virus_vals <- full_join(virus_features_raw, virus_features_dens, by = c('UniversalName', 'Strain')) %>% 
	left_join(latest_names, by = c('UniversalName', 'Strain')) %>% 
	left_join(zoo_status, by = 'LatestSppName') %>% 
	filter(.data$LatestSppName %in% EXAMPLE_VIRUSES) %>% 
	select(one_of('LatestSppName', 'Class', FEATURE_PLAIN, FEATURE_DENS)) %>% 
	rename(raw = !!FEATURE_PLAIN,
				 dens = !!FEATURE_DENS)

## Labels for viruses (on both axes)
wrapped_names <- virus_vals$LatestSppName %>% 
	if_else(. == 'Severe acute respiratory syndrome-related coronavirus',
					'Severe acute respiratory\nsyndrome-related coronavirus',
					.)

virus_labs_raw <- virus_vals$raw
virus_labs_dens <- virus_vals$dens

names(virus_labs_raw) <- wrapped_names
names(virus_labs_dens) <- wrapped_names


## Plot
step4 <- ggplot(human_vals, aes(x = raw)) +
	geom_histogram(aes(y = ..density..), bins = 20, fill = 'grey90') +
	geom_density(colour = 'grey50', size = 1) +
	
	geom_segment(aes(xend = raw, y = dens, yend = Inf), linetype = 3, colour = LINE_COLOUR, 
							 data = virus_vals) +
	geom_segment(aes(xend = Inf, y = dens, yend = dens), linetype = 3, colour = LINE_COLOUR, 
							 data = virus_vals) +
	geom_point(aes(y = dens), colour = '#CC3311', size = 2, data = virus_vals) +
	
	labs(x = FEATURE_LABEL, y = 'Density') +
	
	scale_x_continuous(sec.axis = sec_axis(trans = ~ .,
																				 breaks = virus_labs_raw,
																				 labels = names(virus_labs_raw))) +
	scale_y_continuous(expand = expand_scale(mult = c(0, 0.06)),
										 sec.axis = sec_axis(trans = ~ .,
										 										breaks = virus_labs_dens,
										 										labels = names(virus_labs_dens))) +
	PLOT_THEME +
	theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, 
																			 size = 6, lineheight = 0.7,
																			 face = 'italic'),
				axis.text.y.right = element_text(size = 6, lineheight = 0.7, face = 'italic'))


ggsave2(file.path('Plots', 'Supplement_methods_derived_genome_features.pdf'), step4, width = 6, height = 6)
