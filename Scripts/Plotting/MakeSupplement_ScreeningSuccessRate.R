## Plot a comparison of accumulation curves (i.e. fig 1D, but for alternate models too)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))

set.seed(3105221)

## Data
best_testing <- readRDS(file.path('Plots', 'Intermediates', 'figure1_virus_testing.rds'))

tax_preds <- readRDS("RunData/Taxonomy_LongRun/Taxonomy_LongRun_Bagged_predictions.rds") %>% 
	mutate(Rank = rank(-.data$BagScore)) %>% 
	arrange(.data$Rank)

pn_preds <- readRDS("RunData/PN_LongRun/PN_LongRun_Bagged_predictions.rds") %>% 
	mutate(Rank = rank(-.data$BagScore)) %>% 
	arrange(.data$Rank)


## Calculate accumulation curves
simulate_testing <- function(rank_data) {
  virus_testing <- data.frame()
  total_viruses <- nrow(rank_data)
  pos_viruses <- sum(rank_data$InfectsHumans == "True")
  
  for (i in 1:nrow(rank_data)) {
  	current_row <- rank_data[i, ]
  	viruses_before <- rank_data[rank_data$Rank <= current_row$Rank, ]
  	
  	virus_testing <- rbind(virus_testing, data.frame(
  		Rank = current_row$Rank,
  		prop_screened = nrow(viruses_before)/total_viruses,
  		prop_found = sum(viruses_before$InfectsHumans == "True")/pos_viruses
  	))
  }
  
  virus_testing
}

tax_testing <- simulate_testing(tax_preds)

pn_testing <- simulate_testing(pn_preds)


combined_testing <- list("Taxonomic" = tax_testing, 
												 "Phylogenetic\nneighbourhood" = pn_testing, 
												 "All genome\ncomposition\nfeature sets" = best_testing)


## Simulate random screening
observed_prevalence <- sum(tax_preds$InfectsHumans == "True") / nrow(tax_preds)

screen_randomly <- function(original_labels = tax_preds$InfectsHumans) {
	n <- length(original_labels)
	
	labels <- sample(original_labels)
	
	screened <- vector(length = n)
	found <- vector(length = n)
	
	for (i in 1:n) {
		screened[i] <- i / n
		found[i] <- sum(labels[1:i] == "True") / sum(labels == "True")
	}
	
	data.frame(prop_screened = screened,
						 prop_found = found)
}

random_testing <- replicate(1000, screen_randomly(), simplify = FALSE) %>% 
	bind_rows(.id = "replicate")


## Find point at which 50% of human-infecting viruses are found:
get_point <- function(data, example_points) {
	positive_accumulation_fun <- ecdf(data$prop_found)
	
	data.frame(prop_found = example_points,
						 prop_screened = positive_accumulation_fun(example_points))
}

examples <- lapply(combined_testing, get_point, example_points = 0.5) %>% 
	bind_rows(.id = "model")

combined_testing <- combined_testing %>% 
	bind_rows(.id = "model")


## Reduction in effort:
best_model = "All genome\ncomposition\nfeature sets"

reduction_arrow_data <- examples %>% 
	filter(.data$model == best_model | 
				 	.data$prop_screened == min(.data$prop_screened[.data$model != best_model])) %>% 
	summarise(lower_val = min(.data$prop_screened),
						upper_val = max(.data$prop_screened))

reduction_label_data <- reduction_arrow_data %>% 
	summarise(reduction = .data$upper_val/.data$lower_val,
						label_pos = (.data$upper_val - .data$lower_val)/2 + .data$lower_val)


## Plot
p <- ggplot(combined_testing) +
	geom_step(aes(x = prop_screened, y = prop_found, group = replicate), 
						colour = 'grey90', size = 0.2,
						data = random_testing) +
	
	geom_hline(yintercept = 0.5, linetype = 3, colour = LINE_COLOUR) +
	geom_segment(aes(x = prop_screened, xend = prop_screened, y = -Inf, yend = prop_found), 
							 linetype = 3, colour = LINE_COLOUR, data = examples) +
	
	geom_text(aes(label = sprintf("%2.2f-fold\ndifference\nin effort", reduction), x = label_pos, y = 0.08),
						colour = LINE_COLOUR, size = 2.5, data = reduction_label_data) +
	geom_segment(aes(x = lower_val, xend = upper_val, y = 0.155, yend = 0.155), 
								 arrow = arrow(length = unit(0.3, "lines"), type = "closed", ends = "both"),
								 colour = LINE_COLOUR, data = reduction_arrow_data) +
	
	geom_step(aes(x = prop_screened, y = prop_found, colour = model), size = 1) +
	
	scale_x_continuous(labels = percent_format(), expand = expand_scale()) +
	scale_y_continuous(labels = percent_format(accuracy = 1), expand = expand_scale()) +
	scale_colour_brewer(palette = "Set2", direction = -1) +
	
	labs(x = "Proportion of viruses screened", y = "Proportion of known human-infecting viruses encountered", colour = "Feature set") +
	coord_equal() +
	PLOT_THEME +
	theme(legend.key.height = unit(2, "lines"))


ggsave2("Plots/Supplement_ScreeningSuccessRate.pdf", width = 6, height = 4)


## Mentioned in text:
cat("\nScreening required to find 50% of human-infecting viruses:\n")
print(examples)
