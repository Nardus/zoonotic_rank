## =================================================================================================
## Plot model of host effects on predictions from the "all non-training viruses" case study
## =================================================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(mgcv)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))
source(file.path('Utils', 'gam_plot_utils.R'))

## Load predictions and hosts shown in figure 3:
host_data <- read_rds(file.path('Plots', 'Intermediates', 'figure3_host_data.rds'))

novel_taxonomy <- read_rds(file.path('Plots', 'Intermediates', 'figure3_novel_taxonomy.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Model: do hosts differ in ranks assigned? ---------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data for gam model:
model_data <- host_data %>% 
	left_join(novel_taxonomy, by = c("Name" = "Species")) %>% 
	rename(host_order = .data$order,
				 host_class = .data$class) %>% 
	filter(!is.na(.data$host_corrected))


# Clean labels / specify plotting order
stopifnot(all(model_data$host_class %in% c("Aves", "Mammalia", "Insecta", "Arachnida")))

model_data <- model_data %>% 
	mutate(host_phylum = if_else(.data$host_class %in% c("Aves", "Mammalia"), "Chordata", "Artropoda"),
				 facet = case_when(is.na(.data$host_order) ~ "",
				 									.data$host_order %in% c("Diptera", "Ixodida") ~ "",  # No space in plot for "Vector" label
				 									.data$host_class == "Aves" ~ "Bird",
				 									.data$host_class == "Mammalia" ~ "Mammal",
				 									TRUE ~ .data$host_class),
				 facet = factor(.data$facet, levels = c("Mammal", "Bird", "Vector", "")),
				 
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

model_data <- model_data %>% 
	mutate(host_order = factor(.data$host_order),
				 host_class = factor(.data$host_class),
				 is_human = if_else(.data$host_label == "Human", 1, 0),
				 is_arthropod = if_else(.data$host_phylum == "Artropoda", 1, 0))


set.seed(30032021)
fit_all <- gam(calibrated_score_mean ~ is_human + is_arthropod + s(host_order, bs = "re") + s(host_class, bs = "re") + s(family_label, bs = "re"), 
							 family = "betar",
							 method = "ML",
							 select = TRUE,
							 data = model_data)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot Human effect  -------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_binary <- function(plot_data, x_lab = "x") {
	ggplot(plot_data$effects, aes(x = variable, y = y)) +
		geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
		
		geom_jitter(aes(y = Residual), alpha = 0.4, size = 0.6, width = 0.35, height = 0,
								colour = 'grey80',
								data = plot_data$partialResiduals) +
		
		geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper,
										 colour = IsSignificant, fill = IsSignificant), 
								 stat = 'identity', alpha = 0.6, colour = NA) +
		geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y,
										 colour = IsSignificant, fill = IsSignificant), 
								 stat = 'identity', alpha = 0.6, size = 0.5) +
		
		scale_colour_manual(values = c(No = 'grey60', Yes = '#08519c'), guide = F) +
		scale_fill_manual(values = c(No = 'grey60', Yes = '#08519c'), guide = F) +
		
		scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
		
		labs(x = x_lab, y = 'Effect') +
		PLOT_THEME +
		theme(axis.text.x = element_blank(),
					axis.ticks.x = element_blank())
}

plotData <- get_partial_effects_binary(fit = fit_all, vars = "is_human")

p_human_effect <- plot_binary(plotData, "Human\nsampled")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot host phylum effect --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plotData <- get_partial_effects_binary(fit = fit_all, vars = "is_arthropod")

p_phylum_effect <- plot_binary(plotData, "Arthropod\nsampled")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot host class effect ---------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_re <- function(fit, effect_name, facet_data, x_lab = "x", rotate = FALSE) {
	plot_data <- get_partial_effects(fit = fit, var = effect_name)
	
	plot_data$effects <- plot_data$effects %>% 
		left_join(facet_data, by = effect_name) %>% 
		mutate_at(effect_name, ~ factor(., levels = levels(facet_data[[effect_name]])))
	
	plot_data$partialResiduals <- plot_data$partialResiduals %>% 
		left_join(facet_data, by = effect_name) %>% 
		mutate_at(effect_name, ~ factor(., levels = levels(facet_data[[effect_name]])))
		
	
	p <- ggplot(plot_data$effects, aes_string(x = effect_name, y = "y")) +
		geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
		
		geom_jitter(aes(y = Residual), alpha = 0.4, size = 0.6, width = 0.2, height = 0, 
								colour = 'grey80',
								data = plot_data$partialResiduals) +
		
		geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper,
										 colour = IsSignificant, fill = IsSignificant), 
								 stat = 'identity', alpha = 0.6, colour = NA) +
		geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y,
										 colour = IsSignificant, fill = IsSignificant), 
								 stat = 'identity', alpha = 0.6, size = 0.5) +
		
		scale_colour_manual(values = c(No = 'grey60', Yes = '#08519c'), guide = F) +
		scale_fill_manual(values = c(No = 'grey60', Yes = '#08519c'), guide = F) +
		labs(x = x_lab, y = 'Effect') +
		PLOT_THEME +
		theme(strip.background = element_blank(),
					panel.spacing = unit(1, 'pt'))
	
	if (rotate) {
		p <- p +
			facet_grid(rows = vars(facet), scales = "free", space = "free") +
			coord_flip()
	} else {
		p <- p +
			facet_grid(cols = vars(facet), scales = "free", space = "free")
	}
	
	p
}


facet_data <- model_data %>% 
	mutate(facet = .data$host_phylum) %>% 
	distinct(.data$host_class, .data$facet)

p_class_effect <- plot_re(fit_all, effect_name = "host_class", 
												  facet_data = facet_data, x_lab = "Sampled host or\nvector class") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot host order effect ---------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
facet_data <- model_data %>% 
	distinct(.data$host_order, .data$facet)

p_host_effect <- plot_re(fit_all, effect_name = "host_order", 
												 facet_data = facet_data, x_lab = "Sampled host or vector order") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
				axis.title.x = element_text(margin = margin(t = -2)))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Plot family effect  ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
facet_data <- model_data %>% 
	distinct(.data$family_label, facet = .data$GenomeType)

p_family_effect <- plot_re(fit_all, effect_name = "family_label", 
													 facet_data = facet_data, x_lab = "Family") +
	theme(strip.text.x = element_text(angle = 90, hjust = 0),
				axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
top_row <- plot_grid(p_phylum_effect, p_class_effect, p_host_effect,
										 nrow = 1, rel_widths = c(1, 1.35, 6),
										 labels = c('A', 'B', 'C'),
										 align = "h", axis = "tb")

bottom_row <- plot_grid(p_human_effect, p_family_effect,
												nrow = 1, rel_widths = c(1, 7),
												labels = c('D', 'E'),
												align = "h", axis = "tb")

final_plot <- plot_grid(top_row, bottom_row,
												nrow = 2,
												rel_heights = c(1, 1.1))

ggsave2(file.path('Plots', 'Figure4.pdf'), final_plot, width = 7.5, height = 5)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values mentioned in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Print model summary
cat("\nModel summary:")
print(summary(fit_all))


cat("\n\nViruses from humans have higher predictions: mean =",
		mean(model_data$calibrated_score_mean[model_data$is_human == 1]),
		"sd =",
		sd(model_data$calibrated_score_mean[model_data$is_human == 1]),
		"\n\n")

cat("Viruses from other hosts have: mean =",
		mean(model_data$calibrated_score_mean[model_data$is_human == 0]),
		"sd =",
		sd(model_data$calibrated_score_mean[model_data$is_human == 0]),
		"\n\n")
