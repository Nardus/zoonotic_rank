## =================================================================================================
## General plot settings, used to ensure consistency between plots
## =================================================================================================

# Not available on conda, but exact version does not matter anyway
suppressPackageStartupMessages({
	if (!require("khroma")) {
		install.packages("khroma", repos = "https://cloud.r-project.org", verbose = FALSE)
		library(khroma)
	}
})

library(ggplot2)
library(colorspace)

# Overall theme
PLOT_THEME <- theme_bw() +
	theme(axis.title = element_text(size = 8),
				axis.text = element_text(size = 7),
				strip.text = element_text(size = 6, margin = margin(2.5, 2.5, 2.5, 2.5)),
				legend.title = element_text(size = 6),
				legend.text = element_text(size = 6),
				legend.key.height = unit(0.6, 'lines'),
				legend.key.width = unit(0.6, 'lines'),
				legend.margin = margin(t = 2.5, r = 5.5, b = 2.5, l = 0),
				panel.grid = element_blank(),
				plot.background = element_blank())


# Factor order:
MEASURE_TYPE_ORDER <- c('Amino acid bias', 'Codon bias', 'Dinucleotide bias', 
												'Bridge dinucleotide bias', 'Non-bridge dinucleotide bias', 
												'Mixed')

MEASURE_TYPE_LABELS <- c('Amino\nacid bias', 'Codon\nbias', 'Dinucleotide\nbias', 
												 'Bridge\ndinuc. bias', 'Non-bridge\ndinuc. bias', 
												 'Mixed')

FEATURE_SET_ORDER <- c('Taxonomy', 'Phylogenetic neighbourhood', 
											 'Viral genomic features', 'Similarity to ISGs', 
											 'Similarity to housekeeping genes', 
											 'Similarity to remaining genes',
											 'Mixed')

FEATURE_SET_LABELS <- c('Taxonomy', 'Phylogenetic neighbourhood', 
												'Viral genomic\nfeatures', 'Similarity to\nISGs', 
												'Similarity to\nhousekeeping\ngenes', 
												'Similarity to\nremaining\ngenes',
												'Mixed')


# Sizes
LINE_SIZE <- 0.7      # For data lines (i.e. not outlines)
BAR_LINE_SIZE <- 0.3  # Line width around bars


# General colours
LINE_COLOUR <- 'grey30'  # Outline colour for bars, boxplots, etc


# Colours for specific things:
# - These should repeat in different plots
status_palette <- colour('bright')
discrete_palette <- colour('muted') # Supports more colours than 'bright'

ZOONOTIC_STATUS_COLOURS <- c(status_palette(2), '#AA3377', '#999999')
names(ZOONOTIC_STATUS_COLOURS) <- c('No human infections', 'Zoonotic', 'Human virus', 'Novel')

ZOONOTIC_STATUS_COLOURS_LIGHT <- lighten(ZOONOTIC_STATUS_COLOURS, amount = 0.3)

TF_STATUS_COLOURS <- c(status_palette(2), '#999999')
names(TF_STATUS_COLOURS) <- c('FALSE', 'TRUE', 'NA')


# - From colorbrewer's Set2 (chosen manually to avoid clash with zoonotic colours):
DISCRETE_COLOURS <- c('#66c2a5', '#fc8d62', '#8da0cb', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3')


# - From colorbrewer's Dark2
FEATURE_SET_COLOURS <- c('grey20', 'grey20', 'grey45', '#d95f02', '#1b9e77', '#e6ab02', '#666666')
names(FEATURE_SET_COLOURS) <- FEATURE_SET_ORDER

FEATURE_SET_COLOURS_LB <- FEATURE_SET_COLOURS   # Names in this version contains line breaks
names(FEATURE_SET_COLOURS_LB) <- FEATURE_SET_LABELS


# Priority categories
PRIORITY_COLOURS <- c(Low = '#99d594', Medium = '#fdd568', High = '#fdae61', `Very high` = '#d53e4f')


# Labels for commonly plotted quantities (to ensure consistency):
SCORE_LABEL <- "Predicted probability"

SCORE_LABEL_2LINE <- "Predicted\nprobability"

PRIORITY_LABEL <- "Zoonotic potential"
PRIORITY_LABEL_2LINE <- "Zoonotic\npotential"