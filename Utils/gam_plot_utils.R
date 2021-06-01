##
## Calculation function for gam plotting GAMMs
##

library(parallel)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## ---- To plot model rankings: ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Extract effect names
extract_effect_names <- function(model) {
	modelTerms <- names(model$coefficients) %>% 
		gsub('\\.[[:digit:]]+', '', .) %>% 
		unique()
	
	# ID fixed vs smooth effects
	fixed <- ! grepl('s\\([[:alnum:][:punct:]]+\\)', modelTerms)
	
	# Clean names
	cleanNames <- modelTerms %>% 
		gsub('[s]?\\(([[:alnum:][:punct:]]+)\\)', '\\1', .) %>% 
		gsub(' ', '', .)
	
	# Return
	data.frame(term = modelTerms,
						 termClean = cleanNames,
						 fixedEffect = fixed)
}

get_model_detail <- function(model) {
	modelSummary <- summary(model)
	
	terms <- extract_effect_names(model) %>% 
		filter(termClean != 'Intercept')
	
	dev <- get_relative_deviance(model, modelSummary$formula) %>% 
		select(-Model) %>% 
		mutate(PercentDevianceExplained = PropDevianceExplained * 100)
	
	# Return:
	terms %>% 
		left_join(dev, by = c('termClean' = 'RemovedTerm'))
}


# Summarise the top n models:
summarise_ranked_models <- function(rankedModels, cores = 1) {
	sums <- mclapply(rankedModels$ModelFit, get_model_detail, mc.cores = cores)
	names(sums) <- rankedModels$Formula
	sums <- bind_rows(sums, .id = 'Formula')
	
	rankedModels %>% 
		select(-ModelFit) %>% 
		mutate(Rank = rank(AIC, ties.method = 'random')) %>% 
		left_join(sums, by = 'Formula')
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## ---- To plot partial effects: ---------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Calculate partial residuals:
get_partial_resids <- function(gamFit, terms, seWithMean) {
	predType <- ifelse(seWithMean, 'iterms', 'terms')  # Doesn't have much meaning here, but included for consistency with get_partial_preds
	
	linearTerm <- predict(gamFit, type = predType, terms = terms) %>% 
		rowSums()
	
	partialResids <- residuals(gamFit) + linearTerm     # TODO: unclear if using deviance residuals is the best idea (differs from plot.gam)
	
	partialResids
}

# Make predictions for partial effects:
get_partial_preds <- function(gamFit, newdata, terms, seWithMean) {
	predType <- ifelse(seWithMean, 'iterms', 'terms')
	
	predict(gamFit, newdata = newdata, se.fit = T,
					type = predType, terms = terms) %>% 
		lapply(rowSums) %>% 
		as.data.frame() %>% 
		rename(y = fit, se = se.fit) %>% 
		mutate(ylower = y - 1.96*se,
					 yupper = y + 1.96*se)
}



## Calculate data for a partial effects plot containing interacting variables
#   If var2 is not specified, partial effects for a single variable are calculated
#  'gamFit': a fitted gam model
#  'var1' / 'var2': the variables that form the interaction
#  'seWithMean': models the same option in plot.gam: "If TRUE the component smooths 
#                are shown with confidence intervals that include the uncertainty about 
#                the overall mean [...]"
#  fixedEffect: whether the terms are part of a smooth or are fixed effects
#
get_partial_effects_interaction <- function(gamFit, var1, var2, seWithMean = TRUE, fixedEffect = FALSE) {
  ## Term names: wrap in s():
	if (!is.null(var2)) {
		if (fixedEffect) stop('Non-smooth interactions not implemented')
		
		termnames <- c(paste0('s(', var1, ',', var2, ')'))
		
	} else {
		if (!fixedEffect) {
			termnames <- paste0('s(', var1, ')')
		} else {
			termnames <- var1
		}
	}
		
	
	## Variables not part of the effect / interaction are kept constant:
	modelData <- gamFit$model
	responseIndex <- attr(modelData, 'terms') %>% attr('response')
	responseName <- colnames(modelData)[responseIndex]
	
	otherData <- modelData %>% 
		select(-one_of(responseName, var1, var2))
	
	numericData <- otherData %>% 
		summarise_if(is.numeric, ~ median(.))
	
	factorData <- otherData %>% 
		summarise_if(is.factor, ~ names(which.max(table(.))))		
	
	stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
	
	
	## Calculate partial residuals
	partialDat <- modelData %>% 
		data.frame() %>% 
		select(one_of(var1, var2))
	
	if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
	if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
	
	
	partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
	partialResids <- cbind(partialDat,
												 Residual = partialResids)
	
	
	## Predictions
	# - Make a prediction for each level of the interaction var (so for all interactions that occur)
	newData <- modelData %>% 
		data.frame() %>% 
		select(one_of(var1, var2)) %>% 
		unique()
	
	# - All other data get set to their median (or the most common factor level)
	if (length(numericData) > 0) newData <- cbind(newData, numericData)
	if (length(factorData) > 0) newData <- cbind(newData, factorData)
	
	
	# - Make predictions
	newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
		mutate(IsSignificant = if_else(ylower <= 0 & yupper >= 0, 'No', 'Yes')) %>%    # Check if CI crosses zero
		cbind(newData)
	
	# Add significance to residuals (for plotting):
	partialResids <- newPredictions %>% 
		select(one_of(var1, var2), IsSignificant) %>% 
		right_join(partialResids)
	
	# Return:
	list(effects = newPredictions, partialResiduals = partialResids)
}



## Calculate data for a partial effects plot of a single effect
#   - Note: this simply predicts over the observed values, so may not yield a very
#     smooth plot
get_partial_effects <- function(fit, var, seWithMean = TRUE) {
	get_partial_effects_interaction(fit, var, NULL, seWithMean)
}



## Calculate data for a partial effects plot of a single binary variable
#   - Assumed to be a fixed effect by default
#   - 'removeNegatives': should values for the negative case (corresponding to the intercept/baseline)
#     be returned?
get_partial_effects_binary_single <- function(fit, var, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
	plotData <- get_partial_effects_interaction(fit, var1 = var, NULL, seWithMean, fixedEffect)
	
	# Remove negatives
	if (removeNegatives) {
		plotData$effects <- plotData$effects[plotData$effects[[var]] == 1, ]
		plotData$partialResiduals <- plotData$partialResiduals[plotData$partialResiduals[[var]] == 1, ]
	}
	
	# Add a column containing var as a label
	plotData$effects$variable <- var
	plotData$partialResiduals$variable <- var
	
	# Return
	plotData
}


## Calculate data for partial effects plots of several independent binary variables
# 'vars': a character vector of parameter names
get_partial_effects_binary <- function(fit, vars, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
	allData <- lapply(vars, get_partial_effects_binary_single, fit = fit, 
										seWithMean = seWithMean, 
										fixedEffect = fixedEffect, 
										removeNegatives = removeNegatives)
	
	extract_by_name <- function(x, name) x[[name]]
	effects <- lapply(allData, extract_by_name, 'effects')
	partialResiduals <- lapply(allData, extract_by_name, 'partialResiduals')
	
	effects <- do.call(rbind, effects)
	partialResiduals <- do.call(rbind, partialResiduals)
	
	list(effects = effects, partialResiduals = partialResiduals)
}




## Calculate partial effects for a continious smooth variable
#  Similar to above, but use an even range of values rather than the observed
#  ones to get a smoother plot
get_partial_effects_continuous <- function(gamFit, var, resolution = 1, seWithMean = TRUE) {
	## Term names: wrap in s():
	termnames <- paste0('s(', var, ')')
	
	
	## Data not part of effect kept constant:
	modelData <- gamFit$model
	responseIndex <- attr(modelData, 'terms') %>% attr('response')
	responseName <- colnames(modelData)[responseIndex]
	
	otherData <- modelData %>% 
		select(-one_of(responseName, var))
	
	numericData <- otherData %>% 
		summarise_if(is.numeric, ~ median(.))
	
	factorData <- otherData %>% 
		summarise_if(is.factor, ~ names(which.max(table(.))))		
	
	stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
	
	
	## Calculate partial residuals
	partialDat <- modelData %>% 
		data.frame() %>% 
		select(one_of(var))
	
	if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
	if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
	
	
	partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
	partialResids <- cbind(partialDat,
												 Residual = partialResids)
	
	
	## Predictions
	# - Predictions over a smooth range of values spanning the range of var:
	newData <- seq(min(modelData[, var]), max(modelData[, var]), by = resolution) %>% 
		data.frame()
	
	colnames(newData) <- var
	
	# - All other data get set to their median (or the most common factor level)
	if (length(numericData) > 0) newData <- cbind(newData, numericData)
	if (length(factorData) > 0) newData <- cbind(newData, factorData)
	
	# - Make predictions
	newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
		mutate(NotSignificant = ylower <= 0 & yupper >= 0,
					 IsSignificant = if_else(all(NotSignificant), 'No', 'Yes')) %>%    # Check if CI crosses 0 over entire range
		cbind(newData)
	
	partialResids$IsSignificant <- unique(newPredictions$IsSignificant)
	
	# Return:
	list(effects = newPredictions, partialResiduals = partialResids)
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## ---- Plot partial effects on a phylogeny: ------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

plot_phylo_effect <- function(tree, plotData, showtiplabs = T) {
	# Put family column first to allow joining to tree:
	effects <- plotData$effects %>% 
		select(Family, everything()) %>% 
		rename(x = y)
	
	resids <- plotData$partialResiduals %>% 
		select(Family, everything())
	
	# Plot tree
	p <- ggtree(tree) +
		theme_tree2()
	
	if (showtiplabs) {
		p <- p + 
			geom_tiplab(size = 2.5, colour = 'grey20') +
			xlim_tree(c(0, 1.2))
		
		segmentSize <- 2.8
		jitterHeight <- 0.05
		pointSize <- 1
		
	} else {
		# Assume fig will be small:
		segmentSize <- 1
		jitterHeight <- 0.005
		pointSize <- 0.1
	}
	
	
	p <- facet_plot(p, panel = 'Values', geom = geom_segment, 
									x = 0, xend = 0, y = -Inf, yend = Inf,
									linetype = 2, size = 0.3, colour = 'grey40',
									data = resids)
	
	p <- facet_plot(p, panel = 'Values', 
									aes(x = ylower, xend = yupper, y = y, yend = y, colour = IsSignificant),
									geom = geom_segment,
									size = segmentSize, alpha = 0.5,
									data = effects)
	
	p <- facet_plot(p, panel = 'Values',
									aes(x = x-0.1, xend = x+0.1, y = y, yend = y, colour = IsSignificant),
									geom = geom_segment,
									size = segmentSize,
									data = effects)
	
	p <- facet_plot(p, panel = 'Values', 
									aes(x = Residual, colour = IsSignificant),
									geom = geom_point,
									position = position_jitter(height = jitterHeight),
									alpha = 0.6,
									size = pointSize,
									data = resids)
	
	p +
		scale_colour_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
		scale_fill_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
		labs(x = 'log(Odds zoonotic)', y = 'Virus family') +
		PLOT_THEME +
		theme(strip.text = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank())
}
