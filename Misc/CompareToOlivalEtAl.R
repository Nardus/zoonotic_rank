##
## Zoonosis predictor - calculate prediction accuracy of Olival et al's models
## 

library(stringi)
library(purrr)
library(tibble)
library(tidyr)
library(dplyr)
library(mgcv)
library(ModelMetrics)


# From Olival et al.'s model_reduction.R (no changes made)
rm_terms <- function(mod, terms) {
	fr = as.character(formula(mod))
	lhs = fr[2]
	rhs = fr[3]
	vars_regex = paste0("(", paste(terms, collapse = "|"), ")")
	new_rhs = stri_replace_all_regex(rhs, paste0("\\s*s\\(", vars_regex, "\\,[^\\)]+\\)\\s*\\+?"), "")
	new_rhs = stri_replace_all_fixed(new_rhs, "+, k = 7) ", "")
	new_rhs = stri_replace_all_fixed(new_rhs, "+ +s", "+ s")
	new_formula = paste(lhs, "~", new_rhs)
	new_formula = stri_replace_all_regex(new_formula, "[\\s\\n]+", " ")
	new_formula = stri_replace_all_regex(new_formula, "[+\\s]*$", "")
	return(new_formula)
}

# Derived from Olival et al.'s cv_gam() function in their cross_validation.R script:
#  - Modified to do an 85:15 train:test split to match our model evaluations
#    (we used 70% for training + 15% for calibration, so entire training procedure used 85% 
#     of the data - that's more than each fold model would have in cross-validation)
test_gam <- function(mod, train_prop = 0.85) {
	dat = mod$model
	dat$fold = base::sample(c("train", "test"), size = nrow(dat), replace = TRUE,
													prob = c(train_prop, 1 - train_prop))
	
	terms <- attributes(mod$terms)
	
	training_dat = dat[dat$fold == "train", ]
	
	zero_terms = names(training_dat)[map_lgl(training_dat, ~ length(unique(.x)) == 1)]
	new_formula = rm_terms(mod, zero_terms)
	
	new_mod = gam(
		formula = as.formula(new_formula),
		data = training_dat,
		family = mod$family
	)
	
	preds = as.vector(predict(new_mod, newdata = dat[dat$fold == "test",], type = "response"))
	
	actuals = dat[dat$fold == "test", names(dat)[attributes(terms(mod))$response]]
	
	diffs = actuals - preds
	error_term = mean(abs(diffs))
	
	auc_val <- ModelMetrics::auc(actual = actuals, predicted = preds)
	
	return(
		data_frame(
			n_fit = nrow(dat) - length(preds),
			n_validate = length(preds),
			mean_error = error_term,
			auc = auc_val
		)
	)
}


## Apply
# Get the model fit / data:
download.file("https://zenodo.org/record/807517/files/ecohealthalliance/HP3-v1.0.9.zip", "ExternalData/HP3-v1.0.9.zip")
unzip("ExternalData/HP3-v1.0.9.zip", "ecohealthalliance-HP3-928327a/intermediates/vtraits_strict_models.rds", exdir = "ExternalData")
unzip("ExternalData/HP3-v1.0.9.zip", "ecohealthalliance-HP3-928327a/intermediates/postprocessed_database.rds", exdir = "ExternalData")

set.seed(290321)
virus_traits_strict <- readRDS("ExternalData/ecohealthalliance-HP3-928327a/intermediates/vtraits_strict_models.rds")

db <- readRDS("ExternalData/ecohealthalliance-HP3-928327a/intermediates/postprocessed_database.rds")
hp3_viruses <- db$viruses


## Remove duplicated viruses
#  - Sample one example randomly
name_matches <- read.csv("InternalData/NameMatches_All.csv", stringsAsFactors = FALSE) %>% 
	filter(!.data$Olival == "") %>% 
	mutate(Olival = stringr::str_replace_all(.data$Olival, " ", "_"),
				 LatestSppName = if_else(is.na(.data$SppName_ICTV_MSL2018b), .data$UniversalName, .data$SppName_ICTV_MSL2018b)) %>% 
	distinct(.data$LatestSppName, .data$Olival) %>% 
	group_by(.data$LatestSppName) %>% 
	sample_n(size = 1) %>% 
	ungroup()


hp3_viruses <- hp3_viruses %>% 
	select(.data$IsZoonotic.stringent,
				 .data$vVirusNameCorrected, .data$cb_dist_noHoSa_max.stringent,
				 .data$vCytoReplicTF, .data$Vector, .data$vWOKcitesLn) %>% 
	na.omit()

stopifnot(nrow(hp3_viruses) == nrow(virus_traits_strict$model[[1]]$model))

hp3_viruses <- hp3_viruses %>% 
	filter(.data$vVirusNameCorrected %in% name_matches$Olival)

# Code from Olival et al's "04-fit-models.R":
logp = function(x){   # Fn to take log but make zeros less 10x less than min
	# x[is.na(x)] <- 0
	m = min(x[ x > 0], na.rm = T)
	x = log( x + m )
	return(x)
}

PD_centers = names(hp3_viruses)[stri_detect_regex(names(hp3_viruses), "^(cb|st)_.*stringent$")]

data_set = hp3_viruses %>%
	mutate_each_(funs("logp"), vars = PD_centers)

names(data_set)[names(data_set) %in% PD_centers] <- paste0(PD_centers, "Ln")



## Find best model
best_model <- virus_traits_strict$model[[which(virus_traits_strict$daic == 0)]]
best_model <- update(best_model, 
										 data = data_set,
										 family = "binomial")

# Best model without publication counts:
no_pub_count_model <- update(best_model, . ~ . -s(vWOKcitesLn, bs = "tp", k = 7) - s(vPubMedCitesLn, bs = "tp", k = 7),
														 data = data_set,
														 family = "binomial")

# What if host breadth wasn't avialable?
no_host_breadth_model <- update(no_pub_count_model, . ~ . -s(cb_dist_noHoSa_max.stringentLn, bs = "tp", k = 7),
																data = data_set,
																family = "binomial")


## Repeated training and evaluation
cv_results_best <- replicate(n = 1000, test_gam(best_model), simplify = FALSE)
cv_results_no_pub_counts <- replicate(n = 1000, test_gam(no_pub_count_model), simplify = FALSE)
cv_results_no_host_breadth <- replicate(n = 1000, test_gam(no_host_breadth_model), simplify = FALSE)


## Summarise results
report <- function(cv_results) {
	cv_results <- bind_rows(cv_results)
	
	with(cv_results, {
		cat("\n\tMean accuracy:", mean(1 - mean_error), "| SD:", sd(mean_error), 
				"\n\tMedian AUC:", median(auc), "| SD", sd(auc),
				"\n")
	})
}

cat("\nBest model:")
report(cv_results_best)


cat("\nWithout publication counts:")
report(cv_results_no_pub_counts)

cat("\nWithout host breadth or publication counts:")
report(cv_results_no_host_breadth)
