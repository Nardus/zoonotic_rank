library(dplyr)
library(ggplot2)
library(cowplot)
library(doParallel)
library(caret)
# library(pheatmap)
# library(scales)

source(file.path('Scripts', 'Plotting', 'PlottingConstants.R'))

gene_set_colours <- c('Housekeeping' = FEATURE_SET_COLOURS[['Housekeeping gene mimicry']],
											'ISG' = FEATURE_SET_COLOURS[['ISG mimicry']],
											'Remaining' = FEATURE_SET_COLOURS[['Remaining gene mimicry']])

## Data
human_gene_features <- readRDS('CalculatedData/GenomicFeatures-HumanCombined.rds')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Heatmap ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# NOTE: This figure is unused - across all features, we expect genes to cluster by their ancestral 
#       relationships, but this does not mean they are not distinguishable on the basis of a few 
#       features which may have changed because of their role/when they are expressed

# gene_sets <- human_gene_features %>% 
# 	select(.data$GeneID, .data$GeneSet)
# 
# stopifnot(length(unique(gene_sets$GeneID)) == nrow(gene_sets))  # Expect rows to be individual genes
# 
# 
# # Create a matrix of feature values for each gene
# feature_mat <- human_gene_features %>% 
# 	select(-GeneID, -meanCPM, -TranscriptID, -GeneSet) %>% 
# 	as.matrix() %>% 
# 	t()
# 
# colnames(feature_mat) <- human_gene_features$GeneID
# 
# # Remove genes with only NA values
# all_na <- apply(feature_mat, 2, function(col) all(is.na(col)))
# message('Removing ', sum(all_na), ' genes with no feature values')
# 
# feature_mat <- feature_mat[, !all_na]
# 
# 
# # Scale all features to [0-1], so no feature dominates clustering, and all rows are directly comparable
# # - Note: this also transposes the matrix, but that is needed anyway since we have more vertical 
# #   than horizontal space (so putting genes on y axis)
# scaled_feature_mat <- apply(feature_mat, MARGIN = 1, FUN = rescale)
# rownames(scaled_feature_mat) <- colnames(feature_mat)
# 
# # Plot
# row_annotations <- data.frame(Set = gene_sets$GeneSet)
# rownames(row_annotations) <- gene_sets$GeneID
# 
# pheatmap(scaled_feature_mat, annotation_row = row_annotations,
# 				 annotation_colors = list(Set = gene_set_colours),
# 				 show_rownames = FALSE, show_colnames = FALSE, 
# 				 annotation_names_col = FALSE,
# 				 filename = file.path('Plots', 'Supplement_HumanGeneSetSimilarity_Heatmap.pdf'))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Predictions --------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Train a classifier to see whether the sets of genes can be distinguished based on their 
# genome features
cl <- makeCluster(8)
registerDoParallel(cl)

train_data <- human_gene_features %>% 
	select(-GeneID, -meanCPM, -TranscriptID)


trainSettings <- trainControl(method = 'repeatedcv',
															number = 10, 
															repeats = 3,
															search = 'random',
															returnResamp = 'all',
															sampling = 'down',
															savePredictions = 'final',
															classProbs = TRUE,
															summaryFunction = mnLogLoss)

trainedModel <- train(GeneSet ~ .,
											data = train_data,
											method = 'xgbTree',  # Also try xgbDART and xgbTree
											trControl = trainSettings,
											metric = 'logLoss',
											tuneLength = 100,
											na.action = 'na.pass',
											nthread = 1)


# Plot accuracy across CV folds
final_preds <- trainedModel$pred %>% 
	group_by(.data$obs, .data$Resample) %>% 
	mutate(obs_class_size = n()) %>% 
	group_by(.data$obs, .data$pred, .data$Resample) %>% 
	summarise(Proportion = n()/unique(.data$obs_class_size))


preds_plot <- ggplot(final_preds, aes(x = obs, y = Proportion, colour = pred)) +
	geom_boxplot() +
	scale_colour_manual(values = gene_set_colours) +
	labs(x = 'Actual class', y = 'Proportion of actual class', colour = 'Predicted class') +
	PLOT_THEME



registerDoSEQ()
stopCluster(cl)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Save trained model
saveRDS(trainedModel, 
				file.path('Plots', 'Intermediates', 'Supplement_HumanGeneSetSimilarity_trainedmodel_xgbTree.rds'))


# Save plot
ggsave2(file.path('Plots', 'Supplement_HumanGeneSetSimilarity_xgbTree.pdf'), preds_plot, width = 7, height = 4)
