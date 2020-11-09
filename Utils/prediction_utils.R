## Utility functions: predict from saved model objects
## - To reduce size, only the trained xgboost object and a few relevant components are saved
## - This means we can no longer use caret directly to handle data transformations required 
## 		for prediction (though the code below is based on caret's implementation)

# NOTE: script should be loaded alongside xgboost_utils.R, which provides the 'as_caret_data' function


prepare_prediction_data <- function(fit, newdata, labelcol = "Name") {
	# Get data in the format expected by caret:
	modelFit <- fit$finalModel
	
	caretdata <- newdata[, colnames(newdata) %in% c(modelFit$xNames, labelcol), drop = FALSE]
	caretdata <- as_caret_data(caretdata, labelcol = labelcol, removecols = "")  # labelcol does not matter here
	caretdata <- caretdata[, modelFit$xNames, drop = FALSE]  # Need exactly the same order
	
	caretdata
}

get_prediction <- function(fit, newdata, virusnames, model_type = 'xgbTree') {
	modelFit <- fit$finalModel
	
	# Prediction functions for this model type:
	modelInfo <- getModelInfo(model_type, regex = FALSE)[[1]]
	
	# Predictions:
	probs <- modelInfo$prob(modelFit = modelFit, newdata = newdata)
	preds <- modelInfo$predict(modelFit = modelFit, newdata = newdata)
	
	# Return
	data.frame(Iteration = fit$Iteration,
						 Name = virusnames,
						 Prediction = preds,
						 probs,
						 stringsAsFactors = FALSE)
}
