# Zoonosis predictor: Utility functions to calibrate predictions
require(betacal)

calibrate_preds <- function(predictions, calibration_preds, positive_name = 'True') {
	# 'predictions': a data frame of predictions requiring calibration
	# 'calibration_preds': a data frame of raw predictions for a calibration dataset, used to fit the calibration model
	# 
	# Note: this function applies to a single iteration (use with dplyr::group_map)
	invisible(capture.output(
		calibration_model <- tryCatch({
			beta_calibration(p = calibration_preds$RawScore,
											 y = calibration_preds$InfectsHumans == positive_name,
											 parameters = 'abm')
			
		}, error = function(e) {
			# Fitting 3-param version not always possible, but two-param version not as flexible, so
			# use this only as a fallback (some calibration is better than none...)
			beta_calibration(p = calibration_preds$RawScore,
											 y = calibration_preds$InfectsHumans == positive_name,
											 parameters = 'ab')
		})
	))
	
	predictions %>% 
		mutate(CalibratedScore = beta_predict(.data$RawScore, calibration_model),
					 CalibrationMethod = calibration_model$parameters)
}
