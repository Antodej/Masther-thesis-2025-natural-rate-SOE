#------------------------------------------------------------------------------#
# File:        percentiles.SOE.R
#
# Description: This is a function that output the values necessary in the 
#       posterior distribution table of the thesis. Inputs: posterior.draws
# Antoine DEJEAN (2025)
#------------------------------------------------------------------------------#
percentiles.SOE <- function(posterior.draws) {
  param_names <- colnames(posterior.draws)
  #returning a string "[X-Y]"
  ci_5_95 <- apply(posterior.draws, 2, function(x) {
    q <- quantile(x, probs = c(0.05, 0.95))
    paste0("[", round(q[1], 3), ", ", round(q[2], 3), "]")
  })
  
  means <- apply(posterior.draws, 2, function(x) round(mean(x), 3))
  
  summary_df <- rbind(ci_5_95, means)
  
  colnames(summary_df) <- param_names
  rownames(summary_df) <- c("5th_95th", "Mean")
  
  return(summary_df)}