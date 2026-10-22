# Pull in Data ----
movers    <- read.csv("data/movers_deident.csv")
nonmovers <- read.csv("data/nonmovers_deident.csv")

# ---- Functions ---------------------------------------------------------------
# Root Mean Square Error
rmse <- function(resid) sqrt(mean(resid^2))

# Mean Absolute Error
mae  <- function(resid) mean(abs(resid))

# Single Error Metric Bootstrap
bootstrap_error_single <- function(resid, metric_fn, n_boot = 10000) {
  ## Remove missing values
  resid <- resid[!is.na(resid)]
  ## Storage
  boot_distribution <- numeric(n_boot)
  ## Loop for bootstrap iterations
  for (i in seq_len(n_boot)) {
    ### Resample residuals with replacement
    idx <- sample(length(resid), replace = TRUE)
    ### Calculate metric on bootstrapped residuals
    boot_distribution[i] <- metric_fn(resid[idx])
  }
  
  ## Observed metric from original residuals
  observed <- metric_fn(resid)
  
  ## Output
  list(
    observed          = observed,
    boot_distribution = boot_distribution,
    ci_lower          = quantile(boot_distribution, 0.025),
    ci_upper          = quantile(boot_distribution, 0.975)
  )
}

# Pairwise Bootstrap on Error Metric
bootstrap_error_pair <- function(resid1, resid2, metric_fn, n_boot = 10000) {
  ## Storage
  boot_distribution <- numeric(n_boot)
  ## Loop for iters.
  for (i in seq_len(n_boot)) {
    ### Sample a random residual location - assumes dfs are ordered the same
    idx <- sample(length(resid1), replace = TRUE)
    ### Use the index to calculate the difference in the error metric
    boot_distribution[i] <- metric_fn(resid1[idx]) - metric_fn(resid2[idx])
  }
  
  ## Take raw difference
  observed_diff <- metric_fn(resid1) - metric_fn(resid2)
  
  ## Out
  list(
    observed_diff     = observed_diff,
    boot_distribution = boot_distribution,
    ci_lower          = quantile(boot_distribution, 0.025),
    ci_upper          = quantile(boot_distribution, 0.975)
  )
}

# Convert Bootstrap List to Distribution DF
make_boot_df <- function(rmse_obj, mae_obj, comp_label) {
  rbind(
    data.frame(
      iter = seq_along(rmse_obj$boot_distribution),
      comp = comp_label,
      metric = "RMSE",
      diff = rmse_obj$boot_distribution,
      observed = rmse_obj$observed_diff,
      lower = rmse_obj$ci_lower,
      upper = rmse_obj$ci_upper,
      stringsAsFactors = FALSE
    ),
    data.frame(
      iter = seq_along(mae_obj$boot_distribution),
      comp = comp_label,
      metric = "MAE",
      diff = mae_obj$boot_distribution,
      observed = mae_obj$observed_diff,
      lower = mae_obj$ci_lower,
      upper = mae_obj$ci_upper,
      stringsAsFactors = FALSE
    )
  )
}

# Table for GLMM Results
glmm_summary <- function(model, model_name = NA) {
  # Extract Summary
  s <- summary(model)
  
  # Pull Model Coeffecients
  coef_tab <- as.data.frame(s$coefficients)
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  # Convert Estimate From Log, Report Percent Change
  coef_tab$model_name <- model_name
  coef_tab$res_ratio  <- exp(coef_tab$Estimate)
  coef_tab$perc_change <- 100 * (coef_tab$res_ratio - 1)
  
  return(coef_tab)
}
