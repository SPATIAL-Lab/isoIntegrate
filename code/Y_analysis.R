library(tidymodels)
library(lme4)
library(emmeans)

# Load Data and Functions
source("code/X_setup.R")

# ---- Regressions and Predictions --------------------------------------------
## Nonmovers
nm.lm <- lm(d18O.m ~ d18O.twzcta.mean, data = nonmovers)
nonmovers$residuals <- residuals(nm.lm)

## Movers - d18O Predictions
movers <- movers %>%
  mutate(
    d18O.pred = predict(nm.lm, newdata = .)
  )

# Movers Longest Res
m.lm <- lm(d18O.m ~ d18O.twzcta.mean, 
           data = movers %>%
             group_by(participant_id) %>%
             slice_max(order_by = timelength, n = 1, with_ties = FALSE) %>%
             ungroup())

# ---- Forward Models and Residuals --------------------------------------------
## Longest ----
longest <- movers %>%
  group_by(participant_id) %>%
  slice_max(order_by = timelength, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    participant_id,
    model = "longest",
    d18O.m = d18O.m,
    d18O.modeled = d18O.pred,
    residuals = d18O.m - d18O.pred,
    abs_residual = abs(residuals)
  )

## Simple Average ----
simple <- movers %>%
  group_by(participant_id) %>%
  summarise(
    model = "simple",
    d18O.m = first(d18O.m),
    d18O.modeled = mean(d18O.pred, na.rm = TRUE),
    residuals = first(d18O.m) - d18O.modeled,
    abs_residual = abs(residuals),
    .groups = "drop"
  )

## Time Weighted Average ----
weighted <- movers %>%
  group_by(participant_id) %>%
  summarise(
    model = "weighted",
    d18O.m = first(d18O.m),
    d18O.modeled = sum(d18O.pred * timelength, na.rm = TRUE) /
      sum(timelength, na.rm = TRUE),
    residuals = first(d18O.m) - d18O.modeled,
    abs_residual = abs(residuals),
    .groups = "drop"
  )

resid_long <- bind_rows(longest, simple, weighted)

## Shapiro-Wilks - Normality Check ----
# Non-Mover 
shapiro.test(nonmovers$residuals) 
# Longest
shapiro.test(longest$residuals) 
# Simple 
shapiro.test(simple$residuals)
# Mover Durational Model Residuals
shapiro.test(weighted$residuals)

## Levene's Test ----
lev_df <- bind_rows(
  nonmovers %>%
    transmute(
      residuals = residuals,
      model = "static"
    ),
  resid_long
)

car::leveneTest(residuals ~ model, data = lev_df)
car::leveneTest(residuals ~ model, 
                data = resid_long)
# ---- Error Metrics ----------------------------------------------------------
set.seed(123)
nm_mae <- bootstrap_error_single(nonmovers$residuals, mae, n_boot = 10000)
nm_rmse <- bootstrap_error_single(nonmovers$residuals, rmse, n_boot = 10000)

l_mae <- bootstrap_error_single(longest$residuals, mae, n_boot = 10000)
l_rmse <- bootstrap_error_single(longest$residuals, rmse, n_boot = 10000)

s_mae <- bootstrap_error_single(simple$residuals, mae, n_boot = 10000)
s_rmse <- bootstrap_error_single(simple$residuals, rmse, n_boot = 10000)

w_mae <- bootstrap_error_single(weighted$residuals, mae, n_boot = 10000)
w_rmse <- bootstrap_error_single(weighted$residuals, rmse, n_boot = 10000)

error_table <- data.frame(
  Model = c("Nonmover", "Longest", "Simple", "Weighted"),
  Median_Signed_Error = c(
    round(median(nonmovers$residuals), 3),
    round(median(longest$residuals), 3),
    round(median(simple$residuals), 3),
    round(median(weighted$residuals), 3)),
  Mean_Signed_Error = c(
    round(mean(nonmovers$residuals), 3),
    round(mean(longest$residuals), 3),
    round(mean(simple$residuals), 3),
    round(mean(weighted$residuals), 3)),
  Error_SD = c(
    round(sd(nonmovers$residuals), 3),
    round(sd(longest$residuals), 3),
    round(sd(simple$residuals), 3),
    round(sd(weighted$residuals), 3)),
  Median_Absolute_Error = c(
    round(median(abs(nonmovers$residuals)), 3),
    round(median(abs(longest$residuals)), 3),
    round(median(abs(simple$residuals)), 3),
    round(median(abs(weighted$residuals)), 3)),
  MAE_lower = c(
    round(nm_mae$ci_lower, 3),
    round(l_mae$ci_lower, 3),
    round(s_mae$ci_lower, 3),
    round(w_mae$ci_lower, 3)
  ),
  MAE = c(
    round(mae(nonmovers$residuals), 3),
    round(mae(longest$residuals), 3),
    round(mae(simple$residuals), 3),
    round(mae(weighted$residuals), 3)),
  MAE_upper = c(
    round(nm_mae$ci_upper, 3),
    round(l_mae$ci_upper, 3),
    round(s_mae$ci_upper, 3),
    round(w_mae$ci_upper, 3)
  ),
  RMSE_lower = c(
    round(nm_rmse$ci_lower, 3),
    round(l_rmse$ci_lower, 3),
    round(s_rmse$ci_lower, 3),
    round(w_rmse$ci_lower, 3)
  ),
  RMSE = c(
    round(rmse(nonmovers$residuals), 3),
    round(rmse(longest$residuals), 3),
    round(rmse(simple$residuals), 3),
    round(rmse(weighted$residuals), 3)
  ),
  RMSE_upper = c(
    round(nm_rmse$ci_upper, 3),
    round(l_rmse$ci_upper, 3),
    round(s_rmse$ci_upper, 3),
    round(w_rmse$ci_upper, 3)
  )
)

# Forward Model Pairwise Bootstrap -----
set.seed(124)
# Pairwise: Simple - Longest
rmse_svl <- bootstrap_error_pair(simple$residuals, longest$residuals,
                                 rmse, 10000)
mae_svl  <- bootstrap_error_pair(simple$residuals, longest$residuals,
                                 mae, 10000)

# Pairwise: Time Weighted - Longest
rmse_wvl <- bootstrap_error_pair(weighted$residuals,
                                 longest$residuals, rmse, 10000)

mae_wvl  <- bootstrap_error_pair(weighted$residuals,
                                 longest$residuals, mae, 10000)

# Pairwise: Time Weighted - Simple
rmse_wvs <- bootstrap_error_pair(weighted$residuals, simple$residuals,
                                 rmse, 10000)
mae_wvs  <- bootstrap_error_pair(weighted$residuals, simple$residuals,
                                 mae, 10000)

## Save the bootstrap results
all_boot <- list(
  "nm_rmse" = nm_rmse,
  "nm_mae"  = nm_mae,
  "l_rmse" = l_rmse,
  "l_mae"  = l_mae,
  "s_rmse" = s_rmse,
  "s_mae"  = s_mae,
  "w_rmse" = w_rmse,
  "w_mae"  = w_mae,
  "rmse_svl" = rmse_svl,
  "mae_svl"  = mae_svl,
  "rmse_wvl" = rmse_wvl,
  "mae_wvl"  = mae_wvl,
  "rmse_wvs" = rmse_wvs,
  "mae_wvs"  = mae_wvs)

# Save as DFs
svl_df <- make_boot_df(rmse_svl, mae_svl, "Simp. - Long.")
wvl_df <- make_boot_df(rmse_wvl, mae_wvl, "Time - Long.")
wvs_df <- make_boot_df(rmse_wvs, mae_wvs, "Time - Simp.")

# Bind
bootstrap_bound <- bind_rows(svl_df, wvl_df, wvs_df)
bootstrap_bound <- bootstrap_bound %>%
  mutate(
    comp = factor(comp, levels = c("Simp. - Long.", "Time - Long.",
                                   "Time - Simp."))
  )

# ---- Model Performance Mechanisms -------------------------------------------
## ---- Dataframe Set Up ----
### Mechanism 1: Isotopic Contrast Between Residence Drinking Waters
iso_diff <- movers %>%
  select(participant_id, res_type, d18O.twzcta.mean) %>%
  group_by(participant_id) %>%
  mutate(
    # Difference in d18O in water between residences
    diff_water = d18O.twzcta.mean - lag(d18O.twzcta.mean),
    # Average Difference
    avg_diff = mean(diff_water, na.rm = TRUE),
    abs_avg_diff = abs(avg_diff)
  ) %>%
  distinct(avg_diff, .keep_all = TRUE) %>%
  ungroup() %>%
  # Join with Residuals
  left_join(resid_long, by = "participant_id") %>%
  select(participant_id, avg_diff, abs_avg_diff, model, residuals,
         abs_residual)


### Mechanisms 2 & 3: Timing of Residences & Extent of Residence Contribution
time_summary <- movers %>%
  arrange(participant_id, timelength_beg) %>%
  group_by(participant_id) %>%
  mutate(
    # Proportion of Time
    time_prop = timelength/sum(timelength),
    # Index for residence number
    res_idx = row_number(),
    # Recalc number of residence
    n_res = n(),
    # Index for which residence number has greatest time proportion
    max_prop_idx = which.max(time_prop)
  ) %>%
  summarise(
    # Number of residences
    n_res = first(n_res),
    # Timing of longest-lived residence
    dom_age_beg  = timelength_beg[match(max_prop_idx[1], res_idx)],
    dom_age_end  = timelength_end[match(max_prop_idx[1], res_idx)],
    dom_duration = timelength[match(max_prop_idx[1], res_idx)],
    dom_prop     = time_prop[match(max_prop_idx[1], res_idx)],
    .groups   = "drop"
  ) %>%
  mutate(
    # Age Period Category
    incorp_cat = case_when(
      dom_age_beg <= 9 & dom_age_end >= 13 ~ "Consistent",
      dom_age_end <= 11 ~ "Early",
      dom_age_beg >= 11 ~ "Late",
      TRUE ~ "Intermediate" 
    )
  ) %>%
  left_join(resid_long, by = "participant_id")

# Sensitivity analysis for equal proportion longest-lived residences
## Create longest predictions using latest residence
longest_last <- movers %>%
  arrange(participant_id, timelength_beg) %>%
  group_by(participant_id) %>%
  filter(near(timelength, max(timelength))) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  transmute(
    participant_id,
    model = "longest",
    d18O.m,
    d18O.modeled = d18O.pred,
    residuals = d18O.m - d18O.pred,
    abs_residual = abs(residuals)
  )

## Bind with other models
sens_df <- bind_rows(longest_last, simple, weighted)

## Create new time summary selecting the last residence when equal proportions
time_summary_sens <- movers %>%
  arrange(participant_id, timelength_beg) %>%
  group_by(participant_id) %>%
  mutate(
    # Proportion of Time
    time_prop = timelength/sum(timelength),
    # Index for residence number
    res_idx = row_number(),
    # Recalc number of residence
    n_res = n(),
    # Index for which residence number has greatest time proportion - select
    # latest residence
    max_prop_idx = max(which(time_prop == max(time_prop)))
  ) %>%
  summarise(
    # Number of residences
    n_res = first(n_res),
    # Timing of longest-lived residence
    dom_age_beg  = timelength_beg[match(max_prop_idx[1], res_idx)],
    dom_age_end  = timelength_end[match(max_prop_idx[1], res_idx)],
    dom_duration = timelength[match(max_prop_idx[1], res_idx)],
    dom_prop     = time_prop[match(max_prop_idx[1], res_idx)],
    .groups   = "drop"
  ) %>%
  mutate(
    # Age Period Category
    incorp_cat = case_when(
      dom_age_beg <= 9 & dom_age_end >= 13 ~ "Consistent",
      dom_age_end <= 11 ~ "Early",
      dom_age_beg >= 11 ~ "Late",
      TRUE ~ "Intermediate" 
    )
  ) %>%
  left_join(sens_df, by = "participant_id")


### Mechanism 4: Number of Residences
num_res <- movers %>%
  group_by(participant_id) %>%
  summarize(
    num_res = n(),
    res_group = if_else(num_res > 2, "Multiple Moves", "One Move"),
    .keep_all = TRUE) %>%
  ungroup() %>%
  select(participant_id, num_res, res_group) %>%
  left_join(resid_long, by = "participant_id")

## ---- Mixed Effects Models ------
### Controls ----
resid_long$model <- relevel(
  factor(resid_long$model),
  ref = "weighted"
)

linear_control <- lmer(abs_residual ~ model + (1 | participant_id),
                       data = resid_long)

glmm_control <- lme4::glmer(abs_residual ~ model + (1 | participant_id),
                            family = Gamma(link = "log"),
                            control = glmerControl(optimizer = "nlminbwrap"),
                            data = resid_long)

### Mechanism 1: Isotopic Contrast ----
### Absolute Average Isotopic Difference Between Residences from Modeled Tap Water
iso_diff$model <- relevel(
  factor(iso_diff$model),
  ref = "weighted"
)

iso_model <- lme4::glmer(
  abs_residual ~ abs_avg_diff * model + (1 | participant_id),
  family = Gamma(link = "log"),
  control = glmerControl(optimizer = "nlminbwrap"),
  data = iso_diff
)

### Mechanisms 2 & 3: Residence Timing ----
time_summary$model <- relevel(
  factor(time_summary$model),
  ref = "weighted"
)

time_summary$incorp_cat <- relevel(
  factor(time_summary$incorp_cat),
  ref = "Early"
)

time_summary_sens$model <- relevel(
  factor(time_summary_sens$model),
  ref = "weighted"
)

time_summary_sens$incorp_cat <- relevel(
  factor(time_summary_sens$incorp_cat),
  ref = "Early"
)

#### Enamel Formation Timing at Longest-Lived Residence
incorp_model <- lme4::glmer(
  abs_residual ~ incorp_cat * model + (1 | participant_id),
  family = Gamma(link = "log"),
  control = glmerControl(optimizer = "nlminbwrap"),
  data = time_summary
)

incorp_model_sens <- lme4::glmer(
  abs_residual ~ incorp_cat * model + (1 | participant_id),
  family = Gamma(link = "log"),
  control = glmerControl(optimizer = "nlminbwrap"),
  data = time_summary_sens
)

#### Proportion of Time During Enamel Formation Lived at Longest-Lived Residence
prop_model <- lme4::glmer(
  abs_residual ~ dom_prop * model + (1 | participant_id),
  family = Gamma(link = "log"),
  control = glmerControl(optimizer = "nlminbwrap"),
  data = time_summary
)

### Mechanism 4: Number of Moves ----
num_res$model <- relevel(
  factor(num_res$model),
  ref = "weighted"
)

group_model <- lme4::glmer(
  abs_residual ~ res_group * model + (1 | participant_id),
  family = Gamma(link = "log"),
  control = glmerControl(optimizer = "nlminbwrap"),
  data = num_res
)

## GLMM Tables ----
model_list <- list(
  linear_control = linear_control,
  glmm_control = glmm_control,
  iso_model = iso_model,
  incorp_model = incorp_model,
  incorp_model_sens = incorp_model_sens,
  prop_model = prop_model,
  group_model = group_model
)

#### Metadata -----
model_metadata <- data.frame()

for (model_name in names(model_list)) {
  mod <- model_list[[model_name]]
  
  info <- data.frame(
    model = model_name,
    outcome = "abs_residual",
    family = family(mod)$family,
    link = family(mod)$link,
    formula = deparse(formula(mod)),
    n_obs = nobs(mod),
    n_participants = length(unique(model.frame(mod)$participant_id)),
    AIC = AIC(mod),
    BIC = BIC(mod),
    logLik = as.numeric(logLik(mod)),
    df_resid = df.residual(mod)
  )
  
  model_metadata <- bind_rows(model_metadata, info)
}

### Coefficient Summary ----
model_list$linear_control <- NULL
coef_big <- data.frame()

for (model_name in names(model_list)) {
  
  mod <- model_list[[model_name]]
  
  tab <- glmm_summary(mod, model_name = model_name)
  
  coef_big <- rbind(coef_big, tab)
}

### Coefficient Interpretation ----
## OPTIONAL: Pull coefficient table used in primary analysis
# coef_big <- readRDS("data/coef_big.rds") %>%
#   ## Optional Cleaning of DF
#   select("model_name", "term", "Estimate", "Std. Error", "res_ratio", "perc_change",
#          "t value", "Pr(>|z|)")

# IsoCon Interpretation
sub_df <- coef_big %>% filter(model_name == "iso_model")
## Extract Approach Estimates
### Time Weighted (Reference)
est_1 <- sub_df[sub_df$term == "abs_avg_diff", "Estimate"]
### Longest
est_2 <- sub_df[sub_df$term == "abs_avg_diff:modellongest", "Estimate"]
### Simple
est_3 <- sub_df[sub_df$term == "abs_avg_diff:modelsimple", "Estimate"]

## Create Interpretation DF
IsoCon_int <- data.frame(
  GLMM = "IsoCon",
  ### Approach Name
  approach = c("Time-Weighted", "Longest", "Simple"),
  ### Log Change Per Unit of IsoCon
  log_change = c(
    est_1,
    est_1 + est_2,
    est_1 + est_3
  ) 
)

IsoCon_int <- IsoCon_int %>%
  mutate(
    # Permille Response per Unit Change
    permille_unit_response = exp(log_change),
    # Percent Change Response per Unit Change
    perc_unit_response = 100 * (permille_unit_response - 1),
  )

# EnamForm Interpretaion
sub_df <- coef_big %>% filter(model_name == "incorp_model")
## Extract Approach Estimates
### Intercept
est_intercept <- sub_df[sub_df$term == "(Intercept)", "Estimate"]

### Incorporation Main Effects
est_consistent <- sub_df[sub_df$term == "incorp_catConsistent", "Estimate"]
est_inter      <- sub_df[sub_df$term == "incorp_catIntermediate", "Estimate"]
est_late       <- sub_df[sub_df$term == "incorp_catLate", "Estimate"]

### Model Effects
est_longest <- sub_df[sub_df$term == "modellongest", "Estimate"]
est_simple  <- sub_df[sub_df$term == "modelsimple", "Estimate"]

### Longest Interactions
est_longest_consistent <- sub_df[sub_df$term == "incorp_catConsistent:modellongest", "Estimate"]
est_longest_inter      <- sub_df[sub_df$term == "incorp_catIntermediate:modellongest", "Estimate"]
est_longest_late       <- sub_df[sub_df$term == "incorp_catLate:modellongest", "Estimate"]

### Simple Interactions
est_simple_consistent <- sub_df[sub_df$term == "incorp_catConsistent:modelsimple", "Estimate"]
est_simple_inter      <- sub_df[sub_df$term == "incorp_catIntermediate:modelsimple", "Estimate"]
est_simple_late       <- sub_df[sub_df$term == "incorp_catLate:modelsimple", "Estimate"]

## Create Interpretation DF
EnamForm_int <- data.frame(
  GLMM = "EnamForm",
  ### Approach Name
  approach = rep(c("Time-Weighted", "Longest", "Simple"), each = 4),
  ### Incorporation Category
  incorp_cat = rep(c("Early", "Consistent", "Intermediate", "Late"), times = 3),
  ### Expected Log Error per Incorporation Category
  log_expected_error = c(
    ### Time-Weighted
    est_intercept,
    est_intercept + est_consistent,
    est_intercept + est_inter,
    est_intercept + est_late,
    
    ### Longest
    est_intercept + est_longest,
    est_intercept + est_consistent + est_longest + est_longest_consistent,
    est_intercept + est_inter + est_longest + est_longest_inter,
    est_intercept + est_late + est_longest + est_longest_late,
    
    ### Simple
    est_intercept + est_simple,
    est_intercept + est_consistent + est_simple + est_simple_consistent,
    est_intercept + est_inter + est_simple + est_simple_inter,
    est_intercept + est_late + est_simple + est_simple_late
  ),
  ### Log Change Relative to Time-Weighted Within Each Category
  log_change = c(
    ### Time-Weighted
    0,
    0,
    0,
    0,
    
    ### Longest
    est_longest,
    est_longest + est_longest_consistent,
    est_longest + est_longest_inter,
    est_longest + est_longest_late,
    
    ### Simple
    est_simple,
    est_simple + est_simple_consistent,
    est_simple + est_simple_inter,
    est_simple + est_simple_late
  ) 
)

EnamForm_int <- EnamForm_int %>%
  mutate(
    # Expected Permille Error
    expected_error = exp(log_expected_error),
    # Response Ratio Relative to Time-Weighted Within Each Category
    response_ratio = exp(log_change),
    # Percent Change Relative to Time-Weighted Within Each Category
    perc_response = 100 * (response_ratio - 1)
  )

# EnamForm Sensitivity Interp
sub_df <- coef_big %>% filter(model_name == "incorp_model_sens")
## Extract Approach Estimates
### Intercept
est_intercept <- sub_df[sub_df$term == "(Intercept)", "Estimate"]

### Incorporation Main Effects
est_consistent <- sub_df[sub_df$term == "incorp_catConsistent", "Estimate"]
est_inter      <- sub_df[sub_df$term == "incorp_catIntermediate", "Estimate"]
est_late       <- sub_df[sub_df$term == "incorp_catLate", "Estimate"]

### Model Effects
est_longest <- sub_df[sub_df$term == "modellongest", "Estimate"]
est_simple  <- sub_df[sub_df$term == "modelsimple", "Estimate"]

### Longest Interactions
est_longest_consistent <- sub_df[sub_df$term == "incorp_catConsistent:modellongest", "Estimate"]
est_longest_inter      <- sub_df[sub_df$term == "incorp_catIntermediate:modellongest", "Estimate"]
est_longest_late       <- sub_df[sub_df$term == "incorp_catLate:modellongest", "Estimate"]

### Simple Interactions
est_simple_consistent <- sub_df[sub_df$term == "incorp_catConsistent:modelsimple", "Estimate"]
est_simple_inter      <- sub_df[sub_df$term == "incorp_catIntermediate:modelsimple", "Estimate"]
est_simple_late       <- sub_df[sub_df$term == "incorp_catLate:modelsimple", "Estimate"]

## Create Interpretation DF
EnamForm_sens_int <- data.frame(
  GLMM = "EnamForm_sens",
  ### Approach Name
  approach = rep(c("Time-Weighted", "Longest", "Simple"), each = 4),
  ### Incorporation Category
  incorp_cat = rep(c("Early", "Consistent", "Intermediate", "Late"), times = 3),
  ### Expected Log Error per Incorporation Category
  log_expected_error = c(
    ### Time-Weighted
    est_intercept,
    est_intercept + est_consistent,
    est_intercept + est_inter,
    est_intercept + est_late,
    
    ### Longest
    est_intercept + est_longest,
    est_intercept + est_consistent + est_longest + est_longest_consistent,
    est_intercept + est_inter + est_longest + est_longest_inter,
    est_intercept + est_late + est_longest + est_longest_late,
    
    ### Simple
    est_intercept + est_simple,
    est_intercept + est_consistent + est_simple + est_simple_consistent,
    est_intercept + est_inter + est_simple + est_simple_inter,
    est_intercept + est_late + est_simple + est_simple_late
  ),
  ### Log Change Relative to Time-Weighted Within Each Category
  log_change = c(
    ### Time-Weighted
    0,
    0,
    0,
    0,
    
    ### Longest
    est_longest,
    est_longest + est_longest_consistent,
    est_longest + est_longest_inter,
    est_longest + est_longest_late,
    
    ### Simple
    est_simple,
    est_simple + est_simple_consistent,
    est_simple + est_simple_inter,
    est_simple + est_simple_late
  ) 
)

EnamForm_sens_int <- EnamForm_sens_int %>%
  mutate(
    # Expected Permille Error
    expected_error = exp(log_expected_error),
    # Response Ratio Relative to Time-Weighted Within Each Category
    response_ratio = exp(log_change),
    # Percent Change Relative to Time-Weighted Within Each Category
    perc_response = 100 * (response_ratio - 1)
  )

# PropTime Interpretaion
sub_df <- coef_big %>% filter(model_name == "prop_model")
## Extract Approach Estimates
### Time Weighted (Reference)
est_1 <- sub_df[sub_df$term == "dom_prop", "Estimate"]
### Longest
est_2 <- sub_df[sub_df$term == "dom_prop:modellongest", "Estimate"]
### Simple
est_3 <- sub_df[sub_df$term == "dom_prop:modelsimple", "Estimate"]

## Create Interpretation DF
PropTime_int <- data.frame(
  GLMM = "PropTime",
  ### Approach Name
  approach = c("Time-Weighted", "Longest", "Simple"),
  ### Log Change Over Range of PropTime
  log_change = c(
    est_1,
    est_1 + est_2,
    est_1 + est_3
  ) 
)

PropTime_int <- PropTime_int %>%
  mutate(
    # Permille Response per Unit Change
    permille_unit_response = exp(log_change),
    # Percent Change Response per Unit Change
    perc_unit_response = 100 * (permille_unit_response - 1),
  )

# NumMove Interpretaion
sub_df <- coef_big %>% filter(model_name == "group_model")
## Extract Approach Estimates
### Time Weighted (Reference)
est_1 <- sub_df[sub_df$term == "res_groupOne Move", "Estimate"]
### Longest
est_2 <- sub_df[sub_df$term == "res_groupOne Move:modellongest", "Estimate"]
### Simple
est_3 <- sub_df[sub_df$term == "res_groupOne Move:modelsimple", "Estimate"]

## Create Interpretation DF
NumMove_int <- data.frame(
  GLMM = "NumMove",
  ### Approach Name
  approach = c("Time-Weighted", "Longest", "Simple"),
  ### Log Change for One Move
  log_change = c(
    est_1,
    est_1 + est_2,
    est_1 + est_3
  ) 
)

NumMove_int <- NumMove_int %>%
  mutate(
    # Permille Response per Unit Change
    permille_unit_response = exp(log_change),
    # Percent Change Response per Unit Change
    perc_unit_response = 100 * (permille_unit_response - 1),
  )

# Bind into Tables
interp_df <- bind_rows(IsoCon_int, PropTime_int, NumMove_int)
enamform_int_df <- bind_rows(EnamForm_int, EnamForm_sens_int)

## Estimated Marginal Means Comparisons ----
#IsoCon pairs 
pairs(emtrends(
  iso_model, 
  specs = ~ model,
  var = "abs_avg_diff"),
  type = "response"
)

# EnamForm pairs 
pairs(emmeans(
  incorp_model, 
  specs = ~ model | incorp_cat),
  type = "response"
)

# EnamForm sensitivity pairs 
pairs(emmeans(
  incorp_model_sens, 
  specs = ~ model | incorp_cat),
  type = "response"
)

# PropTime pairs
pairs(emtrends(
  prop_model, 
  specs = ~ model,
  var = "dom_prop"),
  type = "response"
)

# NumMove pairs
pairs(emmeans(
  group_model, 
  specs = ~ model |res_group),
  type = "response"
)