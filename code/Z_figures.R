library(ggplot)

# Load Analytical Code
source("code/Y_analysis.R")

# Figures ---------------------------------------------------------------------
model_levels <- c("longest", "simple", "weighted")

model_pallete <- c(
  longest = "#1464BA",
  simple = "#F2AF29",
  weighted = "#75B855"
)

model_labels <- c(
  longest = "Longest",
  simple = "Simple",
  weighted = "Time Weighted"
)

model_color_scale <- scale_color_manual(
  values = model_pallete,
  labels = model_labels,
  guide = guide_legend(nrow = 3)
)

model_fill_scale <- scale_fill_manual(
  values = model_pallete,
  labels = model_labels,
  guide = guide_legend(nrow = 3)
)

base_plot_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.8, "lines"),
    legend.position = "right",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

## Regression Figures ----
nm_labels <- nonmovers %>%
  summarise(
    r2 = summary(nm.lm)$adj.r.squared,
    intercept = coef(nm.lm)[1],
    slope = coef(nm.lm)[2],
    x = -18,
    y = -3,
    label = sprintf(
      'atop(delta*phantom(",")^18*O[C] == %.2f %s %.2f %%.%% delta*phantom(",")^18*O[tw], italic(R)^2 == %.3f)',
      intercept,
      ifelse(slope >= 0, "+", "-"),
      abs(slope),
      r2
    ),
    .groups = "drop"
  )

m_labels <- movers %>%
  summarise(
    r2 = summary(m.lm)$adj.r.squared,
    intercept = coef(m.lm)[1],
    slope = coef(m.lm)[2],
    x = -18,
    y = -3,
    label = sprintf(
      'atop(delta*phantom(",")^18*O[C] == %.2f %s %.2f %%.%% delta*phantom(",")^18*O[tw], italic(R)^2 == %.3f)',
      intercept,
      ifelse(slope >= 0, "+", "-"),
      abs(slope),
      r2
    ),
    .groups = "drop"
  )

# Plot w/ regression per predictor and source  
nm_regression <- ggplot(nonmovers, aes(x = d18O.twzcta.mean, y = d18O.m)) + 
  geom_point(shape = 21, color = "black", fill = "#CA2B30", 
             alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", se = FALSE, 
              color = "black") + 
  labs(
    x = bquote("Tap Water " * delta * phantom(",")^{18} * O * " (‰, VSMOW)"),
    y = bquote("Measured Enamel " * delta* phantom(",")^{18} * O[C] * " (‰, VPDB)")
  ) +
  xlim(-19, -1) +
  scale_y_continuous(limits= c(-2, -14), breaks = c(-5, -7.5, -10, -12.5)) +
  geom_text(
    data = nm_labels,
    aes(x = x, y = y, label = label),
    hjust = 0.15,
    vjust = 0.4,
    parse = TRUE,
    inherit.aes = FALSE,
    size = 2.25
  ) + 
  base_plot_theme +  
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

m_regression <- ggplot(movers, aes(x = d18O.twzcta.mean, y = d18O.m)) + 
  geom_point(shape = 21, color = "black", fill = "#1464BA", 
             alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", se = FALSE, 
              color = "black") + 
  labs(
    x = bquote("Tap Water " * delta * phantom(",")^{18} * O * " (‰, VSMOW)"),
    y = bquote("Measured Enamel " * delta* phantom(",")^{18} * O[C] * " (‰, VPDB)")
  ) +
  xlim(-19, -1) +
  scale_y_continuous(limits= c(-2, -14), breaks = c(-5, -7.5, -10, -12.5)) +
  geom_text(
    data = m_labels,
    aes(x = x, y = y, label = label),
    hjust = 0.15,
    vjust = 0.4,
    parse = TRUE,
    inherit.aes = FALSE,
    size = 2.25
  ) + 
  base_plot_theme +  
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## Error Metric Bootstrap Figure ----
### OPTIONAL: Load initial pairwise bootstrap 
#### bootstrap_bound <- readRDS("data/pairwise_results.rds")

# CI band for plotting
band_df <- bootstrap_bound %>%
  distinct(comp, metric, observed, lower, upper)

# Plot and Facet
error_bootstrap <- ggplot(bootstrap_bound, aes(x = diff, fill = comp)) +
  geom_rect(
    data = band_df,
    aes(xmin = lower, xmax = upper,
        ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.15,
    fill = "steelblue"
  ) +
  geom_histogram(bins = 20, boundary = 0) +
  scale_fill_manual(
    values = c(
      "Simp. - Long." = "gray80",
      "Time - Long." = "gray60",
      "Time - Simp." = "gray40"
    )
  ) +
  geom_vline(
    data = band_df,
    aes(xintercept = observed),
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_vline(aes(xintercept = 0),
             linetype = "solid",
             linewidth = 0.2) +
  scale_y_continuous(breaks = c(0, 2000, 4000)) + 
  facet_grid(comp ~ metric) +
  labs(
    x = "Error Difference (\u2030)",
    y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7)
  )

# ---- Mixed Effect Figures ---------------------------------------------------
### Isotopic Contrast ----
iso_new <- expand.grid(
  abs_avg_diff = seq(
    min(iso_diff$abs_avg_diff),
    max(iso_diff$abs_avg_diff),
    length.out = 100
  ),
  model = levels(factor(iso_diff$model))
)

iso_new$participant_id <- iso_diff$participant_id[1]

iso_new$fit <- predict(
  iso_model, iso_new,
  type = "response",
  re.form = NA
)

iso_new$model <- factor(iso_new$model, levels = model_levels)
iso_diff$model <- factor(iso_diff$model, levels = model_levels)

iso_model_plot <- ggplot() + 
  geom_point(data = iso_diff, aes(x = abs_avg_diff, y = abs_residual, 
                                  color = model), alpha = 0.45) +
  geom_line(data = iso_new, aes(x = abs_avg_diff, y = fit, color = model),
            linewidth = 0.7) +
  model_color_scale +
  labs(
    title = "Average Isotopic Contrast",
    x = expression("|" * Delta * delta^{18} * "O"["tw"] * "| (\u2030, VSMOW)" ),
    y = "|e|",
    color = "Approach"
  ) +
  base_plot_theme + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

### Incorporation Timing  ----
incorp_new <- expand.grid(
  incorp_cat = unique(time_summary$incorp_cat),
  model = unique(time_summary$model)
)

incorp_new$participant_id <- time_summary$participant_id[1]

incorp_new$fit <- predict(
  incorp_model,
  newdata = incorp_new,
  type = "response",
  re.form = NA
)

time_summary$incorp_cat <- factor(time_summary$incorp_cat,
                                  levels = c("Early", "Intermediate", "Late", "Consistent")
)

incorp_new$incorp_cat <- factor(incorp_new$incorp_cat,
                                levels = c("Early", "Intermediate", "Late", "Consistent")
)

time_summary$model <- factor(time_summary$model, levels = model_levels)
incorp_new$model <- factor(incorp_new$model, levels = model_levels)

group_label <- time_summary %>% 
  filter(model == "longest") %>%
  group_by(incorp_cat) %>%
  summarise(
    n = n(),
    y = 4.5
  )

incorp_model_plot <- ggplot() +
  geom_boxplot(
    data = time_summary,
    aes(x = incorp_cat, y = abs_residual, fill = model),
    position = position_dodge(width = 0.8),
    alpha = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    data = time_summary,
    aes(x = incorp_cat, y = abs_residual, color = model),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.8
    ),
    alpha = 0.20,
    size = 1.3
  ) +
  geom_point(
    data = incorp_new,
    aes(x = incorp_cat, y = fit, color = model),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  geom_text(
    data = group_label, 
    aes(label = paste("n = ", n), 
        x = incorp_cat, y = y), 
    vjust = -0.2
  ) +
  model_color_scale +
  model_fill_scale +
  labs(
    title = "Enamel Formation Timing at Longest-Lived Residence",
    x = NULL,
    y = "|e|",
    fill = "Approach",
    color = "GLMM Predicted |e|"
  ) +
  base_plot_theme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
         fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")) + 
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = 10, hjust = 1),
        legend.text  = element_text(size = 8),
        legend.spacing.y = grid::unit(4, "pt"),
        legend.key.height = grid::unit(10, "pt"),
        legend.margin = margin(0)
  )

### Proportion of Enamel Formation Occurring at Dominant Residence ----
prop_new <- expand.grid(
  dom_prop = seq(
    min(time_summary$dom_prop),
    max(time_summary$dom_prop),
    length = 100
  ),
  model = unique(time_summary$model)
)

prop_new$participant_id <- time_summary$participant_id[1]

prop_new$fit <- predict(
  prop_model,
  newdata = prop_new,
  type = "response",
  re.form = NA
)

prop_new$model <- factor(prop_new$model, levels = model_levels)

prop_model_plot <- ggplot() +
  geom_point(data = time_summary, aes(x = dom_prop, y = abs_residual, 
                                      color = model), alpha = 0.45) +
  geom_line(
    data = prop_new,
    aes(x = dom_prop, y = fit, color = model, group = model),
    linewidth = 1.1
  ) +
  model_color_scale +
  scale_x_continuous(labels = scales::label_percent())+
  labs(
    title = "Proportion of Time at Longest-Lived Residence",
    x = NULL,
    y = "|e|",
    color = "Approach"
  ) +
  coord_cartesian(xlim = c(0.2, 0.9), ylim = c(0, 2)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  base_plot_theme + theme(legend.position = "bottom")

### Number of Moves ----
group_new <- expand.grid(
  res_group = unique(num_res$res_group),
  model = unique(num_res$model)
)

group_new$participant_id <- num_res$participant_id[1]

group_new$fit <- predict(
  group_model,
  newdata = group_new,
  type = "response",
  re.form = NA
)

group_label <- num_res %>% 
  filter(model == "longest") %>%
  group_by(res_group) %>%
  summarise(
    n = n(),
    y = 4.5
  )

num_res$model <- factor(num_res$model, levels = model_levels)
group_new$model <- factor(group_new$model, levels = model_levels)

group_model_plot <- ggplot() +
  geom_boxplot(
    data = num_res,
    aes(x = res_group, y = abs_residual, fill = model),
    position = position_dodge(width = 0.8),
    alpha = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    data = num_res,
    aes(x = res_group, y = abs_residual, color = model),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.8
    ),
    alpha = 0.20,
    size = 1.3
  ) +
  geom_point(
    data = group_new,
    aes(x = res_group, y = fit, color = model),
    position = position_dodge(width = 0.8),
    size = 3
  ) +
  geom_text(
    data = group_label, 
    aes(label = paste("n = ", n), 
        x = res_group, y = y), 
    vjust = -0.2
  ) +
  model_color_scale +
  model_fill_scale +
  labs(
    title = "Number of Moves",
    x = NULL,
    y = NULL,
    fill = "Approach",
    color = "GLMM Predicted |e|"
  ) +
  base_plot_theme 

# ---- Other Figures ----------------------------------------------------------
## Residual Density -----
ref_df <- nonmovers %>% 
  transmute(
    ref_model = "Static", 
    residuals = residuals)

resid_df <- bind_rows(
  longest %>% transmute(
    model = "Longest",
    residuals = residuals),
  simple %>% transmute(
    model = "Simple",
    residuals = residuals),
  weighted %>% transmute(
    model = "Weighted",
    residuals = residuals)) %>%
  mutate(
    model = factor(model,
                   levels = c("Longest", "Simple", "Weighted")
    ))

# Plot Facet
resid_density <- ggplot() +
  geom_density(data = ref_df, 
               aes(x = residuals, fill = "Static"), 
               color = "black", 
               linewidth = 0.5,
               alpha = 0.8) +
  geom_density(data = resid_df, 
               aes(x = residuals, fill = model),
               color = "black", 
               linewidth = 0.5,
               alpha = 0.6) +
  facet_wrap(~ model, nrow = 3) +
  scale_fill_manual(
    name = "Model",
    breaks =c("Static", "Longest", "Simple", "Weighted"),
    values = c(
      Static = "#CA2B30",
      Longest = "#1464BA",
      Simple = "#F2AF29",
      Weighted = "#75B855"
    ),
    labels = c(
      Static = "Non-Mobile",
      Longest = "Longest",
      Simple = "Simple",
      Weighted = "Time Weighted"
    )
  ) +
  facet_wrap(~ model, nrow = 1) +
  coord_flip(xlim = c(-5, 5)) +
  labs(
    x = "Residuals",
    y = "Density",
    fill = NULL
  ) +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) 

## QQ Plot ----- 
resid_df <- bind_rows(
  nonmovers %>% 
    transmute(
      model = "Non-Mobile Reference", 
      residuals = residuals),
  longest %>% transmute(
    model = "Longest",
    residuals = residuals),
  simple %>% transmute(
    model = "Simple",
    residuals = residuals),
  weighted %>% transmute(
    model = "Time Weighted",
    residuals = residuals)) %>%
  mutate(
    model = factor(model,
                   levels = c("Non-Mobile Reference", "Longest", "Simple", 
                              "Time Weighted")
    ))

resid_qq <- ggplot(resid_df, aes(sample = residuals)) +
  stat_qq(size = 1) +
  stat_qq_line(color = "red", linewidth = 0.5) +
  facet_wrap(~ model, nrow = 2) +
  labs(
    x = "Expected Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )