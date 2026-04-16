library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(patchwork)
library(gridExtra)
library(grid)
library(ggplotify)

options(contrasts = c("contr.sum", "contr.poly"))

## READ MEAN + SE CSVs
pc_var <- read.csv("place_cell_variable_trial_mean_over_runs.csv")
dc_var <- read.csv("distance_cell_variable_trial_mean_over_runs.csv")

## CHANGE "TARGET QUADRANT" AND "OPPOSITE QUADRANT" PMs TO "POSSIBLE QUADRANT" FOR VARIABLE PLATFORM CASE
pc_var <- pc_var %>%
  mutate(
    possible_quadrants_pc = target_quadrant_pc + opposite_quadrant_pc,
    possible_quadrants_pc_se = sqrt(target_quadrant_pc_se^2 + opposite_quadrant_pc_se^2)
  )

dc_var <- dc_var %>%
  mutate(
    possible_quadrants_dc = target_quadrant_dc + opposite_quadrant_dc,
    possible_quadrants_dc_se = sqrt(target_quadrant_dc_se^2 + opposite_quadrant_dc_se^2)
  )

## MERGE BOTH DATASETS
df_var_plot <- merge(
  pc_var, dc_var,
  by = c("day", "trial", "trial_global"),
  all = TRUE
)

## AVERAGE THE TRIALS PMs PER DAY FOR VISUALIZATION PLOTS
day_summary_var <- df_var_plot %>%
  group_by(day) %>%
  summarise(
    latency_pc = mean(latency_pc, na.rm = TRUE),
    latency_pc_se = mean(latency_pc_se, na.rm = TRUE),
    latency_dc = mean(latency_dc, na.rm = TRUE),
    latency_dc_se = mean(latency_dc_se, na.rm = TRUE),
    
    dist_pc = mean(dist_pc, na.rm = TRUE),
    dist_pc_se = mean(dist_pc_se, na.rm = TRUE),
    dist_dc = mean(dist_dc, na.rm = TRUE),
    dist_dc_se = mean(dist_dc_se, na.rm = TRUE),
    
    possible_quadrants_pc = mean(possible_quadrants_pc, na.rm = TRUE),
    possible_quadrants_pc_se = mean(possible_quadrants_pc_se, na.rm = TRUE),
    possible_quadrants_dc = mean(possible_quadrants_dc, na.rm = TRUE),
    possible_quadrants_dc_se = mean(possible_quadrants_dc_se, na.rm = TRUE),
    
    wall_zone_pc = mean(wall_zone_pc, na.rm = TRUE),
    wall_zone_pc_se = mean(wall_zone_pc_se, na.rm = TRUE),
    wall_zone_dc = mean(wall_zone_dc, na.rm = TRUE),
    wall_zone_dc_se = mean(wall_zone_dc_se, na.rm = TRUE),
    .groups = "drop"
  )

## HELPER FUNCTION FOR PLOTS
plot_pm_day_var <- function(df, pc_mean, pc_se, dc_mean, dc_se, ylab_text, title_text, panel_letter) {
  
  plot_df <- data.frame(
    day = rep(df$day, 2),
    mean = c(df[[pc_mean]], df[[dc_mean]]),
    se = c(df[[pc_se]], df[[dc_se]]),
    model = rep(c("Place-cell", "Distance-cell"), each = nrow(df))
  )
  
  ymax <- max(plot_df$mean + plot_df$se, na.rm = TRUE)
  ymin <- min(plot_df$mean - plot_df$se, na.rm = TRUE)
  
  ggplot(plot_df, aes(x = day, y = mean, color = model, group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15, linewidth = 0.4) +
    annotate(
      "text",
      x = 0.6,
      y = ymax * 1.08,
      label = panel_letter,
      fontface = "bold",
      size = 5,
      hjust = 0
    ) +
    scale_x_continuous(breaks = 1:8) +
    scale_color_manual(values = c("Distance-cell" = "red1", "Place-cell" = "navyblue")) +
    labs(
      x = "Day",
      y = ylab_text,
      title = title_text,
      color = "Model"
    ) +
    coord_cartesian(ylim = c(ymin, ymax * 1.12), clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

## PLOTS
p1_var_day <- plot_pm_day_var(
  day_summary_var,
  "latency_pc", "latency_pc_se",
  "latency_dc", "latency_dc_se",
  "Latency (s)", "Latency", "A"
)

p2_var_day <- plot_pm_day_var(
  day_summary_var,
  "dist_pc", "dist_pc_se",
  "dist_dc", "dist_dc_se",
  "Distance travelled", "Distance", "B"
)

p3_var_day <- plot_pm_day_var(
  day_summary_var,
  "possible_quadrants_pc", "possible_quadrants_pc_se",
  "possible_quadrants_dc", "possible_quadrants_dc_se",
  "Possible quadrants (%)", "Time in possible quadrants (%)", "C"
)

p4_var_day <- plot_pm_day_var(
  day_summary_var,
  "wall_zone_pc", "wall_zone_pc_se",
  "wall_zone_dc", "wall_zone_dc_se",
  "Wall zone (%)", "Wall zone time (%)", "D"
)

print(p1_var_day)
print(p2_var_day)
print(p3_var_day)
print(p4_var_day)

plot_figure_day_var <- wrap_plots(
  p1_var_day, p2_var_day, p3_var_day, p4_var_day,
  ncol = 2,
  guides = "collect"
) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  )

plot_figure_day_var <- plot_figure_day_var +
  plot_annotation(
    title = "Variable Platform PMs",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

print(plot_figure_day_var)

ggsave("variable_platform_figure.png", plot_figure_day_var, width = 12, height = 8, dpi = 300)

## READ RUN-LEVEL CSVs FOR MIXED MODELS
pc_run <- read.csv("place_cell_variable_run_level.csv")
dc_run <- read.csv("distance_cell_variable_run_level.csv")

pc_run$model <- "Place-cell"
dc_run$model <- "Distance-cell"

# CHANGE "TARGET QUADRANT" AND "OPPOSITE QUADRANT" PMs TO "POSSIBLE QUADRANT" FOR VARIABLE PLATFORM CASE
pc_run <- pc_run %>%
  mutate(
    possible_quadrants = target_quadrant + opposite_quadrant,
    wall = wall_zone
  )

dc_run <- dc_run %>%
  mutate(
    possible_quadrants = target_quadrant + opposite_quadrant,
    wall = wall_zone
  )

## MERGE DATASETS
df_var <- bind_rows(pc_run, dc_run)

## FACTORIZE THE EFFECTS 
df_var$model <- factor(df_var$model)
df_var$day <- factor(df_var$day)
df_var$trial <- factor(df_var$trial)
df_var$run <- factor(df_var$run)

## MIXED MODELS
model_latency_var <- lmer(
  latency ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_var
)

model_distance_var <- lmer(
  distance ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_var
)

model_possible_var <- lmer(
  possible_quadrants ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_var
)

model_wall_var <- lmer(
  wall ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_var
)

## SIGNIFICANCE STARS
p_to_star <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("n.s.")
}

## TYPE 3 ANOVA TO TEST WHETHER THE MODEL'S SIGNIFICANCE
extract_fixed_effects <- function(model_obj) {
  tab <- as.data.frame(anova(model_obj, type = 3))
  tab$Effect <- rownames(tab)
  rownames(tab) <- NULL
  
  effects_keep <- c("model", "day", "trial", "model:day", "model:trial")
  
  out <- sapply(effects_keep, function(eff) {
    pval <- tab$`Pr(>F)`[tab$Effect == eff]
    if (length(pval) == 0) return("")
    p_to_star(pval[1])
  })
  
  return(out)
}

lat_res_var <- extract_fixed_effects(model_latency_var)
dis_res_var <- extract_fixed_effects(model_distance_var)
pos_res_var <- extract_fixed_effects(model_possible_var)
wal_res_var <- extract_fixed_effects(model_wall_var)

summary_table_var <- data.frame(
  Effect = c("Model", "Day", "Trial", "Model × Day", "Model × Trial"),
  Latency = lat_res_var,
  Distance = dis_res_var,
  Possible_Quadrants = pos_res_var,
  Wall_Zone = wal_res_var,
  row.names = NULL
)

print(summary_table_var)

write.csv(summary_table_var, "variable_platform_mixed_model_summary.csv", row.names = FALSE)

## BUILD TABLE
table_theme <- ttheme_minimal(
  base_size = 11,
  padding = unit(c(2.5, 4), "mm"),
  core = list(
    fg_params = list(fontsize = 11, fontface = 1),
    bg_params = list(
      fill = rep(c("#F2F2F2", "#EAEAEA"), length.out = nrow(summary_table_var)),
      col = "white"
    )
  ),
  colhead = list(
    fg_params = list(fontsize = 12, fontface = "bold"),
    bg_params = list(fill = "#D9D9D9", col = "white")
  )
)

table_grob <- tableGrob(
  summary_table_var,
  rows = NULL,
  theme = table_theme
)

table_grob$widths <- unit(
  c(1.6, 1.0, 1.0, 1.8, 1.1),
  "in"
)

table_grob$heights <- rep(unit(0.42, "in"), length(table_grob$heights))

title_box <- grobTree(
  rectGrob(gp = gpar(fill = "#E6E6E6", col = "grey80")),
  textGrob(
    "Variable platform: mixed-effects model (fixed = model, day, trial; random = run)",
    gp = gpar(fontsize = 11, fontface = "bold")
  )
)

compact_table_var <- arrangeGrob(
  title_box,
  table_grob,
  ncol = 1,
  heights = c(0.18, 1)
)

table_only_var <- as.ggplot(compact_table_var)

print(table_only_var)

ggsave("variable_platform_table.png", table_only_var, width = 8, height = 4, dpi = 300)
