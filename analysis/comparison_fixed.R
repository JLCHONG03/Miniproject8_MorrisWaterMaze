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
pc_fixed <- read.csv("place_cell_trial_mean_over_runs.csv")
dc_fixed <- read.csv("distance_cell_trial_mean_over_runs.csv")

## MERGE BOTH DATASETS
df_fixed_plot <- merge(
  pc_fixed, dc_fixed,
  by = c("day", "trial", "trial_global"),
  all = TRUE
)

## AVERAGE THE TRIALS PMs PER DAY FOR VISUALIZATION PLOTS
day_summary_fixed <- df_fixed_plot %>%
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
    
    target_quadrant_pc = mean(target_quadrant_pc, na.rm = TRUE),
    target_quadrant_pc_se = mean(target_quadrant_pc_se, na.rm = TRUE),
    target_quadrant_dc = mean(target_quadrant_dc, na.rm = TRUE),
    target_quadrant_dc_se = mean(target_quadrant_dc_se, na.rm = TRUE),
    
    opposite_quadrant_pc = mean(opposite_quadrant_pc, na.rm = TRUE),
    opposite_quadrant_pc_se = mean(opposite_quadrant_pc_se, na.rm = TRUE),
    opposite_quadrant_dc = mean(opposite_quadrant_dc, na.rm = TRUE),
    opposite_quadrant_dc_se = mean(opposite_quadrant_dc_se, na.rm = TRUE),
    
    wall_zone_pc = mean(wall_zone_pc, na.rm = TRUE),
    wall_zone_pc_se = mean(wall_zone_pc_se, na.rm = TRUE),
    wall_zone_dc = mean(wall_zone_dc, na.rm = TRUE),
    wall_zone_dc_se = mean(wall_zone_dc_se, na.rm = TRUE),
    .groups = "drop"
  )

## HELPER FUNCTION FOR PLOTS
plot_pm_day_fixed <- function(df, pc_mean, pc_se, dc_mean, dc_se, ylab_text, title_text, panel_letter) {
  
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
p1_fixed_day <- plot_pm_day_fixed(
  day_summary_fixed,
  "latency_pc", "latency_pc_se",
  "latency_dc", "latency_dc_se",
  "Latency (s)", "Latency", "A"
)

p2_fixed_day <- plot_pm_day_fixed(
  day_summary_fixed,
  "dist_pc", "dist_pc_se",
  "dist_dc", "dist_dc_se",
  "Distance travelled", "Distance", "B"
)

p3_fixed_day <- plot_pm_day_fixed(
  day_summary_fixed,
  "target_quadrant_pc", "target_quadrant_pc_se",
  "target_quadrant_dc", "target_quadrant_dc_se",
  "Target quadrant (%)", "Target quadrant time (%)", "C"
)

p4_fixed_day <- plot_pm_day_fixed(
  day_summary_fixed,
  "opposite_quadrant_pc", "opposite_quadrant_pc_se",
  "opposite_quadrant_dc", "opposite_quadrant_dc_se",
  "Opposite quadrant (%)", "Opposite quadrant time (%)", "D"
)

p5_fixed_day <- plot_pm_day_fixed(
  day_summary_fixed,
  "wall_zone_pc", "wall_zone_pc_se",
  "wall_zone_dc", "wall_zone_dc_se",
  "Wall zone (%)", "Wall zone time (%)", "E"
)

print(p1_fixed_day)
print(p2_fixed_day)
print(p3_fixed_day)
print(p4_fixed_day)
print(p5_fixed_day)

plot_figure_day <- wrap_plots(
  p1_fixed_day, p2_fixed_day, p3_fixed_day, p4_fixed_day, p5_fixed_day,
  ncol = 3,
  guides = "collect"
) &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.size = unit(0.9, "cm")
  )

plot_figure_day <- plot_figure_day +
  plot_annotation(
    title = "Fixed Platform PMs",
    theme = theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5)
    )
  )

print(plot_figure_day)

ggsave("fixed_platform_figure.png", plot_figure_day, width = 16, height = 9, dpi = 300)

## READ RUN-LEVEL CSVs FOR MIXED MODELS
pc_run <- read.csv("place_cell_fixed_run_level.csv")
dc_run <- read.csv("distance_cell_fixed_run_level.csv")

pc_run$model <- "Place-cell"
dc_run$model <- "Distance-cell"

## RENAME PMs COLUMNS TO COMMON NAME BEFORE MERGING
pc_run <- pc_run %>%
  rename(
    target = target_quadrant,
    opposite = opposite_quadrant,
    wall = wall_zone
  )

dc_run <- dc_run %>%
  rename(
    target = target_quadrant,
    opposite = opposite_quadrant,
    wall = wall_zone
  )

## MERGE DATASETS
df_fixed <- bind_rows(pc_run, dc_run)

df_fixed$model <- factor(df_fixed$model)
df_fixed$day <- factor(df_fixed$day)
df_fixed$trial <- factor(df_fixed$trial)
df_fixed$run <- factor(df_fixed$run)

## MIXED MODELS
model_latency_fixed <- lmer(
  latency ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_fixed
)

model_target_fixed <- lmer(
  target ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_fixed
)

model_wall_fixed <- lmer(
  wall ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_fixed
)

model_distance_fixed <- lmer(
  distance ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_fixed
)

model_opposite_fixed <- lmer(
  opposite ~ model + day + trial + model:day + model:trial + (1 | run),
  data = df_fixed
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

lat_res_fixed <- extract_fixed_effects(model_latency_fixed)
dis_res_fixed <- extract_fixed_effects(model_distance_fixed)
tar_res_fixed <- extract_fixed_effects(model_target_fixed)
opp_res_fixed <- extract_fixed_effects(model_opposite_fixed)
wal_res_fixed <- extract_fixed_effects(model_wall_fixed)

summary_table_fixed <- data.frame(
  Effect = c("Model", "Day", "Trial", "Model × Day", "Model × Trial"),
  Latency = lat_res_fixed,
  Distance = dis_res_fixed,
  Target_Quadrant = tar_res_fixed,
  Opposite_Quadrant = opp_res_fixed,
  Wall_Zone = wal_res_fixed,
  row.names = NULL
)

print(summary_table_fixed)

##BUILD TABLE
table_theme <- ttheme_minimal(
  base_size = 11,
  padding = unit(c(2.5, 4), "mm"),
  core = list(
    fg_params = list(fontsize = 11, fontface = 1),
    bg_params = list(
      fill = rep(c("#F2F2F2", "#EAEAEA"), length.out = nrow(summary_table_fixed)),
      col = "white"
    )
  ),
  colhead = list(
    fg_params = list(fontsize = 12, fontface = "bold"),
    bg_params = list(fill = "#D9D9D9", col = "white")
  )
)

table_grob <- tableGrob(
  summary_table_fixed,
  rows = NULL,
  theme = table_theme
)

table_grob$widths <- unit(
  c(1.4, 1.0, 1.0, 1.5, 1.7, 1.1),
  "in"
)

table_grob$heights <- rep(unit(0.42, "in"), length(table_grob$heights))

title_box <- grobTree(
  rectGrob(gp = gpar(fill = "#E6E6E6", col = "grey80")),
  textGrob(
    "Fixed platform: mixed-effects model (fixed = model, day, trial; random = run)",
    gp = gpar(fontsize = 11, fontface = "bold")
  )
)

compact_table <- arrangeGrob(
  title_box,
  table_grob,
  ncol = 1,
  heights = c(0.18, 1)
)

table_only <- as.ggplot(compact_table)

print(table_only)

ggsave("fixed_platform_table.png", table_only, width = 8, height = 4, dpi = 300)
