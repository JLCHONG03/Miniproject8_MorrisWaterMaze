library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(gridExtra)
library(grid)
library(ggplotify)

options(contrasts = c("contr.sum", "contr.poly"))

## FILES (VARIABLE VERSION)
var_summary_files <- c(
  "combined_place1p0_wall0p0_var.csv",
  "combined_place0p8_wall0p2_var.csv",
  "combined_place0p6_wall0p4_var.csv",
  "combined_place0p5_wall0p5_var.csv",
  "combined_place0p4_wall0p6_var.csv",
  "combined_place0p2_wall0p8_var.csv",
  "combined_place0p0_wall1p0_var.csv"
)

var_run_files <- c(
  "combined_place1p0_wall0p0_var_run_level.csv",
  "combined_place0p8_wall0p2_var_run_level.csv",
  "combined_place0p6_wall0p4_var_run_level.csv",
  "combined_place0p5_wall0p5_var_run_level.csv",
  "combined_place0p4_wall0p6_var_run_level.csv",
  "combined_place0p2_wall0p8_var_run_level.csv",
  "combined_place0p0_wall1p0_var_run_level.csv"
)

## EXTRACT MODEL RATIO
extract_ratios <- function(file) {
  name <- tools::file_path_sans_ext(basename(file))
  
  place_txt <- stringr::str_extract(name, "(?<=place)[0-9p]+")
  wall_txt  <- stringr::str_extract(name, "(?<=wall)[0-9p]+")
  
  place_val <- as.numeric(gsub("p", ".", place_txt))
  wall_val  <- as.numeric(gsub("p", ".", wall_txt))
  
  group_label <- dplyr::case_when(
    place_val == 1 & wall_val == 0 ~ "place-only",
    place_val == 0 & wall_val == 1 ~ "distance-only",
    TRUE ~ paste0("place:", place_val)
  )
  
  list(
    wall_ratio = wall_val,
    place_ratio = place_val,
    group = group_label
  )
}

group_levels <- c(
  "place-only",
  "place:0.8",
  "place:0.6",
  "place:0.5",
  "place:0.4",
  "place:0.2",
  "distance-only"
)

## HELPERS FUNCTION TO READ CSV FILE
read_summary_file <- function(file) {
  rr <- extract_ratios(file)
  df <- read.csv(file)
  
  df %>%
    mutate(
      wall_ratio = rr$wall_ratio,
      place_ratio = rr$place_ratio,
      group = rr$group
    )
}

read_run_file <- function(file) {
  rr <- extract_ratios(file)
  df <- read.csv(file)
  
  df %>%
    mutate(
      wall_ratio = rr$wall_ratio,
      place_ratio = rr$place_ratio,
      group = rr$group
    )
}

## HELPERS FUNCTION TO PLOT SIGNIFICANCE STAR
p_to_star <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("n.s.")
}

## HELPERS FUNCTION TO RUN TYPE 3 ANOVA TO TEST WHETHER THE MODEL'S SIGNIFICANCE
extract_fixed_effects <- function(model_obj) {
  tab <- as.data.frame(anova(model_obj, type = 3))
  tab$Effect <- rownames(tab)
  rownames(tab) <- NULL
  
  effects_keep <- c("group", "day", "trial", "group:day", "group:trial")
  
  out <- sapply(effects_keep, function(eff) {
    pval <- tab$`Pr(>F)`[tab$Effect == eff]
    if (length(pval) == 0) return("")
    p_to_star(pval[1])
  })
  
  return(out)
}

## HELPERS FUNCTION FOR VISUALIZATION PLOT
plot_pm <- function(df, mean_col, se_col, ylab, title, letter) {
  ymax <- max(df[[mean_col]] + df[[se_col]], na.rm = TRUE)
  ymin <- min(df[[mean_col]] - df[[se_col]], na.rm = TRUE)
  
  ggplot(df, aes(x = day, y = .data[[mean_col]], color = group, group = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_errorbar(
      aes(ymin = .data[[mean_col]] - .data[[se_col]],
          ymax = .data[[mean_col]] + .data[[se_col]]),
      width = 0.15
    ) +
    annotate("text", x = 0.6, y = ymax * 1.08,
             label = letter, fontface = "bold", size = 5) +
    scale_x_continuous(breaks = 1:8) +
    labs(x = "Day", y = ylab, title = title, color = "Group") +
    coord_cartesian(ylim = c(ymin, ymax * 1.12), clip = "off") +
    theme_classic() +
    theme(legend.position = "none")
}

## READ MEAN + SE CSVs FILES
pm_var <- bind_rows(lapply(var_summary_files, read_summary_file))

## CHANGE "TARGET QUADRANT" AND "OPPOSITE QUADRANT" PMs TO "POSSIBLE QUADRANT" FOR VARIABLE PLATFORM CASE
pm_var <- pm_var %>%
  mutate(
    possible = target_quadrant + opposite_quadrant,
    possible_se = sqrt(target_quadrant_se^2 + opposite_quadrant_se^2)
  )

pm_var_day <- pm_var %>%
  group_by(group, place_ratio, wall_ratio, day) %>%
  summarise(
    latency = mean(latency),
    latency_se = mean(latency_se),
    
    distance = mean(dist),
    distance_se = mean(dist_se),
    
    possible = mean(possible),
    possible_se = mean(possible_se),
    
    wall = mean(wall_zone),
    wall_se = mean(wall_zone_se),
    .groups = "drop"
  )

pm_var_day$group <- factor(pm_var_day$group, levels = group_levels)

## PLOTS
pA <- plot_pm(pm_var_day, "latency", "latency_se", "Latency (s)", "Latency", "A")
pB <- plot_pm(pm_var_day, "distance", "distance_se", "Distance", "Distance", "B")
pC <- plot_pm(pm_var_day, "possible", "possible_se",
              "Possible quadrants (%)", "Time in possible quadrants (%)", "C")
pD <- plot_pm(pm_var_day, "wall", "wall_se",
              "Wall zone (%)", "Wall zone time (%)", "D")

plot_var <- (pA + pB) / (pC + pD) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  )

plot_var <- plot_var +
  plot_annotation(
    title = "Variable Platforms PMs",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

print(plot_var)

ggsave("variable_platform_comparison.png", plot_var, width = 12, height = 8)

## READ RUN-LEVEL CSVs FILES
run_var <- bind_rows(lapply(var_run_files, read_run_file)) %>%
  mutate(
    possible = target_quadrant + opposite_quadrant
  )

run_var$group <- factor(run_var$group, levels = group_levels)
run_var$day <- factor(run_var$day)
run_var$trial <- factor(run_var$trial)
run_var$run <- factor(run_var$run)

run_var$run_id <- interaction(run_var$group, run_var$run, drop = TRUE)


## MIXED MODELS
model_latency_var <- lmer(latency ~ group + day + trial + group:day + group:trial + (1 | run_id), data = run_var)

model_distance_var <- lmer(distance ~ group + day + trial + group:day + group:trial + (1 | run_id), data = run_var)

model_possible_var <- lmer(possible ~ group + day + trial + group:day + group:trial + (1 | run_id), data = run_var)

model_wall_var <- lmer(wall_zone ~ group + day + trial + group:day + group:trial + (1 | run_id), data = run_var)

## SIGNIFICANCE TABLE
summary_table_var <- data.frame(
  Effect = c("Group", "Day", "Trial", "Group × Day", "Group × Trial"),
  Latency = extract_fixed_effects(model_latency_var),
  Distance = extract_fixed_effects(model_distance_var),
  Possible = extract_fixed_effects(model_possible_var),
  Wall = extract_fixed_effects(model_wall_var)
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
  c(1.6, 1.0, 1.0, 1.8, 1.2),
  "in"
)

table_grob$heights <- rep(unit(0.42, "in"), length(table_grob$heights))

title_box <- grobTree(
  rectGrob(gp = gpar(fill = "#E6E6E6", col = "grey80")),
  textGrob(
    "Variable platform: mixed-effects model (fixed = group, day, trial; random = run)",
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

ggsave(
  "variable_platform_table.png",
  table_only_var,
  width = 8,
  height = 4,
  dpi = 300
)
