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

## FILES (FIXED VERSION)
fixed_summary_files <- c(
  "combined_place1p0_wall0p0_fixed.csv",
  "combined_place0p8_wall0p2_fixed.csv",
  "combined_place0p6_wall0p4_fixed.csv",
  "combined_place0p5_wall0p5_fixed.csv",
  "combined_place0p4_wall0p6_fixed.csv",
  "combined_place0p2_wall0p8_fixed.csv",
  "combined_place0p0_wall1p0_fixed.csv"
)

fixed_run_files <- c(
  "combined_place1p0_wall0p0_fixed_run_level.csv",
  "combined_place0p8_wall0p2_fixed_run_level.csv",
  "combined_place0p6_wall0p4_fixed_run_level.csv",
  "combined_place0p5_wall0p5_fixed_run_level.csv",
  "combined_place0p4_wall0p6_fixed_run_level.csv",
  "combined_place0p2_wall0p8_fixed_run_level.csv",
  "combined_place0p0_wall1p0_fixed_run_level.csv"
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
plot_pm_fixed <- function(df, mean_col, se_col, ylab, title, letter) {
  ymax <- max(df[[mean_col]] + df[[se_col]], na.rm = TRUE)
  ymin <- min(df[[mean_col]] - df[[se_col]], na.rm = TRUE)
  
  ggplot(df, aes(x = day, y = .data[[mean_col]], color = group, group = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_errorbar(
      aes(
        ymin = .data[[mean_col]] - .data[[se_col]],
        ymax = .data[[mean_col]] + .data[[se_col]]
      ),
      width = 0.15,
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 0.6,
      y = ymax * 1.08,
      label = letter,
      fontface = "bold",
      size = 5,
      hjust = 0
    ) +
    scale_x_continuous(breaks = 1:8) +
    labs(x = "Day", y = ylab, title = title, color = "Group") +
    coord_cartesian(ylim = c(ymin, ymax * 1.12), clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13),
      plot.margin = margin(15, 10, 10, 15)
    )
}


## READ MEAN + SE CSVs FILES
pm_fixed <- bind_rows(lapply(fixed_summary_files, read_summary_file))

pm_fixed_day <- pm_fixed %>%
  group_by(group, place_ratio, wall_ratio, day) %>%
  summarise(
    latency = mean(latency, na.rm = TRUE),
    latency_se = mean(latency_se, na.rm = TRUE),
    
    distance = mean(dist, na.rm = TRUE),
    distance_se = mean(dist_se, na.rm = TRUE),
    
    target = mean(target_quadrant, na.rm = TRUE),
    target_se = mean(target_quadrant_se, na.rm = TRUE),
    
    opposite = mean(opposite_quadrant, na.rm = TRUE),
    opposite_se = mean(opposite_quadrant_se, na.rm = TRUE),
    
    wall = mean(wall_zone, na.rm = TRUE),
    wall_se = mean(wall_zone_se, na.rm = TRUE),
    
    .groups = "drop"
  )

pm_fixed_day$group <- factor(pm_fixed_day$group, levels = group_levels)

## PLOTS 
pA <- plot_pm_fixed(
  pm_fixed_day, "latency", "latency_se",
  "Latency (s)", "Latency", "A"
)

pB <- plot_pm_fixed(
  pm_fixed_day, "distance", "distance_se",
  "Distance travelled", "Distance", "B"
)

pC <- plot_pm_fixed(
  pm_fixed_day, "target", "target_se",
  "Target quadrant (%)", "Target quadrant time (%)", "C"
)

pD <- plot_pm_fixed(
  pm_fixed_day, "opposite", "opposite_se",
  "Opposite quadrant (%)", "Opposite quadrant time (%)", "D"
)

pE <- plot_pm_fixed(
  pm_fixed_day, "wall", "wall_se",
  "Wall zone (%)", "Wall zone time (%)", "E"
)

plot_fixed <- wrap_plots(
  pA, pB, pC, pD, pE,
  ncol = 3,
  guides = "collect"
) &
  theme(legend.position = "right")

plot_fixed <- plot_fixed +
  plot_annotation(
    title = "Fixed Platform PMs",
    theme = theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5)
    )
  )

print(plot_fixed)

ggsave(
  "fixed_platform_comparison.png",
  plot_fixed,
  width = 16,
  height = 9,
  dpi = 300
)


## READ RUN-LEVEL CSVs FILES
run_fixed <- bind_rows(lapply(fixed_run_files, read_run_file))

run_fixed$group <- factor(run_fixed$group, levels = group_levels)
run_fixed$day <- factor(run_fixed$day)
run_fixed$trial <- factor(run_fixed$trial)
run_fixed$run <- factor(run_fixed$run)

if ("pset" %in% names(run_fixed)) {
  run_fixed$pset <- factor(run_fixed$pset)
  run_fixed$run_id <- interaction(run_fixed$group, run_fixed$pset, run_fixed$run, drop = TRUE)
} else {
  run_fixed$run_id <- interaction(run_fixed$group, run_fixed$run, drop = TRUE)
}


## MIXED MODELS
model_latency_fixed <- lmer(
  latency ~ group + day + trial + group:day + group:trial + (1 | run_id),
  data = run_fixed
)

model_distance_fixed <- lmer(
  distance ~ group + day + trial + group:day + group:trial + (1 | run_id),
  data = run_fixed
)

model_target_fixed <- lmer(
  target_quadrant ~ group + day + trial + group:day + group:trial + (1 | run_id),
  data = run_fixed
)

model_opposite_fixed <- lmer(
  opposite_quadrant ~ group + day + trial + group:day + group:trial + (1 | run_id),
  data = run_fixed
)

model_wall_fixed <- lmer(
  wall_zone ~ group + day + trial + group:day + group:trial + (1 | run_id),
  data = run_fixed
)


## SIGNIFICANCE TABLE
lat_res_fixed <- extract_fixed_effects(model_latency_fixed)
dis_res_fixed <- extract_fixed_effects(model_distance_fixed)
tar_res_fixed <- extract_fixed_effects(model_target_fixed)
opp_res_fixed <- extract_fixed_effects(model_opposite_fixed)
wal_res_fixed <- extract_fixed_effects(model_wall_fixed)

summary_table_fixed <- data.frame(
  Effect = c("Group", "Day", "Trial", "Group × Day", "Group × Trial"),
  Latency = lat_res_fixed,
  Distance = dis_res_fixed,
  Target_Quadrant = tar_res_fixed,
  Opposite_Quadrant = opp_res_fixed,
  Wall_Zone = wal_res_fixed,
  row.names = NULL
)

print(summary_table_fixed)

write.csv(
  summary_table_fixed,
  "fixed_platform_mixed_model_summary.csv",
  row.names = FALSE
)


## BUILD TABLE
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
  c(1.6, 1.0, 1.0, 1.5, 1.7, 1.1),
  "in"
)

table_grob$heights <- rep(unit(0.42, "in"), length(table_grob$heights))

title_box <- grobTree(
  rectGrob(gp = gpar(fill = "#E6E6E6", col = "grey80")),
  textGrob(
    "Fixed platform: mixed-effects model (fixed = group, day, trial; random = run)",
    gp = gpar(fontsize = 11, fontface = "bold")
  )
)

compact_table_fixed <- arrangeGrob(
  title_box,
  table_grob,
  ncol = 1,
  heights = c(0.18, 1)
)

table_only_fixed <- as.ggplot(compact_table_fixed)

print(table_only_fixed)

ggsave(
  "fixed_platform_table.png",
  table_only_fixed,
  width = 8,
  height = 4,
  dpi = 300
)

