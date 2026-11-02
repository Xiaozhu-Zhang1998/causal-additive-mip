library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------------
# Method order and colors
# ----------------------------------------------------------

methods <- c(
  "CAM", "NPVAR", "MIP", "NoTears", "GLOBE",
  "TOPIC", "EqVar", "RESIT", "GSGES"
)

colors <- c(
  "CAM"     = "#1F77B4",
  "NPVAR"   = "#FF7F0E",
  "MIP"     = "#2CA02C",
  "NoTears" = "#D62728",
  "GLOBE"   = "#9467BD",
  "TOPIC"   = "#8C564B",
  "EqVar"   = "#E377C2",
  "RESIT"   = "#7F7F7F",
  "GSGES"   = "#BCBD22"
)

# ----------------------------------------------------------
# Read results
# ----------------------------------------------------------

df_results <- read_csv(
  "./causal_chamber/RS.csv",
  show_col_types = FALSE
) %>%
  mutate(
    method = recode(
      method,
      "Globe" = "GLOBE"
    ),
    method = factor(
      method,
      levels = methods
    )
  )

# ----------------------------------------------------------
# Aggregate over the 10 trials
# ----------------------------------------------------------

se <- function(x) {
  x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}

agg_results <- df_results %>%
  group_by(n, method) %>%
  summarise(
    mean_shd = mean(shd, na.rm = TRUE),
    se_shd   = se(shd),

    mean_sid = mean(sid, na.rm = TRUE),
    se_sid   = se(sid),

    mean_f1 = mean(directed_f1, na.rm = TRUE),
    se_f1   = se(directed_f1),

    mean_tdr = mean(tdr, na.rm = TRUE),
    se_tdr   = se(tdr),

    mean_fdr = mean(fdr, na.rm = TRUE),
    se_fdr   = se(fdr),

    .groups = "drop"
  ) %>%
  mutate(
    line_width = 0.8,
    point_size = 1.5
  )

# ----------------------------------------------------------
# Common theme
# ----------------------------------------------------------

theme_cc <- theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# ----------------------------------------------------------
# Function for ordinary metric panels
# ----------------------------------------------------------

make_metric_plot <- function(
  mean_var,
  se_var,
  ylab,
  bounded = FALSE
) {

  p <- ggplot(
    agg_results,
    aes(
      x = n,
      y = .data[[mean_var]],
      color = method,
      group = method
    )
  ) +
    geom_ribbon(
      aes(
        ymin = .data[[mean_var]] - .data[[se_var]],
        ymax = .data[[mean_var]] + .data[[se_var]],
        fill = method
      ),
      alpha = 0.10,
      color = NA
    ) +
    geom_line(
      aes(linewidth = line_width)
    ) +
    geom_point(
      aes(size = point_size)
    ) +
    scale_color_manual(
      values = colors
    ) +
    scale_fill_manual(
      values = colors,
      guide = "none"
    ) +
    scale_linewidth_identity() +
    scale_size_identity() +
    scale_x_continuous(
      breaks = sort(unique(df_results$n))
    ) +
    labs(
      x = expression(
        paste("Sample size ", italic(n))
      ),
      y = ylab
    ) +
    theme_cc

  if (bounded) {
    p <- p +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2)
      )
  }

  p
}

# ----------------------------------------------------------
# Five ordinary panels
# ----------------------------------------------------------

p_shd <- make_metric_plot(
  "mean_shd",
  "se_shd",
  "SHD"
)

p_sid <- make_metric_plot(
  "mean_sid",
  "se_sid",
  "SID"
)

p_f1 <- make_metric_plot(
  "mean_f1",
  "se_f1",
  "F1 score",
  bounded = TRUE
)

p_tdr <- make_metric_plot(
  "mean_tdr",
  "se_tdr",
  "TDR",
  bounded = TRUE
)

p_fdr <- make_metric_plot(
  "mean_fdr",
  "se_fdr",
  "FDR",
  bounded = TRUE
)

# ----------------------------------------------------------
# FDR-TDR trajectory
# ----------------------------------------------------------

p_traj <- ggplot(
  agg_results,
  aes(
    x = mean_fdr,
    y = mean_tdr,
    color = method,
    group = method
  )
) +
  geom_line(
    aes(linewidth = line_width)
  ) +
  geom_point(
    aes(size = point_size)
  ) +
  scale_color_manual(
    values = colors
  ) +
  scale_linewidth_identity() +
  scale_size_identity() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "FDR",
    y = "TDR"
  ) +
  theme_cc

# ----------------------------------------------------------
# Combine 2 x 3
# ----------------------------------------------------------

p <- (
  p_shd | p_sid | p_f1
) / (
  p_tdr | p_fdr | p_traj
) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

p


(p_shd | p_sid | p_traj) +
  plot_layout(guides = "collect") &
  guides(
    fill = "none",
    color = guide_legend(nrow = 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
