library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# ----------------------------------------------------------
# Read results from all t-distribution folders
# ----------------------------------------------------------

base_dir <- "./rawdata_tdisn"

folders <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

folders <- folders[
  str_detect(
    basename(folders),
    "^rawdata_p_10_\\d+$"
  )
]

df_all <- map_dfr(
  folders,
  ~ read_csv(
    file.path(.x, "RS.csv"),
    show_col_types = FALSE
  )
)

# ----------------------------------------------------------
# Select representative degrees of freedom
# ----------------------------------------------------------

df_keep <- c(10, 15, 20, 25)

# ----------------------------------------------------------
# Method order and labels
# ----------------------------------------------------------

method_order <- c(
  "gsges", "resit", "notears", "globe", "topic",
  "eqvar", "npvar", "cam", "mip"
)

method_labels <- c(
  "gsges"   = "GSGES",
  "resit"   = "RESIT",
  "notears" = "NoTears",
  "globe"   = "GLOBE",
  "topic"   = "TOPIC",
  "eqvar"   = "EqVar",
  "npvar"   = "NPVAR",
  "cam"     = "CAM",
  "mip"     = "MIP"
)

metric_labels <- c(
  "shd" = "SHD",
  "sid" = "SID",
  "f1"  = "F1"
)

# Same four panel colors as the earlier simulation figure
panel_colors <- c(
  "#F8766D",
  "#7CAE00",
  "#00BFC4",
  "#C77CFF"
)

# ----------------------------------------------------------
# Prepare plotting data
# ----------------------------------------------------------

plot_df <- df_all %>%
  filter(df %in% df_keep) %>%
  mutate(
    method = factor(
      method,
      levels = method_order,
      labels = method_labels[method_order]
    ),
    df = factor(
      df,
      levels = df_keep
    )
  ) %>%
  pivot_longer(
    cols = c(shd, sid, f1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("shd", "sid", "f1")
    )
  )

# ----------------------------------------------------------
# Facet labels
# ----------------------------------------------------------

df_labels <- c(
  "10"  = "italic(t)[10]",
  "15" = "italic(t)[15]",
  "20" = "italic(t)[20]",
  "25" = "italic(t)[25]"
)

# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------

p <- ggplot(
  plot_df,
  aes(
    x = method,
    y = value,
    color = df
  )
) +
  geom_boxplot(
    width = 0.7,
  ) +
  facet_grid(
    metric ~ df,
    scales = "free_y",
    axes = "all_x",
    axis.labels = "all_x",
    labeller = labeller(
      metric = as_labeller(metric_labels),
      df = as_labeller(df_labels, label_parsed)
    )
  ) +
  scale_color_manual(
    values = panel_colors
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 8),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9),
    legend.position = "none"
  )

p

# ----------------------------------------------------------
# Save
# ----------------------------------------------------------

