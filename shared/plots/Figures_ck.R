library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

methods <- c("gsges", "resit", "notears", "globe", "topic", "eqvar", "npvar", "cam", "mip")

method_labels <- c(
  "gsges" = "GSGES",
  "resit" = "RESIT",
  "notears" = "NoTears",
  "globe" = "GLOBE",
  "topic" = "TOPIC",
  "eqvar" = "EqVar",
  "npvar" = "NPVAR",
  "cam" = "CAM",
  "mip" = "MIP"
)

metric_labels <- c(
  "shd" = "SHD",
  "sid" = "SID",
  "f1" = "Orientation F1"
)

panel_colors <- c(
  "#F8766D",
  "#7CAE00",
  "#00BFC4",
  "#C77CFF"
)


# ============================================================
# Figure 2: n and tau
# Fix p = 15, edge density = 0.30
# ============================================================

conditions <- tibble(
  c = c("c2", "c1", "c4", "c3"),
  condition = c(
    "italic(n) == 300 * ',' ~ tau == 0.1",
    "italic(n) == 500 * ',' ~ tau == 0.1",
    "italic(n) == 300 * ',' ~ tau == 0.2",
    "italic(n) == 500 * ',' ~ tau == 0.2"
  )
)

metric_labels <- c(
  "shd" = "SHD",
  "sid" = "SID",
  "f1" = "F1"
)

df1 <- map2_dfr(
  conditions$c,
  conditions$condition,
  function(ck, label) {
    read_csv(
      file.path(
        paste0("rawdata_", ck),
        paste0("rawdata_p_15_030_", ck),
        "RS.csv"
      ),
      show_col_types = FALSE
    ) %>%
      mutate(condition = label)
  }
) %>%
  filter(method %in% methods) %>%
  mutate(
    method = factor(method, levels = methods),
    condition = factor(condition, levels = conditions$condition)
  ) %>%
  pivot_longer(
    cols = c(shd, sid, f1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("shd", "sid", "f1")))


fig1 <- ggplot(
  df1,
  aes(x = method, y = value, color = condition)
) +
  geom_boxplot(width = 0.7) +
  facet_grid(
    metric ~ condition,
    scales = "free_y",
    labeller = labeller(
      metric = as_labeller(metric_labels),
      condition = label_parsed
    )
  ) +
  scale_color_manual(
    values = setNames(panel_colors, conditions$condition),
    guide = "none"
  ) +
  scale_x_discrete(labels = method_labels) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 10)
  )

fig1



# ============================================================
# Figure 3: graph structure / complexity
# Fix c3: n = 500, tau = 0.2
# ============================================================

graph_settings <- tibble(
  folder = c(
    "rawdata_p_10_050_c3",
    "rawdata_p_15_030_c3",
    "rawdata_p_20_020_c3",
    "rawdata_p_25_020_c3"
  ),
  setting = c(
    "italic(p) == 10 * ',' ~ rho == 0.5",
    "italic(p) == 15 * ',' ~ rho == 0.3",
    "italic(p) == 20 * ',' ~ rho == 0.2",
    "italic(p) == 25 * ',' ~ rho == 0.2"
  )
)

df2 <- map2_dfr(
  graph_settings$folder,
  graph_settings$setting,
  function(folder, label) {
    read_csv(
      file.path("rawdata_c3", folder, "RS.csv"),
      show_col_types = FALSE
    ) %>%
      mutate(setting = label)
  }
) %>%
  filter(method %in% methods) %>%
  mutate(
    method = factor(method, levels = methods),
    setting = factor(setting, levels = graph_settings$setting)
  ) %>%
  pivot_longer(
    cols = c(shd, sid, f1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("shd", "sid", "f1")))


library(ggh4x)

fig2 <- ggplot(
  df2,
  aes(x = method, y = value, color = setting)
) +
  geom_boxplot(width = 0.7) +
  ggh4x::facet_grid2(
    metric ~ setting,
    scales = "free_y",
    independent = "y",
    labeller = labeller(
      metric = as_labeller(metric_labels),
      setting = label_parsed
    )
  ) +
  scale_color_manual(
    values = setNames(panel_colors, graph_settings$setting),
    guide = "none"
  ) +
  scale_x_discrete(labels = method_labels) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 10)
  )

fig2
