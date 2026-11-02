library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggh4x)

methods <- c(
  "gsges", "resit", "notears", "globe", "topic",
  "eqvar", "npvar", "cam", "mip", "mip_rad"
)

method_labels <- c(
  "gsges"="GSGES",
  "resit"="RESIT",
  "notears"="NoTears",
  "globe"="GLOBE",
  "topic"="TOPIC",
  "eqvar"="EqVar",
  "npvar"="NPVAR",
  "cam"="CAM",
  "mip"="MIP Spl",
  "mip_rad"="MIP Rad"
)

metric_labels <- c(
  "shd"="SHD",
  "sid"="SID",
  "f1"="F1"
)

panel_colors <- c(
  "#F8766D",
  "#7CAE00",
  "#00BFC4",
  "#C77CFF"
)

graph_settings <- tibble(
  folder = c(
    "rawdata_p_10_030",
    "rawdata_p_10_040",
    "rawdata_p_15_015",
    "rawdata_p_15_030"
  ),
  setting = c(
    "italic(p) == 10 * ',' ~ rho == 0.30",
    "italic(p) == 10 * ',' ~ rho == 0.40",
    "italic(p) == 15 * ',' ~ rho == 0.15",
    "italic(p) == 15 * ',' ~ rho == 0.30"
  )
)

df <- map2_dfr(
  graph_settings$folder,
  graph_settings$setting,
  function(folder, label) {
    read_csv(
      file.path("rawdata_newbasis", folder, "RS.csv"),
      show_col_types=FALSE
    ) %>%
      mutate(setting=label)
  }
) %>%
  filter(method %in% methods) %>%
  mutate(
    method=factor(method, levels=methods),
    setting=factor(setting, levels=graph_settings$setting)
  ) %>%
  pivot_longer(
    cols=c(shd, sid, f1),
    names_to="metric",
    values_to="value"
  ) %>%
  mutate(metric=factor(metric, levels=c("shd", "sid", "f1")))

fig_newbasis <- ggplot(
  df,
  aes(x=method, y=value, color=setting)
) +
  geom_boxplot(width=.7) +
  ggh4x::facet_grid2(
    metric ~ setting,
    scales="free_y",
    independent="y",
    labeller=labeller(
      metric=as_labeller(metric_labels),
      setting=label_parsed
    )
  ) +
  scale_color_manual(
    values=setNames(panel_colors, graph_settings$setting),
    guide="none"
  ) +
  scale_x_discrete(labels=method_labels) +
  labs(x="", y="") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 10)
  )

fig_newbasis
