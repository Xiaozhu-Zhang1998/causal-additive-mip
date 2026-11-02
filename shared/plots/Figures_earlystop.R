library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

df <- read_csv("./rawdata_earlystop/RS.csv", show_col_types=FALSE) %>%
  mutate(tau_f=factor(tau, levels=sort(unique(tau)))) %>%
  filter(tau <= 0.1)

theme_es <- theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

p_shd <- ggplot(df, aes(x=tau_f, y=shd)) +
  geom_boxplot(width=.65, outlier.size=0) +
  labs(
    x=expression("Early stopping threshold " * tau^early),
    y="SHD"
  ) +
  theme_es

p_time <- ggplot(df, aes(x=tau_f, y=time)) +
  geom_boxplot(width=.65, outlier.size=0) +
  labs(
    x=expression("Early stopping threshold " * tau^early),
    y="Running time (seconds)"
  ) +
  theme_es

fig <- p_shd | p_time
fig
