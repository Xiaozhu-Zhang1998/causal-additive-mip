library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

df <- read_csv("rawdata_gamma/RS.csv", show_col_types = FALSE) %>%
  mutate(gamma = factor(gamma, levels = sort(unique(gamma))))

theme_gamma <- theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

p_shd <- ggplot(df, aes(x = gamma, y = shd)) +
  geom_boxplot(width = 0.7) +
  labs(
    x = expression("Tuning parameter " * gamma[n]^{-1}),
    y = "SHD"
  ) +
  theme_gamma

p_sid <- ggplot(df, aes(x = gamma, y = sid)) +
  geom_boxplot(width = 0.7) +
  labs(
    x = expression("Tuning parameter " * gamma[n]^{-1}),
    y = "SID"
  ) +
  theme_gamma

p_f1 <- ggplot(df, aes(x = gamma, y = f1)) +
  geom_boxplot(width = 0.7) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = expression("Tuning parameter " * gamma[n]^{-1}),
    y = "F1"
  ) +
  theme_gamma

p_time <- ggplot(df, aes(x = gamma, y = time/8)) +
  geom_boxplot(width = 0.7) +
  scale_y_log10() +
  labs(
    x = expression("Tuning parameter " * gamma[n]^{-1}),
    y = "Runtime (secs, log scale)"
  ) +
  theme_gamma

p <- (p_shd | p_sid) /
  (p_f1  | p_time)

p
