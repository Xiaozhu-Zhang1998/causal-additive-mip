library(dplyr)
library(ggplot2)
library(patchwork)
method_order <- c("GSGES","RESIT","NoTears","GLOBE","TOPIC", "EqVar","NPVAR","CAM","MIP")

metrics_file <- file.path(result_root,'Expression/influence_scores.csv')
if(!file.exists(metrics_file)) metrics_file <- file.path(result_root,'Expression/subsamples/RS.csv')
RS <- readr::read_csv(metrics_file,show_col_types=FALSE) %>%
  mutate(
    method = recode(method,
                    gsges="GSGES", resit="RESIT", notears="NoTears",
                    globe="GLOBE", topic="TOPIC", eqvar="EqVar",
                    npvar="NPVAR", cam="CAM", mip="MIP"
    ),
    method = factor(method, levels = method_order)
  )

p_overall <- ggplot(RS, aes(method, net_influence_score)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8) +
  labs(x = NULL, y = "Net influence score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cond <- ggplot(RS, aes(method, influence_score)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8, na.rm = TRUE) +
  labs(x = NULL, y = "Influence score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig_microarray <- p_overall | p_cond
