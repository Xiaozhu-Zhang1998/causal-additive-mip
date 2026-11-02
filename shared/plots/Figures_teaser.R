library(dplyr)
library(ggplot2)
library(patchwork)



colnames = c("n", "m", "var", "method", "rec", "SHD", "func")
read_teaser_results <- function(path) {
  tab <- read.csv(path)
  if (!all(colnames %in% names(tab))) stop('Missing required teaser result columns')
  tab[, colnames, drop=FALSE]
}
tab = rbind(
  read_teaser_results("./rawdata_teaser/0.5x^3/RS.csv"),
  read_teaser_results("./rawdata_teaser/atan(x^2)/RS.csv"),
  read_teaser_results("./rawdata_teaser/abs(x)/RS.csv"),
  read_teaser_results("./rawdata_teaser/exp(x)/RS.csv")
) 


p1 = tab %>%
  filter(var == "unequal") %>%
  mutate(rec = case_when(
    rec %in% c("0", "FALSE", "False") ~ 0,
    rec %in% c("1", "TRUE", "True") ~ 1
  )) %>%
  group_by(n, var, method, func) %>% 
  summarise(rec = mean(rec)) %>%
  mutate(var = recode(var,
                      "equal" = "Equal variances",
                      "unequal" = "Unequal variances"),
         method = factor(method,
                         levels = c("CAM", "NPVAR", "linear", "Mips_unequal", "Mips_equal", "Mips"),
                         labels = c("CAM-IncEdge", "NPVAR", "MIP-linear", 
                                    "MIP-nonlinear (our approach)", "MIP-nonlinear (our approach)", "MIP-nonlinear (our approach)")),
         func = factor(func, 
                       levels = c("sin(x)", "atan(x^2)", "exp(x)", "0.5x^3", "abs(x)"),
                       labels = c("h(x)==sin(x)", "h(x)==arctan(x^2)", "h(x)==exp(x)", "h(x)==0.5*x^3", "h(x)==abs(x)"))
  ) %>%
  ggplot(aes(x = n, y = rec, col = method, shape = method)) +
  geom_point(size = 1.2) +
  facet_grid(cols = vars(func), labeller = labeller(func = label_parsed)) +
  scale_shape_manual(
    values = c("CAM-IncEdge" = 7, "NPVAR" = 9, "MIP-linear" = 8, "MIP-nonlinear (our approach)" = 19),
  ) +
  labs(x = "Sample size n", y = "P (correct recovery)", col = "Method", shape = "Method") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10), 
        plot.title = element_text(size = 10),  # Title size
        axis.title.x = element_text(size = 10),  # x-axis label size
        axis.title.y = element_text(size = 10)   # y-axis label size
  )


p2 = tab %>%
  filter(var == "unequal") %>%
  group_by(n, var, method, func) %>%
  summarise(SHD = mean(SHD)) %>%
  mutate(var = recode(var,
                      "equal" = "Equal variances",
                      "unequal" = "Unequal variances"),
         method = factor(method,
                         levels = c("CAM", "NPVAR", "linear", "Mips_unequal", "Mips_equal", "Mips"),
                         labels = c("CAM-IncEdge", "NPVAR", "MIP-linear", 
                                    "MIP-nonlinear (our approach)", "MIP-nonlinear (our approach)", "MIP-nonlinear (our approach)")),
         func = factor(func, 
                       levels = c("sin(x)", "atan(x^2)", "exp(x)", "0.5x^3", "abs(x)"),
                       labels = c("h(x)==sin(x)", "h(x)==arctan(x^2)", "h(x)==exp(x)", "h(x)==0.5*x^3", "h(x)==abs(x)"))
  ) %>%
  ggplot(aes(x = n, y = SHD, col = method, shape = method)) +
  geom_point(size = 1.2) +
  facet_grid(cols = vars(func), labeller = labeller(func = label_parsed)) +
  scale_shape_manual(
    values = c("CAM-IncEdge" = 7, "NPVAR" = 9, "MIP-linear" = 8, "MIP-nonlinear (our approach)" = 19),
  ) +
  labs(x = "Sample size n", y = "SHD", col = "Method", shape = "Method") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10), 
        plot.title = element_text(size = 10),  # Title size
        axis.title.x = element_text(size = 10),  # x-axis label size
        axis.title.y = element_text(size = 10)   # y-axis label size
  )

(p1 / p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
