library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

graphs <- c("dsep", "asia", "bowling", "inssmall", "insurance")

graph_info <- tibble(
  graph = graphs,
  p = c(6, 8, 9, 15, 27),
  s = c(6, 8, 11, 25, 52)
)

df <- map_dfr(graphs, \(g) {
  read_csv(file.path("rawdata_vanillacomp", g, "RS.csv"), show_col_types=FALSE) %>%
    mutate(graph=g)
})

df_mip <- df %>%
  filter(method %in% c("MIP_OA", "MIP_vanilla")) %>%
  select(graph, l, method, objval, time, shd) %>%
  pivot_wider(names_from=method, values_from=c(objval, time, shd))


summary_df <- df_mip %>%
  mutate(
    objval_diff = objval_MIP_vanilla - objval_MIP_OA,
    shd_diff = shd_MIP_vanilla - shd_MIP_OA 
  ) %>%
  group_by(graph) %>%
  summarise(
    objval_diff = mean(objval_diff, na.rm=TRUE),
    shd_diff = mean(shd_diff, na.rm=TRUE),
    cp_time = mean(time_MIP_OA, na.rm=TRUE),
    vanilla_time = mean(time_MIP_vanilla, na.rm=TRUE),
    .groups="drop"
  ) %>%
  left_join(graph_info, by="graph") %>%
  mutate(dataset=paste0(graph, ".", p, ".", s)) %>%
  select(dataset, graph, objval_diff, shd_diff, cp_time, vanilla_time)

summary_df <- summary_df %>% arrange(match(graph, graphs))
summary_df
