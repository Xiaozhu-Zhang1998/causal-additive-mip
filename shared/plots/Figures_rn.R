library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# ----------------------------------------------------------
# Read results
# ----------------------------------------------------------

df <- read_csv("./rawdata_rn/RS.csv", show_col_types=FALSE)

n <- 500
zeta_levels <- sort(unique(df$zeta))

# R_n = int(n^zeta)
rn_df <- tibble(
  zeta = zeta_levels,
  Rn = as.integer(n^zeta_levels)
)

# ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

zeta_labels <- ifelse(
  zeta_levels == 0,
  "Linear",
  format(zeta_levels, trim=TRUE)
)

df <- df %>%
  mutate(zeta_f=factor(zeta, levels=zeta_levels, labels=zeta_labels))

rn_df <- rn_df %>%
  mutate(zeta_f=factor(zeta, levels=zeta_levels, labels=zeta_labels))

print(rn_df)

# ----------------------------------------------------------
# Plot function
# ----------------------------------------------------------

make_plot <- function(metric, ylab) {
  d <- df %>% select(zeta_f, value=all_of(metric))
  
  y_min <- min(d$value, na.rm=TRUE)
  y_max <- max(d$value, na.rm=TRUE)
  r_min <- min(rn_df$Rn)
  r_max <- max(rn_df$Rn)
  
  # Map R_n to approximately 10%-90% of the left-axis range
  a <- 0.8 * (y_max - y_min) / (r_max - r_min)
  b <- y_min + 0.1 * (y_max - y_min) - a * r_min
  
  ggplot(d, aes(x=zeta_f, y=value)) +
    geom_boxplot(width=.65, outlier.size=0.1) +
    geom_point(
      data=rn_df,
      aes(x=zeta_f, y=a * Rn + b),
      inherit.aes=FALSE,
      shape=18, size=3, color="red"
    ) +
    scale_y_continuous(
      name=ylab,
      sec.axis=sec_axis(
        ~ (. - b) / a,
        breaks=seq(0, 30, 5),
        name=expression("Number of basis " * R[n])
      )
    ) +
    labs(x=expression(zeta)) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 10),
      plot.title = element_text(size = 10)
    )
}

# ----------------------------------------------------------
# Three panels
# ----------------------------------------------------------
p_shd  <- make_plot("shd", "SHD")
p_sid  <- make_plot("sid", "SID")
p_f1   <- make_plot("f1", "F1")

fig <- p_sid | p_f1 
fig

fig
