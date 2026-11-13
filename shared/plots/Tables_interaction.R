# Mean directed SHD and standard errors for the 20 trials in each setting.
write_interaction_table <- function(root, tau, table_name) {
  saved <- file.path(root, "results", "interaction")
  data <- file.path(root, "data", "simulated", "interaction")
  truth <- read.csv(file.path(data, "ground_truth.csv"))
  truth <- subset(truth, kind == "edge" & noise_tau == tau)
  mip <- read.csv(file.path(saved, "mip_summary.csv"))
  mip <- subset(mip, bic == "edf" & noise_tau == tau)
  methods <- c("MIP-main", "MIP-interaction", "GLOBE", "TOPIC", "NOTEARS", "GSGES")
  baseline_names <- tolower(methods[3:6])
  labels <- c(
    paste0("One interaction, $\\eta=", c(0.5, 1, 2), "$"),
    paste0("Three interactions, $\\eta=", c(0.5, 1, 2), "$"),
    paste0("Linear mains, three interactions, $\\eta=", c(0.5, 2), "$"),
    "Additive", "Pure interaction", "Three parents, $b=0.2$",
    "Three parents, $b=0.5$", "Cosine mains, three interactions, $\\eta=1$"
  )

  rows <- list()
  for (setting in 0:12) {
    scores <- matrix(NA_real_, 20, length(baseline_names))
    folder <- file.path(saved, paste0("tau_", tau), sprintf("setting_%02d", setting))
    for (trial in 0:19) {
      edges <- truth[truth$setting == setting & truth$run == trial, ]
      if (!nrow(edges)) stop("Missing true edges for setting ", setting, ", trial ", trial)
      actual <- matrix(0, 8, 8)
      actual[cbind(edges$parent1 + 1L, edges$child + 1L)] <- 1
      for (i in seq_along(baseline_names)) {
        path <- file.path(folder, sprintf("dag_%s_%d.csv", baseline_names[i], trial))
        estimate <- as.matrix(read.csv(path, header = FALSE))
        if (!identical(dim(estimate), c(8L, 8L)) || !all(estimate %in% c(0, 1))) {
          stop("Invalid adjacency matrix: ", path)
        }
        scores[trial + 1L, i] <- sum(actual != estimate)
      }
    }
    m <- mip[mip$setting == setting, ]
    m <- m[match(c("MIP-main", "MIP-interactions"), m$method), ]
    stopifnot(all(m$completed_runs == 20))
    rows[[setting + 1L]] <- data.frame(
      setting = setting, noise_tau = tau, method = methods,
      mean_shd = c(m$mean_shd, colMeans(scores)),
      se_shd = c(m$se_shd, apply(scores, 2, sd) / sqrt(nrow(scores)))
    )
  }

  out <- file.path(root, "outputs", "tables")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(do.call(rbind, rows), file.path(out, paste0(table_name, ".csv")), row.names = FALSE)
  formatted <- vapply(seq_along(rows), function(i) {
    row <- rows[[i]]
    paste(c(labels[i], sprintf("%.2f (%.2f)", row$mean_shd, row$se_shd)), collapse = " & ")
  }, character(1))
  end_row <- " \\\\"
  tex <- c(
    "\\begin{table}", "\\centering", "\\begin{tabular}{lrrrrrr}", "\\hline",
    paste0(paste(c("Setting", methods), collapse = " & "), end_row), "\\hline",
    paste0(formatted, end_row), "\\hline", "\\end{tabular}",
    sprintf("\\caption{Mean directed SHD (standard error) over 20 replications with $p=8$, $n=300$, and $\\tau=%.1f$. Both MIP estimators use $\\mathrm{BIC}_{\\mathrm{EDF}}$ for final selection.}", tau),
    "\\end{table}"
  )
  writeLines(tex, file.path(out, paste0(table_name, ".tex")))
  cat("Wrote", table_name, "to", out, "\n")
}
