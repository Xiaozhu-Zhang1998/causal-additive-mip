# Summarize directed SHD or SID over the 20 trials in each interaction setting.
write_interaction_table <- function(root, noise_levels, table_name, metric = "shd") {
  metric <- match.arg(metric, c("shd", "sid"))
  if (metric == "sid" && !requireNamespace("SID", quietly = TRUE)) {
    stop("Install the R package SID to render Table 4.")
  }

  saved <- file.path(root, "results", "interaction")
  data <- file.path(root, "data", "simulated", "interaction")
  truth <- subset(read.csv(file.path(data, "ground_truth.csv")), kind == "edge")
  mip <- subset(read.csv(file.path(saved, "mip_summary.csv")), bic == "edf")
  methods <- c("MIP-main", "MIP-interaction", "GLOBE", "TOPIC", "NOTEARS", "GSGES")
  baseline_names <- tolower(methods[3:6])
  labels <- c(
    paste0("One interaction, $\\eta=", c(0.5, 1, 2), "$"),
    paste0("Three interactions, $\\eta=", c(0.5, 1, 2), "$"),
    paste0("Linear mains, three interactions, $\\eta=", c(0.5, 2), "$"),
    "Additive", "Pure interaction", "Three parents, $b=0.2$",
    "Three parents, $b=0.5$", "Cosine mains, three interactions, $\\eta=1$"
  )
  mean_col <- paste0("mean_", metric)
  se_col <- paste0("se_", metric)

  rows <- list()
  for (tau in noise_levels) {
    for (setting in 0:12) {
      scores <- matrix(NA_real_, 20, length(baseline_names))
      folder <- file.path(saved, paste0("tau_", tau), sprintf("setting_%02d", setting))
      for (trial in 0:19) {
        edges <- truth[truth$noise_tau == tau & truth$setting == setting & truth$run == trial, ]
        if (!nrow(edges)) stop("Missing true edges for tau ", tau, ", setting ", setting, ", trial ", trial)
        actual <- matrix(0, 8, 8)
        actual[cbind(edges$parent1 + 1L, edges$child + 1L)] <- 1
        for (i in seq_along(baseline_names)) {
          path <- file.path(folder, sprintf("dag_%s_%d.csv", baseline_names[i], trial))
          estimate <- as.matrix(read.csv(path, header = FALSE))
          if (!identical(dim(estimate), c(8L, 8L)) || !all(estimate %in% c(0, 1))) {
            stop("Invalid adjacency matrix: ", path)
          }
          scores[trial + 1L, i] <- if (metric == "shd") {
            sum(actual != estimate)
          } else {
            SID::structIntervDist(actual, estimate)$sid
          }
        }
      }
      m <- mip[mip$noise_tau == tau & mip$setting == setting, ]
      m <- m[match(c("MIP-main", "MIP-interactions"), m$method), ]
      stopifnot(all(m$completed_runs == 20))
      row <- data.frame(
        setting = setting, noise_tau = tau, method = methods,
        mean = c(m[[mean_col]], colMeans(scores)),
        se = c(m[[se_col]], apply(scores, 2, sd) / sqrt(nrow(scores)))
      )
      names(row)[4:5] <- c(mean_col, se_col)
      rows[[length(rows) + 1L]] <- row
    }
  }

  out <- file.path(root, "outputs", "tables")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  result <- do.call(rbind, rows)
  write.csv(result, file.path(out, paste0(table_name, ".csv")), row.names = FALSE)
  end_row <- " \\\\"
  formatted <- unlist(lapply(noise_levels, function(tau) {
    header <- if (length(noise_levels) > 1L) {
      paste0(sprintf("\\multicolumn{7}{l}{$\\tau=%.1f$}", tau), end_row)
    }
    block <- vapply(0:12, function(setting) {
      row <- result[result$noise_tau == tau & result$setting == setting, ]
      cells <- sprintf("%.2f (%.2f)", row[[mean_col]], row[[se_col]])
      paste0(paste(c(labels[setting + 1L], cells), collapse = " & "), end_row)
    }, character(1))
    c(header, block)
  }), use.names = FALSE)
  metric_label <- if (metric == "shd") "directed SHD" else "SID"
  tau_label <- paste(sprintf("%.1f", noise_levels), collapse = ", ")
  layout <- if (metric == "sid") {
    c("\\footnotesize", "\\setlength{\\tabcolsep}{2pt}", "\\begin{tabular}{p{0.28\\linewidth}rrrrrr}")
  } else {
    "\\begin{tabular}{lrrrrrr}"
  }
  tex <- c(
    "\\begin{table}", "\\centering", layout, "\\hline",
    paste0(paste(c("Setting", methods), collapse = " & "), end_row), "\\hline",
    formatted, "\\hline", "\\end{tabular}",
    sprintf("\\caption{Mean %s (standard error) over 20 replications with $p=8$, $n=300$, and $\\tau=%s$. Both MIP estimators use $\\mathrm{BIC}_{\\mathrm{EDF}}$ for final selection.}", metric_label, tau_label),
    "\\end{table}"
  )
  writeLines(tex, file.path(out, paste0(table_name, ".tex")))
  cat("Wrote", table_name, "to", out, "\n")
}
