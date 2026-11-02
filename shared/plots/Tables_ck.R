library(dplyr)
library(purrr)
library(readr)
library(stringr)

methods <- c(
  "gsges", "resit", "notears", "globe", "topic",
  "eqvar", "npvar", "cam", "mip"
)

method_labels <- c(
  "gsges" = "GSGES",
  "resit" = "RESIT",
  "notears" = "NoTears",
  "globe" = "Globe",
  "topic" = "TOPIC",
  "eqvar" = "EqVar",
  "npvar" = "NPVAR",
  "cam" = "CAM",
  "mip" = "MIP"
)

blocks <- list(
  c("gsges", "resit", "notears"),
  c("globe", "topic", "eqvar")
)

graph_order <- tibble(
  p = c(10, 10, 10, 15, 15, 20, 20, 25),
  edge_density = c(.30, .40, .50, .15, .30, .15, .20, .20)
)

settings <- list(
  c1 = c(n = 500, tau = 0.1),
  c2 = c(n = 300, tau = 0.1),
  c3 = c(n = 500, tau = 0.2),
  c4 = c(n = 300, tau = 0.2)
)


# ----------------------------------------------------------
# Read all RS.csv files for one c_k
# ----------------------------------------------------------

read_ck <- function(ck) {
  base_dir <- paste0("rawdata_", ck)
  
  folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  folders <- folders[
    str_detect(
      basename(folders),
      paste0("^rawdata_p_\\d+_\\d+_", ck, "$")
    )
  ]
  
  map_dfr(
    folders,
    ~ read_csv(file.path(.x, "RS.csv"), show_col_types = FALSE)
  )
}


# ----------------------------------------------------------
# Compute means and adjust parallelizable runtimes
# ----------------------------------------------------------

summarize_ck <- function(ck) {
  read_ck(ck) %>%
    mutate(
      time = case_when(
        method == "mip" ~ time / if_else(p == 10 & (abs(edge_density-.3)<1e-8 | abs(edge_density-.4)<1e-8), 6, 12),
        method == "notears" ~ time / 3,
        method == "gsges" ~ time / 3,
        method == "bootstrap" ~ time / 20,
        TRUE ~ time
      )
    ) %>%
    group_by(p, edge_density, method) %>%
    summarise(
      shd = mean(shd, na.rm = TRUE),
      sid = mean(sid, na.rm = TRUE),
      f1 = mean(f1, na.rm = TRUE),
      time = mean(time, na.rm = TRUE),
      .groups = "drop"
    )
}


# ----------------------------------------------------------
# Formatting
# ----------------------------------------------------------

fmt <- function(x, digits = 1) {
  if (length(x) == 0 || is.na(x) || is.nan(x)) return("--")
  sprintf(paste0("%.", digits, "f"), x)
}


get_value <- function(df, p0, rho0, method0, metric) {
  x <- df %>%
    filter(
      p == p0,
      abs(edge_density - rho0) < 1e-8,
      method == method0
    ) %>%
    pull(all_of(metric))
  
  if (length(x) == 0) return(NA_real_)
  x[1]
}


get_method_values <- function(df, p0, rho0, method) {
  c(
    fmt(get_value(df, p0, rho0, method, "shd"), 1),
    fmt(get_value(df, p0, rho0, method, "sid"), 1),
    fmt(get_value(df, p0, rho0, method, "f1"), 2),
    fmt(get_value(df, p0, rho0, method, "time"), 1)
  )
}


# ----------------------------------------------------------
# Generate LaTeX table
# ----------------------------------------------------------

make_latex_table <- function(ck) {
  
  df <- summarize_ck(ck)
  
  n0 <- settings[[ck]]["n"]
  tau0 <- settings[[ck]]["tau"]
  
  # Vertical separators are inserted between blocks.
  # Each separator has 15pt space on both sides.
  col_spec <- paste0(
    "c",
    "@{\\hspace{5pt}\\vrule\\hspace{5pt}}",
    "cccc",
    "@{\\hspace{5pt}\\vrule\\hspace{5pt}}",
    "cccc",
    "@{\\hspace{5pt}\\vrule\\hspace{5pt}}",
    "cccc",
    "@{\\hspace{5pt}\\vrule\\hspace{5pt}}",
    "c"
  )
  
  out <- c(
    "\\begingroup",
    "\\small",
    "\\setlength{\\tabcolsep}{2.5pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    paste0("\\begin{longtable}{", col_spec, "}"),
    "\\toprule"
  )
  
  
  # ========================================================
  # Blocks 1 and 2:
  # GSGES | RESIT | NoTears
  # Globe | TOPIC | EqVar
  # ========================================================
  
  for (b in seq_along(blocks)) {
    
    m1 <- blocks[[b]][1]
    m2 <- blocks[[b]][2]
    m3 <- blocks[[b]][3]
    
    out <- c(
      out,
      sprintf(
        paste0(
          "& \\multicolumn{4}{c}{%s}",
          " & \\multicolumn{4}{c}{%s}",
          " & \\multicolumn{4}{c}{%s} & \\\\"
        ),
        method_labels[m1],
        method_labels[m2],
        method_labels[m3]
      ),
      
      # Untrimmed cmidrule so the line spans the full block
      "\\cmidrule{2-5}\\cmidrule{6-9}\\cmidrule{10-13}",
      
      paste0(
        "$(p,\\,\\rho)$",
        " & SHD & SID & F1 & Time",
        " & SHD & SID & F1 & Time",
        " & SHD & SID & F1 & Time & \\\\"
      ),
      
      "\\midrule"
    )
    
    for (i in seq_len(nrow(graph_order))) {
      
      p0 <- graph_order$p[i]
      rho0 <- graph_order$edge_density[i]
      
      v1 <- get_method_values(df, p0, rho0, m1)
      v2 <- get_method_values(df, p0, rho0, m2)
      v3 <- get_method_values(df, p0, rho0, m3)
      
      out <- c(
        out,
        sprintf(
          paste0(
            "$(%d,\\,%.2f)$",
            " & %s & %s & %s & %s",
            " & %s & %s & %s & %s",
            " & %s & %s & %s & %s & \\\\"
          ),
          p0, rho0,
          v1[1], v1[2], v1[3], v1[4],
          v2[1], v2[2], v2[3], v2[4],
          v3[1], v3[2], v3[3], v3[4]
        )
      )
    }
    
    out <- c(out, "\\midrule")
  }
  
  
  # ========================================================
  # Block 3:
  # NPVAR | CAM | MIP | Bootstrap
  # ========================================================
  
  out <- c(
    out,
    
    paste0(
      "& \\multicolumn{4}{c}{NPVAR}",
      " & \\multicolumn{4}{c}{CAM}",
      " & \\multicolumn{4}{c}{MIP}",
      " & Bootstrap \\\\"
    ),
    
    # Again, no l/r trimming
    paste0(
      "\\cmidrule{2-5}",
      "\\cmidrule{6-9}",
      "\\cmidrule{10-13}",
      "\\cmidrule{14-14}"
    ),
    
    paste0(
      "$(p,\\,\\rho)$",
      " & SHD & SID & F1 & Time",
      " & SHD & SID & F1 & Time",
      " & SHD & SID & F1 & Time",
      " & Time \\\\"
    ),
    
    "\\midrule"
  )
  
  for (i in seq_len(nrow(graph_order))) {
    
    p0 <- graph_order$p[i]
    rho0 <- graph_order$edge_density[i]
    
    v1 <- get_method_values(df, p0, rho0, "npvar")
    v2 <- get_method_values(df, p0, rho0, "cam")
    v3 <- get_method_values(df, p0, rho0, "mip")
    
    boot_time <- fmt(
      get_value(df, p0, rho0, "bootstrap", "time"),
      1
    )
    
    out <- c(
      out,
      sprintf(
        paste0(
          "$(%d,\\,%.2f)$",
          " & %s & %s & %s & %s",
          " & %s & %s & %s & %s",
          " & %s & %s & %s & %s",
          " & %s \\\\"
        ),
        p0, rho0,
        v1[1], v1[2], v1[3], v1[4],
        v2[1], v2[2], v2[3], v2[4],
        v3[1], v3[2], v3[3], v3[4],
        boot_time
      )
    )
  }
  
  
  # ========================================================
  # Caption at bottom
  # ========================================================
  
  out <- c(
    out,
    "\\bottomrule",
    
    sprintf(
      paste0(
        "\\caption{\\small Full simulation results for $n=%d$ and $\\tau=%.1f$. ",
        "Each entry reports the mean over 50 independent simulation trials. ",
        "Runtime is reported in seconds. ",
        "For parallelizable computations in \\textit{MIP}, \\textit{NoTears}, \\textit{GSGES} ",
        "and bootstrap, runtime is reported per run.}"
      ),
      n0, tau0
    ),
    
    sprintf("\\label{tab:compare-%s}", ck),
    
    "\\end{longtable}",
    "\\endgroup"
  )
  
  writeLines(out, file.path(OUTPUT_DIR, paste0("table_", ck, ".tex")))
  
  invisible(out)
}


# ----------------------------------------------------------
# Generate c1-c4 tables
# ----------------------------------------------------------

for (ck in SELECTED_CK) {
  make_latex_table(ck)
}
