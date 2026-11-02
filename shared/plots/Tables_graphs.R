library(dplyr)
library(purrr)
library(readr)

graphs <- c(
  "dsep", "asia", "bowling", "rain", "inssmall",
  "cloud", "funnel", "galaxy", "insurance", "factors"
)

graph_labels <- c(
  "dsep" = "Dsep",
  "asia" = "Asia",
  "bowling" = "Bowling",
  "rain" = "Rain",
  "inssmall" = "InsuranceSmall",
  "cloud" = "Cloud",
  "funnel" = "Funnel",
  "galaxy" = "Galaxy",
  "insurance" = "Insurance",
  "factors" = "Factors"
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


# ----------------------------------------------------------
# Read all RS.csv files
# ----------------------------------------------------------

read_graph_results <- function() {
  map_dfr(
    graphs,
    ~ read_csv(
      file.path(paste0("rawdata_graphs/rawdata_graph_", .x), "RS.csv"),
      show_col_types = FALSE
    )
  )
}


# ----------------------------------------------------------
# Compute graph size p and true edge count s*
# ----------------------------------------------------------

get_graph_info <- function(graph) { c(p=unname(c(dsep=6,asia=8,bowling=9,rain=14,inssmall=15,cloud=16,funnel=18,galaxy=20,insurance=27,factors=27)[graph]),s=unname(c(dsep=6,asia=8,bowling=11,rain=18,inssmall=25,cloud=19,funnel=18,galaxy=22,insurance=52,factors=68)[graph])) }


# ----------------------------------------------------------
# Compute means and adjust parallelizable runtimes
# ----------------------------------------------------------

summarize_graph_results <- function() {
  
  read_graph_results() %>%
    mutate(
      time = case_when(
        method == "mip" ~ time / 12,
        method == "notears" ~ time / 3,
        method == "gsges" ~ time / 3,
        method == "bootstrap" ~ time / 20,
        TRUE ~ time
      )
    ) %>%
    group_by(graph, method) %>%
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


get_value <- function(df, graph0, method0, metric) {
  
  x <- df %>%
    filter(
      graph == graph0,
      method == method0
    ) %>%
    pull(all_of(metric))
  
  if (length(x) == 0) return(NA_real_)
  x[1]
}


get_method_values <- function(df, graph0, method) {
  
  c(
    fmt(get_value(df, graph0, method, "shd"), 1),
    fmt(get_value(df, graph0, method, "sid"), 1),
    fmt(get_value(df, graph0, method, "f1"), 2),
    fmt(get_value(df, graph0, method, "time"), 1)
  )
}


# ----------------------------------------------------------
# Generate LaTeX table
# ----------------------------------------------------------

make_graph_latex_table <- function() {
  
  df <- summarize_graph_results()
  
  col_spec <- paste0(
    "l",
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
  # Blocks 1 and 2
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
      
      "\\cmidrule{2-5}\\cmidrule{6-9}\\cmidrule{10-13}",
      
      paste0(
        "Network.$p$.$s^\\star$",
        " & SHD & SID & F1 & Time",
        " & SHD & SID & F1 & Time",
        " & SHD & SID & F1 & Time & \\\\"
      ),
      
      "\\midrule"
    )
    
    for (g in graphs) {
      
      info <- get_graph_info(g)
      
      row_name <- sprintf(
        "%s.%d.%d",
        graph_labels[g],
        as.integer(info["p"]),
        as.integer(info["s"])
      )
      
      v1 <- get_method_values(df, g, m1)
      v2 <- get_method_values(df, g, m2)
      v3 <- get_method_values(df, g, m3)
      
      out <- c(
        out,
        
        sprintf(
          paste0(
            "%s",
            " & %s & %s & %s & %s",
            " & %s & %s & %s & %s",
            " & %s & %s & %s & %s & \\\\"
          ),
          row_name,
          v1[1], v1[2], v1[3], v1[4],
          v2[1], v2[2], v2[3], v2[4],
          v3[1], v3[2], v3[3], v3[4]
        )
      )
    }
    
    out <- c(out, "\\midrule")
  }
  
  
  # ========================================================
  # Block 3: NPVAR | CAM | MIP | Bootstrap
  # ========================================================
  
  out <- c(
    out,
    
    paste0(
      "& \\multicolumn{4}{c}{NPVAR}",
      " & \\multicolumn{4}{c}{CAM}",
      " & \\multicolumn{4}{c}{MIP}",
      " & Bootstrap \\\\"
    ),
    
    paste0(
      "\\cmidrule{2-5}",
      "\\cmidrule{6-9}",
      "\\cmidrule{10-13}",
      "\\cmidrule{14-14}"
    ),
    
    paste0(
      "Network.$p$.$s^\\star$",
      " & SHD & SID & F1 & Time",
      " & SHD & SID & F1 & Time",
      " & SHD & SID & F1 & Time",
      " & Time \\\\"
    ),
    
    "\\midrule"
  )
  
  for (g in graphs) {
    
    info <- get_graph_info(g)
    
    row_name <- sprintf(
      "%s.%d.%d",
      graph_labels[g],
      as.integer(info["p"]),
      as.integer(info["s"])
    )
    
    v1 <- get_method_values(df, g, "npvar")
    v2 <- get_method_values(df, g, "cam")
    v3 <- get_method_values(df, g, "mip")
    
    boot_time <- fmt(
      get_value(df, g, "bootstrap", "time"),
      1
    )
    
    out <- c(
      out,
      
      sprintf(
        paste0(
          "%s",
          " & %s & %s & %s & %s",
          " & %s & %s & %s & %s",
          " & %s & %s & %s & %s",
          " & %s \\\\"
        ),
        row_name,
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
    
    paste0(
      "\\caption{\\small Full simulation results on the benchmark graph structures. ",
      "Each entry reports the mean over 50 independent simulation trials. ",
      "Runtime is reported in seconds; for parallelizable computations in ",
      "\\textit{MIP}, \\textit{NoTears}, \\textit{GSGES}, and bootstrap, the reported runtime ",
      "corresponds to a single run.}"
    ),
    
    "\\label{tab:graph-simulation}",
    
    "\\end{longtable}",
    "\\endgroup"
  )
  
  writeLines(
    out,
    file.path(OUTPUT_DIR,"table_graphs.tex")
  )
  
  invisible(out)
}


# ----------------------------------------------------------
# Generate table
# ----------------------------------------------------------

make_graph_latex_table()
