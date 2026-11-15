# Render Table 4 from saved interaction results at both noise levels.
cmd <- commandArgs(FALSE)
this <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
root <- normalizePath(file.path(dirname(this), "..", ".."))
source(file.path(root, "shared", "plots", "Tables_interaction.R"))
write_interaction_table(root, c(0.1, 0.2), "Table04", metric = "sid")
