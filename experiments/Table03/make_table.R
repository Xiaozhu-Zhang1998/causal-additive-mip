# Render Table03 from saved interaction results.
cmd <- commandArgs(FALSE)
this <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
root <- normalizePath(file.path(dirname(this), "..", ".."))
source(file.path(root, "shared", "plots", "Tables_interaction.R"))
write_interaction_table(root, 0.2, "Table03")
