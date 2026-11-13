# Table08: numerical summaries and LaTeX export from saved RS files.
cmd<-commandArgs(FALSE);this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
ROOT<-normalizePath(dirname(this));while(!file.exists(file.path(ROOT,'shared/settings.py'))) ROOT<-dirname(ROOT)
OUTPUT_DIR<-file.path(ROOT,'outputs/tables');dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
SELECTED_CK<-character()
setwd(Sys.getenv('JASA_RESULTS_ROOT',file.path(ROOT,'results')))
source(file.path(ROOT,'shared/plots/Tables_graphs.R'),local=TRUE)

cat('Wrote tables to',OUTPUT_DIR,'\n')

write.csv(summarize_graph_results(),file.path(OUTPUT_DIR,"Table08_means.csv"),row.names=FALSE)
