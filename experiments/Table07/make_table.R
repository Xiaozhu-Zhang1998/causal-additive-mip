# Table07: numerical summaries and LaTeX export from saved RS files.
cmd<-commandArgs(FALSE);this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
ROOT<-normalizePath(dirname(this));while(!file.exists(file.path(ROOT,'shared/settings.py'))) ROOT<-dirname(ROOT)
OUTPUT_DIR<-file.path(ROOT,'outputs/tables');dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
SELECTED_CK<-'c3'
setwd(Sys.getenv('JASA_RESULTS_ROOT',file.path(ROOT,'results')))
source(file.path(ROOT,'shared/plots/Tables_ck.R'),local=TRUE)

cat('Wrote tables to',OUTPUT_DIR,'\n')

write.csv(summarize_ck("c3"),file.path(OUTPUT_DIR,"Table07_means.csv"),row.names=FALSE)
