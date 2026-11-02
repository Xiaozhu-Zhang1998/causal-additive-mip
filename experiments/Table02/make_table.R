# Table02: numerical summaries and LaTeX export from saved RS files.
cmd<-commandArgs(FALSE);this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
ROOT<-normalizePath(dirname(this));while(!file.exists(file.path(ROOT,'shared/settings.py'))) ROOT<-dirname(ROOT)
OUTPUT_DIR<-file.path(ROOT,'outputs/tables');dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
SELECTED_CK<-character()
setwd(Sys.getenv('JASA_RESULTS_ROOT',file.path(ROOT,'results')))
source(file.path(ROOT,'shared/plots/Figures_vanillacomp.R'),local=TRUE)
write.csv(summary_df,file.path(OUTPUT_DIR,'Table02.csv'),row.names=FALSE)
cols<-c('Dataset','Objective difference','SHD difference','CP time','Vanilla time')
lines<-c(paste(cols,collapse=' & '),apply(as.data.frame(summary_df[,c('dataset','objval_diff','shd_diff','cp_time','vanilla_time')]),1,function(x) paste(x,collapse=' & ')))
writeLines(paste0(lines,' \\\\'),file.path(OUTPUT_DIR,'Table02_rows.tex'))

cat('Wrote tables to',OUTPUT_DIR,'\n')

obj<-ifelse(abs(summary_df$objval_diff)<1e-4,'$<10^{-4}$',sprintf('%.3f',summary_df$objval_diff))
shd<-sprintf('%.1f',ifelse(abs(summary_df$shd_diff)<.05,0,summary_df$shd_diff))
body<-sprintf('%s & %s & %s & %.3f & %.3f',summary_df$dataset,obj,shd,summary_df$cp_time,summary_df$vanilla_time)
tex<-c('\\begin{table}','\\centering','\\begin{tabular}{lrrrr}','\\hline',
       'Dataset & ObjVal diff & SHD diff & CP time & Vanilla time',body,'\\hline','\\end{tabular}',
       '\\caption{Cutting-plane versus vanilla MIP; differences are vanilla minus cutting-plane. Times are in seconds.}',
       '\\end{table}')
tex[5:(5+length(body))]<-paste0(tex[5:(5+length(body))],' \\\\')
writeLines(tex,file.path(OUTPUT_DIR,'Table02.tex'))
