# Plot Figure 8 from the saved numerical results.
args_full<-commandArgs(FALSE)
this<-sub('^--file=','',args_full[grepl('^--file=',args_full)][1])
ROOT<-normalizePath(dirname(this));while(!file.exists(file.path(ROOT,'shared/settings.py'))) ROOT<-dirname(ROOT)
OUTPUT_DIR<-file.path(ROOT,'outputs/figures');dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
result_root<-Sys.getenv('JASA_RESULTS_ROOT',file.path(ROOT,'results'))
setwd(result_root)
source(file.path(ROOT,'shared/plots/Figures_gamma.R'),local=TRUE)
figure<-p
ggplot2::ggsave(file.path(OUTPUT_DIR,'Fig08.pdf'),figure,width=11,height=5)
cat(file.path(OUTPUT_DIR,'Fig08.pdf'),'\n')
