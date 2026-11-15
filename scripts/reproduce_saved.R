# Render Figures 1-10 and Tables 1-9 from saved numerical results.
cmd<-commandArgs(FALSE);this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
root<-normalizePath(file.path(dirname(this),'..'))
scripts<-c(paste0('experiments/Fig',sprintf('%02d',1:10),'/plot.R'),
 paste0('experiments/Table',sprintf('%02d',1:9),'/make_table.R'))
for(script in scripts) {
 status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(file.path(root,script))))
 if(status!=0) stop(paste('Rendering failed:',script))
}
cat('Saved-result figures/tables completed.\n')
