# Fig03: terminal R entry point, using the same settings resolver as Python.
cmd<-commandArgs(FALSE)
this<-sub('^--file=','',cmd[grepl('^--file=',cmd)][1])
root<-normalizePath(dirname(this))
while(!file.exists(file.path(root,'shared/settings.py'))) root<-dirname(root)
python<-Sys.getenv('JASA_PYTHON','python3')
status<-system2(python,c(shQuote(file.path(root,'scripts/run_r_adapter.py')),'er','r_baselines',
                       'Fig03', vapply(commandArgs(TRUE),shQuote,character(1))))
quit(status=status)
