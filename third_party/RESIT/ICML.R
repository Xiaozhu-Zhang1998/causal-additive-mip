# Copyright (c) 2010-2012 Jonas Peters. Original COPYING was not supplied.
# Retained fit/test interface only; the unused ICML driver is omitted.
fit_and_test_independence <- function(x,y,z,alpha,model,parsModel = list(),indtest, parsIndtest)
    # fits x using y and tests against z
{
    y <- as.matrix(y)
    z <- as.matrix(z)
    
    ####
    # fit x using y
    ####
    mod_fit <- train_model(model,y,x,parsModel)
    r2 <- mod_fit$resid
    
    Tquan <- indtestAll(indtest,z,r2,alpha,parsIndtest)
    return(Tquan)
}
