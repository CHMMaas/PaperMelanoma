# helper functions for calibration plots
km <- function(sel, S.km, horizon){
  S.km.sel <- S.km[sel,]
  sf <- survfit(S.km.sel~1)
  year1 <- max(sf$time[sf$time<=horizon])
  1-sf$surv[sf$time==year1]
}
km.lower <- function(sel, S.km, horizon){
  S.km.sel <- S.km[sel,]
  sf <- survfit(S.km.sel~1)
  year1 <- max(sf$time[sf$time<=horizon])
  1-sf$lower[sf$time==year1]
}
km.upper <- function(sel, S.km, horizon){
  S.km.sel <- S.km[sel,]
  sf <- survfit(S.km.sel~1)
  year1 <- max(sf$time[sf$time<=horizon])
  1-sf$upper[sf$time==year1]
}
fun.event <- function(lp, h0)
{
  h <- h0*exp(lp)
  p <- 1-exp(-h)
  return(p)
}
missings.table <- function(data){
  missing.values.table <- data.frame(m=rep(0, ncol(data)))
  rownames(missing.values.table) <- colnames(data)
  for (variable in colnames(data)){
    missing.values.table[variable,] <- sum(is.na(data[,variable]))
  }
  print(missing.values.table[order(missing.values.table$m, decreasing=TRUE), , drop=FALSE])
}
