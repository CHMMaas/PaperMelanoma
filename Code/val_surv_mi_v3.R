# # Input for val.surv.mi
# p<-as.matrix(prediction)   # matrix with predicted probabilities for imputation i in columns (complete case analyis: one column)
# y<-Surv(time,status)  # Time to event outcome as Surv object (time,status)
# g<-4 # number of risk groups; default=5
# main<-"Plot label" # Default = ""
# save_plots<-FALSE # save plots to file? default=FALSE
# time<-30 #time point at which to evaluate the predicted probabilities; default=NULL (not entered), the maximum time point will be taken
#
# val.surv.mi(p=p,y=y,g=g,main=main,save_plots = save_plots,time=time)

Rubin.combine<-function(est,se){
  m<-length(est)
  est.mi<-mean(est)
  var.w<-mean(se^2)
  var.b<-0
  if (m>1) var.b<-sum((est - est.mi)^2)/(m - 1)
  se.mi<-sqrt(var.w+(1+1/m)*var.b)
  return(list(est=est.mi,se=se.mi))
}

val.surv.mi<-function(p,y,g=5,main="",file_path="",save_plots=FALSE,time=NULL,
                      lim=c(0,1),dist=TRUE,CI=FALSE, df=3){
  if (!is.null(time)){
    y[y[,1]>time,2]<-0
    y[y[,1]>time,1]<-time
  }

  lp<-log(-log(1-p))

  n<-length(y)
  m.imp.val<-ncol(lp)

  cindex<-rep(0,m.imp.val)
  cindex.se<-rep(0,m.imp.val)
  slope<-rep(0,m.imp.val)
  slope.se<-rep(0,m.imp.val)
  int<-rep(0,m.imp.val)
  int.se<-rep(0,m.imp.val)
  cindex<-rep(0,m.imp.val)
  cindex.se<-rep(0,m.imp.val)

  p.groups<-array(rep(0,g*m.imp.val),dim=c(m.imp.val,g),dimnames=list(1:m.imp.val,1:g))
  y.groups<-array(rep(0,2*g*m.imp.val),dim=c(m.imp.val,g,2),dimnames=list(1:m.imp.val,1:g,c("obs","se")))
  lp.range<-min(lp)+0:100*(max(lp)-min(lp))/100
  lp.sm<-array(rep(0,101*m.imp.val),dim=c(101,m.imp.val),dimnames=list(1:101,1:m.imp.val))
  lp.sm.se<-lp.sm

  for (i in 1:m.imp.val){
    lp.val<-lp[,i]

    f.val<-coxph(y~lp.val)
    f.val.offset<-coxph(y~offset(lp.val))
    f.val.rcs<-cph(y~rcs(lp.val,df),x=TRUE,y=TRUE) # TODO: make df adjustable

    surv.sm<-Predict(f.val.rcs,lp.val=lp.range,conf.int = 0.95, conf.type = "simultaneous",time=max(y[,1]))
    # print(surv.sm)
    # print(c(sum(surv.sm$upper==0), sum(surv.sm$lower==0), sum(surv.sm$upper==1), sum(surv.sm$lower==1)))
    lp.sm[,i]<-log(-log(surv.sm$yhat))
    lp.sm.se[,i]<-(log(-log(surv.sm$upper))-log(-log(surv.sm$lower)))/(qnorm(.975)-qnorm(.025))
    #cindex[i]<-f.val$stats["C"]

    rc<-rcorr.cens(-lp.val,y)
    cindex[i]<-rc["C Index"]
    cindex.se[i]<-rc["S.D."]/2

    slope[i]<-f.val$coefficients[[1]]
    slope.se[i]<-sqrt(vcov(f.val)[[1,1]])

    sf<-survfit(f.val.offset,conf.type="log-log")
    log.H<-log(-log(tail(sf$surv,1)))
    log.H.upper<-log(-log(tail(sf$upper,1)))
    int[i]<-log.H-mean(f.val.offset$linear.predictors)
    int.se[i]<-(log.H-log.H.upper)/qnorm(.975)

    p.val<-p[,i]
    quants<-quantile(p.val,(1:(g-1))/g)
    cuts<-cut(p.val,breaks=c(0,quants,1))
    p.groups[i,]<-tapply(p.val,cuts,mean)
    for (j in 1:g){
      sub<-(cuts==levels(cuts)[j])
      if (sum(y[sub,2])>0){
        # sf<-survfit(y~1,subset=sub)
        # y.groups[i,j,1]<-tail(sf$cumhaz,1)
        # y.groups[i,j,2]<-tail(sf$std.chaz,1)} else {y.groups[i,j,]<- -Inf}
        sf<-survfit(y~1,subset=sub,conf.type="log-log")
        y.groups[i,j,1]<-log(-log(tail(sf$surv,1)))
        y.groups[i,j,2]<-(log(-log(tail(sf$surv,1)))-log(-log(tail(sf$upper,1))))/qnorm(.975)} else {y.groups[i,j,]<- -Inf}
    }
  }

  p.mi<-colMeans(p.groups)
  obs.mi<-rep(0,g)
  obs.mi.lower<-rep(0,g)
  obs.mi.upper<-rep(0,g)
  for (j in 1:g)
  {
    RC<-Rubin.combine(y.groups[,j,1],y.groups[,j,2])
    obs.mi[j]<-1-exp(-exp(RC$est))
    obs.mi.lower[j]<-1-exp(-exp(RC$est+qnorm(.025)*RC$se))
    obs.mi.upper[j]<-1-exp(-exp(RC$est+qnorm(.975)*RC$se))
  }

  p.sm.mi<-1-exp(-exp(lp.range))
  obs.sm.mi<-rep(0,101)
  obs.sm.mi.lower<-rep(0,101)
  obs.sm.mi.upper<-rep(0,101)
  for (j in 1:101){
    RC<-Rubin.combine(lp.sm[j,],lp.sm.se[j,])
    obs.sm.mi[j]<-1-exp(-exp(RC$est))
    obs.sm.mi.lower[j]<-1-exp(-exp(RC$est+qnorm(.025)*RC$se))
    obs.sm.mi.upper[j]<-1-exp(-exp(RC$est+qnorm(.975)*RC$se))
  }

  if (save_plots){
    png(file=file_path)
  }
  par(mar = c(5,5,2,1))

  par(mar = c(5,5,2,1))
  plot(lim,lim,type='l',xlab="Predicted probability",ylab="Observed frequency",
       main=main, lwd=1, bty='n', cex.axis=2, cex.lab=2, cex.main=2)

  if (CI){
    polygon(x=c(p.sm.mi,rev(p.sm.mi)),y=c(obs.sm.mi.lower,rev(obs.sm.mi.upper)),border = NA,col="Lightgray")
    lines(p.sm.mi,obs.sm.mi,lwd=2,col="Darkgray")
  }

  lines(lim,lim)
  abline(v=quants,col="darkgrey",lwd=1,lty=2)

  if (dist){
    line.bins <- 0.0
    length.seg <- 1
    dist.label <- 0.04
    dist.label2 <- 0.03
    d0lab <- 0
    d1lab <- 1
    cex.d01 <- 0.7
    x <- rowMeans(p)
    bins <- seq(0, min(1,max(lim[2])), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0	<-table(cut(x,bins))
    j0	<-f0 > 0
    bins0 <-(bins[-101])[j0]
    f0	<-f0[j0]
    maxf <-max(f0)
    f0	<-(0.1*f0)/maxf
    #segments(bins0,line.bins,bins0,length.seg*-f0+line.bins)
    segments(bins0,line.bins,bins0,length.seg*f0+line.bins, col="grey")
  }

  # quantiles
  segments(p.mi,obs.mi.lower,p.mi,obs.mi.upper)
  points(p.mi,obs.mi,pch=20)

  int.mi<-Rubin.combine(int,int.se)
  slope.mi<-Rubin.combine(slope,slope.se)
  cindex.mi<-Rubin.combine(cindex,cindex.se)

  legend(lim[1], lim[2], c(paste("n =",format(n,big.mark=",")),
                           paste("a =",format(round(int.mi$est,2),nsmall=2)),
                           paste("b =",format(round(slope.mi$est,2),nsmall=2)),
                           paste("c =",format(round(cindex.mi$est,2),nsmall=2))),
         box.col="white",  bg = "white",cex=1)

  if (save_plots){
    dev.off()
  }

  return(list(main=main,
              n=n,quants=quants,
              p.mi=p.mi,
              obs.mi=obs.mi,
              obs.mi.lower=obs.mi.lower,
              obs.mi.upper=obs.mi.upper,
              int=int.mi$est,
              int.lower=int.mi$est+qnorm(.025)*int.mi$se,
              int.upper=int.mi$est+qnorm(.975)*int.mi$se,
              slope=slope.mi$est,
              slope.lower=slope.mi$est+qnorm(.025)*slope.mi$se,
              slope.upper=slope.mi$est+qnorm(.975)*slope.mi$se,
              cindex=cindex.mi$est,
              cindex.lower=cindex.mi$est+qnorm(.025)*cindex.mi$se,
              cindex.upper=cindex.mi$est+qnorm(.975)*cindex.mi$se
              ))

}

# n <- 100
# m <- 5 # number of imputations
# p <- matrix(runif(n*m, 0, 1), n, m)
# time_until_event <- runif(n, 0, 5)
# status <- rbinom(n, 1, 0.5)
# y <- Surv(time_until_event, status)
# g <- 4
# main <- "Plot label"
# save_plots <- FALSE
# time <- 30
# val.surv.mi(p=p,y=y,g=g,main=main,save_plots=save_plots,time=time)
