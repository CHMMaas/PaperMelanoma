h0.MSM <- 0.1438848
MSM.cal.fact <- 1.0973015

rscript <- function(status, age, ulceration, location, breslow, tumor_burden) {
  h0.Rec <- 0.2534864
  coef.Rec <- c(1.0264235, 0.3514440, 0.4914014, 0.2621760, 0.3803103,
                0.7062802, 0.8212041, 0.1866868, -0.3838754)
  names(coef.Rec) <- c("SNstatus=Positive", "Age.SN", "Ulceration=Yes",
                       "Loc_CAT=leg", "Loc_CAT=trunk", "Loc_CAT=headneck",
                       "Breslow", "Rdamcrit", "SNstatus=Positive * Breslow")
  c.Breslow <- 0.7154091
  center.Rec <- 2.031039
  lp <- coef.Rec["SNstatus=Positive"]*status+
    coef.Rec["Age.SN"]*log(age)+
    coef.Rec["Ulceration=Yes"]*ulceration+
    coef.Rec["Loc_CAT=leg"]*I(location == 1)+
    coef.Rec["Loc_CAT=trunk"]*I(location == 2)+
    coef.Rec["Loc_CAT=headneck"]*I(location == 3)+
    coef.Rec["Breslow"]*(log(breslow)-c.Breslow)+
    coef.Rec["Rdamcrit"]*log(tumor_burden)+
    coef.Rec["SNstatus=Positive * Breslow"]*(log(breslow)-c.Breslow)*status-center.Rec
  p.Rec <- 1-exp(-h0.Rec*exp(lp))
  return(as.numeric(p.Rec))
}
rscript(0, 55, 0, 0, 7, 1)
