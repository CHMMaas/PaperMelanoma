# clear all
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) dev.off()
cat("\014")

#libraries
library(Hmisc)
library(survivalROC)
library(mice)
library(survAUC) # Uno's C-index not by Uno
library(survC1)  # Uno's C-index by Uno
library(dplyr)

# load data
load("Z:/Project Melanoom/PaperMelanoma/Data/training/training.data.RData")
load("Z:/Project Melanoom/PaperMelanoma/Data/validation/validation.data.RData")
file.path <- "Z:/Project Melanoom/PaperMelanoma/Results/"

# calculate Harrell's C and Uno's C
lp.train <- c()
lp.val <- c()
Harrell.Rec.train <- c()
Uno.Rec.train <- c()
Uno.Rec.train.correct <- c()
Harrell.Rec.val <- c()
Uno.Rec.val <- c()
Uno.Rec.val.correct <- c()
for (i in 1:m){
  lp.train.i <- predict(f.mi.BS.Rec.5,
                      newdata=mice::complete(imputed.training.data, i),
                      type="lp")
  lp.train <- cbind(lp.train, lp.train.i)
  Harrell.Rec.train <- c(Harrell.Rec.train,
                         Hmisc::rcorr.cens(-lp.train.i, S=S.Rec.5.training)["C Index"])
  Uno.Rec.train <- c(Uno.Rec.train, survAUC::UnoC(Surv.rsp=S.Rec.5.training,
                                                  Surv.rsp.new=S.Rec.5.training,
                                                  lpnew=lp.train.i, time=5))

  # right function for Uno?
  df.unoC.train <- cbind(S.Rec.5.training[,1],
                         S.Rec.5.training[,2],
                         lp.train.i)
  unoC.train <- Est.Cval(mydata=df.unoC.train, tau=5, nofit=TRUE)
  Uno.Rec.train.correct <- c(Uno.Rec.train.correct, unoC.train$Dhat)

  lp.val.i <- predict(f.mi.BS.Rec.5,
                    newdata=mice::complete(imputed.validation.data, i),
                    type="lp")
  lp.val <- cbind(lp.val, lp.val.i)
  Harrell.Rec.val <- c(Harrell.Rec.val,
                       Hmisc::rcorr.cens(-lp.val.i, S=S.Rec.5.validation)["C Index"])
  Uno.Rec.val <- c(Uno.Rec.val, survAUC::UnoC(Surv.rsp=S.Rec.5.validation,
                                              Surv.rsp.new=S.Rec.5.validation,
                                              lpnew=lp.val.i, time=5))

  # right function for Uno?
  df.unoC.val <- cbind(S.Rec.5.validation[,1],
                   S.Rec.5.validation[,2],
                   lp.val.i)
  unoC.val <- Est.Cval(mydata=df.unoC.val, tau=5, nofit=TRUE)
  Uno.Rec.val.correct <- c(Uno.Rec.val.correct, unoC.val$Dhat)
}
write.table(round(data.frame(Training=c(mean(Harrell.Rec.train),
                                 mean(Uno.Rec.train),
                                 mean(Uno.Rec.train.correct)),
                       Validation=c(mean(Harrell.Rec.val),
                             mean(Uno.Rec.val),
                             mean(Uno.Rec.val.correct))), 2),
            file="Z:/Project Melanoom/PaperMelanoma/Results/Harrell.Uno.txt",
            sep=",")

# # test
# imp.data.i <- complete(imputed.validation.data, 1)
# df.unoC.val <- cbind(S.Rec.5.validation[,1],
#                      S.Rec.5.validation[,2],
#                      as.numeric(imp.data.i[,"SNstatus"])-1,
#                      imp.data.i["Age.SN"],
#                      as.numeric(imp.data.i[,"Ulceration"]),
#                      as.numeric(imp.data.i[,"Loc_CAT"]),
#                      log(imp.data.i$Breslow)-c.Breslow,
#                      log(imp.data.i$Rdamcrit),
#                      I(log(imp.data.i$Breslow)-c.Breslow)*(as.numeric(imp.data.i[,"SNstatus"])-1))
# unoC.i <- Inf.Cval(mydata=df.unoC.val, tau=5, itr=1, seed=1)
# (unoC.i$Dhat+1)/2       # 0.8621638
# Uno.Rec.val.correct[1]  # 0.8586629

# linear predictor
lp.train <- c()
lp.val <- c()
for (i in 1:m){
  lp.train.i <- predict(f.mi.BS.Rec.5,
                        newdata=mice::complete(imputed.training.data, i),
                        type="lp")
  lp.train <- cbind(lp.train, lp.train.i)

  lp.val.i <- predict(f.mi.BS.Rec.5,
                      newdata=mice::complete(imputed.validation.data, i),
                      type="lp")
  lp.val <- cbind(lp.val, lp.val.i)
}

# ROC curve
labels <- data.frame(outcomes=c("Rec", "MSM"),
                     long.outcomes=c("Recurrence-free survival", "Melanoma-specific survival"),
                     types=c("training", "validation"),
                     long.types=c("(Europe)", "(Australia)"))
for (type in c("training", "validation")){
  if (type=="training"){
    lp <- lp.train
  } else{
    lp <- lp.val
  }
  for (outcome in labels$outcomes){
    sens.df <- c()
    spec.df <- c()
    for (i in 1:m){
      if (outcome=="Rec"){
        y <- eval(parse(text=paste0("mice::complete(imputed.", type, ".data, i)$Rec.composite")))
      } else{
        y <- eval(parse(text=paste0("mice::complete(imputed.", type, ".data, i)$MSM==\"Death by melanoma\"")))
      }

      roc.i <- pROC::roc(response=as.numeric(y),
                         predictor=lp[,i],
                         direction="<")
      if (i==1){
        sens.df <- data.frame(thresholds=roc.i$thresholds,
                              sensitivities=roc.i$sensitivities) %>%
          filter_all(all_vars(is.finite(.)))
        spec.df <- data.frame(thresholds=roc.i$thresholds,
                               specificities=roc.i$specificities) %>%
          filter_all(all_vars(is.finite(.)))
      } else{
        sens.df <- sens.df %>% left_join(sens.df,
                                         as.numeric(roc.i$sensitivities),
                                         by="thresholds")
        spec.df <- spec.df %>% left_join(spec.df,
                                         as.numeric(roc.i$specificities),
                                         by="thresholds")
      }
    }
    sens.avg <- rowMeans(sens.df[, -1])
    spec.avg <- rowMeans(spec.df[, -1])

    assign(paste0("sens.avg.", outcome, ".", type), sens.avg)
    assign(paste0("spec.avg.", outcome, ".", type), spec.avg)
  }
}

for (type in c("training", "validation")){
  for (outcome in labels$outcomes){
    png(filename=paste0(file.path, "/Figures/ROC.", outcome, ".", type, ".MI.png"), res=300, height=1500, width=1500)
    par(pty="s", mar=rep(4, 4))
    sens.outcome <- eval(parse(text=paste0("sens.avg.", outcome, ".", type)))
    spec.outcome <- eval(parse(text=paste0("spec.avg.", outcome, ".", type)))

    sens.q <- quantile(sens.outcome)
    spec.q <- quantile(spec.outcome)
    FPR.q <- 1-spec.q

    plot(x=1-spec.outcome, y=sens.outcome, lty=1, type="l",
         main=paste(labels[labels$outcomes==outcome, "long.outcomes"],
                    labels[labels$types==type, "long.types"]),
         xlab="1-Specificity", ylab="Sensitivty",
         asp=1, xlim=c(0, 1), ylim=c(0, 1))
    abline(a=0, b=1, col="grey")
    points(x=rev(FPR.q), y=sens.q, pch=16)
    text(x=rev(FPR.q)[1:4], y=sens.q[1:4],
         label=paste0("(", sprintf("%.2f", rev(FPR.q)[1:4]), "; ",
                      sprintf("%.2f", sens.q[1:4]), ")"), cex=0.5, pos=4)
    text(x=1-FPR.q[5], y=sens.q[5],
         label=paste0("(", sprintf("%.2f", FPR.q[5]), "; ",
                      sprintf("%.2f", sens.q[5]), ")"), cex=0.5, pos=2)
    dev.off()
  }
}
png(file=paste0(file.path, "/Review/ROC.curves.png"), width=16, height=16, units="cm", res=300)
par(mar=rep(0, 4))
layout(matrix(1:4, ncol=2, byrow=TRUE))
Centers.val <- rep(1:4, 3)
for (name in c("ROC.Rec.training.MI", "ROC.MSM.training.MI",
               "ROC.Rec.validation.MI", "ROC.MSM.validation.MI")){
  plot(NA, xlim=0:1, ylim=0:1, xaxt="n", yaxt="n", bty="n")
  img <- png::readPNG(paste0(file.path, "Figures/", name, ".png"))
  rasterImage(img, 0, 0, 1, 1)
}
grDevices::dev.off()
