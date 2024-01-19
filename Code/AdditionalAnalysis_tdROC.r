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

# time-dependent ROC for a single imputation
labels <- data.frame(outcomes=c("Rec", "MSM"),
                     long.outcomes=c("Recurrence-free survival", "Melanoma-specific survival"),
                     types=c("training", "validation"),
                     long.types=c("(Europe)", "(Australia)"))
for (data.type in labels$types){
  imputed.data <- eval(parse(text=paste0("imputed.", data.type, ".data")))
  single.imputation <- complete(imputed.data, 1)
  lp.train <- predict(f.mi.BS.Rec.5,
                      newdata=single.imputation,
                      type="lp")

  for (outcome in labels$outcomes){
    S <- eval(parse(text=paste0("S.", outcome, ".", data.type)))
    time.ROC <- timeROC::timeROC(T=S[,1],
                                 delta=S[,2],
                                 marker=lp.train,
                                 times=5,
                                 cause=1,
                                 iid=FALSE)
    png(file=paste0("Z:/Project Melanoom/PaperMelanoma/Results/Figures/tdROC.", outcome, ".", data.type, ".png"),
        width=1800, height=1800, units="px", res=300)
    plot(y=time.ROC$TP[,2], x=time.ROC$FP[,1], type="l",
         xlab="False positive fraction (1-specificity)",
         ylab="True positive fraction (sensitivity)",
         main=paste(labels[labels$outcomes==outcome, "long.outcomes"], labels[labels$types==data.type, "long.types"]))
    abline(a=0, b=1, col="grey")
    text(x=0.1, y=1, label=paste("AUC:", sprintf("%.2f", time.ROC$AUC[2])))
    dev.off()
  }
}

# combine calibration plots into one
png(file=paste0("Z:/Project Melanoom/PaperMelanoma/Results/tdROC.png"),
    width=16, height=16, units="cm", res=300)
par(mar=rep(0, 4))
layout(matrix(1:4, ncol=2, byrow=TRUE))
Centers.val <- rep(1:4, 3)
for (data.type in labels$types){
  for (outcome in labels$outcomes){
    plot(NA, xlim=0:1, ylim=0:1, xaxt="n", yaxt="n", bty="n")
    img <- png::readPNG(paste0("Z:/Project Melanoom/PaperMelanoma/Results/Figures/tdROC.", outcome, ".", data.type, ".png"))
    rasterImage(img, 0, 0, 1, 1)
  }
}
grDevices::dev.off()
