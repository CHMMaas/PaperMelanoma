#libraries
library(Hmisc)
library(survAUC)
library(pROC)
library(mice)

# load data
load("Z:/Project Melanoom/PaperMelanoma/Data/training/training.data.RData")
load("Z:/Project Melanoom/PaperMelanoma/Data/validation/validation.data.RData")
file.path <- "Z:/Project Melanoom/PaperMelanoma/Results/"

# calculate Harrell's C and Uno's C
Harrell.Rec.train <- c()
Uno.Rec.train <- c()
Harrell.Rec.val <- c()
Uno.Rec.val <- c()
for (i in 1:m){
  lp.train <- predict(f.mi.BS.Rec.5, newdata=mice::complete(imputed.training.data, i), type="lp")
  Harrell.Rec.train <- c(Harrell.Rec.train, Hmisc::rcorr.cens(-lp.train, S=S.Rec.5.training)["C Index"])
  Uno.Rec.train <- c(Uno.Rec.train, survAUC::UnoC(Surv.rsp=S.Rec.5.training, Surv.rsp.new=S.Rec.5.training, lpnew=lp.train, time=5))

  lp.val <- predict(f.mi.BS.Rec.5, newdata=mice::complete(imputed.validation.data, i), type="lp")
  Harrell.Rec.val <- c(Harrell.Rec.val, Hmisc::rcorr.cens(-lp.val, S=S.Rec.5.validation)["C Index"])
  Uno.Rec.val <- c(Uno.Rec.val, survAUC::UnoC(Surv.rsp=S.Rec.5.validation, Surv.rsp.new=S.Rec.5.validation, lpnew=lp.val, time=5))
}
write.table(round(data.frame(Harrell=c(mean(Harrell.Rec.train),
                                 mean(Uno.Rec.train)),
                       Uno=c(mean(Harrell.Rec.val),
                             mean(Uno.Rec.val))), 2),
            file="Z:/Project Melanoom/PaperMelanoma/Results/Harrell.Uno.txt",
            sep=",")

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

      roc.i <- pROC::roc(y, lp)
      sens.df <- cbind(sens.df, as.numeric(roc.i$sensitivities))
      spec.df <- cbind(spec.df, as.numeric(roc.i$specificities))
    }
    sens.avg <- rowMeans(sens.df)
    spec.avg <- rowMeans(spec.df)

    assign(paste0("sens.avg.", outcome, ".", type), sens.avg)
    assign(paste0("spec.avg.", outcome, ".", type), spec.avg)
  }
}

for (type in c("training", "validation")){
  for (outcome in labels$outcomes){
    png(filename=paste0(file.path, "/Figures/ROC.", outcome, ".", type, ".MI.png"), res=300, height=1500, width=1500)
    par(pty="s", mar=rep(3, 4))
    spec.outcome <- eval(parse(text=paste0("spec.avg.", outcome, ".", type)))
    sens.outcome <- eval(parse(text=paste0("sens.avg.", outcome, ".", type)))
    plot(x=spec.outcome, y=sens.outcome, lty=1, type="l",
         main=paste(labels[labels$outcomes==outcome, "long.outcomes"],
                    labels[labels$types==type, "long.types"]),
         xlab="Specificity", ylab="Sensitivty",
         asp=1,
         xlim=c(1, 0),
         ylim=c(0, 1))
    abline(a=1, b=-1, col="grey")
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
