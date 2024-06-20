# Make the Kaplan-Meier plots using
load("Z:/Project Melanoom/PaperMelanoma/Results/Kaplan.Meier.curve.RFS.Rdata")

load("Z:/Project Melanoom/PaperMelanoma/Results/Kaplan.Meier.curve.MSS.Rdata")

survminer::ggsurvplot(KM.curve.RFS, conf.int=TRUE)

survminer::ggsurvplot(KM.curve.MSS, conf.int=TRUE)
