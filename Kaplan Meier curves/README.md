# Kaplan-Meier curve of recurrence-free survival
load("Z:/.../Kaplan.Meier.curve.RFS.Rdata")

survminer::ggsurvplot(KM.curve.RFS, conf.int=TRUE)


# Kaplan-Meier curve of melanoma-specific survival 
load("Z:/.../Kaplan.Meier.curve.MSS.Rdata")

survminer::ggsurvplot(KM.curve.MSS, conf.int=TRUE)
