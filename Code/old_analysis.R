# De belangrijkste primaire factoren die we mee zouden willen nemen zijn (waarvan de laatste 2 iets minder dan de eerste 6):
# -	Leeftijd (continu)
# -	Gender (male/female)
# -	Locatie melanoom (o.a. arm, hoofd/hals)
# -	Breslow dikte (continu)
# -	Histologie (o.a. NM, ALM)
# -	Ulceratie (yes/no). (Er zijn slechts 128 unknown = 4.1%, ik was in de war met andere factor, dus eigenlijk weinig unknown!).
# -	Clark level. Dit was altijd een zeer sterke prognostische factor in allerlei studies, mogelijk wordt deze factor in de toekomst iets minder van belang maar we denken dat hij misschien toch wel meegenomen zou moeten worden.
# -	Hoeveelheid verwijderde sentinel nodes (=Tot_nr_SNs)
#
# Daarnaast mogelijk interessant:
# -	Hoeveelheid SN fields (totaal aantal locaties waar de SN is afgenomen, er zijn er een aantal met 2 of zelfs 3, zoals beide oksels en/of ??n lies)? En/of locatie SN field. Maar wij denken dat deze beide factoren te maken hebben met de locatie van het melanoom (bijv. een melanoom op de rug geeft de kans op meerdere SN fields, terwijl een melanoom op het been altijd maar ??n SN field heeft namelijk in de lies aan dezelfde kant). Locatie melanoom staat reeds hierboven genoemd om mee te nemen.
# -	Mitosis; een nieuwere sterke prognostische factor. Het woord nieuw zegt het al, veel missing data (no n=38, yes n=109, unknown/missing n=2957) (Zelfde geldt voor Regression, Sattelites, Invasion: ook nieuwere factoren. Maar idem, veel missing data. Dus ik denk dat deze allen zullen afvallen).
#
# Nog extra vragen:
# -	Voor de overleving; er zijn een aantal pati?nten die na hun recidief ziekte een aanvullende therapie hebben ondergaan (o.a. lymfkliertoilet, bestraling). Moet hier nog voor worden gecorrigeerd voor de overlevingsanalyses?


# Activation of required libraries  ### INSTALL THESE PACKAGES FIRST (Packages menu)
library(foreign) # Reading external (spps) data
library(rms) # Harrell's regression library
library(mice) # Multipe imputation
library(mitools)

# Path to data
# setwd("C:/Users/438007/Dropbox/Melanoma")
# setwd("C:/Users/David/Dropbox/Melanoma")
# setwd("C:/Users/dvanklaveren/Dropbox/Melanoma")

# Read data (note the forward slash, use.value.labels=FALSE to avoid strange labels)
dat.orig<-read.spss("./Data/BRWA SN negatives prediction new inclusioncriteria 09-05-2017.sav",to.data.frame=TRUE)
head(dat.orig)

#number of patients in original file
nrow(dat.orig)

#sink("output.txt",append = TRUE)
#cat("number of patients in original file \n")
#nrow(dat.orig)
#sink()

labels(dat.orig)[[2]]
#  [1] "EXCLUSION"
#  [2] "Exclusion_reason"
#  [3] "SNstatus"
#  [4] "Center"
#  [5] "Sex"
#  [6] "DOB"
#  [7] "Age"
#  [8] "Diagdate"
#  [9] "Ulceration"
# [10] "Loc_detail"
# [11] "LOC_D0"
# [12] "LOC_D1"
# [13] "Loc_CAT"
# [14] "Loc_CAT2"
# [15] "Side"
# [16] "Histology"
# [17] "Histo_comm"
# [18] "Breslow"
# [19] "Breslow_CAT"
# [20] "Clark"
# [21] "Mitosis"
# [22] "Reg"
# [23] "Sat"
# [24] "Invasion"
# [25] "SENTINELNODE______________"
# [26] "SNdate"
# [27] "SN_Field_1"
# [28] "SN_FI0"
# [29] "SN_FI1"
# [30] "Site1PAloc"
# [31] "Site1Side"
# [32] "nr_SNs_Field_1"
# [33] "Snpos_Field_1"
# [34] "Site1Dewar"
# [35] "Site1tumordiam"
# [36] "Site1Rdam"
# [37] "Site1Starz"
# [38] "SN_field_2"
# [39] "nr_SNs_Field_2"
# [40] "Snpos_Field_2"
# [41] "Site2tunordiam"
# [42] "Site2Dewar"
# [43] "Site2Rdam"
# [44] "Site2Starz"
# [45] "SNfield3"
# [46] "nr_SNs_Field_3"
# [47] "Snpos_field_3"
# [48] "Tot_nr_fields"
# [49] "Tot_nr_SNs"
# [50] "Tot_SNs_pos"
# [51] "nr_nonSNs"
# [52] "nonSN_pos"
# [53] "nonSN_Field_1"
# [54] "Dewar"
# [55] "subcap"
# [56] "AnnOncol."
# [57] "ann_onco"
# [58] "Rdamcrit"
# [59] "Starzd"
# [60] "Starzn"
# [61] "Rdam_cat"
# [62] "Rdamcat_simple"
# [63] "CLNDcat"
# [64] "CLNDdate"
# [65] "totalnr.nodes"
# [66] "nr.pos.nodes"
# [67] "add_posLN_CAT"
# [68] "add_posLN_CAT1"
# [69] "CLNDcomm"
# [70] "RECURRENCE___________"
# [71] "Recurrence"
# [72] "Date_1st_Rec"
# [73] "Site1stRec"
# [74] "Rec_Type"
# [75] "Rec_CAT"
# [76] "Rec_comm"
# [77] "REC_C0"
# [78] "REC_C1"
# [79] "FOLLOWUP_______________"
# [80] "TLND"
# [81] "AdjTher"
# [82] "Trials"
# [83] "Last_FU"
# [84] "Status"
# [85] "Date_today"
# [86] "COMMENT_FU"
# [87] "dead"
# [88] "SURVIVAL_______________"
# [89] "FU_time_months"
# [90] "Time_1st_recur_mo"
# [91] "DFS_recurrence"
# [92] "DFS"
# [93] "filter_."


vars<-c(
	"Center",
	"Sex",
	"Age",
	"DOB",
	"Ulceration",
	"Loc_CAT2",
	"Histology",
	"Breslow",
	"Clark",
	"Mitosis",
	"Tot_nr_SNs",
	"nr_SNs_Field_1",
	"nr_SNs_Field_2",
	"nr_SNs_Field_3",
	"SNdate",
	"Recurrence",
	"Date_1st_Rec",
	"Last_FU",
	"dead",
	"Status"
)

dat<-dat.orig[dat.orig$EXCLUSION=="Inclusion",vars]
nrow(dat)
# 36 patients excluded

### dead missing: everything is missing
dat[is.na(dat$dead),]
dat<-dat[!is.na(dat$dead),]
nrow(dat)
# 4 patients excluded

# Delete Clark i (in situ excluderen)
dat<-dat[dat$Clark!="i",]
nrow(dat)
# 12 patients excluded


#number of patients that meet the inclusion criteria
nrow(dat)
describe(dat)

# Disease specific death
table(dat$Status,dat$dead,useNA="ifany")
dat$DOD<-1*(dat$Status=="DOD")

# Define multiple fields
table(dat$nr_SNs_Field_1,dat$nr_SNs_Field_2,dat$nr_SNs_Field_3,useNA="ifany")
dat$multiple.fields<-ifelse(!is.na(dat$nr_SNs_Field_2)&dat$nr_SNs_Field_2>0,1,0)

# Age at SN is better than age at diagnosis
dat$SNdate[dat$SNdate==0]<-NA # Set SNdate=0 to NA
dat$Age.SN<-floor((dat$SNdate-dat$DOB)/(24*60*60)/365.25)
table(dat$Age,dat$Age.SN)


# Define times since diagnosis
dat$FU.time<-(dat$Last_FU-dat$SNdate)/(24*60*60)/365.25
dat[is.na(dat$FU.time),]
median(dat$FU.time[!is.na(dat$FU.time)])  # 4.5
dat$FU.time[is.na(dat$FU.time)]<-median(dat$FU.time[!is.na(dat$FU.time)])
dat$DOD[is.na(dat$FU.time)]<-0

# 39 missing reccurrence times
dat$Rec.time<-(dat$Date_1st_Rec-dat$SNdate)/(24*60*60)/365.25
dat$Rec.time[dat$Recurrence=="no recurrence"|is.na(dat$SNdate)]<-dat$FU.time[dat$Recurrence=="no recurrence"|is.na(dat$SNdate)]
dat[is.na(dat$Rec.time),]
median(dat$Rec.time[dat$Recurrence=="recurred"&!is.na(dat$Date_1st_Rec)])  # 1.68
dat$Rec.time[dat$Recurrence=="recurred"&is.na(dat$Date_1st_Rec)]<-median(dat$Rec.time[dat$Recurrence=="recurred"&!is.na(dat$Date_1st_Rec)])


# Replace "Unknown" by NA
dat$Sex[dat$Sex=="Unknown"]<-NA
dat$Ulceration[dat$Ulceration=="unknown"]<-NA
dat$Loc_CAT2[dat$Loc_CAT2=="unknown"]<-NA
dat$Histology[dat$Histology=="unknown"]<-NA
dat$Clark[dat$Clark=="unknown"]<-NA
dat$Mitosis[dat$Mitosis=="Unknown"]<-NA

# Delete unused levels
dat<-droplevels(dat)
#levels(dat$Clark)[1]<-"i/ii"

# Variable description
head(dat)
describe(dat)

# Save distributions of all variables, necessary for Harrell summaries
dd<-datadist(dat)
options(datadist='dd')


#######################
### Multiple imputation
#######################

vars.mi<-c(
	"Center",
	"Sex",
	"Age.SN",
	"Ulceration",
	"Loc_CAT2",
	"Histology",
	"Breslow",
	"Clark",
	"Tot_nr_SNs",
	"multiple.fields",
	"Recurrence",
	"dead",
	"DOD"
)

dat.mi<-dat[,vars.mi]
dat.mi$logRec.time<-log(dat$Rec.time+1/365.25)
dat.mi$logFU.time<-log(dat$FU.time+1/365.25)

# Save distributions of all variables, necessary for Harrell summaries
dd<-datadist(dat.mi)
options(datadist='dd')
options(digits=8)

mi<-mice(dat.mi,seed=13)
# save (mi,file="mi.RData",compress=TRUE)
# load('mi.RData')


### Hulpfuncties
fun.event<-function(lp)
{
	h<-h0*exp(lp)
	p<-100*(1-exp(-h))
	return(round(p,ifelse(p<10,1,0)))
}

fun.event.DOD<-function(lp)
{
	h<-h0.DOD*exp(lp)
	p<-100*(1-exp(-h))
	return(round(p,ifelse(p<10,1,0)))
}


km<-function(sel){
	S.km.sel<-S.km[sel,]
	sf<-survfit(S.km.sel~1)
	year1<-max(sf$time[sf$time<=horizon])
	1-sf$surv[sf$time==year1]
	}

km.lower<-function(sel){
	S.km.sel<-S.km[sel,]
	sf<-survfit(S.km.sel~1)
	year1<-max(sf$time[sf$time<=horizon])
	1-sf$lower[sf$time==year1]
	}

km.upper<-function(sel){
	S.km.sel<-S.km[sel,]
	sf<-survfit(S.km.sel~1)
	year1<-max(sf$time[sf$time<=horizon])
	1-sf$upper[sf$time==year1]
	}
###


#############
### Modelling recurrence
#############

## Recurrence
S<-Surv(dat$Rec.time,dat$Recurrence=="recurred")
SF<-survfit(S~1,data=dat)
horizon <- 5
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])


## Center-specific survival
plot(survfit(S~Center,data=dat))
plot(survfit(S~1,data=dat))


#### All variables
form<-S ~
	Sex+
	rcs(Age.SN,4)+
	Ulceration+
	Loc_CAT2+
	Histology+
	rcs(Breslow,4)+
      rcs(Tot_nr_SNs,4)+
	multiple.fields

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum1 <- summary(f.mi)
test1 <- coef(f.mi)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

f.mi.boot<-f.mi

#### Backward selection
form<-S ~
###	Sex+
#1 	rcs(Age.SN,4)+
	Ulceration+
	Loc_CAT2+
#3	Histology+
	rcs(Breslow,4)
#2    rcs(Tot_nr_SNs,4)+
#4	multiple.fields


f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum2 <- summary(f.mi)
test2 <- coef(f.mi)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))



#### Non-linearity simplified
dat.complete<-complete(mi,1)
describe(dat.complete$Breslow)
f<-cph(S ~ rcs(Breslow,4),data=dat)
f
plot(Predict(f,Breslow))
f<-cph(S ~ log(Breslow),data=dat)
f
plot(Predict(f,Breslow))

#### Non-linearity simplified
form<-S ~
	Ulceration+
	Loc_CAT2+
	log(Breslow)

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum3 <- summary(f.mi)
test3 <- coef(f.mi)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))


#### Disease specific mortality
## Unknown should be 0
dat$DOD[is.na(dat$DOD)]<-0
S.DOD<-Surv(dat$FU.time,dat$DOD)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),S.DOD)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))








### Proportionality
z<-cox.zph(f.mi$fits[[1]])
z
plot(z)


### Limit FU to 5 year?
S.5<-S
S.5[S[,1]>5,1]<-5
S.5[S[,1]>5,2]<-0

# CHECK difference in recurrence probability
SF<-survfit(S~1)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])
SF<-survfit(S.5~1)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])

form.5<-update(form, S.5  ~ . )
f.mi.5<-fit.mult.impute(form.5,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi.5
anova(f.mi.5)
# CHECKIT
sum4 <- summary(f.mi.5)
test4 <- coef(f.mi.5,Clark="i/ii")

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi.5,newdata=complete(mi,i)),f.mi.5$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

z<-cox.zph(f.mi.5$fits[[1]])
z
plot(z)


# complete case
f<-cph(form,data=dat.mi)
f
anova(f)



###### BOOTSTRAP
opt<-NULL
slope<-NULL
intercept<-NULL
set.seed(0)
for (i in 1:5)
{
	v<-validate(f.mi.boot$fits[[i]],bw=TRUE,rule='p',sls=0.05,B=100,pr=FALSE,type='individual')
	opt<-c(opt,v["Dxy","optimism"]/2)
	slope<-c(slope,v["Slope","test"])
}
opt
mean(opt)
mean(slope)


### Additional value Mitosis?
Mit<-dat$Mitosis
lp<-f.mi$linear
f<-cph(S~offset(lp)+Mit)
f
summary(f,Mit=c("No"))



### Fit disease specific mortality to recurrence lp
lp<-f.mi$linear
f.DOD<-cph(S.DOD~lp,y=TRUE,x=TRUE)
f.DOD

####### Nomogram
#### Limits of Breslow
table(dat$Breslow)
Breslow.class<-c((2:10)/10,1+(1:5)/5,2.5,3:7)
nom.0<-nomogram(f.mi,maxscale=10,Breslow=Breslow.class)

nom<-nomogram(f.mi,maxscale=10,Breslow=Breslow.class,lp=FALSE)
attributes(nom)$names[attributes(nom)$names=="Loc_CAT2"]<-"Location"
plot(nom,total.sep.page=TRUE,col.grid = gray(c(0.8, 0.95)))
nom.0
nom

# Slope and intercept of recurrence and DOD
rc<-(1/attributes(nom.0)$info$sc)
int<-attributes(nom.0)$info$Intercept
int
log(0.2)*f.mi$coef["Breslow"]-f.mi$center # check

rc.DOD<-f.DOD$coef[1]*rc
int.DOD<-f.DOD$coef[1]*int

points<-0:15
lp.points<-int+rc*points
lp.points.DOD<-int.DOD+rc.DOD*points


#### Absolute risk predictions
f.basehaz<-basehaz(f.mi,TRUE)
f.DOD.basehaz<-basehaz(f.DOD,TRUE)


#### Set the time horizon to 5 years
horizon<-5

h0<-f.basehaz$hazard[f.basehaz$time==max(f.basehaz$time[f.basehaz$time<=horizon])]
mean(fun.event(lp))
#### Kaplan Meier horizon-year recurrence
SF<-survfit(S~1,data=dat)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])

h0.DOD<-f.DOD.basehaz$hazard[f.DOD.basehaz$time==max(f.DOD.basehaz$time[f.DOD.basehaz$time<=horizon])]
lp.DOD<-f.DOD$linear
mean(fun.event.DOD(lp.DOD))
#### Kaplan Meier horizon-year recurrence
SF<-survfit(S.DOD~1,data=dat)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])

dat.print<-data.frame('Predicted Recurrence'=fun.event(lp.points),'Predicted DOD'=fun.event.DOD(lp.points.DOD))
rownames(dat.print)<-0:15
print(dat.print)

# Risk score distribution
dat.complete<-complete(mi,1)
Breslow.trunc<-dat.complete$Breslow
Breslow.trunc[dat.complete$Breslow<.2]<-.2
Breslow.trunc[dat.complete$Breslow>7]<-7
lp.score<-predict(f.mi,newdata=data.frame(Ulceration=dat.complete$Ulceration,Loc_CAT2=dat.complete$Loc_CAT2,Breslow=Breslow.trunc))
mean(lp.score)
score<-round((lp.score-int)/rc,0)
h<-hist(score,plot=FALSE,breaks=0:16,right=FALSE)
quantile(score,c(.25,.5,.75))
quantile(score,c(1/3,2/3))
table(score)/length(score)
data.frame(table(score)/length(score))
# Check mean event rates of risk distribution
sum(h$density*fun.event(lp.points))
sum(h$density*fun.event.DOD(lp.points.DOD))

h.plot<-h
h.plot$density<-100*h$density
par(mar = c(5,5,2,5))
plot(points,fun.event(lp.points),xlab="Risk score",ylab="Predicted 5-year risk (%)",type="b",lwd=1,xlim=c(-.5,15.5),ylim=c(0,50),pch=15)
points(points,fun.event.DOD(lp.points.DOD),type="b",lwd=1,pch=16)
#grid(nx=NA,ny=NULL,col="light grey")
abline(h=(0:13)*5,col="light grey",lty=2)
legend(x=0,y=50,legend=c("Recurrence","Disease-specific mortality"),pch=c(15,16),cex=.9,bg='white')
par(new = TRUE)
plot(h.plot,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,25),main=NA,col='white')
axis(side = 4)
mtext(side = 4, line = 3, 'Risk score distribution (%)')

par(new = TRUE)
h.q<-h.plot
h.q$density[h.plot$breaks>6]<-0
plot(h.q,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,25),main=NA,col="#CCFF33")

par(new = TRUE)
h.q<-h.plot
h.q$density[h.plot$breaks<=6|h.plot$breaks>9]<-0
plot(h.q,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,25),main=NA,col= "#FFCC33")

par(new = TRUE)
h.q<-h.plot
h.q$density[h.plot$breaks<=9]<-0
plot(h.q,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,25),main=NA,col="#FF6633")

points(points+.5,fun.event(lp.points)/2,xlim=c(-.5,15.5),type="b",ylim=c(0,66),pch=15)
points(points+.5,fun.event.DOD(lp.points.DOD)/2,xlim=c(-.5,15.5),type="b",ylim=c(0,66),pch=16)


#### Set the time horizon to 1 year
horizon<-1

h0<-f.basehaz$hazard[f.basehaz$time==max(f.basehaz$time[f.basehaz$time<=horizon])]
mean(fun.event(lp))
#### Kaplan Meier horizon-year recurrence
SF<-survfit(S~1,data=dat)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])

h0.DOD<-f.DOD.basehaz$hazard[f.DOD.basehaz$time==max(f.DOD.basehaz$time[f.DOD.basehaz$time<=horizon])]
lp.DOD<-f.DOD$linear
mean(fun.event.DOD(lp.DOD))
#### Kaplan Meier horizon-year recurrence
SF<-survfit(S.DOD~1,data=dat)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=horizon])])


dat.print<-data.frame('Predicted Recurrence'=fun.event(lp.points),'Predicted DOD'=fun.event.DOD(lp.points.DOD))
rownames(dat.print)<-0:15
print(dat.print)


# Risk score distribution
dat.complete<-complete(mi,1)
Breslow.trunc<-dat.complete$Breslow
Breslow.trunc[dat.complete$Breslow<.2]<-.2
Breslow.trunc[dat.complete$Breslow>7]<-7
lp.score<-predict(f.mi,newdata=data.frame(Ulceration=dat.complete$Ulceration,Loc_CAT2=dat.complete$Loc_CAT2,Breslow=Breslow.trunc))
mean(lp.score)
score<-round((lp.score-int)/rc,0)
h<-hist(score,plot=FALSE,breaks=0:16,right=FALSE)
quantile(score,c(.25,.5,.75))
quantile(score,c(1/3,2/3))
table(score)/length(score)
# Check mean event rates of risk distribution
sum(h$density*fun.event(lp.points))
sum(h$density*fun.event.DOD(lp.points.DOD))

h.plot<-h
h.plot$density<-100*h$density
par(mar = c(5,5,2,5))
plot(points,fun.event(lp.points),xlab="Risk score",ylab="Predicted 1-year risk (%)",type="b",lwd=1,xlim=c(-.5,15.5),ylim=c(0,15),pch=15)
points(points,fun.event.DOD(lp.points.DOD),type="b",lwd=1,pch=16)
#grid(nx=NA,ny=NULL,col="light grey")
abline(h=(0:13)*2.5,col="light grey",lty=2)
legend(x=0,y=15,legend=c("Recurrence","Disease-specific mortality"),pch=c(15,16),cex=.9,bg='white')
par(new = TRUE)
plot(h.plot,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,30),main=NA,col='white')
axis(side = 4)
mtext(side = 4, line = 3, 'Risk score distribution (%)')

par(new = TRUE)
h.q<-h.plot
h.q$density[h.plot$breaks>6]<-0
plot(h.q,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,30),main=NA,col="#CCFF33")

par(new = TRUE)
h.q<-h.plot
h.q$density[h.plot$breaks<=6|h.plot$breaks>9]<-0
plot(h.q,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,30),main=NA,col= "#FFCC33")

par(new = TRUE)
h.q<-h.plot
h.q$density[h.plot$breaks<=9]<-0
plot(h.q,freq=FALSE,axes=FALSE, xlab=NA, ylab=NA,ylim=c(0,30),main=NA,col="#FF6633")

points(points+.5,fun.event(lp.points)*2,xlim=c(-.5,15.5),type="b",ylim=c(0,15),pch=15)
points(points+.5,fun.event.DOD(lp.points.DOD)*2,xlim=c(-.5,15.5),type="b",ylim=c(0,15),pch=16)














### Cross validation recurrence 5-year

Centers<-unique(dat$Center)
discrimination<-list()
discrimination.DOD<-list()
cuts<-5
horizon<-5

for (j in 1:4){
	S.j<-S[dat.mi$Center==Centers[j],]
	S.notj<-S[dat.mi$Center!=Centers[j],]
	S.DOD.j<-S.DOD[dat.mi$Center==Centers[j],]
	S.DOD.notj<-S.DOD[dat.mi$Center!=Centers[j],]

	form.notj<-update(form, S.notj  ~ . )
	f.mi.notj<-fit.mult.impute(form.notj,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE,sub=dat.mi$Center!=Centers[j])

	lp.notj<-f.mi.notj$linear
	f.DOD.notj<-cph(S.DOD.notj~lp.notj,y=TRUE,x=TRUE)

	f.basehaz.notj<-basehaz(f.mi.notj,TRUE)
	h0<-f.basehaz.notj$hazard[f.basehaz.notj$time==max(f.basehaz.notj$time[f.basehaz.notj$time<=horizon])]

	f.DOD.basehaz.notj<-basehaz(f.DOD.notj,TRUE)
	h0.DOD<-f.DOD.basehaz.notj$hazard[f.DOD.basehaz.notj$time==max(f.DOD.basehaz.notj$time[f.DOD.basehaz.notj$time<=horizon])]

	cindex<-NULL
	cindex.DOD<-NULL
	lp<-NULL
	lp.DOD<-NULL
	for (i in 1:5){
		lp.j<-predict(f.mi.notj,newdata=complete(mi,i)[dat.mi$Center==Centers[j],])
		rc<-rcorr.cens(-lp.j,S.j)
		cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
		lp<-cbind(lp,lp.j)

		lp.DOD.j<-predict(f.DOD.notj,newdata=lp.j)
		rc<-rcorr.cens(-lp.DOD.j,S.DOD.j)
		cindex.DOD<-rbind(cindex.DOD,c(rc["C Index"],rc["S.D."]/2))
		lp.DOD<-cbind(lp.DOD,lp.DOD.j)
	}

	discrimination[[j]]<-summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))
	discrimination.DOD[[j]]<-summary(MIcombine(as.list(cindex.DOD[,1]),as.list(cindex.DOD[,2]^2)))

	lp.avg<-rowMeans(lp)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.j[,2],cut.all,length)

	S.km<-S.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.6)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)

	lp.avg<-rowMeans(lp.DOD)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.DOD.j[,2],cut.all,length)

	S.km<-S.DOD.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event.DOD(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.4)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)


}

names(discrimination)<-Centers
discrimination
names(discrimination.DOD)<-Centers
discrimination.DOD




### Cross validation recurrence 1-year

Centers<-unique(dat$Center)
discrimination<-list()
discrimination.DOD<-list()
cuts<-5
horizon<-1

for (j in 1:4){
	S.j<-S[dat.mi$Center==Centers[j],]
	S.notj<-S[dat.mi$Center!=Centers[j],]
	S.DOD.j<-S.DOD[dat.mi$Center==Centers[j],]
	S.DOD.notj<-S.DOD[dat.mi$Center!=Centers[j],]

	form.notj<-update(form, S.notj  ~ . )
	f.mi.notj<-fit.mult.impute(form.notj,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE,sub=dat.mi$Center!=Centers[j])

	lp.notj<-f.mi.notj$linear
	f.DOD.notj<-cph(S.DOD.notj~lp.notj,y=TRUE,x=TRUE)

	f.basehaz.notj<-basehaz(f.mi.notj,TRUE)
	h0<-f.basehaz.notj$hazard[f.basehaz.notj$time==max(f.basehaz.notj$time[f.basehaz.notj$time<=horizon])]

	f.DOD.basehaz.notj<-basehaz(f.DOD.notj,TRUE)
	h0.DOD<-f.DOD.basehaz.notj$hazard[f.DOD.basehaz.notj$time==max(f.DOD.basehaz.notj$time[f.DOD.basehaz.notj$time<=horizon])]

	cindex<-NULL
	cindex.DOD<-NULL
	lp<-NULL
	lp.DOD<-NULL
	for (i in 1:5){
		lp.j<-predict(f.mi.notj,newdata=complete(mi,i)[dat.mi$Center==Centers[j],])
		rc<-rcorr.cens(-lp.j,S.j)
		cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
		lp<-cbind(lp,lp.j)

		lp.DOD.j<-predict(f.DOD.notj,newdata=lp.j)
		rc<-rcorr.cens(-lp.DOD.j,S.DOD.j)
		cindex.DOD<-rbind(cindex.DOD,c(rc["C Index"],rc["S.D."]/2))
		lp.DOD<-cbind(lp.DOD,lp.DOD.j)
	}

	discrimination[[j]]<-summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))
	discrimination.DOD[[j]]<-summary(MIcombine(as.list(cindex.DOD[,1]),as.list(cindex.DOD[,2]^2)))

	lp.avg<-rowMeans(lp)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.j[,2],cut.all,length)

	S.km<-S.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.25)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)

	lp.avg<-rowMeans(lp.DOD)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.DOD.j[,2],cut.all,length)

	S.km<-S.DOD.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event.DOD(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.1)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)


}

names(discrimination)<-Centers
discrimination
names(discrimination.DOD)<-Centers
discrimination.DOD












### Cross validation recurrence 1-year

Centers<-unique(dat$Center)
discrimination<-list()
discrimination.DOD<-list()
cuts<-5
horizon<-1

for (j in 1:4){
	S.j<-S[dat.mi$Center==Centers[j],]
	S.notj<-S[dat.mi$Center!=Centers[j],]
	form.notj<-update(form, S.notj  ~ . )
	f.mi.notj<-fit.mult.impute(form.notj,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE,sub=dat.mi$Center!=Centers[j])

	f.basehaz.notj<-basehaz(f.mi.notj,TRUE)
	h0<-f.basehaz.notj$hazard[f.basehaz.notj$time==max(f.basehaz.notj$time[f.basehaz.notj$time<=horizon])]

	cindex<-NULL
	lp<-NULL
	for (i in 1:5){
		lp.j<-predict(f.mi.notj,newdata=complete(mi,i)[dat.mi$Center==Centers[j],])
		rc<-rcorr.cens(-lp.j,S.j)
		cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
		lp<-cbind(lp,lp.j)
		}

	discrimination[[j]]<-summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

	lp.avg<-rowMeans(lp)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.j[,2],cut.all,length)

	S.km<-S.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.25)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)
}

names(discrimination)<-Centers
discrimination

























### Disease specific mortality

## Unknown should be 0
dat$DOD[is.na(dat$DOD)]<-0

## Mortality
S.DOD<-Surv(dat$FU.time,dat$DOD)


lp<-f.mi$linear
#f.DOD<-coxph(S.DOD~lp)
f.DOD<-coxph(S.DOD~offset(lp))





#############
### Modelling all cause mortality
#############

## Mortality
S<-Surv(dat$FU.time,dat$dead=="Yes")

## Center-specific survival
plot(survfit(S~Center,data=dat))
plot(survfit(S~1,data=dat))


#### All variables
form<-S ~
	Sex+
	rcs(Age.SN,4)+
	Ulceration+
	Loc_CAT2+
	Histology+
	rcs(Breslow,4)+
	Clark+
      rcs(Tot_nr_SNs,4)+
	multiple.fields

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum5 <- summary(f.mi)
test5 <- coef(f.mi,Clark="i/ii")

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

f.mi.boot<-f.mi

#### Backward selection, leave Clark in
form<-S ~
	Sex+
 	rcs(Age.SN,4)+
	Ulceration+
	Loc_CAT2+
#3	Histology+
	rcs(Breslow,4)+
	Clark
#1	rcs(Tot_nr_SNs,4)
#2	multiple.fields


f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum6 <- summary(f.mi)
test6 <- coef(f.mi,Clark="i/ii")

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))



#### Non-linearity simplified
describe(complete(mi,1)$Breslow)
f<-cph(S ~ rcs(Breslow,4),data=dat)
f
anova(f)
plot(Predict(f,Breslow))
f<-cph(S ~ log(Breslow),data=dat)
f
plot(Predict(f,Breslow))


f<-cph(S ~ rcs(Age.SN,4),data=dat)
f
anova(f)
plot(Predict(f,Age.SN))
f<-cph(S ~ ifelse(Age.SN>60,Age.SN,60),data=dat)
f
plot(Predict(f,Age.SN))


#### Non-linearity simplified
form<-S ~
	Sex+
 	ifelse(Age.SN>60,Age.SN,60)+
	Ulceration+
	Loc_CAT2+
	log(Breslow)+
	Clark

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum7 <- summary(f.mi)
test7 <- coef(f.mi,Clark="i/ii")

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))



#### Non-linearity simplified  WITHOUT CLARK
form<-S ~
	Sex+
 	ifelse(Age.SN>60,Age.SN,60)+
	Ulceration+
	Loc_CAT2+
	log(Breslow)

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum8 <- summary(f.mi)
test8 <- coef(f.mi)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))


### Proportionality
z<-cox.zph(f.mi$fits[[1]])
z
plot(z)


### Limit FU to 5 year?
S.5<-S
S.5[S[,1]>5,1]<-5
S.5[S[,1]>5,2]<-0
form.5<-update(form, S.5  ~ . )
f.mi.5<-fit.mult.impute(form.5,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi.5
anova(f.mi.5)
# CHECKIT
sum9 <- summary(f.mi)
test9 <- coef(f.mi.5)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi.5,newdata=complete(mi,i)),f.mi.5$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

z<-cox.zph(f.mi.5$fits[[1]])
z
plot(z)


# complete case
f<-cph(form,data=dat.mi)
f
anova(f)



###### BOOTSTRAP
opt<-NULL
slope<-NULL
intercept<-NULL
set.seed(0)
for (i in 1:5)
{
	v<-validate(f.mi.boot$fits[[i]],bw=TRUE,rule='p',sls=0.05,B=100,pr=FALSE,type='individual')
	opt<-c(opt,v["Dxy","optimism"]/2)
	slope<-c(slope,v["Slope","test"])
}
opt
mean(opt)
mean(slope)


### Additional value Mitosis?
Mit<-dat$Mitosis
lp<-f.mi$linear
f<-cph(S~offset(lp)+Mit)
f
summary(f,Mit=c("No"))


####### Nomogram
f.basehaz<-basehaz(f.mi,TRUE)
h0<-f.basehaz$hazard[f.basehaz$time==max(f.basehaz$time[f.basehaz$time<=5])]
mean(fun.event(lp))
#### Kaplan Meier 5-year recurrence
SF<-survfit(S~1,data=dat)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=5])])

#### Limits of Breslow
table(dat$Breslow)
Breslow.class<-c((2:10)/10,1+(1:5)/5,2.5,3:7)

nom.0<-nomogram(f.mi,maxscale=10,Breslow=Breslow.class)
# Baseline hazard
rc<-(1/attributes(nom.0)$info$sc)
int<-attributes(nom.0)$info$Intercept

nom<-nomogram(f.mi,maxscale=10,Breslow=Breslow.class,lp=FALSE,fun=fun.event,fun.at=fun.event(int+rc*(0:10)*2))
plot(nom)
nom.0
nom

points<-0:24
lp.points<-int+rc*points
print(data.frame(points,lp.points,fun.event(lp.points)))

points<-2*(0:10)
lp.points<-int+rc*points
print(data.frame('Total points'=points,'Predicted Value'=fun.event(lp.points)))




### Cross validation

Centers<-unique(dat$Center)
discrimination<-list()
cuts<-5


for (j in 1:4){
	S.j<-S[dat.mi$Center==Centers[j],]
	S.notj<-S[dat.mi$Center!=Centers[j],]
	form.notj<-update(form, S.notj  ~ . )
	f.mi.notj<-fit.mult.impute(form.notj,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE,sub=dat.mi$Center!=Centers[j])

	f.basehaz.notj<-basehaz(f.mi.notj,TRUE)
	h0<-f.basehaz.notj$hazard[f.basehaz.notj$time==max(f.basehaz.notj$time[f.basehaz.notj$time<=5])]

	cindex<-NULL
	lp<-NULL
	for (i in 1:5){
		lp.j<-predict(f.mi.notj,newdata=complete(mi,i)[dat.mi$Center==Centers[j],])
		rc<-rcorr.cens(-lp.j,S.j)
		cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
		lp<-cbind(lp,lp.j)
		}

	discrimination[[j]]<-summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

	lp.avg<-rowMeans(lp)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.j[,2],cut.all,length)

	S.km<-S.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.5)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)
}

names(discrimination)<-Centers
discrimination





#############
### Modelling disease specific mortality
#############

## Unknown should be 0
dat$DOD[is.na(dat$DOD)]<-0

## Mortality
S<-Surv(dat$FU.time,dat$DOD)

## Center-specific survival
plot(survfit(S~Center,data=dat))
plot(survfit(S~1,data=dat))


#### All variables
form<-S ~
	Sex+
	rcs(Age.SN,4)+
	Ulceration+
	Loc_CAT2+
	Histology+
	rcs(Breslow,4)+
	Clark+
      rcs(Tot_nr_SNs,4)+
	multiple.fields

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum10 <- summary(f.mi)
test10 <- coef(f.mi,Clark="i/ii")

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

f.mi.boot<-f.mi

#### Backward selection
form<-S ~
####	Sex+
#4 	rcs(Age.SN,4)+
	Ulceration+
	Loc_CAT2+
#2	Histology+
	rcs(Breslow,4)
#5	Clark
#1	rcs(Tot_nr_SNs,4)
#3	multiple.fields


f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum11 <- summary(f.mi)
test11 <- coef(f.mi)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))



#### Non-linearity simplified
describe(complete(mi,1)$Breslow)
f<-cph(S ~ rcs(Breslow,4),data=dat)
f
anova(f)
plot(Predict(f,Breslow))
f<-cph(S ~ log(Breslow),data=dat)
f
plot(Predict(f,Breslow))



#### Non-linearity simplified
form<-S ~
	Ulceration+
	Loc_CAT2+
	log(Breslow)

f.mi<-fit.mult.impute(form,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi
anova(f.mi)
# CHECKIT
sum12 <- summary(f.mi)
test12 <- coef(f.mi)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi,newdata=complete(mi,i)),f.mi$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))




### Proportionality
z<-cox.zph(f.mi$fits[[1]])
z
plot(z)


### Limit FU to 5 year?
S.5<-S
S.5[S[,1]>5,1]<-5
S.5[S[,1]>5,2]<-0
form.5<-update(form, S.5  ~ . )
f.mi.5<-fit.mult.impute(form.5,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE)
f.mi.5
anova(f.mi.5)
# CHECKIT
sum13 <- summary(f.mi.5)
test13 <- coef(f.mi.5)

cindex<-NULL
for (i in 1:5){
	rc<-rcorr.cens(-predict(f.mi.5,newdata=complete(mi,i)),f.mi.5$fits[[i]]$y)
	cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
	}
cindex
summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

z<-cox.zph(f.mi.5$fits[[1]])
z
plot(z)


# complete case
f<-cph(form,data=dat.mi)
f
anova(f)



###### BOOTSTRAP
opt<-NULL
slope<-NULL
intercept<-NULL
set.seed(0)
for (i in 1:5)
{
	v<-validate(f.mi.boot$fits[[i]],bw=TRUE,rule='p',sls=0.01,B=100,pr=FALSE,type='individual')
	opt<-c(opt,v["Dxy","optimism"]/2)
	slope<-c(slope,v["Slope","test"])
}
opt
mean(opt)
mean(slope)


### Additional value Mitosis?
Mit<-dat$Mitosis
lp<-f.mi$linear
f<-cph(S~offset(lp)+Mit)
f
summary(f,Mit=c("No"))


####### Nomogram
f.basehaz<-basehaz(f.mi,TRUE)
h0<-f.basehaz$hazard[f.basehaz$time==max(f.basehaz$time[f.basehaz$time<=5])]
mean(fun.event(lp))
#### Kaplan Meier 5-year recurrence
SF<-survfit(S~1,data=dat)
100*(1-SF$surv[SF$time==max(SF$time[SF$time<=5])])

#### Limits of Breslow
table(dat$Breslow)
Breslow.class<-c((2:10)/10,1+(1:5)/5,2.5,3:7)

nom.0<-nomogram(f.mi,maxscale=10,Breslow=Breslow.class)
# Baseline hazard
rc<-(1/attributes(nom.0)$info$sc)
int<-attributes(nom.0)$info$Intercept

nom<-nomogram(f.mi,maxscale=10,Breslow=Breslow.class,lp=FALSE,fun=fun.event,fun.at=fun.event(int+rc*(0:10)*2))
plot(nom)
nom.0
nom

points<-0:22
lp.points<-int+rc*points
print(data.frame(points,lp.points,fun.event(lp.points)))

points<-2*(0:10)
lp.points<-int+rc*points
print(data.frame('Total points'=points,'Predicted Value'=fun.event(lp.points)))




### Cross validation

Centers<-unique(dat$Center)
discrimination<-list()
cuts<-5


for (j in 1:4){
	S.j<-S[dat.mi$Center==Centers[j],]
	S.notj<-S[dat.mi$Center!=Centers[j],]
	form.notj<-update(form, S.notj  ~ . )
	f.mi.notj<-fit.mult.impute(form.notj,cph,xtrans=mi,data=dat.mi,n.impute=5,pr=FALSE,fit.reps=TRUE,y=TRUE,x=TRUE,se.fit=TRUE,sub=dat.mi$Center!=Centers[j])

	f.basehaz.notj<-basehaz(f.mi.notj,TRUE)
	h0<-f.basehaz.notj$hazard[f.basehaz.notj$time==max(f.basehaz.notj$time[f.basehaz.notj$time<=5])]

	cindex<-NULL
	lp<-NULL
	for (i in 1:5){
		lp.j<-predict(f.mi.notj,newdata=complete(mi,i)[dat.mi$Center==Centers[j],])
		rc<-rcorr.cens(-lp.j,S.j)
		cindex<-rbind(cindex,c(rc["C Index"],rc["S.D."]/2))
		lp<-cbind(lp,lp.j)
		}

	discrimination[[j]]<-summary(MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))

	lp.avg<-rowMeans(lp)
	cut.all<-cut2(lp.avg,g=cuts)
	tapply(S.j[,2],cut.all,length)

	S.km<-S.j
	CI.q<-tapply(1:length(lp.avg),cut.all,km)
	CI.q.l<-tapply(1:length(lp.avg),cut.all,km.upper)
	CI.q.u<-tapply(1:length(lp.avg),cut.all,km.lower)
	pred.q<-tapply(fun.event(lp.avg)/100,cut.all,mean)

	### Calibration plot
	lim<-c(0,.5)
win.graph(width = 6, height = 6)
	plot(pred.q,CI.q,xlab='Predicted probability',ylab='Observed frequency',xlim=lim,ylim=lim,cex=1,cex.axis=1,,cex.lab=1,lwd=1.5,main=paste("Calibration in",Centers[j]))
	lines(c(-1,1),c(-1,1),lwd=1.5)
	segments(pred.q,CI.q.l,pred.q,CI.q.u)
}

names(discrimination)<-Centers
discrimination



test1["Ulceration=Yes"]
test2["Ulceration=Yes"]
test3["Ulceration=Yes"]
test4["Ulceration=Yes"]
test5["Ulceration=Yes"]
test6["Ulceration=Yes"]
test7["Ulceration=Yes"]
test8["Ulceration=Yes"]
test9["Ulceration=Yes"]
test10["Ulceration=Yes"]
test11["Ulceration=Yes"]
test12["Ulceration=Yes"]
test13["Ulceration=Yes"]

sum1
sum2
sum3
sum4
sum5
sum6
sum7
sum8
sum9
sum10
sum11
sum12
sum13
