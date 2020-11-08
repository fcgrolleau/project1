library(survival)
library(rms)
library(mfp)
library(mice)
require(devtools) 
# install_git("https://github.com/BavoDC/CalibrationCurves") 
library(CalibrationCurves)
library(gplots)
library(survminer)
library(gtools)


load("rcts.RData")

rcts$sexe<-as.factor(rcts$sexe)
rcts$heart_failure<-as.factor(rcts$heart_failure)
rcts$hypertension<-as.factor(rcts$hypertension)
rcts$respiratory_disease<-as.factor(rcts$respiratory_disease)

rcts$aides<-as.factor(rcts$aides)
rcts$organ_graft<-as.factor(rcts$organ_graft)
rcts$immunosupressive_drug<-as.factor(rcts$immunosupressive_drug)
rcts$hemopathy<-as.factor(rcts$hemopathy)
rcts$type.eer<-as.factor(rcts$type.eer)
rcts$reprisediurese<-as.factor(rcts$reprisediurese)

rcts$diabetes[rcts$diabetes=="Aucun diabète OU diabète\ntraité par régime uniquement"]<-"NON"
rcts$diabetes[rcts$diabetes=="Diabète sans complication"]<-"NON"
rcts$diabetes[rcts$diabetes=="Diabète avec ou sans complications"]<-"OUI"
rcts$diabetes[rcts$diabetes=="Diabète compliqué (lésions\norganiques liées au diabète)"]<-"OUI"
rcts$diabetes<-as.factor(rcts$diabetes)

rcts$cirrhosis[rcts$cirrhosis=="Intermédiaire ou sévère"]<-"OUI"
rcts$cirrhosis[rcts$cirrhosis=="Modérée (sans hypertension\nportale, en incluant les\nhépatites chroniques)"]<-"OUI"
rcts$cirrhosis<-as.factor(rcts$cirrhosis)

rcts$cancer[rcts$cancer=="Oui, AVEC métastases"]<-"OUI"
rcts$cancer[rcts$cancer=="Oui, SANS métastases"]<-"OUI"
rcts$cancer<-as.factor(rcts$cancer)
rcts$anurie<-as.factor(rcts$kidney_sofa==4)

rcts$delaihospirea<-as.numeric(difftime(rcts$datehospirea, rcts$datehospi, units = "days"))

rcts$event<-"columninitiation"

for (i in 1:nrow(rcts)) {
        if(rcts$etatsortierea[i]=="Vivant") {
                rcts$event[i]<-"SORTIEm"
        }
}

for (i in 1:nrow(rcts)) {
        if(rcts$etatsortierea[i]=="Décédé") {
                rcts$event[i]<-"DCm"
        }
}

threshold_rrt_day<-2

for (i in 1:nrow(rcts)) {
        if(is.na(rcts$date.eer)[i]==FALSE & as.numeric(difftime(rcts$date.eer[i], rcts$date.rando[i], units = "days"))<=threshold_rrt_day) {
                rcts$event[i]<-"EERm"
        }
}


rcts$event<-as.factor(rcts$event)
with(rcts[rcts$bras=="STRATEGIE D ATTENTE",], table(etude, event))
rcts2<-rcts[,-which(names(rcts) %in% c("naiss", "datehospi", "datehospirea", "datehospi","kidney_sofa" ,"type.eer", "lastrrt", "vasofree"))]

# model developpement on pooled dataset

imp <- mice(rcts2, maxit=0)
predM <- imp$predictorMatrix
predM[,c("bras", "date.rando", "hospitalmortality", "ventilfree", "etude", "date.eer", "eer.initiation", "rrtfree")]<-0
meth <- imp$method
set.seed(123)
imp_pooled_tardif<-mice(rcts2[rcts2$bras=="STRATEGIE D ATTENTE",], m=5, seed = 1, predictorMatrix = predM, method = meth)
imp_pooled_precoce<-mice(rcts2[rcts2$bras=="STRATEGIE PRECOCE",], m=5, seed = 1, predictorMatrix = predM, method = meth)
save.image("imputedpooleddata2.RData")
load("imputedpooleddata2.RData")

## Selection in each imputed (We'll have to do AIC-like backward slection on pooled Wald tests with cut-off p=0.157)

var.test<-"age + sexe + pot_e + I(creat_e/creat_b) + anurie + sofa_e + weight_e + heart_failure + hypertension + diabetes + cirrhosis + immunosupressive_drug + urea_e + ph_e"
var.test<-strsplit(var.test,split=' + ', fixed=TRUE)[[1]]
uni.mod<-vector("list", length(var.test))

for (i in 1:length(var.test)) {
uni.mod[[i]]<-summary(pool(with(imp_pooled_tardif, glm(reformulate(var.test[i], "I(event=='EERm')"), family="binomial"))))
}

uni.mod<-data.frame(t(sapply(uni.mod, function(x) c(as.character(x[2,1]), x[2,2],x[2,3], x[2,6]))))
colnames(uni.mod)<-c("term","coef","std.error", "pval")
uni.mod[,2]<-as.numeric(uni.mod[,2])
uni.mod[,3]<-as.numeric(uni.mod[,3])
uni.mod[,4]<-as.numeric(uni.mod[,4])

#HR for increase in 0.01 in ph_e
uni.mod$lb<-uni.mod[,2]-qnorm(.975)*uni.mod[,3]
uni.mod$ub<-uni.mod[,2]+qnorm(.975)*uni.mod[,3]
uni.mod[uni.mod$term=="ph_e",-c(1,4)]<-uni.mod[uni.mod$term=="ph_e",-c(1,4)]/100
uni.mod
uni.table<-cbind(uni.mod[,1], round(exp(uni.mod[,c(2,5,6)]),3), round(uni.mod[,4],3))
colnames(uni.table)<-c("Variable", "HR", "lb", "ub", "p-value")
uni.table
write.csv(uni.table, "uni.table.csv")

fullmodel<-summary(pool(with(imp_pooled_tardif,
                             glm(I(event=="EERm")
                                 ~ age + sexe + pot_e + I(creat_e/creat_b)+anurie + sofa_e + weight_e + heart_failure + hypertension + diabetes + cirrhosis + immunosupressive_drug +
                                         urea_e + ph_e, family="binomial"))))


fullmodel
full.table<-cbind(
#as.character(fullmodel[-c(1,nrow(fullmodel)),1]),
exp(fullmodel[-c(1,nrow(fullmodel)),2]),
exp(fullmodel[-c(1,nrow(fullmodel)),2]-qnorm(.975)*fullmodel[-c(1,nrow(fullmodel)),3]),
exp(fullmodel[-c(1,nrow(fullmodel)),2]+qnorm(.975)*fullmodel[-c(1,nrow(fullmodel)),3]),
fullmodel[-c(1,nrow(fullmodel)),6]
)

#HR for increase in 0.01 in ph_e
full.table<-rbind(full.table,
c( #as.character(fullmodel[nrow(fullmodel),1]),
  exp(c(fullmodel[nrow(fullmodel),2],
  fullmodel[nrow(fullmodel),2]-qnorm(.975)*fullmodel[nrow(fullmodel),3],
  fullmodel[nrow(fullmodel),2]+qnorm(.975)*fullmodel[nrow(fullmodel),3])/100),
  fullmodel[nrow(fullmodel),6]))

colnames(full.table)<-c("HR", "lb", "ub", "p-value")
full.table<-cbind(uni.table[,1], round(full.table,3))
full.table
write.csv(full.table, "full.table.csv")

imp1bis<-complete(imp_pooled_tardif,1)
library(gam)
par(mfrow=c(3,2))
plot(gam(I(event=="EERm")~s(age), data=imp1bis))
plot(gam(I(event=="EERm")~s(pot_e), data=imp1bis))
plot(gam(I(event=="EERm")~s(ph_e), data=imp1bis))
plot(gam(I(event=="EERm")~s(sofa_e), data=imp1bis))
plot(gam(I(event=="EERm")~s(urea_e), data=imp1bis))
plot(gam(I(event=="EERm")~s(weight_e), data=imp1bis))

imp1bis$creatratio<-imp1bis$creat_e/imp1bis$creat_b
plot(gam(I(event=="EERm")~s(creatratio), data=imp1bis))
plot(gam(I(event=="EERm")~creatratio, data=imp1bis))
plot(gam(I(event=="EERm")~poly(creatratio,2), data=imp1bis))

selmod <- as.list(1: imp_pooled_tardif$m)
for (i in 1: imp_pooled_tardif$m) {
        data.i <- complete(imp_pooled_tardif, i)
        data.i$rrt4 <- as.numeric(data.i$event=="EERm")
        data.i$h_sofagreater3<-factor(data.i$hemodynamic_sofa>=3)
        data.i$r_sofagreater3<-factor(data.i$respi_sofa>=3)
        data.i$p_sofagreater2<-factor(data.i$platelet_sofa>=2)
        data.i$creatratiosquared<-(data.i$creat_e/data.i$creat_b)^2
        data.i$phsquared<-data.i$ph_e^2
        data.i$potsquared<-data.i$pot_e^2
        data.i$sofasquared<-data.i$sofa_e^2
        data.i$new.pot_e<-ifelse(data.i$pot_e>4,data.i$pot_e,4)
        data.i$new.ph_e<-ifelse(data.i$ph_e>7.1,data.i$ph_e,7.1)
        selmod[[i]] <- fastbw(lrm(I(event=="EERm")~ age + sexe + pot_e + I(creat_e/creat_b)+anurie + sofa_e + weight_e + heart_failure + hypertension + diabetes + cirrhosis + immunosupressive_drug +
        urea_e + ph_e, data=data.i),rule="p",type="individual",sls=0.157)
}

selvar <- lapply(selmod,function(x){x$factors.kept})
matselvar <- matrix(0, nrow=5, ncol=length(fullmodel$term)-1)
colnames(matselvar) <- as.character(fullmodel$term)[-1]

for(i in 1: imp_pooled_tardif$m) { matselvar[i, selvar[[i]]] <- 1 }

apply(matselvar,2,sum)*20 # Percent selection of each variable

       
# final model
finalmodel<-summary(pool(with(imp_pooled_tardif, glm(I(event=="EERm") ~  pot_e + sofa_e + weight_e + immunosupressive_drug +
                                                             urea_e + ph_e, family="binomial"))))
finalmodel
final.table<-cbind(
        #as.character(fullmodel[-c(1,nrow(fullmodel)),1]),
        exp(finalmodel[-c(1,nrow(finalmodel)),2]),
        exp(finalmodel[-c(1,nrow(finalmodel)),2]-qnorm(.975)*finalmodel[-c(1,nrow(finalmodel)),3]),
        exp(finalmodel[-c(1,nrow(finalmodel)),2]+qnorm(.975)*finalmodel[-c(1,nrow(finalmodel)),3]),
        finalmodel[-c(1,nrow(finalmodel)),6]
)

#HR for increase in 0.01 in ph_e
final.table<-rbind(final.table,
                  c( #as.character(fullmodel[nrow(fullmodel),1]),
                          exp(c(finalmodel[nrow(finalmodel),2],
                                finalmodel[nrow(finalmodel),2]-qnorm(.975)*finalmodel[nrow(finalmodel),3],
                                finalmodel[nrow(finalmodel),2]+qnorm(.975)*finalmodel[nrow(finalmodel),3])/100),
                          finalmodel[nrow(finalmodel),6]))

final.table<-data.frame(as.character(finalmodel[-1,1]), round(final.table,3), stringsAsFactors=TRUE)
colnames(final.table)<-c("variable", "HR", "lb", "ub", "p-value")
final.table
write.csv(final.table, "final.table.csv")

### Validation
validationbw <- as.list(1: imp_pooled_tardif$m)
for (i in 1:imp_pooled_tardif$m) {
        data.i <- complete(imp_pooled_tardif, i)
        data.i$rrt4 <- as.numeric(data.i$event=="EERm")
        data.i$h_sofagreater3<-factor(data.i$hemodynamic_sofa>=3)
        data.i$r_sofagreater3<-factor(data.i$respi_sofa>=3)
        data.i$p_sofagreater2<-factor(data.i$platelet_sofa>=2)
        data.i$creatratiosquared<-(data.i$creat_e/data.i$creat_b)^2
        data.i$phsquared<-data.i$ph_e^2
        data.i$potsquared<-data.i$pot_e^2
        data.i$sofasquared<-data.i$sofa_e^2
        data.i$new.pot_e<-ifelse(data.i$pot_e>4,data.i$pot_e,4)
        data.i$new.ph_e<-ifelse(data.i$ph_e>7.1,data.i$ph_e,7.1)
        validationbw[[i]] <- validate(lrm(I(event=="EERm")~ age + sexe + pot_e + I(creat_e/creat_b)+anurie + sofa_e + weight_e + heart_failure + hypertension + diabetes + cirrhosis + immunosupressive_drug + urea_e + ph_e, data=data.i, x=T, y=T),
                                      bw=TRUE, rule='p', sls=0.157, type= 'individual',B=200)
}

cindex<-matrix(0,2,imp_pooled_tardif$m)
rownames(cindex)<-c("apparent c-index", "corrected c-index")
for (i in 1:imp_pooled_tardif$m) {
        cindex[1,i]<-.5*(validationbw[[i]][1,1]+1)
        cindex[2,i]<-.5*(validationbw[[i]][1,5]+1)
}
cindex
apply(cindex, 1, mean)

calcurves<-matrix(0,2,imp_pooled_tardif$m)
rownames(calcurves)<-c("corrected intercept", "corrected slope")
for (i in 1:imp_pooled_tardif$m) {
        calcurves[1,i]<-validationbw[[i]][3,5]
        calcurves[2,i]<-validationbw[[i]][4,5]
}
calcurves
apply(calcurves, 1, mean)

###

calibrationcurves<-as.list(1:imp_pooled_tardif$m)
for (i in 1: imp_pooled_tardif$m) {
        data.i <- complete(imp_pooled_tardif, i)
        data.i$rrt4 <- as.numeric(data.i$event=="EERm")
        data.i$h_sofagreater3<-factor(data.i$hemodynamic_sofa>=3)
        data.i$r_sofagreater3<-factor(data.i$respi_sofa>=3)
        data.i$p_sofagreater2<-factor(data.i$platelet_sofa>=2)
        data.i$creatratiosquared<-(data.i$creat_e/data.i$creat_b)^2
        data.i$phsquared<- (data.i$ph_e)*(data.i$ph_e)
        data.i$sofasquared<- (data.i$sofa_e)*(data.i$sofa_e)
        calibrationcurves[[i]] <- calibrate(lrm(I(event=="EERm")~ age + sexe + pot_e + I(creat_e/creat_b)+anurie + sofa_e + weight_e + heart_failure + hypertension + diabetes + cirrhosis + immunosupressive_drug + urea_e + ph_e, data=data.i, x=T, y=T),
                                      bw=TRUE, rule='p', sls=0.157, type= 'individual',B=200)
        
        # calibrationcurves[[i]]<-calibrate(lrm(I(event=="EERm") ~ pot_e + sofa_e + weight_e + immunosupressive_drug +
        #                                              urea_e + ph_e, data=data.i, x=T,y=T), B=200)
}

par(mfrow=c(2,3))
for (i in 1:imp_pooled_tardif$m) {
        plot(calibrationcurves[[i]])
}

## Internal Calibration plot (contains all usefull information)
par(mfrow=c(2,3))
expit <- function(x){exp(x)/(1+exp(x))}

for (i in 1:imp_pooled_tardif$m) {
        fit1 <- glm(I(event=="EERm")~ pot_e + sofa_e + weight_e + immunosupressive_drug +
                            urea_e + ph_e, family="binomial", data=complete(imp_pooled_tardif,i))
        Pred1 <- predict(fit1)
        Probs1<-expit(Pred1)
        val.prob.ci.2(Probs1, complete(imp_pooled_tardif,i)$event=="EERm", col.ideal="blue", smooth = FALSE, statloc=FALSE, legendloc=FALSE)
        temp<-Probs1
        lines(calibrationcurves[[i]][,c(1,2)], lty=2, col="black")
        lines(calibrationcurves[[i]][,c(1,3)], lty=1, col="black")
        legend(.7,.25, legend=c("Ideal", "Apparent", "Bias-corrected"),
               col=c("blue", "black", "black"), lty=c(1,2,1), cex=.8, seg.len=1.6, bg="transparent",x.intersp=.5, y.intersp=.8,
               box.lty=2, box.lwd=0)

        n <- 5
        plotCI(tapply(temp,quantcut(temp,seq(0,1,by=1/n)),mean), tapply(complete(imp_pooled_tardif,i)$event=="EERm", quantcut(temp,seq(0,1,by=1/n)),mean), ui=tapply(complete(imp_pooled_tardif,i)$event=="EERm", quantcut(temp,seq(0,1,by=1/n)),function(x){binom.test(sum(x==1), length(x))$conf.int[2]}), li=tapply(complete(imp_pooled_tardif,i)$event=="EERm", quantcut(temp,seq(0,1,by=1/n)),function(x){binom.test(sum(x==1), length(x))$conf.int[1]}), add=T, pch=18, gap=0, sfrac=0.005, col="red")
        seqpos<-seq(0,1,by=.05)
        textpos<-rev(seqpos[(length(seqpos)-6):length(seqpos)])
        space<-.03
        text(0,textpos[1], cex=.8, pos=4, paste("Calibration"))
        text(0,textpos[2], cex=.8, pos=4, paste("   biais-corrected intercept:", format(round(calcurves[1,i],2), nsmall=2) ))
        text(0,textpos[3], cex=.8, pos=4, paste("   biais-corrected slope:", format(round(calcurves[2,i],2), nsmall=2) ))
        text(0,textpos[4]-space, cex=.8, pos=4, paste("Discrimination"))
        text(0,textpos[5]-space, cex=.8, pos=4, paste("   apparent c-statistic:", format(round(cindex[1,i],2), nsmall=2) ))
        text(0,textpos[6]-space, cex=.8, pos=4, paste("   biais-corrected c-statistic:", format(round(cindex[2,i],2), nsmall=2) ))
        box()
        }
#dev.copy2pdf(file="calibrationplotq5.pdf")

### Treatment effect across quantiles

## model predictions for akiki dataset
finalmodelcoefs<-finalmodel[,"estimate"]

# predictions tardif
tempprobs_tardif<-as.list(1:imp_pooled_tardif$m)
for (i in 1:imp_pooled_tardif$m) {
        X <- model.matrix(~ pot_e + sofa_e + weight_e + immunosupressive_drug +
                                  urea_e + ph_e, data= complete(imp_pooled_tardif,i))
        tempprobs_tardif[[i]]<-expit(X %*% finalmodelcoefs)
        tempprobs_tardif<-sapply(tempprobs_tardif, c)
}
probs_tardif<-apply(tempprobs_tardif, 1, mean)

# predictions precoce
tempprobs_precoce<-as.list(1:imp_pooled_precoce$m)
for (i in 1:imp_pooled_tardif$m) {
        X <- model.matrix(~  pot_e + sofa_e + weight_e + immunosupressive_drug +
                                  urea_e + ph_e, data= complete(imp_pooled_precoce,i))
        tempprobs_precoce[[i]]<-expit(X %*% finalmodelcoefs)
        tempprobs_precoce<-sapply(tempprobs_precoce, c)
}
probs_precoce<-apply(tempprobs_precoce, 1, mean)

# put prababilities in the first completed dataset (non imputations of the primary outcome)
pooled_po_dataset<-rbind(
        cbind(complete(imp_pooled_tardif,1), probs=probs_tardif),
        cbind(complete(imp_pooled_precoce,1), probs=probs_precoce)
)

overallcox<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset)
HR_overall<-c(
        exp(coef(overallcox)),
        exp(confint(overallcox))
)

overallcox_akiki<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset, subset=etude=="akiki")
HR_overall_akiki<-c(
  exp(coef(overallcox_akiki)),
  exp(confint(overallcox_akiki))
)

overallcox_idealicu<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset, subset=etude=="idealicu")
HR_overall_idealicu<-c(
  exp(coef(overallcox_idealicu)),
  exp(confint(overallcox_idealicu))
)

survfit.obj<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset)
splots<-ggsurvplot(survfit.obj, data=pooled_po_dataset,
                        ggtheme = theme_survminer(),
                        title    = "Pooled AKIKI and IDEAL-ICU",
                        legend.title = "Strategy",
                        pval=FALSE, pval.method = F,
                        risk.table = FALSE, 
                        break.time.by = 7,
                        tables.theme = theme_cleantable(), 
                        tables.y.text = F)

splots

survfit.obj.akiki<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="akiki", ])
splots<-ggsurvplot(survfit.obj.akiki, data=pooled_po_dataset[pooled_po_dataset$etude=="akiki", ],
                   ggtheme = theme_survminer(),
                   title    = "AKIKI",
                   legend.title = "Strategy",
                   pval=FALSE, pval.method = F,
                   risk.table = FALSE, 
                   break.time.by = 7,
                   tables.theme = theme_cleantable(), 
                   tables.y.text = F)

splots

survfit.obj.idealicu<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu", ])
splots<-ggsurvplot(survfit.obj.idealicu, data=pooled_po_dataset[pooled_po_dataset$etude=="idealicu", ],
                   ggtheme = theme_survminer(),
                   title    = "IDEAL-ICU",
                   legend.title = "Strategy",
                   pval=FALSE, pval.method = F,
                   risk.table = FALSE, 
                   break.time.by = 7,
                   tables.theme = theme_cleantable(), 
                   tables.y.text = F)

splots

ARD_overalltable<-c(
1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="A"),1)[2],
1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="B"),1)[2],
(1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="B"),1)[2])-(1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="A"),1)[2]),
survfit.obj$n[1], survfit.obj$n[2] )
names(ARD_overalltable)<-c("death rate in late", "death rate in early", "overall ARD", "n late", "n early")
ARD_overalltable
ardoverallvar<-function(x) {(x[1]*(1-x[1])/x[4]) + (x[2]*(1-x[2])/x[5])}

ARD_overall<-c(
ARD_overalltable[3],
ARD_overalltable[3]+c(-1,1)*qnorm(.975)*sqrt(ardoverallvar(ARD_overalltable))
)

### akiki subanalysis ###
ARD_overalltable_akiki<-c(
  1-tail(weightedsurvivaltables(survfit.obj.akiki$time, survfit.obj.akiki$surv, group="A"),1)[2],
  1-tail(weightedsurvivaltables(survfit.obj.akiki$time, survfit.obj.akiki$surv, group="B"),1)[2],
  (1-tail(weightedsurvivaltables(survfit.obj.akiki$time, survfit.obj.akiki$surv, group="B"),1)[2])-(1-tail(weightedsurvivaltables(survfit.obj.akiki$time, survfit.obj.akiki$surv, group="A"),1)[2]),
  survfit.obj.akiki$n[1], survfit.obj.akiki$n[2] )
names(ARD_overalltable_akiki)<-c("death rate in late", "death rate in early", "overall ARD", "n late", "n early")
ARD_overalltable_akiki
ardoverallvar<-function(x) {(x[1]*(1-x[1])/x[4]) + (x[2]*(1-x[2])/x[5])}

ARD_overall_akiki<-c(
  ARD_overalltable_akiki[3],
  ARD_overalltable_akiki[3]+c(-1,1)*qnorm(.975)*sqrt(ardoverallvar(ARD_overalltable_akiki))
)

### ideal-icu subanalysis ###
ARD_overalltable_idealicu<-c(
  1-tail(weightedsurvivaltables(survfit.obj.idealicu$time, survfit.obj.idealicu$surv, group="A"),1)[2],
  1-tail(weightedsurvivaltables(survfit.obj.idealicu$time, survfit.obj.idealicu$surv, group="B"),1)[2],
  (1-tail(weightedsurvivaltables(survfit.obj.idealicu$time, survfit.obj.idealicu$surv, group="B"),1)[2])-(1-tail(weightedsurvivaltables(survfit.obj.idealicu$time, survfit.obj.idealicu$surv, group="A"),1)[2]),
  survfit.obj.idealicu$n[1], survfit.obj.idealicu$n[2] )
names(ARD_overalltable_idealicu)<-c("death rate in late", "death rate in early", "overall ARD", "n late", "n early")
ARD_overalltable_idealicu
ardoverallvar<-function(x) {(x[1]*(1-x[1])/x[4]) + (x[2]*(1-x[2])/x[5])}

ARD_overall_idealicu<-c(
  ARD_overalltable_idealicu[3],
  ARD_overalltable_idealicu[3]+c(-1,1)*qnorm(.975)*sqrt(ardoverallvar(ARD_overalltable_idealicu))
)
###

RR_overall<-ARD_overalltable[2]/ARD_overalltable[1]
seRR_overall<-sqrt(
        ((ARD_overalltable[4]-ARD_overalltable[4]*ARD_overalltable[1])/ARD_overalltable[4]*ARD_overalltable[1])/ARD_overalltable[4] +
        ((ARD_overalltable[5]-ARD_overalltable[5]*ARD_overalltable[2])/ARD_overalltable[5]*ARD_overalltable[2])/ARD_overalltable[5]        
)

RR_overall<-c(RR_overall, RR_overall+c(-1,1)*qnorm(.975)*seRR_overall)

HRs<-list()
par(mfrow=c(2,4))
for (i in 3:10) {
        n<-i # number of quantiles
        pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        HRs[[i]]<-exp(sapply(levels(pooled_po_dataset$quantile), function(x) {mod<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset, subset=quantile==x)
        c(coef(mod), confint(mod))}))
        HRs
        plotCI(1:n, HRs[[i]][1,], ui=HRs[[i]][3,], li=HRs[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="Hazard Ratio")
        abline(h=1)
}

q<-5

pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/q))
quantilemean<-tapply(pooled_po_dataset$probs, pooled_po_dataset$quantile, function(x) summary(x)[4])
quantilecuts<-as.numeric(substr(names(quantilemean), start=2, stop=5)[-1])

########## HR penalized smoothing ##########
library(mgcv)
nfittp = gam(derniere.nouvelles.censureJ60~bras+s(probs,by=bras, bs='ps'), 
             data=pooled_po_dataset, weights = etat.censureJ60, 
             family = 'cox.ph', method = "REML")

df.precoce  = expand.grid(bras = 'STRATEGIE PRECOCE', probs = seq(0.02468, .55, by=.005))
df.tardif = expand.grid(bras = 'STRATEGIE D ATTENTE', probs = seq(0.02468, .55, by=.005))

y0.precoce_ = predict(nfittp,newdata=df.precoce,type="lpmatrix")
y0.tardif_ = predict(nfittp,newdata=df.tardif,type="lpmatrix")

y0.precoce = predict(nfittp,newdata=df.precoce,type="link", se.fit=TRUE)
y0.tardif = predict(nfittp,newdata=df.tardif,type="link", se.fit=TRUE)

dplot = data.frame(
        probs = seq(0.02468, .55, by=.005),
        b.precoce = y0.precoce$fit - y0.tardif$fit,
        s.precoce = sapply(1:length(seq(0.02468, .55, by=.005)), function(i){
                l = matrix(y0.precoce_[i,] - y0.tardif_[i,])
                t(l) %*% vcov(nfittp) %*% l
        })
)


par(mfrow=c(1,2))
plot(seq(0.02468, .55, by=.005), exp(y0.tardif$fit), type='l', log='y', lwd=2, xlab="Probs", ylab="Relative Risk")
points(seq(0.02468, .55, by=.005), exp(y0.tardif$fit + qnorm(.975) * y0.tardif$se.fit), type='l', col=1, lwd=2, lty=2)
points(seq(0.02468, .55, by=.005), exp(y0.tardif$fit - qnorm(.975) * y0.tardif$se.fit), type='l', col=1, lwd=2, lty=2)
points(seq(0.02468, .55, by=.005), exp(y0.precoce$fit), type='l', col='red', lwd=2)
points(seq(0.02468, .55, by=.005), exp(y0.precoce$fit + qnorm(.975) * y0.precoce$se.fit), type='l', col=2, lwd=2, lty=2)
points(seq(0.02468, .55, by=.005), exp(y0.precoce$fit - qnorm(.975) * y0.precoce$se.fit), type='l', col=2, lwd=2, lty=2)
legend(.7, .8, c("Late", "Early"), lty=1, lwd=2, col=1:2, bty='n')
abline(h=1)


plot(seq(0.02468, .55, by=.005), exp(dplot$b.precoce), type = 'l', ylab = 'Hazard Ratio', xlab = 'Probs', lwd=2, ylim=c(.5,1.3))
points(seq(0.02468, .55, by=.005), exp(dplot$b.precoce - qnorm(.975) * sqrt(dplot$s.precoce)), type= 'l', lty = 2, lwd=2)
points(seq(0.02468, .55, by=.005), exp(dplot$b.precoce + qnorm(.975) * sqrt(dplot$s.precoce)), type= 'l', lty = 2, lwd=2)
abline(h=1)

polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c(exp(dplot$b.precoce + qnorm(.975) * sqrt(dplot$s.precoce)),
                                                                          rev(exp(dplot$b.precoce - qnorm(.975) * sqrt(dplot$s.precoce)))),
        col= rgb(1, 0, 0,0.5))


#### Penalized splines (doesn't work for now..) struggle to extract relevant basis from pspline functions
cent<-.3
probs.ps.p <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE PRECOCE", probs, cent))
probs.ps.t <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE D ATTENTE", probs, cent))

degf<-3
nt<-trunc(2.5*degf)+1
fit.ps.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + pspline(probs.ps.p, nterm=nt, df=degf, penalty=TRUE) + pspline(probs.ps.t, nterm=nt, df=degf, penalty=TRUE), data = pooled_po_dataset, x=TRUE)

ypred<-predict(fit.ps.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.ps.p=seq(0.02468, .55, by=.005), probs.ps.t=cent), se=TRUE)
ypred2<-predict(fit.ps.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.ps.p=cent, probs.ps.t=seq(0.02468, .55, by=.005)), se=TRUE)

yy <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2<- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')
matplot(seq(0.02468, .55, by=.005), exp(matrix(cbind(yy,yy2), ncol=6)), type='l', lty=c(1,2,2), col=c(1,1,1,2,2,2),
        lwd=2, log='y',
        xlab="Probs", ylab="Relative risk")


y0.precoce_<-data.frame(bras = 1, cbind(as.matrix(pspline(seq(0.02468, .55, by=.005), nterm=nt, df=degf, penalty=FALSE)), as.matrix(pspline(rep(cent,106), nterm=nt, df=degf, penalty=FALSE))))
y0.tardif_<-data.frame(bras = 0, cbind(as.matrix(pspline(rep(cent,106), nterm=nt, df=degf, penalty=FALSE)), as.matrix(pspline(seq(0.02468, .55, by=.005), nterm=nt, df=degf, penalty=FALSE))))

plot(seq(0.02468, .55, by=.005),exp(
        as.matrix(y0.precoce_) %*% coef(fit.ps.int)), type="l")
lines(seq(0.02468, .55, by=.005), exp(as.matrix(y0.tardif_) %*% coef(fit.ps.int)), col="red")


l<-sapply(1:length(seq(0.02468, .55, by=.005)), function(i){
        matrix(y0.precoce_[i,] - y0.tardif_[i,]) })

l<-matrix(as.numeric(l), nrow(l), 106)

fits<-as.numeric(t(l) %*% coef(fit.ps.int))
se_s<-sapply(1:length(seq(0.02468, .55, by=.005)), function(i){
        l = as.numeric(y0.precoce_[i,] - y0.tardif_[i,])
        t(l) %*% vcov(fit.ps.int) %*% l
})

yy <- fits + outer(sqrt(se_s), c(0, -qnorm(.975), qnorm(.975)), '*')

matplot(seq(0.02468, .55, by=.005), exp(matrix(yy, ncol=3)), type='l', lty=c(1,2,2),
        lwd=2, col=1, log='y',
        xlab="Probs", ylab="Relative risk")
abline(a=0, b=0)


##### RCS one knot at the median (-> changed for two knots)
probs.int.p <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset$probs-.25

dfs<-3
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)
library(splines)
fit.sep<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.int.p, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset, x=TRUE)

fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.c, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset, x=TRUE)

ypred<-predict(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0), se=TRUE)
ypred2<-predict(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25), se=TRUE)
yy <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2<- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')
matplot(seq(0.02468, .55, by=.005), exp(matrix(cbind(yy,yy2), ncol=6)), type='l', lty=c(1,2,2), col=c(1,1,1,2,2,2),
        lwd=2, log='y',
        xlab="Probs", ylab="Relative risk")


y0.precoce_<-data.frame(bras = 1, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
y0.tardif_<-data.frame(bras = 0, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg))


l<-sapply(1:length(seq(0.02468, .55, by=.005)), function(i){
        matrix(y0.precoce_[i,] - y0.tardif_[i,]) })

l<-matrix(as.numeric(l), nrow(l), 106)

fits<-as.numeric(t(l) %*% coef(fit.int))
se_s<-sapply(1:length(seq(0.02468, .55, by=.005)), function(i){
        l = as.numeric(y0.precoce_[i,] - y0.tardif_[i,])
        t(l) %*% vcov(fit.int) %*% l
})

yy <- fits + outer(sqrt(se_s), c(0, -qnorm(.975), qnorm(.975)), '*')

matplot(seq(0.02468, .55, by=.005), exp(matrix(yy, ncol=3)), type='l', lty=c(1,2,2),
        lwd=2, col=1, log='y',
        xlab="Probs", ylab="Relative risk")
abline(a=0, b=0)


## HR plot
dev.off()
plot(1, xlim=c(0, .6), ylim = c(.4,2.5), yaxt="n", bty="n", type="n", log="y", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Hazard Ratio", main="60 Days Survival")
axis(2, at=c(1,seq(0.5,2.5,by=0.5)),las=1) #, labels=NA)
lines(c(0,.6), y=c(1,1), lwd=2)
segments(0,HR_overall[1], .6, HR_overall[1], lty=2, lwd=2)
plotCI(quantilemean, HRs[[q]][1,], ui=HRs[[q]][3,], li=HRs[[q]][2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", add=TRUE, las=1, bty="n")
#plotCI(quantilemean, HRs[[q]][1,], ui=HRs[[q]][3,], li=HRs[[q]][2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Hazard Ratio", main="60 Days Survival", xlim=c(0, .6), ylim = c(.4,2.3))
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.45, quantilecuts[i], 2.2, lty=3)
}
segments(0,.45, 0, 2.2, lty=3)
text(quantilemean, 2.3, paste("Q", 1:q, sep=""), cex=.8)
arrows(.6, c(1.05,1/1.05), .6, c(1.05+.4, 1/(1.05+.4)), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=1.05, cex=1, adj=0)
mtext(text="favors early", side=4,line=-0, at=2-1.05, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), exp(yy[,1]), type = 'l', lwd=1.2,col="#6A6599FF")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c(exp(yy[,2]), rev(exp(yy[,3]))),
        col= rgb(0, .1,.5,0.15), border=NA)
#dev.copy2pdf(file="hrplot.pdf")

#points(seq(0.02468, .55, by=.005), exp(dplot$b.precoce), type = 'l', lwd=1.2, col="#6A6599FF")
#polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c(exp(dplot$b.precoce + qnorm(.975) * sqrt(dplot$s.precoce)),
#                                                                          rev(exp(dplot$b.precoce - qnorm(.975) * sqrt(dplot$s.precoce)))),
#                                                                col= rgb(0, .1,.5,0.15), border=NA)




pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/q))
survfit.obj<-list()
splots<-list()
for (j in 1:q) {
        survfit.obj[[j]]<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$quantile==colnames(HRs[[q]])[j],])
        splots[[j]]<-ggsurvplot(survfit.obj[[j]], data=pooled_po_dataset[pooled_po_dataset$quantile==colnames(HRs[[q]])[j],],
                                ggtheme = theme_survminer() + theme(plot.title = element_text(hjust = 0.5)),
                                title    = paste("Q",j, " ",colnames(HRs[[q]])[j], sep=""),
                                font.title=12,
                                legend.title = "",
                                legend.labs = c("Delayed strategy", "Early strategy"),
                                pval=FALSE, pval.method = F,
                                risk.table = TRUE,
                                risk.table.fontsize = 5,
                                break.time.by = 7,
                                tables.theme = theme_cleantable(), 
                                tables.y.text = F)

        splots[[j]]$table <- splots[[j]]$table + labs(title = "", subtitle = "No. at risk")
        splots[[j]]<- splots[[j]] + labs(x  = "Days", y = "Overall Survival")       
}

kmplots_pooled<-arrange_ggsurvplots(splots, ncol = q, nrow = 1, print = TRUE)
ggsave("km_pooled.pdf", kmplots_pooled, width = 26, height = 7)

### akiki sub analysis
survfit.obj.akiki<-list()
splots<-list()
for (j in 1:q) {
  survfit.obj.akiki[[j]]<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="akiki" & pooled_po_dataset$quantile==colnames(HRs[[q]])[j],])
        splots[[j]]<-ggsurvplot(survfit.obj.akiki[[j]], data=pooled_po_dataset[pooled_po_dataset$etude=="akiki" & pooled_po_dataset$quantile==colnames(HRs[[q]])[j],],
                                ggtheme = theme_survminer() + theme(plot.title = element_text(hjust = 0.5)),
                                title    = paste("Q",j, " ",colnames(HRs[[q]])[j], sep=""),
                                font.title=12,
                                legend.title = "",
                                legend.labs = c("Delayed strategy", "Early strategy"),
                                pval=FALSE, pval.method = F,
                                risk.table = TRUE,
                                risk.table.fontsize = 5,
                                break.time.by = 7,
                                tables.theme = theme_cleantable(), 
                                tables.y.text = F)
        
        splots[[j]]$table <- splots[[j]]$table + labs(title = "", subtitle = "No. at risk")
        splots[[j]]<- splots[[j]] + labs(x  = "Days", y = "Overall Survival")      
}

kmplots_akiki<-arrange_ggsurvplots(splots, ncol = q, nrow = 1, print = TRUE)
ggsave("km_akiki.pdf", kmplots_akiki, width = 26, height = 7)


### ideal-icu sub analysis
survfit.obj.idealicu<-list()
splots<-list()
for (j in 1:q) {
  survfit.obj.idealicu[[j]]<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu" & pooled_po_dataset$quantile==colnames(HRs[[q]])[j],])
        splots[[j]]<-ggsurvplot(survfit.obj.idealicu[[j]], data=pooled_po_dataset[pooled_po_dataset$etude=="idealicu" & pooled_po_dataset$quantile==colnames(HRs[[q]])[j],],
                                ggtheme = theme_survminer() + theme(plot.title = element_text(hjust = 0.5)),
                                title    = paste("Q",j, " ",colnames(HRs[[q]])[j], sep=""),
                                font.title=12,
                                legend.title = "",
                                legend.labs = c("Delayed strategy", "Early strategy"),
                                pval=FALSE, pval.method = F,
                                risk.table = TRUE,
                                risk.table.fontsize = 5,
                                break.time.by = 7,
                                tables.theme = theme_cleantable(), 
                                tables.y.text = F)
        
        splots[[j]]$table <- splots[[j]]$table + labs(title = "", subtitle = "No. at risk")
        splots[[j]]<- splots[[j]] + labs(x  = "Days", y = "Overall Survival")       
}

kmplots_idealicu<-arrange_ggsurvplots(splots, ncol = q, nrow = 1, print = TRUE)
ggsave("km_idealicu.pdf", kmplots_idealicu, width = 26, height = 7)


### akiki sub analysis

HRs_akiki<-list()
par(mfrow=c(2,4))
for (i in 3:10) {
  n<-i # number of quantiles
  pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
  
  HRs_akiki[[i]]<-exp(sapply(levels(pooled_po_dataset$quantile), function(x) {mod<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="akiki",], subset=quantile==x)
  c(coef(mod), confint(mod))}))
  HRs_akiki
  plotCI(1:n, HRs_akiki[[i]][1,], ui=HRs_akiki[[i]][3,], li=HRs_akiki[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="Hazard Ratio")
  abline(h=1)
}


### idealicu sub analysis

HRs_idealicu<-list()
par(mfrow=c(2,4))
for (i in 3:10) {
  n<-i # number of quantiles
  pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
  
  HRs_idealicu[[i]]<-exp(sapply(levels(pooled_po_dataset$quantile), function(x) {mod<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], subset=quantile==x)
  c(coef(mod), confint(mod))}))
  HRs_idealicu
  plotCI(1:n, HRs_idealicu[[i]][1,], ui=HRs_idealicu[[i]][3,], li=HRs_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="Hazard Ratio")
  abline(h=1)
}

q<-5


### investigate the first and last quantile
round(prop.table(with(pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], table(eer.initiation, quantile)),2)*100,2)
tapply(as.numeric(with(pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], difftime(date.eer, date.rando, units = "days"))), pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$quantile, function(x) 24*summary (x)[3]) 
tapply(as.numeric(with(pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], difftime(date.eer, date.rando, units = "days"))), pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$quantile, function(x) 24*summary (x)[3])

tapply(with(pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], sofa_e), pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$quantile, function(x) summary (x)[3])
tapply(with(pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], sofa_e), pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$quantile, function(x) summary (x)[3])

## ARD at day 60

ARDd60_val<-list()
for (j in 1:q) {
ARDd60_val[[j]]<-c((1-tail(weightedsurvivaltables(survfit.obj[[j]]$time, survfit.obj[[j]]$surv, group="B"),1)[2]),
                         (1-tail(weightedsurvivaltables(survfit.obj[[j]]$time, survfit.obj[[j]]$surv, group="A"),1)[2]),
                         (1-tail(weightedsurvivaltables(survfit.obj[[j]]$time, survfit.obj[[j]]$surv, group="B"),1)[2])-(1-tail(weightedsurvivaltables(survfit.obj[[j]]$time, survfit.obj[[j]]$surv, group="A"),1)[2])
)
}
ARDd60_val
ARDd60_table<-cbind(t(sapply(ARDd60_val,c)),table(pooled_po_dataset$quantile, pooled_po_dataset$bras))
colnames(ARDd60_table)[1:3]<-c("PROP STRATEGIE PRECOCE", "PROP STRATEGIE D ATTENTE", "ARD")
ARDd60_table

ARDd60_list<-as.list(as.data.frame(t(ARDd60_table)))
ardvar<-function(x) {(x[1]*(1-x[1])/x[5]) + (x[2]*(1-x[2])/x[4])}

ARDd60_ci<-rbind(
ARDd60_table[,"ARD"],
ARDd60_table[,"ARD"]-qnorm(.975)*sqrt(sapply(ARDd60_list, ardvar)),
ARDd60_table[,"ARD"]+qnorm(.975)*sqrt(sapply(ARDd60_list, ardvar))
)


## ARD at day 60 akiki sub analysis
q<-5
pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/q))

ARDd60_val_akiki<-list()
for (j in 1:q) {
  ARDd60_val_akiki[[j]]<-c((1-tail(weightedsurvivaltables(survfit.obj.akiki[[j]]$time, survfit.obj.akiki[[j]]$surv, group="B"),1)[2]),
                     (1-tail(weightedsurvivaltables(survfit.obj.akiki[[j]]$time, survfit.obj.akiki[[j]]$surv, group="A"),1)[2]),
                     (1-tail(weightedsurvivaltables(survfit.obj.akiki[[j]]$time, survfit.obj.akiki[[j]]$surv, group="B"),1)[2])-(1-tail(weightedsurvivaltables(survfit.obj.akiki[[j]]$time, survfit.obj.akiki[[j]]$surv, group="A"),1)[2])
  )
}
ARDd60_val_akiki
ARDd60_table_akiki<-cbind(t(sapply(ARDd60_val_akiki,c)),table(pooled_po_dataset$quantile, pooled_po_dataset$bras))
colnames(ARDd60_table_akiki)[1:3]<-c("PROP STRATEGIE PRECOCE", "PROP STRATEGIE D ATTENTE", "ARD")
ARDd60_table_akiki

ARDd60_list_akiki<-as.list(as.data.frame(t(ARDd60_table_akiki)))
ardvar<-function(x) {(x[1]*(1-x[1])/x[5]) + (x[2]*(1-x[2])/x[4])}

ARDd60_ci_akiki<-rbind(
  ARDd60_table_akiki[,"ARD"],
  ARDd60_table_akiki[,"ARD"]-qnorm(.975)*sqrt(sapply(ARDd60_list_akiki, ardvar)),
  ARDd60_table_akiki[,"ARD"]+qnorm(.975)*sqrt(sapply(ARDd60_list_akiki, ardvar))
)

## ARD at day 60 idealicu sub analysis

ARDd60_val_idealicu<-list()
for (j in 1:q) {
  ARDd60_val_idealicu[[j]]<-c((1-tail(weightedsurvivaltables(survfit.obj.idealicu[[j]]$time, survfit.obj.idealicu[[j]]$surv, group="B"),1)[2]),
                           (1-tail(weightedsurvivaltables(survfit.obj.idealicu[[j]]$time, survfit.obj.idealicu[[j]]$surv, group="A"),1)[2]),
                           (1-tail(weightedsurvivaltables(survfit.obj.idealicu[[j]]$time, survfit.obj.idealicu[[j]]$surv, group="B"),1)[2])-(1-tail(weightedsurvivaltables(survfit.obj.idealicu[[j]]$time, survfit.obj.idealicu[[j]]$surv, group="A"),1)[2])
  )
}
ARDd60_val_idealicu
ARDd60_table_idealicu<-cbind(t(sapply(ARDd60_val_idealicu,c)),table(pooled_po_dataset$quantile, pooled_po_dataset$bras))
colnames(ARDd60_table_idealicu)[1:3]<-c("PROP STRATEGIE PRECOCE", "PROP STRATEGIE D ATTENTE", "ARD")
ARDd60_table_idealicu

ARDd60_list_idealicu<-as.list(as.data.frame(t(ARDd60_table_idealicu)))
ardvar<-function(x) {(x[1]*(1-x[1])/x[5]) + (x[2]*(1-x[2])/x[4])}

ARDd60_ci_idealicu<-rbind(
  ARDd60_table_idealicu[,"ARD"],
  ARDd60_table_idealicu[,"ARD"]-qnorm(.975)*sqrt(sapply(ARDd60_list_idealicu, ardvar)),
  ARDd60_table_idealicu[,"ARD"]+qnorm(.975)*sqrt(sapply(ARDd60_list_idealicu, ardvar))
)

##### RCS two knots for ARD 
probs.int.p <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset$probs-.25

dfs<-4
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)

a1<-ns(probs.c, knots = kn, Boundary.knots=rg)[,1]
a2<-ns(probs.c, knots = kn, Boundary.knots=rg)[,2]
a3<-ns(probs.c, knots = kn, Boundary.knots=rg)[,3]
b1<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,1]
b2<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,2]
b3<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,3]



#fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+b1+b2, data = pooled_po_dataset, x=TRUE)
fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+a3+b1+b2+b3, data = pooled_po_dataset, x=TRUE)

y0.precoce_<-data.frame(bras = 1, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
y0.tardif_<-data.frame(bras = 0, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg))

y0.precoce_bis<-data.frame(bras="STRATEGIE PRECOCE", y0.precoce_[,-1])
#colnames(y0.precoce_bis)<-c("bras","a1","a2","b1","b2")
colnames(y0.precoce_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
precoced60d<-1-survfit(fit.int, newdata=y0.precoce_bis)$surv[56,]

y0.tardif_bis<-data.frame(bras="STRATEGIE D ATTENTE", y0.tardif_[,-1])
#colnames(y0.tardif_bis)<-c("bras","a1","a2","b1","b2")
colnames(y0.tardif_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
tardifd60d<-1-survfit(fit.int, newdata=y0.tardif_bis)$surv[56,]

se_s<-ardbootfunction(dataset=pooled_po_dataset, nboot = 1000)

plot(seq(0.02468, .55, by=.005), precoced60d-tardifd60d, type="l", ylim=c(-.2,.2))
abline(a=0,b=0)
lines(seq(0.02468, .55, by=.005), (precoced60d-tardifd60d)-qnorm(.975)*se_s, lty=2)
lines(seq(0.02468, .55, by=.005), (precoced60d-tardifd60d)+qnorm(.975)*se_s, lty=2)


### ARD plot
dev.off()
plotCI(quantilemean, ARDd60_ci[1,], ui=ARDd60_ci[3,], li=ARDd60_ci[2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Absolute Risk Difference", main="Death at Day 60", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
lines(c(0,.6), y=c(0,0))
text(quantilemean, 0.3, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-.35, quantilecuts[i], .27, lty=3)
}
segments(0,-.35, 0, .27, lty=3)
segments(0,ARD_overall[1], .6, ARD_overall[1], lty=2, lwd=2)
arrows(.6, c(0.02,-0.02), .6, c(0.02+.15, -0.02-.15), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=0.02, cex=1, adj=0)
mtext(text="favors early", side=4,line=-0, at=-0.02, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), precoced60d-tardifd60d, type = 'l', lwd=1.2,col="#6A6599FF")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c((precoced60d-tardifd60d)+qnorm(.975)*se_s,
                                                                           rev((precoced60d-tardifd60d)-qnorm(.975)*se_s)),
        col= rgb(0, .1,.5,0.15), border=NA)
#dev.copy2pdf(file="ardplotq5.pdf")




##### AKIKI subanalysis : RCS two knots for ARD
probs.int.p <- with(pooled_po_dataset[pooled_po_dataset$etude=="akiki",], ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset[pooled_po_dataset$etude=="akiki",], ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]$probs-.25

dfs<-4
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)

a1<-ns(probs.c, knots = kn, Boundary.knots=rg)[,1]
a2<-ns(probs.c, knots = kn, Boundary.knots=rg)[,2]
a3<-ns(probs.c, knots = kn, Boundary.knots=rg)[,3]
b1<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,1]
b2<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,2]
b3<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,3]



#fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+b1+b2, data = pooled_po_dataset[pooled_po_dataset$etude=="akiki",], x=TRUE)
fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+a3+b1+b2+b3, data = pooled_po_dataset[pooled_po_dataset$etude=="akiki",], x=TRUE)

y0.precoce_<-data.frame(bras = 1, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
y0.tardif_<-data.frame(bras = 0, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg))

y0.precoce_bis<-data.frame(bras="STRATEGIE PRECOCE", y0.precoce_[,-1])
#colnames(y0.precoce_bis)<-c("bras","a1","a2","b1","b2")
colnames(y0.precoce_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
precoced60d<-1-tail(survfit(fit.int, newdata=y0.precoce_bis)$surv,1)

y0.tardif_bis<-data.frame(bras="STRATEGIE D ATTENTE", y0.tardif_[,-1])
#colnames(y0.tardif_bis)<-c("bras","a1","a2","b1","b2")
colnames(y0.tardif_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
tardifd60d<-1-tail(survfit(fit.int, newdata=y0.tardif_bis)$surv,1)

# caution: function wroks only with two knots at position 0.13143840 & 0.32789328
# for more knots or different knots positions function adjustment is needed

se_s<-ardbootfunction(dataset=pooled_po_dataset[pooled_po_dataset$etude=="akiki",], nboot = 1000)

plot(seq(0.02468, .55, by=.005), precoced60d-tardifd60d, type="l", ylim=c(-.2,.2))
abline(a=0,b=0)
lines(seq(0.02468, .55, by=.005), (precoced60d-tardifd60d)-qnorm(.975)*se_s, lty=2)
lines(seq(0.02468, .55, by=.005), (precoced60d-tardifd60d)+qnorm(.975)*se_s, lty=2)


### ARD plot : AKIKI subanalysis : RCS two knots for ARD
dev.off()
plotCI(quantilemean, ARDd60_ci_akiki[1,], ui=ARDd60_ci_akiki[3,], li=ARDd60_ci_akiki[2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Absolute Risk Difference", main="Death at Day 60 in AKIKI", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
lines(c(0,.6), y=c(0,0))
text(quantilemean, 0.3, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
  segments(quantilecuts[i],-.35, quantilecuts[i], .27, lty=3)
}
segments(0,-.35, 0, .27, lty=3)
segments(0,ARD_overall_akiki[1], .6, ARD_overall_akiki[1], lty=2, lwd=2)
arrows(.6, c(0.02,-0.02), .6, c(0.02+.15, -0.02-.15), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=0.02, cex=1, adj=0)
mtext(text="favors early", side=4,line=-0, at=-0.02, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), precoced60d-tardifd60d, type = 'l', lwd=1.2,col="#6A6599FF")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c((precoced60d-tardifd60d)+qnorm(.975)*se_s,
                                                                          rev((precoced60d-tardifd60d)-qnorm(.975)*se_s)),
        col= rgb(0, .1,.5,0.15), border=NA)


##### IDEALICU subanalysis : RCS two knots for ARD
probs.int.p <- with(pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]$probs-.25

dfs<-4
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)

a1<-ns(probs.c, knots = kn, Boundary.knots=rg)[,1]
a2<-ns(probs.c, knots = kn, Boundary.knots=rg)[,2]
a3<-ns(probs.c, knots = kn, Boundary.knots=rg)[,3]
b1<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,1]
b2<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,2]
b3<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,3]



#fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+b1+b2, data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], x=TRUE)
fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+a3+b1+b2+b3, data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], x=TRUE)

y0.precoce_<-data.frame(bras = 1, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
y0.tardif_<-data.frame(bras = 0, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg))

y0.precoce_bis<-data.frame(bras="STRATEGIE PRECOCE", y0.precoce_[,-1])
#colnames(y0.precoce_bis)<-c("bras","a1","a2","b1","b2")
colnames(y0.precoce_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
precoced60d<-1-tail(survfit(fit.int, newdata=y0.precoce_bis)$surv,1)

y0.tardif_bis<-data.frame(bras="STRATEGIE D ATTENTE", y0.tardif_[,-1])
#colnames(y0.tardif_bis)<-c("bras","a1","a2","b1","b2")
colnames(y0.tardif_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
tardifd60d<-1-tail(survfit(fit.int, newdata=y0.tardif_bis)$surv,1)

# caution: function wroks only with two knots at position 0.13143840 & 0.32789328
# for more knots or different knots positions function adjustment is needed

se_s<-ardbootfunction(dataset=pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], nboot = 1000)

plot(seq(0.02468, .55, by=.005), precoced60d-tardifd60d, type="l", ylim=c(-.2,.2))
abline(a=0,b=0)
lines(seq(0.02468, .55, by=.005), (precoced60d-tardifd60d)-qnorm(.975)*se_s, lty=2)
lines(seq(0.02468, .55, by=.005), (precoced60d-tardifd60d)+qnorm(.975)*se_s, lty=2)


### ARD plot : IDEALICU subanalysis : RCS two knots for ARD
dev.off()
plotCI(quantilemean, ARDd60_ci_idealicu[1,], ui=ARDd60_ci_idealicu[3,], li=ARDd60_ci_idealicu[2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Absolute Risk Difference", main="Death at Day 60 in IDEAL-ICU", xlim=c(0, .6), ylim=c(-.35,.35), las=1, bty="n")
lines(c(0,.6), y=c(0,0))
text(quantilemean, 0.35, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
  segments(quantilecuts[i],-.35, quantilecuts[i], .33, lty=3)
}
segments(0,-.35, 0, .27, lty=3)
segments(0,ARD_overall_idealicu[1], .6, ARD_overall_idealicu[1], lty=2, lwd=2)
arrows(.6, c(0.02,-0.02), .6, c(0.02+.15, -0.02-.15), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=0.02, cex=1, adj=0)
mtext(text="favors early", side=4,line=-0, at=-0.02, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), precoced60d-tardifd60d, type = 'l', lwd=1.2,col="#6A6599FF")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c((precoced60d-tardifd60d)+qnorm(.975)*se_s,
                                                                          rev((precoced60d-tardifd60d)-qnorm(.975)*se_s)),
        col= rgb(0, .1,.5,0.15), border=NA)


## RR at day 60

RRd60_table<-ARDd60_table
colnames(RRd60_table)[3]<-"RR"
RRd60_table[,3]<-ARDd60_table[,1]/ARDd60_table[,2]
RRd60_table<-as.data.frame(RRd60_table)
RRd60_table$DC_STRATEGIEPRECOCE<-RRd60_table[,1]*RRd60_table[,5]
RRd60_table$DC_STRATEGIEDATTENTE<-RRd60_table[,2]*RRd60_table[,4]

rrvar<-function(x) {(((x[5]-x[6])/x[6])/x[5]) + (((x[4]-x[7])/x[7])/x[4]) }

RRd60_list<-as.list(as.data.frame(t(RRd60_table)))

RRd60_ci<-
exp(
rbind(
log(RRd60_table[,"RR"]),
log(RRd60_table[,"RR"])-qnorm(.975)*sqrt(sapply(RRd60_list, rrvar)),
log(RRd60_table[,"RR"])+qnorm(.975)*sqrt(sapply(RRd60_list, rrvar))
)
)


plotCI(quantilemean, RRd60_ci[1,], ui=RRd60_ci[3,], li=RRd60_ci[2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="Relative Risk", main="Death at Day 60", xlim=c(0, .6), ylim=c(.5,1.9), las=1, bty="n")
lines(c(0,.6), y=c(1,1))
text(quantilemean, 1.9, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
segments(quantilecuts[i],.55, quantilecuts[i], 1.8, lty=3)
}
segments(0,RR_overall[1], .6, RR_overall[1], lty=2, lwd=2)
segments(0,.55, 0, 1.8, lty=3)
arrows(.6, c(1.1,2-1.1), .6, c(1.4,2-1.4), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=1.1, cex=1, adj=0)
mtext(text="favors early", side=4,line=-0, at=.9, cex=1, adj=1)

## death rate at day 60
ARDd60_table
se_deathrate_early<-function(x) { sqrt((x[1]*(1-x[1]))/x[5]) }
se_deathrate_late<-function(x) { sqrt((x[2]*(1-x[2]))/x[4]) }

deathrate_early<-rbind(
ARDd60_table[,1],
ARDd60_table[,1]-qnorm(.975)*sapply(as.list(as.data.frame(t(ARDd60_table))), se_deathrate_early),
ARDd60_table[,1]+qnorm(.975)*sapply(as.list(as.data.frame(t(ARDd60_table))), se_deathrate_early)
)
rownames(deathrate_early)<-c("death rate", "li", "ui")
deathrate_early

deathrate_late<-rbind(
ARDd60_table[,2],
ARDd60_table[,2]-qnorm(.975)*sapply(as.list(as.data.frame(t(ARDd60_table))), se_deathrate_late),
ARDd60_table[,2]+qnorm(.975)*sapply(as.list(as.data.frame(t(ARDd60_table))), se_deathrate_late)
)
rownames(deathrate_late)<-c("death rate", "li", "ui")
deathrate_late

###

##### Smoothing for event rate plot
# pooled
probs.int.p <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset, ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset$probs-.25

dfs<-3
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)
library(splines)
fit.sep<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.int.p, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset, x=TRUE)

fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.c, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset, x=TRUE)

ypred<-predict(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0), se=TRUE)
ypred2<-predict(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25), se=TRUE)

ypred$fit<-as.numeric(1-tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0))$surv,1))
ypred$se<-as.numeric(tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0))$std.err,1))

ypred2$fit<-as.numeric(1-tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25))$surv,1))
ypred2$se<-as.numeric(tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25))$std.err,1))

yy <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2<- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')

matplot(seq(0.02468, .55, by=.005), matrix(cbind(yy,yy2), ncol=6), type='l', lty=c(1,2,2), col=c(1,1,1,2,2,2), lwd=2, xlab="Probs", ylab="Relative risk")

# akiki
probs.int.p <- with(pooled_po_dataset[pooled_po_dataset$etude=="akiki",], ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset[pooled_po_dataset$etude=="akiki",], ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]$probs-.25

dfs<-3
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)
library(splines)
fit.sep<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.int.p, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset[pooled_po_dataset$etude=="akiki",], x=TRUE)

fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.c, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset[pooled_po_dataset$etude=="akiki",], x=TRUE)

ypred<-predict(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0), se=TRUE)
ypred2<-predict(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25), se=TRUE)

ypred$fit<-as.numeric(1-tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0))$surv,1))
ypred$se<-as.numeric(tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0))$std.err,1))

ypred2$fit<-as.numeric(1-tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25))$surv,1))
ypred2$se<-as.numeric(tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25))$std.err,1))

yy.akiki <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2.akiki <- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')

matplot(seq(0.02468, .55, by=.005), matrix(cbind(yy.akiki,yy2.akiki), ncol=6), type='l', lty=c(1,2,2), col=c(1,1,1,2,2,2), lwd=2, xlab="Probs", ylab="Relative risk")

# idealicu
probs.int.p <- with(pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]$probs-.25

dfs<-5
#kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
rg <-range(probs.c)
library(splines)
fit.sep<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.int.p, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], x=TRUE)

fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + ns(probs.c, knots = kn, Boundary.knots=rg) + ns(probs.int.t, knots = kn, Boundary.knots=rg), data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], x=TRUE)

#ypred<-predict(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0), se=TRUE)
#ypred2<-predict(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25), se=TRUE)

ypred<-data.frame(fit=rep(NA,106), se=rep(NA,106))
ypred2<-data.frame(fit=rep(NA,106), se=rep(NA,106))

ypred$fit<-as.numeric(1-tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0))$surv,1))
ypred$se<-as.numeric(tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=0))$std.err,1))

ypred2$fit<-as.numeric(1-tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25))$surv,1))
ypred2$se<-as.numeric(tail(survfit(fit.int, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq(0.02468, .55, by=.005)-.25, probs.int.t=seq(0.02468, .55, by=.005)-.25))$std.err,1))

yy.idealicu <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2.idealicu <- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')

matplot(seq(0.02468, .55, by=.005), matrix(cbind(yy.idealicu,yy2.idealicu), ncol=6), type='l', lty=c(1,2,2), col=c(1,1,1,2,2,2), lwd=2, xlab="Probs", ylab="Relative risk")


## accurate event rate tables in pooled studies and subanalysis
coxph.obj<-list()
deathrate_early<-list()
deathrate_late<-list()

for (j in 1:q) {
    coxph.obj[[j]]<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$quantile==colnames(HRs[[q]])[j],])
    deathrate_early[[j]]<-c(
    1-tail(survfit(coxph.obj[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$surv,1), 
    1-tail(survfit(coxph.obj[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$upper,1), 
    1-tail(survfit(coxph.obj[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$lower,1)
    )
    
    deathrate_late[[j]]<-c(
    1-tail(survfit(coxph.obj[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$surv,1), 
    1-tail(survfit(coxph.obj[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$upper,1), 
    1-tail(survfit(coxph.obj[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$lower,1)
    )
    }

deathrate_early<-sapply(deathrate_early, c)
colnames(deathrate_early)<-colnames(HRs[[q]])
rownames(deathrate_early)<-c("death rate", "li", "ui")

deathrate_late<-sapply(deathrate_late, c)
colnames(deathrate_late)<-colnames(HRs[[q]])
rownames(deathrate_late)<-c("death rate", "li", "ui")


coxph.obj.akiki<-list()
deathrate_early.akiki<-list()
deathrate_late.akiki<-list()

for (j in 1:q) {
  coxph.obj.akiki[[j]]<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="akiki" & pooled_po_dataset$quantile==colnames(HRs[[q]])[j],])
  deathrate_early.akiki[[j]]<-c(
    1-tail(survfit(coxph.obj.akiki[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$surv,1), 
    1-tail(survfit(coxph.obj.akiki[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$upper,1), 
    1-tail(survfit(coxph.obj.akiki[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$lower,1)
  )
  
  deathrate_late.akiki[[j]]<-c(
    1-tail(survfit(coxph.obj.akiki[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$surv,1), 
    1-tail(survfit(coxph.obj.akiki[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$upper,1), 
    1-tail(survfit(coxph.obj.akiki[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$lower,1)
  )
}

deathrate_early.akiki<-sapply(deathrate_early.akiki, c)
colnames(deathrate_early.akiki)<-colnames(HRs[[q]])
rownames(deathrate_early.akiki)<-c("death rate", "li", "ui")

deathrate_late.akiki<-sapply(deathrate_late.akiki, c)
colnames(deathrate_late.akiki)<-colnames(HRs[[q]])
rownames(deathrate_late.akiki)<-c("death rate", "li", "ui")

coxph.obj.idealicu<-list()
deathrate_early.idealicu<-list()
deathrate_late.idealicu<-list()

for (j in 1:q) {
  coxph.obj.idealicu[[j]]<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = pooled_po_dataset[pooled_po_dataset$etude=="idealicu" & pooled_po_dataset$quantile==colnames(HRs[[q]])[j],])
  deathrate_early.idealicu[[j]]<-c(
    1-tail(survfit(coxph.obj.idealicu[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$surv,1), 
    1-tail(survfit(coxph.obj.idealicu[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$upper,1), 
    1-tail(survfit(coxph.obj.idealicu[[j]], newdata=data.frame(bras="STRATEGIE PRECOCE"))$lower,1)
  )
  
  deathrate_late.idealicu[[j]]<-c(
    1-tail(survfit(coxph.obj.idealicu[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$surv,1), 
    1-tail(survfit(coxph.obj.idealicu[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$upper,1), 
    1-tail(survfit(coxph.obj.idealicu[[j]], newdata=data.frame(bras="STRATEGIE D ATTENTE"))$lower,1)
  )
}

deathrate_early.idealicu<-sapply(deathrate_early.idealicu, c)
colnames(deathrate_early.idealicu)<-colnames(HRs[[q]])
rownames(deathrate_early.idealicu)<-c("death rate", "li", "ui")

deathrate_late.idealicu<-sapply(deathrate_late.idealicu, c)
colnames(deathrate_late.idealicu)<-colnames(HRs[[q]])
rownames(deathrate_late.idealicu)<-c("death rate", "li", "ui")

## event rate plot        
# pooled
plotCI(quantilemean+.005, deathrate_early[1,], ui=deathrate_early[3,], li=deathrate_early[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="60 Days Death Rate",xlim=c(0, .6),ylim=c(.3, .83), las=1, bty="n")
plotCI(quantilemean-.005, deathrate_late[1,], ui=deathrate_late[3,], li=deathrate_late[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="Death Rate", add=TRUE)

segments(0,.3, 0, .8, lty=3)
for (i in 1:length(quantilecuts)) {
  segments(quantilecuts[i],.3, quantilecuts[i], .8, lty=3)
}
text(quantilemean, 0.82, paste("Q", 1:q, sep=""), cex=.8)

lines(seq(0.02468, .55, by=.005), yy[,1], col="#00a1d5")
lines(seq(0.02468, .55, by=.005), yy2[,1], col="#b24745")

legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

# akiki
plotCI(quantilemean+.005, deathrate_early.akiki[1,], ui=deathrate_early.akiki[3,], li=deathrate_early.akiki[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="60 Days Death Rate",xlim=c(0, .6),ylim=c(.2, .83), las=1, bty="n")
plotCI(quantilemean-.005, deathrate_late.akiki[1,], ui=deathrate_late.akiki[3,], li=deathrate_late.akiki[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="Death Rate", add=TRUE)

segments(0,.3, 0, .8, lty=3)
for (i in 1:length(quantilecuts)) {
  segments(quantilecuts[i],.3, quantilecuts[i], .8, lty=3)
}
text(quantilemean, 0.82, paste("Q", 1:q, sep=""), cex=.8)

lines(seq(0.02468, .55, by=.005), yy.akiki[,1], col="#00a1d5")
lines(seq(0.02468, .55, by=.005), yy2.akiki[,1], col="#b24745")

legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

# idealicu
plotCI(quantilemean+.005, deathrate_early.idealicu[1,], ui=deathrate_early.idealicu[3,], li=deathrate_early.idealicu[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="60 Days Death Rate",xlim=c(0, .6),ylim=c(.2, .83), las=1, bty="n")
plotCI(quantilemean-.005, deathrate_late.idealicu[1,], ui=deathrate_late.idealicu[3,], li=deathrate_late.idealicu[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="Predicted Probability of RRT Initiation Within 48 hours", ylab="Death Rate", add=TRUE)

segments(0,.3, 0, .8, lty=3)
for (i in 1:length(quantilecuts)) {
  segments(quantilecuts[i],.3, quantilecuts[i], .8, lty=3)
}
text(quantilemean, 0.82, paste("Q", 1:q, sep=""), cex=.8)

lines(seq(0.02468, .55, by=.005), yy.idealicu[,1], col="#00a1d5")
lines(seq(0.02468, .55, by=.005), yy2.idealicu[,1], col="#b24745")

legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

