## ardboot function

ardbootfunction<-function(dataset, nboot=1000) {

ardmat<-matrix(999, nrow=nboot, ncol = length(seq(0.02468, .55, by=.005)))
                
probs.int.p <- with(dataset, ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
probs.int.t <- with(dataset, ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
probs.c<-dataset$probs-.25

set.seed(1000)
for (i in 1:nboot) {
        bootindex<-c(sample(which(dataset$bras!="STRATEGIE PRECOCE"), replace=TRUE), sample(which(dataset$bras=="STRATEGIE PRECOCE"), replace=TRUE))
        zdata<-dataset[bootindex,]
        zdata$probs.int.p <- with(zdata, ifelse(bras=="STRATEGIE PRECOCE", probs-0.25, 0))
        zdata$probs.int.t <- with(zdata, ifelse(bras=="STRATEGIE D ATTENTE", probs-0.25, 0))
        zdata$probs.c<-zdata$probs-.25
        
        dfs<-4
        #kn <- quantile(probs.c, seq(0,1,length=dfs)[-c(1,dfs)])
        kn<-c(0.13143840, 0.32789328)-.25 # 2nd & 5th quantile mean for 6 quantiles
        rg <-range(probs.c)
        
        zdata$a1<-ns(probs.c, knots = kn, Boundary.knots=rg)[,1]
        zdata$a2<-ns(probs.c, knots = kn, Boundary.knots=rg)[,2]
        zdata$a3<-ns(probs.c, knots = kn, Boundary.knots=rg)[,3]
        zdata$b1<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,1]
        zdata$b2<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,2]
        zdata$b3<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,3]
        
        
        
        fit.int<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2+a3+b1+b2+b3, data = zdata, x=TRUE)
        
        temp_y0.precoce_<-data.frame(bras = 1, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
        temp_y0.tardif_<-data.frame(bras = 0, ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg), ns(seq(0.02468, .55, by=.005)-.25, knots = kn, Boundary.knots=rg))
        
        
        
        temp_y0.precoce_bis<-data.frame(bras="STRATEGIE PRECOCE", temp_y0.precoce_[,-1])
        colnames(temp_y0.precoce_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
        temp_precoced60d<-1-tail(survfit(fit.int, newdata=temp_y0.precoce_bis)$surv,1)
        
        temp_y0.tardif_bis<-data.frame(bras="STRATEGIE D ATTENTE", temp_y0.tardif_[,-1])
        colnames(temp_y0.tardif_bis)<-c("bras","a1","a2","a3","b1","b2", "b3")
        temp_tardifd60d<-1-tail(survfit(fit.int, newdata=temp_y0.tardif_bis)$surv,1)
        
        ardmat[i,]<-temp_precoced60d-temp_tardifd60d
}
apply(ardmat, 2, function(x) sd(x))
}


#ardbootfunction(pooled_po_dataset, nboot = 100)
#plot(ardbootfunction(pooled_po_dataset[pooled_po_dataset$etude=="idealicu",], nboot = 100), type="l")
