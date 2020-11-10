library(boot)
#icu free days

negzero<-function(x) ifelse(x<0, 0, x)
pooled_po_dataset$icufree<-trunc(with(pooled_po_dataset, negzero(28-difftime(datesortierea, date.rando, unit="days"))))        
pooled_po_dataset$icufree[pooled_po_dataset$derniere.nouvelles.censureJ60<=28 & pooled_po_dataset$etat.censureJ60==1]<-0


diff_icufree<-list()
for (i in 5:5) {
        n<-i # number of quantiles
        #pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        #outcomecolumn<-"icufree"
        
        diff_icufree[[i]]<-sapply(levels(pooled_po_dataset$quantile), function(x) {welsh<-t.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$icufree,
                                                                                                 pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$icufree, 
                                                                                                 var.equal = FALSE)
        wilcox<-wilcox.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$icufree,
                            pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$icufree)
        nprecoce<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",])
        res<-boot(pooled_po_dataset[pooled_po_dataset$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(pooled_po_dataset[pooled_po_dataset$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="icufree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_icufree[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_icufree[[i]][1,], ui=diff_icufree[[i]][3,], li=diff_icufree[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", main="AKIKI & IDEAL-ICU", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_icufree[[i]][1,], ui=diff_icufree[[i]][8,], li=diff_icufree[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_icufree[[i]][4,],2), nsmall=2))
}

diff_ventilfree<-list()
for (i in 5:5) {
        n<-i # number of quantiles
        #pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        diff_ventilfree[[i]]<-sapply(levels(pooled_po_dataset$quantile), function(x) {welsh<-t.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                                                                                                 pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$ventilfree, 
                                                                                                 var.equal = FALSE)
        wilcox<-wilcox.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                            pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$ventilfree)
        nprecoce<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",])
        res<-boot(pooled_po_dataset[pooled_po_dataset$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(pooled_po_dataset[pooled_po_dataset$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="ventilfree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_ventilfree[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_ventilfree[[i]][1,], ui=diff_ventilfree[[i]][3,], li=diff_ventilfree[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="ventil free days (late-early)", main="AKIKI & IDEAL-ICU", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_ventilfree[[i]][1,], ui=diff_ventilfree[[i]][8,], li=diff_ventilfree[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_ventilfree[[i]][4,],2), nsmall=2))
}

diff_rrtfree<-list()
for (i in 5:5) {
        n<-i # number of quantiles
        #pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        diff_rrtfree[[i]]<-sapply(levels(pooled_po_dataset$quantile), function(x) {welsh<-t.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                                                                                                    pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$rrtfree, 
                                                                                                    var.equal = FALSE)
        wilcox<-wilcox.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                            pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$rrtfree)
        nprecoce<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",])
        res<-boot(pooled_po_dataset[pooled_po_dataset$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(pooled_po_dataset[pooled_po_dataset$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="rrtfree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_rrtfree[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_rrtfree[[i]][1,], ui=diff_rrtfree[[i]][3,], li=diff_rrtfree[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="rrt free days (late-early)", main="AKIKI & IDEAL-ICU", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_rrtfree[[i]][1,], ui=diff_rrtfree[[i]][8,], li=diff_rrtfree[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_rrtfree[[i]][4,],2), nsmall=2))
}

### akiki subanalysis

diff_icufree_akiki<-list()
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]
for (i in 5:5) {
        n<-i # number of quantiles
        #zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        #outcomecolumn<-"icufree"
        
        diff_icufree_akiki[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                                                                                     zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree, 
                                                                                                 var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        res<-boot(zdata[zdata$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(zdata[zdata$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="icufree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_icufree_akiki[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_icufree_akiki[[i]][1,], ui=diff_icufree_akiki[[i]][3,], li=diff_icufree_akiki[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", main="AKIKI", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_icufree_akiki[[i]][1,], ui=diff_icufree_akiki[[i]][8,], li=diff_icufree_akiki[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_icufree_akiki[[i]][4,],2), nsmall=2))
}

diff_ventilfree_akiki<-list()
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]
for (i in 5:5) {
        n<-i # number of quantiles
        #zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))

        diff_ventilfree_akiki[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                                                                                     zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree, 
                                                                                     var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        res<-boot(zdata[zdata$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(zdata[zdata$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="ventilfree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_ventilfree_akiki[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_ventilfree_akiki[[i]][1,], ui=diff_ventilfree_akiki[[i]][3,], li=diff_ventilfree_akiki[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="ventil free days (late-early)", main="AKIKI", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_ventilfree_akiki[[i]][1,], ui=diff_ventilfree_akiki[[i]][8,], li=diff_ventilfree_akiki[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_ventilfree_akiki[[i]][4,],2), nsmall=2))
}

diff_rrtfree_akiki<-list()
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]
for (i in 5:5) {
        n<-i # number of quantiles
        #zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_rrtfree_akiki[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                                                                                              zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$rrtfree, 
                                                                                              var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$rrtfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        res<-boot(zdata[zdata$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(zdata[zdata$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="rrtfree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_rrtfree_akiki[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_rrtfree_akiki[[i]][1,], ui=diff_rrtfree_akiki[[i]][3,], li=diff_rrtfree_akiki[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="rrt free days (late-early)", main="AKIKI", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_rrtfree_akiki[[i]][1,], ui=diff_rrtfree_akiki[[i]][8,], li=diff_rrtfree_akiki[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_rrtfree_akiki[[i]][4,],2), nsmall=2))
}

### idealicu subanalysis

diff_icufree_idealicu<-list()
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]
for (i in 5:5) {
        n<-i # number of quantiles
        #zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        #outcomecolumn<-"icufree"
        
        diff_icufree_idealicu[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                                                                                           zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree, 
                                                                                           var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        res<-boot(zdata[zdata$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(zdata[zdata$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="icufree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_icufree_idealicu[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_icufree_idealicu[[i]][1,], ui=diff_icufree_idealicu[[i]][3,], li=diff_icufree_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", main="IDEAL-ICU", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_icufree_idealicu[[i]][1,], ui=diff_icufree_idealicu[[i]][8,], li=diff_icufree_idealicu[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_icufree_idealicu[[i]][4,],2), nsmall=2))
}

diff_ventilfree_idealicu<-list()
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]
for (i in 5:5) {
        n<-i # number of quantiles
        #zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_ventilfree_idealicu[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                                                                                              zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree, 
                                                                                              var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        res<-boot(zdata[zdata$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(zdata[zdata$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="ventilfree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_ventilfree_idealicu[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_ventilfree_idealicu[[i]][1,], ui=diff_ventilfree_idealicu[[i]][3,], li=diff_ventilfree_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="ventil free days (late-early)", main="IDEAL-ICU", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_ventilfree_idealicu[[i]][1,], ui=diff_ventilfree_idealicu[[i]][8,], li=diff_ventilfree_idealicu[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_ventilfree_idealicu[[i]][4,],2), nsmall=2))
}

diff_rrtfree_idealicu<-list()
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]
for (i in 5:5) {
        n<-i # number of quantiles
        #zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_rrtfree_idealicu[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                                                                                           zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$rrtfree, 
                                                                                           var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$rrtfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        res<-boot(zdata[zdata$quantile==x,], meandiff_bootfunction, R=1999, strata=as.numeric(zdata[zdata$quantile==x,]$bras=="STRATEGIE D ATTENTE"), outcomecolumn="rrtfree")
        bootcibca<-boot.ci(res)$bca[4:5]
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif, bootcibca)
        })
        rownames(diff_rrtfree_idealicu[[i]])<-c("meandif", "li_welsh", "ui_welsh", "wilcox_pval","nprecoce", "ntardif", "li_boot", "ui_boot")
        plotCI(1:n, diff_rrtfree_idealicu[[i]][1,], ui=diff_rrtfree_idealicu[[i]][3,], li=diff_rrtfree_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="rrt free days (late-early)", main="IDEAL-ICU", ylim=c(-6,6), las=1)
        plotCI((1:n)+.1, diff_rrtfree_idealicu[[i]][1,], ui=diff_rrtfree_idealicu[[i]][8,], li=diff_rrtfree_idealicu[[i]][7,], pch=18, gap=0, sfrac=0.005, col="green", add=TRUE)
        abline(h=0)
        legend("topright", legend=c("Welsh CI", "Bootstrap CI"), col=c("red", "green"), pch=18, lty=1)
        text(1:n, -6, format(round(diff_rrtfree_idealicu[[i]][4,],2), nsmall=2))
}


### plot

dev.off()
dev.new(width=10,height=10,noRStudioGD = TRUE)
par(mfrow=c(3,3), mar = c(5.1, 5.1, 4.1, 2.1))

#RRT free days
# pooled
plot(quantilemean, xlab="", ylab="RRT-Free Days at Day 28\n (Mean Difference)",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="AKIKI & IDEAL-ICU")
lines(c(0,.6), y=c(0,0), lwd=2)
plotCI(quantilemean, diff_rrtfree[[5]][1,], ui=diff_rrtfree[[5]][8,], li=diff_rrtfree[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}
text(quantilemean, 10.3, paste("Q", 1:q, sep=""), cex=.8)

#legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
#       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

# akiki
plot(quantilemean, xlab="", ylab="",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="AKIKI")
lines(c(0,.6), y=c(0,0), lwd=2)
plotCI(quantilemean, diff_rrtfree_akiki[[5]][1,], ui=diff_rrtfree_akiki[[5]][8,], li=diff_rrtfree_akiki[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}
text(quantilemean, 10.3, paste("Q", 1:q, sep=""), cex=.8)

# idealicu
plotCI(quantilemean, diff_rrtfree_idealicu[[5]][1,], ui=diff_rrtfree_idealicu[[5]][8,], li=diff_rrtfree_idealicu[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="", ylab="",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="IDEAL-ICU")

axis(2, at=seq(-10,10, by=2), las=1)
lines(c(0,.6), y=c(0,0), lwd=2)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}
text(quantilemean, 10.3, paste("Q", 1:q, sep=""), cex=.8)

arrows(.6, c(.5,-.5), .6, c(.5+6.5, -.5-6.5), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=.5, cex=2/3, adj=0)
mtext(text="favors early", side=4,line=-0, at=-.5, cex=2/3, adj=1)

# ICU free days
# pooled
plot(quantilemean, xlab="", ylab="ICU-Free Days at Day 28\n (Mean Difference)",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="")
lines(c(0,.6), y=c(0,0), lwd=2)

plotCI(quantilemean, diff_icufree[[5]][1,], ui=diff_icufree[[5]][8,], li=diff_icufree[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}


#legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
#       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

# akiki
plot(quantilemean, xlab="", ylab="",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="")
lines(c(0,.6), y=c(0,0), lwd=2)
plotCI(quantilemean, diff_icufree_akiki[[5]][1,], ui=diff_icufree_akiki[[5]][8,], li=diff_icufree_akiki[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black",add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}


# idealicu
plot(quantilemean, xlab="", ylab="",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="")
lines(c(0,.6), y=c(0,0), lwd=2)

plotCI(quantilemean, diff_icufree_idealicu[[5]][1,], ui=diff_icufree_idealicu[[5]][8,], li=diff_icufree_idealicu[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}


arrows(.6, c(.5,-.5), .6, c(.5+6.5, -.5-6.5), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=.5, cex=2/3, adj=0)
mtext(text="favors early", side=4,line=-0, at=-.5, cex=2/3, adj=1)


# ventilator free days
# pooled
plot(quantilemean, xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Ventilator-Free Days at Day 28\n (Mean Difference)",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="")
lines(c(0,.6), y=c(0,0), lwd=2)
plotCI(quantilemean, diff_ventilfree[[5]][1,], ui=diff_ventilfree[[5]][8,], li=diff_ventilfree[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black",add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}


# akiki
plot(quantilemean, xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="")
lines(c(0,.6), y=c(0,0), lwd=2)
plotCI(quantilemean, diff_ventilfree_akiki[[5]][1,], ui=diff_ventilfree_akiki[[5]][8,], li=diff_ventilfree_akiki[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE)

axis(2, at=seq(-10,10, by=2), las=1)

segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}


# idealicu
plotCI(quantilemean, diff_ventilfree_idealicu[[5]][1,], ui=diff_ventilfree_idealicu[[5]][8,], li=diff_ventilfree_idealicu[[5]][7,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="",xlim=c(0, .6),ylim=c(-10, 10), yaxt="n", bty="n", main="")

axis(2, at=seq(-10,10, by=2), las=1)
lines(c(0,.6), y=c(0,0), lwd=2)
segments(0, -10, 0, 10, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-10, quantilecuts[i], 10, lty=3)
}


arrows(.6, c(.5,-.5), .6, c(.5+6.5, -.5-6.5), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=.5, cex=2/3, adj=0)
mtext(text="favors early", side=4,line=-0, at=-.5, cex=2/3, adj=1)

dev.copy2pdf(file="freedaysmeandif.pdf")




