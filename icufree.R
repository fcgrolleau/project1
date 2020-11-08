#icu free days

negzero<-function(x) ifelse(x<0, 0, x)
pooled_po_dataset$icufree<-trunc(with(pooled_po_dataset, negzero(28-difftime(datesortierea, date.rando, unit="days"))))        
pooled_po_dataset$icufree[pooled_po_dataset$derniere.nouvelles.censureJ60<=28 & pooled_po_dataset$etat.censureJ60==1]<-0


diff_icufree<-list()
par(mfrow=c(2,4))
for (i in 3:10) {
        n<-i # number of quantiles
        pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        diff_icufree[[i]]<-sapply(levels(pooled_po_dataset$quantile), function(x) {welsh<-t.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$icufree,
                                                                                                    pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$icufree, 
                                                                                                    var.equal = FALSE)
        wilcox<-wilcox.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$icufree,
                            pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$icufree)
        nprecoce<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_icufree[[i]][1,], ui=diff_icufree[[i]][3,], li=diff_icufree[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_icufree[[i]][4,],2), nsmall=2))
}


diff_ventilfree<-list()
par(mfrow=c(2,4))
for (i in 3:10) {
        n<-i # number of quantiles
        pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        diff_ventilfree[[i]]<-sapply(levels(pooled_po_dataset$quantile), function(x) {welsh<-t.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                                                                                                 pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$ventilfree, 
                                                                                                 var.equal = FALSE)
        wilcox<-wilcox.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                            pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$ventilfree)
        nprecoce<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_ventilfree[[i]][1,], ui=diff_ventilfree[[i]][3,], li=diff_ventilfree[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="rrt free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_ventilfree[[i]][4,],2), nsmall=2))
}


diff_rrtfree<-list()
par(mfrow=c(2,4))
for (i in 3:10) {
        n<-i # number of quantiles
        pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/n))
        
        diff_rrtfree[[i]]<-sapply(levels(pooled_po_dataset$quantile), function(x) {welsh<-t.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                                                                                                    pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$rrtfree, 
                                                                                                    var.equal = FALSE)
        wilcox<-wilcox.test(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                            pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",]$rrtfree)
        nprecoce<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(pooled_po_dataset[pooled_po_dataset$quantile==x & pooled_po_dataset$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_rrtfree[[i]][1,], ui=diff_rrtfree[[i]][3,], li=diff_rrtfree[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="ventilator free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_rrtfree[[i]][4,],2), nsmall=2))
}


### akiki subanalysis

diff_icufree_akiki<-list()
par(mfrow=c(2,4))
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]
for (i in 3:10) {
        n<-i # number of quantiles
        zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_icufree_akiki[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                                                                                     zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree, 
                                                                                                 var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_icufree_akiki[[i]][1,], ui=diff_icufree_akiki[[i]][3,], li=diff_icufree_akiki[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_icufree_akiki[[i]][4,],2), nsmall=2))
}

diff_ventilfree_akiki<-list()
par(mfrow=c(2,4))
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]
for (i in 3:10) {
        n<-i # number of quantiles
        zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_ventilfree_akiki[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                                                                                           zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree, 
                                                                                           var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_ventilfree_akiki[[i]][1,], ui=diff_ventilfree_akiki[[i]][3,], li=diff_ventilfree_akiki[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_ventilfree_akiki[[i]][4,],2), nsmall=2))
}

### ideal-icu subanalysis

diff_icufree_idealicu<-list()
par(mfrow=c(2,4))
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]
for (i in 3:10) {
        n<-i # number of quantiles
        zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_icufree_idealicu[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                                                                                           zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree, 
                                                                                           var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$icufree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$icufree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_icufree_idealicu[[i]][1,], ui=diff_icufree_idealicu[[i]][3,], li=diff_icufree_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_icufree_idealicu[[i]][4,],2), nsmall=2))
}

diff_ventilfree_idealicu<-list()
par(mfrow=c(2,4))
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]
for (i in 3:10) {
        n<-i # number of quantiles
        zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_ventilfree_idealicu[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                                                                                              zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree, 
                                                                                              var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$ventilfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$ventilfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_ventilfree_idealicu[[i]][1,], ui=diff_ventilfree_idealicu[[i]][3,], li=diff_ventilfree_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="icu free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_ventilfree_idealicu[[i]][4,],2), nsmall=2))
}

diff_rrtfree_idealicu<-list()
par(mfrow=c(2,4))
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]
for (i in 3:10) {
        n<-i # number of quantiles
        zdata$quantile<-quantcut(zdata$probs, seq(0,1,by=1/n))
        
        diff_rrtfree_idealicu[[i]]<-sapply(levels(zdata$quantile), function(x) {welsh<-t.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                                                                                          zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$rrtfree, 
                                                                                          var.equal = FALSE)
        wilcox<-wilcox.test(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",]$rrtfree,
                            zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",]$rrtfree)
        nprecoce<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE PRECOCE",])
        ntardif<-nrow(zdata[zdata$quantile==x & zdata$bras=="STRATEGIE D ATTENTE",])
        c(as.numeric(diff(rev(welsh$estimate))), welsh$conf.int, wilcox$p.value, nprecoce, ntardif)
        })
        
        plotCI(1:n, diff_rrtfree_idealicu[[i]][1,], ui=diff_rrtfree_idealicu[[i]][3,], li=diff_rrtfree_idealicu[[i]][2,], pch=18, gap=0, sfrac=0.005, col="red", xlab="Quantile of Predicted Probability", ylab="Ventilator free days (late-early)", ylim=c(-6,6))
        abline(h=0)
        text(1:n, -6, format(round(diff_rrtfree_idealicu[[i]][4,],2), nsmall=2))
}


q<-5
pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/q))        

tapply(pooled_po_dataset$probs, pooled_po_dataset$quantile, function(x) summary(x)[4])
pooled_po_dataset$quantilemean<-as.numeric(quantilemean[match(pooled_po_dataset$quantile, names(quantilemean), nomatch = 0)])

library(grid)
grob <- grobTree(textGrob("N evaluated", just = "left", gp=gpar(fontsize=10)))
grob1 <- grobTree(textGrob("Early Strategy", just = "left", gp=gpar(fontsize=8)))
grob2 <- grobTree(textGrob("Delayed strategy", just = "left", gp=gpar(fontsize=8)))
grob3 <- grobTree(textGrob("P-value", just = "left", gp=gpar(fontsize=10)))

boxplotsgg<-
ggplot() + geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = ventilfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = ventilfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "Ventilator-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
nprecoce[[i]] <-  grobTree(textGrob(diff_ventilfree[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
ntardif[[i]] <-  grobTree(textGrob(diff_ventilfree[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_ventilfree[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
boxplotsgg <- boxplotsgg +
        annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
        annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
        annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
        annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("vfd.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()

### icu free days plot

boxplotsgg<-
        ggplot() + geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = icufree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = icufree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "ICU-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_icufree[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_icufree[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_icufree[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("icufd.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


### rrt free days plot

boxplotsgg<-
        ggplot() +
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = rrtfree, fill="#00a1d5"), width=0.01) +
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = rrtfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "RRT-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_rrtfree[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_rrtfree[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_rrtfree[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("rrtfd.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


### akiki subanalysis free-days outcomes plots 
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]

boxplotsgg<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = ventilfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = ventilfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "Ventilator-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_ventilfree_akiki[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_ventilfree_akiki[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_ventilfree_akiki[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("vfd_akiki.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


boxplotsgg<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = icufree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = icufree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "ICU-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_icufree_akiki[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_icufree_akiki[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_icufree_akiki[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("icufd_akiki.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


boxplotsgg<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = rrtfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = rrtfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "RRT-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_rrtfree_akiki[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_rrtfree_akiki[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_rrtfree_akiki[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("rrtfd_akiki.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


### idealicu subanalysis free-days outcomes plots 
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]

boxplotsgg<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = ventilfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = ventilfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "Ventilator-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_ventilfree_idealicu[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_ventilfree_idealicu[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_ventilfree_idealicu[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("vfd_idealicu.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


boxplotsgg<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = icufree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = icufree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "ICU-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_icufree_idealicu[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_icufree_idealicu[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_icufree_idealicu[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("icufd_idealicu.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()


boxplotsgg<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = rrtfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = rrtfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "RRT-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme( legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,1,13,4), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_rrtfree_idealicu[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        ntardif[[i]] <-  grobTree(textGrob(diff_rrtfree_idealicu[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(format(round(diff_rrtfree_idealicu[[q]][4,i],2),nsmall=2), just = "right", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg <- boxplotsgg +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10.5, ymax=-10.5)
}

boxplotsgg
ggsave("rrtfd_idealicu.pdf", plot = boxplotsgg, height = 12, width = 11)
dev.off()




