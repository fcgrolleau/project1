
q<-5
pooled_po_dataset$quantile<-quantcut(pooled_po_dataset$probs, seq(0,1,by=1/q))        

tapply(pooled_po_dataset$probs, pooled_po_dataset$quantile, function(x) summary(x)[4])
pooled_po_dataset$quantilemean<-as.numeric(quantilemean[match(pooled_po_dataset$quantile, names(quantilemean), nomatch = 0)])

library(grid)
grob <- grobTree(textGrob("N evaluated", just = "left", gp=gpar(fontsize=15)))
grob1 <- grobTree(textGrob("Early Strategy", just = "left", gp=gpar(fontsize=10)))
grob2 <- grobTree(textGrob("Delayed strategy", just = "left", gp=gpar(fontsize=10)))
grob3 <- grobTree(textGrob("P-value", just = "left", gp=gpar(fontsize=15)))

boxplotsgg1.1<-
ggplot() + geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = ventilfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = ventilfree, fill="#b24745"), width=0.01) +  
        labs(x = "", y = "Ventilator-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grobTree(textGrob("AKIKI & IDEAL-ICU", just = "center", gp=gpar(fontsize=20, fontface="bold"))), xmin=.24,xmax=.24,ymin=32, ymax=32) +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
quantname[[i]] <-  grobTree(textGrob(paste("Q",i, sep=""), just = "bottom", gp=gpar(fontsize=8)))
#nprecoce[[i]] <-  grobTree(textGrob(diff_ventilfree[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
#ntardif[[i]] <-  grobTree(textGrob(diff_ventilfree[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_ventilfree[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg1.1 <- boxplotsgg1.1 +
        annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=28.5, ymax=28.5) +
        #annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
        #annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
        annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg1.1
#ggsave("vfd.pdf", plot = boxplotsgg1.1, height = 12, width = 11)
dev.off()

### icu free days plot

boxplotsgg2.1<-
        ggplot() + geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = icufree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = icufree, fill="#b24745"), width=0.01) +  
        labs(x = "", y = "ICU-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        #quantname[[i]] <-  grobTree(textGrob(paste("Q",i, sep=""), just = "center", gp=gpar(fontsize=8)))
        #nprecoce[[i]] <-  grobTree(textGrob(diff_icufree[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        #ntardif[[i]] <-  grobTree(textGrob(diff_icufree[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_icufree[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg2.1 <- boxplotsgg2.1 +
                #annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                #annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                #annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg2.1
#ggsave("icufd.pdf", plot = boxplotsgg2.1, height = 12, width = 11)
dev.off()


### rrt free days plot

boxplotsgg3.1<-
        ggplot() +
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = rrtfree, fill="#00a1d5"), width=0.01) +
        geom_boxplot(data=pooled_po_dataset[pooled_po_dataset$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = rrtfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "RRT-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grob, xmin=-.07,xmax=-.07,ymin=-7.5, ymax=-7.5) + annotation_custom(grobTree(textGrob("AKIKI & IDEAL-ICU", just = "center", gp=gpar(fontsize=15))), xmin=.25,xmax=.25,ymin=-7.5, ymax=-7.5) +
        annotation_custom(grob1, xmin=-.07,xmax=-.07,ymin=-10, ymax=-10) +
        annotation_custom(grob2, xmin=-.07,xmax=-.07,ymin=-11.5, ymax=-11.5) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        #quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_rrtfree[[q]][5,i], just = "right", gp=gpar(fontsize=10)))
        ntardif[[i]] <-  grobTree(textGrob(diff_rrtfree[[q]][6,i], just = "right", gp=gpar(fontsize=10)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_rrtfree[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg3.1 <- boxplotsgg3.1 +
                #annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10, ymax=-10) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-11.5, ymax=-11.5) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg3.1
#ggsave("rrtfd.pdf", plot = boxplotsgg3.1, height = 12, width = 11)
dev.off()


### akiki subanalysis free-days outcomes plots 
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="akiki",]

boxplotsgg1.2<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = ventilfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = ventilfree, fill="#b24745"), width=0.01) +  
        labs(x = "", y = "Ventilator-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grobTree(textGrob("AKIKI", just = "center", gp=gpar(fontsize=20, fontface="bold"))), xmin=.24,xmax=.24,ymin=32, ymax=32) +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i, sep=""), just = "bottom", gp=gpar(fontsize=8)))
        #nprecoce[[i]] <-  grobTree(textGrob(diff_ventilfree_akiki[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        #ntardif[[i]] <-  grobTree(textGrob(diff_ventilfree_akiki[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_ventilfree_akiki[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg1.2 <- boxplotsgg1.2 +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=28.5, ymax=28.5) +
                #annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                #annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg1.2
#ggsave("vfd_akiki.pdf", plot = boxplotsgg1.2, height = 12, width = 11)
dev.off()


boxplotsgg2.2<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = icufree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = icufree, fill="#b24745"), width=0.01) +  
        labs(x = "", y = "ICU-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        #quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        #nprecoce[[i]] <-  grobTree(textGrob(diff_icufree_akiki[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        #ntardif[[i]] <-  grobTree(textGrob(diff_icufree_akiki[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_icufree_akiki[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg2.2 <- boxplotsgg2.2 +
                #annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                #annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                #annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg2.2
#ggsave("icufd_akiki.pdf", plot = boxplotsgg2.2, height = 12, width = 11)
dev.off()


boxplotsgg3.2<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = rrtfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = rrtfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "RRT-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grobTree(textGrob("AKIKI", just = "center", gp=gpar(fontsize=15))), xmin=.25,xmax=.25,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        #quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_rrtfree_akiki[[q]][5,i], just = "right", gp=gpar(fontsize=10)))
        ntardif[[i]] <-  grobTree(textGrob(diff_rrtfree_akiki[[q]][6,i], just = "right", gp=gpar(fontsize=10)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_rrtfree_akiki[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg3.2 <- boxplotsgg3.2 +
                #annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10, ymax=-10) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-11.5, ymax=-11.5) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg3.2
#ggsave("rrtfd_akiki.pdf", plot = boxplotsgg3.2, height = 12, width = 11)
dev.off()


### idealicu subanalysis free-days outcomes plots 
zdata<-pooled_po_dataset[pooled_po_dataset$etude=="idealicu",]

boxplotsgg1.3<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = ventilfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = ventilfree, fill="#b24745"), width=0.01) +  
        labs(x = "", y = "Ventilator-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(1.1,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.4,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        annotation_custom(grobTree(textGrob("IDEAL-ICU", just = "center", gp=gpar(fontsize=20, fontface="bold"))), xmin=.24,xmax=.24,ymin=32, ymax=32) +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        quantname[[i]] <-  grobTree(textGrob(paste("Q",i, sep=""), just = "bottom", gp=gpar(fontsize=8)))
        #nprecoce[[i]] <-  grobTree(textGrob(diff_ventilfree_idealicu[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        #ntardif[[i]] <-  grobTree(textGrob(diff_ventilfree_idealicu[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_ventilfree_idealicu[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg1.3 <- boxplotsgg1.3 +
                annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=28.5, ymax=28.5) +
                #annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                #annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg1.3
#ggsave("vfd_idealicu.pdf", plot = boxplotsgg1.3, height = 12, width = 11)
dev.off()


boxplotsgg2.3<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = icufree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = icufree, fill="#b24745"), width=0.01) +  
        labs(x = "", y = "ICU-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        #quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        #nprecoce[[i]] <-  grobTree(textGrob(diff_icufree_idealicu[[q]][5,i], just = "right", gp=gpar(fontsize=8)))
        #ntardif[[i]] <-  grobTree(textGrob(diff_icufree_idealicu[[q]][6,i], just = "right", gp=gpar(fontsize=8)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_icufree_idealicu[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg2.3 <- boxplotsgg2.3 +
                #annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                #annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-7.5, ymax=-7.5) +
                #annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-9, ymax=-9) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg2.3
#ggsave("icufd_idealicu.pdf", plot = boxplotsgg2.3, height = 12, width = 11)
dev.off()


boxplotsgg3.3<-
        ggplot() + geom_boxplot(data=zdata[zdata$bras=="STRATEGIE PRECOCE",], aes(x = quantilemean-.007, group=quantilemean, y = rrtfree, fill="#00a1d5"), width=0.01) + 
        geom_boxplot(data=zdata[zdata$bras=="STRATEGIE D ATTENTE",], aes(x = quantilemean+.007, group=quantilemean, y = rrtfree, fill="#b24745"), width=0.01) +  
        labs(x = "Predicted Probability of RRT Initiation Within 48 Hours", y = "RRT-Free Days") +
        scale_fill_manual(name="", labels=c("Early strategy","Delayed strategy"), values=c("#00a1d5","#b24745")) +
        scale_x_continuous(breaks=seq(0, .6, .1), limits = c(-.02,.5), expand = c(0,0)) +
        scale_y_continuous(breaks=seq(0, 28, 7), limits = c(-2,28), expand = c(0,0)) +
        theme_classic() + theme(legend.position = c(2,.8), legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"), plot.margin = unit(c(1,7,7,1), "lines"), axis.title.x = element_text(vjust=-2.5) ) + coord_cartesian(clip = "off") +
        #annotation_custom(grob, xmin=-.06,xmax=-.06,ymin=-6, ymax=-6) +
        annotation_custom(grobTree(textGrob("IDEAL-ICU", just = "center", gp=gpar(fontsize=15))), xmin=.25,xmax=.25,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob1, xmin=-.04,xmax=-.04,ymin=-7.5, ymax=-7.5) +
        #annotation_custom(grob2, xmin=-.04,xmax=-.04,ymin=-9, ymax=-9) +
        #annotation_custom(grob3, xmin=-.06,xmax=-.06,ymin=-10.5, ymax=-10.5) +
        geom_vline(xintercept = c(0,quantilecuts), linetype="dotted", color = "black", size=0.5)

quantname<-list()
nprecoce<-list()
ntardif<-list()
wilpavl<-list()

for (i in 1:q) {
        #quantname[[i]] <-  grobTree(textGrob(paste("Q",i," ", names(quantilemean)[i], sep=""), just = "center", gp=gpar(fontsize=8)))
        nprecoce[[i]] <-  grobTree(textGrob(diff_rrtfree_idealicu[[q]][5,i], just = "right", gp=gpar(fontsize=10)))
        ntardif[[i]] <-  grobTree(textGrob(diff_rrtfree_idealicu[[q]][6,i], just = "right", gp=gpar(fontsize=10)))
        wilpavl[[i]] <-  grobTree(textGrob(paste("P = ", format(round(diff_rrtfree_idealicu[[q]][4,i],2),nsmall=2)), just = "center", gp=gpar(fontsize=8)))
}       

for (i in 1:q) {        
        boxplotsgg3.3 <- boxplotsgg3.3 +
                #annotation_custom(quantname[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1) +
                annotation_custom(nprecoce[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-10, ymax=-10) +
                annotation_custom(ntardif[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-11.5, ymax=-11.5) +
                annotation_custom(wilpavl[[i]], xmin=quantilemean[i], xmax=quantilemean[i], ymin=-1, ymax=-1)
}

boxplotsgg3.3
#ggsave("rrtfd_idealicu.pdf", plot = boxplotsgg3.3, height = 12, width = 11)
dev.off()


require(gridExtra)
dev.new(width=20,height=20,noRStudioGD = TRUE)

grid.arrange(nrow=3, ncol=3, 
             boxplotsgg1.1, boxplotsgg1.2, boxplotsgg1.3, 
             boxplotsgg2.1, boxplotsgg2.2, boxplotsgg2.3,
             boxplotsgg3.1, boxplotsgg3.2, boxplotsgg3.3)
dev.copy2pdf(file="boxplotstest.pdf")





