save(list=c("quantilemean", "quantilecuts", 
            "deathrate_early", "deathrate_late",
            "deathrate_early.akiki", "deathrate_late.akiki",
            "deathrate_early.idealicu", "deathrate_late.idealicu",
            "HR_overall", "HR_overall_akiki", "HR_overall_idealicu",
            "HRs","HRs_akiki", "HRs_idealicu",
            "ARD_overall", "ARD_overall_akiki", "ARD_overall_idealicu",
            "ARDd60_ci", "ARDd60_ci_akiki", "ARDd60_ci_idealicu",
            
            "yy.er", "yy2.er",
            "yy.er.akiki", "yy2.er.akiki",
            "yy.er.idealicu", "yy2.er.idealicu",
            "yy.hr", "yy.hr.akiki", "yy.hr.idealicu",
            "precoced60d_pooled", "tardifd60d_pooled", "se_s_pooled",
            "precoced60d_akiki", "tardifd60d_akiki", "se_s_akiki",
            "precoced60d_idealicu", "tardifd60d_idealicu", "se_s_idealicu",
            
            "finalmodel", "t_finalmodel_vcov",
            
            "fit.int.hr", "fit.int.hr.akiki", "fit.int.hr.idealicu",
            "fit.int.ard", "fit.int.ard.akiki", "fit.int.ard.idealicu",
            "fit.int.er", "fit.int.er.akiki", "fit.int.er.idealicu"
            ), file="multiscalefigdata.Rdata")

load("multiscalefigdata.Rdata")

library(gplots)

dev.off()
dev.new(width=10,height=10,noRStudioGD = TRUE)
par(mfrow=c(3,3))
q<-5

# pooled
plotCI(quantilemean+.005, deathrate_early[1,], ui=deathrate_early[3,], li=deathrate_early[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="", ylab="60 Days Death Rate",xlim=c(0, .6),ylim=c(.1, .83), yaxt="n", bty="n", main="AKIKI & IDEAL-ICU")
plotCI(quantilemean-.005, deathrate_late[1,], ui=deathrate_late[3,], li=deathrate_late[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="", ylab="", add=TRUE)

axis(2, at=seq(.1,.8, by=.1), las=1)

segments(0,.1, 0, .8, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.1, quantilecuts[i], .8, lty=3)
}
text(quantilemean, 0.83, paste("Q", 1:q, sep=""), cex=.8)

lines(seq(0.02468, .55, by=.005), yy.er[,1], col="#00a1d5")
lines(seq(0.02468, .55, by=.005), yy2.er[,1], col="#b24745")

#legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
#       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

# akiki
plotCI(quantilemean+.005, deathrate_early.akiki[1,], ui=deathrate_early.akiki[3,], li=deathrate_early.akiki[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="", ylab="",xlim=c(0, .6),ylim=c(.1, .83), yaxt="n", bty="n", main="AKIKI")
plotCI(quantilemean-.005, deathrate_late.akiki[1,], ui=deathrate_late.akiki[3,], li=deathrate_late.akiki[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="", ylab="", add=TRUE)

axis(2, at=seq(.1,.8, by=.1), las=1)

segments(0,.1, 0, .8, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.1, quantilecuts[i], .8, lty=3)
}
text(quantilemean, 0.83, paste("Q", 1:q, sep=""), cex=.8)

lines(seq(0.02468, .55, by=.005), yy.er.akiki[,1], col="#00a1d5")
lines(seq(0.02468, .55, by=.005), yy2.er.akiki[,1], col="#b24745")

#legend(.5, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
#       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)

# idealicu
plotCI(quantilemean+.005, deathrate_early.idealicu[1,], ui=deathrate_early.idealicu[3,], li=deathrate_early.idealicu[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="", ylab="",xlim=c(0, .6),ylim=c(.1, .83), yaxt="n", bty="n", main="IDEAL-ICU")
plotCI(quantilemean-.005, deathrate_late.idealicu[1,], ui=deathrate_late.idealicu[3,], li=deathrate_late.idealicu[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="", ylab="", add=TRUE)
axis(2, at=seq(.1,.8, by=.1), las=1)

segments(0,.1, 0, .8, lty=3)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.1, quantilecuts[i], .8, lty=3)
}
text(quantilemean, 0.83, paste("Q", 1:q, sep=""), cex=.8)

lines(seq(0.02468, .55, by=.005), yy.er.idealicu[,1], col="#00a1d5")
lines(seq(0.02468, .55, by=.005), yy2.er.idealicu[,1], col="#b24745")

legend(.4, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
       lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)


### Hazard Ratio Plots ###
## HR plot : akiki & ideal-icu
plot(1, xlim=c(0, .6), ylim = c(.3,3.7), yaxt="n", bty="n", type="n", log="y", xlab="", ylab="Hazard Ratio")
axis(2, at=seq(.3,3.8, by=.7),las=1)
lines(c(0,.6), y=c(1,1), lwd=2)
segments(0,HR_overall[1], .6, HR_overall[1], lty=2, lwd=2)
plotCI(quantilemean, HRs[[q]][1,], ui=HRs[[q]][3,], li=HRs[[q]][2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE, las=1, bty="n")
#plotCI(quantilemean, HRs[[q]][1,], ui=HRs[[q]][3,], li=HRs[[q]][2,], pch=18, gap=0, sfrac=0.003, col="blue4", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Hazard Ratio", main="60 Days Survival", xlim=c(0, .6), ylim = c(.4,2.3))
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.3, quantilecuts[i], 3.7, lty=3)
}
segments(0,.3, 0, 3.7, lty=3)
#text(quantilemean, 3.8, paste("Q", 1:q, sep=""), cex=.8)
#arrows(.6, c(1.05,1/1.05), .6, c(1.05+.4, 1/(1.05+.4)), length = 0.07, xpd=TRUE)
#mtext(text="favors delayed", side=4,line=0, at=1.05, cex=1, adj=0)
#mtext(text="favors early", side=4,line=-0, at=2-1.05, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), exp(yy.hr[,1]), type = 'l', lwd=1.2,col="#79af97")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c(exp(yy.hr[,2]), rev(exp(yy.hr[,3]))),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)

## HR plot : akiki
plot(1, xlim=c(0, .6), ylim = c(.3,3.7), yaxt="n", bty="n", type="n", log="y", xlab="", ylab="")
axis(2, at=seq(.3,3.8, by=.7),las=1) #, labels=NA)
lines(c(0,.6), y=c(1,1), lwd=2)
segments(0,HR_overall_akiki[1], .6, HR_overall_akiki[1], lty=2, lwd=2)
plotCI(quantilemean, HRs_akiki[[q]][1,], ui=HRs_akiki[[q]][3,], li=HRs_akiki[[q]][2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE, las=1, bty="n")
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.3, quantilecuts[i], 3.7, lty=3)
}
segments(0,.3, 0, 3.7, lty=3)
#text(quantilemean, 2.3, paste("Q", 1:q, sep=""), cex=.8)
#arrows(.6, c(1.05,1/1.05), .6, c(1.05+1.25, 1/(1.05+1.25)), length = 0.07, xpd=TRUE)
#mtext(text="favors delayed", side=4,line=0, at=1.05, cex=2/3, adj=0)
#mtext(text="favors early", side=4,line=-0, at=2-1.05, cex=2/3, adj=1)

points(seq(0.02468, .55, by=.005), exp(yy.hr.akiki[,1]), type = 'l', lwd=1.2,col="#79af97")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c(exp(yy.hr.akiki[,2]), rev(exp(yy.hr.akiki[,3]))),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)

## HR plot : idealicu
plot(1, xlim=c(0, .6), ylim = c(.3,3.7), yaxt="n", bty="n", type="n", log="y", xlab="", ylab="")
axis(2, at=seq(.3,3.8, by=.7),las=1) #, labels=NA)
lines(c(0,.6), y=c(1,1), lwd=2)
segments(0,HR_overall_idealicu[1], .6, HR_overall_idealicu[1], lty=2, lwd=2)
plotCI(quantilemean, HRs_idealicu[[q]][1,], ui=HRs_idealicu[[q]][3,], li=HRs_idealicu[[q]][2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", add=TRUE, las=1, bty="n")
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],.3, quantilecuts[i], 3.7, lty=3)
}
segments(0,.3, 0, 3.7, lty=3)
#text(quantilemean, 2.3, paste("Q", 1:q, sep=""), cex=.8)
arrows(.6, c(1.05,1/1.05), .6, c(1.05+1.25, 1/(1.05+1.25)), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=1.05, cex=2/3, adj=0)
mtext(text="favors early", side=4,line=-0, at=2-1.05, cex=2/3, adj=1)

points(seq(0.02468, .55, by=.005), exp(yy.hr.idealicu[,1]), type = 'l', lwd=1.2,col="#79af97")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c(exp(yy.hr.idealicu[,2]), rev(exp(yy.hr.idealicu[,3]))),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)

### ARD plot
### ARD plot : AKIKI & IDEALICU
plot(quantilemean, xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Absolute Risk Reduction", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
segments(0,ARD_overall[1], .6, ARD_overall[1], lty=2, lwd=2)
lines(c(0,.6), y=c(0,0))
plotCI(quantilemean, ARDd60_ci[1,], ui=ARDd60_ci[3,], li=ARDd60_ci[2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="Absolute Risk Reduction", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n", add=TRUE)
#text(quantilemean, 0.3, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-.35, quantilecuts[i], .3, lty=3)
}
segments(0,-.35, 0, .3, lty=3)
#arrows(.6, c(0.02,-0.02), .6, c(0.02+.15, -0.02-.15), length = 0.07, xpd=TRUE)
#mtext(text="favors delayed", side=4,line=0, at=0.02, cex=1, adj=0)
#mtext(text="favors early", side=4,line=-0, at=-0.02, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), precoced60d_pooled-tardifd60d_pooled, type = 'l', lwd=1.2,col="#79af97")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c((precoced60d_pooled-tardifd60d_pooled)+qnorm(.975)*se_s_pooled,
                                                                          rev((precoced60d_pooled-tardifd60d_pooled)-qnorm(.975)*se_s_pooled)),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)

### ARD plot : AKIKI subanalysis : RCS two knots for ARD
plotCI(quantilemean, ARDd60_ci_akiki[1,], ui=ARDd60_ci_akiki[3,], li=ARDd60_ci_akiki[2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
lines(c(0,.6), y=c(0,0))
#text(quantilemean, 0.3, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-.35, quantilecuts[i], .3, lty=3)
}
segments(0,-.35, 0, .3, lty=3)
segments(0,ARD_overall_akiki[1], .6, ARD_overall_akiki[1], lty=2, lwd=2)
#arrows(.6, c(0.02,-0.02), .6, c(0.02+.15, -0.02-.15), length = 0.07, xpd=TRUE)
#mtext(text="favors delayed", side=4,line=0, at=0.02, cex=1, adj=0)
#mtext(text="favors early", side=4,line=-0, at=-0.02, cex=1, adj=1)

points(seq(0.02468, .55, by=.005), precoced60d_akiki-tardifd60d_akiki, type = 'l', lwd=1.2,col="#79af97")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c((precoced60d_akiki-tardifd60d_akiki)+qnorm(.975)*se_s_akiki,
                                                                          rev((precoced60d_akiki-tardifd60d_akiki)-qnorm(.975)*se_s_akiki)),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)

## ARD plot : IDEALICU subanalysis : RCS two knots for ARD
plotCI(quantilemean, ARDd60_ci_idealicu[1,], ui=ARDd60_ci_idealicu[3,], li=ARDd60_ci_idealicu[2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="Predicted Probability of RRT Initiation Within 48 Hours", ylab="", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
lines(c(0,.6), y=c(0,0))
#text(quantilemean, 0.3, paste("Q", 1:q, sep=""), cex=.8)
for (i in 1:length(quantilecuts)) {
        segments(quantilecuts[i],-.35, quantilecuts[i], .3, lty=3)
}
segments(0,-.35, 0, .3, lty=3)
segments(0,ARD_overall_idealicu[1], .6, ARD_overall_idealicu[1], lty=2, lwd=2)
arrows(.6, c(0.02,-0.02), .6, c(0.02+.2, -0.02-.2), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=0.02, cex=2/3, adj=0)
mtext(text="favors early", side=4,line=-0, at=-0.02, cex=2/3, adj=1)

points(seq(0.02468, .55, by=.005), precoced60d_idealicu-tardifd60d_idealicu, type = 'l', lwd=1.2,col="#79af97")
polygon(c(seq(0.02468, .55, by=.005), rev(seq(0.02468, .55, by=.005))), c((precoced60d_idealicu-tardifd60d_idealicu)+qnorm(.975)*se_s_idealicu,
                                                                          rev((precoced60d_idealicu-tardifd60d_idealicu)-qnorm(.975)*se_s_idealicu)),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)

#dev.copy2pdf(file="multiscalefig.pdf")
#

