library(shiny)
shinyServer(function(input, output) {
  #install.packages("survival")
  #install.packages("splines")
  #install.packages("gplots")
  library(survival)
  library(splines)
  library(gplots)
  # load all necessary data
  load(url("https://github.com/fcgrolleau/Functions-and-explanatory-analysis/raw/master/multiscalefigdata.Rdata"))
  
  expit <- function(x){exp(x)/(1+exp(x))}
  coefs<-finalmodel$estimate
  x_vect<-reactive({as.numeric(c(1, input$potassium, input$sofa, input$weight, input$drug, input$urea, input$ph))})
  pred<-reactive({ expit(x_vect() %*% coefs) })
  
  pred_se<-reactive({ sqrt(t(x_vect()) %*% t_finalmodel_vcov %*% x_vect()) })
  lci_pred<-reactive({ expit(x_vect() %*% coefs - qnorm(.975)*pred_se()) })
  uci_pred<-reactive({ expit(x_vect() %*% coefs + qnorm(.975)*pred_se()) })
  
  output$prediction<-renderText({ paste( sep="",
                                        "The predicted probablity of RRT initiation within 48 hours is ",
                                        100*round(pred(),2), "%") })
  
  output$prediction2<-renderText({ paste(sep = "", 
    100*round(pred(),2), "% ",
    "95%CI [",
    100*round(lci_pred(),2),"%",
    " - ",
    100*round(uci_pred(),2),"%",
    "]") })
  
  
###  
  
  # Plot the data ####
  output$myImage <- renderImage({
    # A temp file to save the output.
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext = '.png')
    
    # Generate the PNG
    png(outfile, 
        width = 600*8, 
        height = 600*8,
        res = 72*8)
    
### compute predicted probablity of RRT initiation in the reactive image
    x_vect<-as.numeric(c(1,input$potassium,input$sofa, input$weight, input$drug, input$urea, input$ph))
    pred<-expit(x_vect() %*% coefs)
    trunc_pred<-ifelse(pred>.6, .62, pred)
    
    
### compute other stats directly in the reactive image
    dr_early<-as.numeric(1-tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=pred-.25, probs.int.t=0))$surv,1))
    dr_early_se<-as.numeric(tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=pred-.25, probs.int.t=0))$std.err,1))
    
    dr_late<-as.numeric(1-tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=pred-.25, probs.int.t=pred()-.25))$surv,1))
    dr_late_se<-as.numeric(tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=pred-.25, probs.int.t=pred()-.25))$std.err,1))
    
    dr_early_ci <- dr_early + outer(dr_early_se, c(0,-qnorm(.975), qnorm(.975)), '*')
    dr_late_ci<- dr_late + outer(dr_late_se, c(0,-qnorm(.975), qnorm(.975)), '*')
    
### actual plot
    
    par(mfrow=c(3,3))
    
    # pooled
    plotCI(quantilemean+.005, deathrate_early[1,], ui=deathrate_early[3,], li=deathrate_early[2,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="", ylab="60 Days Death Rate",xlim=c(0, .6),ylim=c(.1, .83), yaxt="n", bty="n", main="AKIKI & IDEAL-ICU")
    plotCI(quantilemean-.005, deathrate_late[1,], ui=deathrate_late[3,], li=deathrate_late[2,], pch=4,cex=.8, gap=0, sfrac=0.003, col="#b24745", barcol="grey", xlab="", ylab="", add=TRUE)
    
    segments(trunc_pred,0, trunc_pred, .8, col="green")
    mtext(text="Early", side=1,line=3, at=-.2, adj=0)
    mtext(text="Delayed", side=1,line=5, at=-.2, adj=0)
    
    mtext(text=paste0(round(dr_early_ci[1],2)*100, "%", "  95% CI [", round(dr_early_ci[2],2)*100,"%-", round(dr_early_ci[3],2)*100,"%","]"), side=1,line=3, at=0.03, adj=0) 
    mtext(text=paste0(round(dr_late_ci[1],2)*100, "%", "  95% CI [", round(dr_late_ci[2],2)*100,"%-", round(dr_late_ci[3],2)*100,"%","]"), side=1,line=5, at=0.03, adj=0)
    
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
    
    segments(trunc_pred,0, trunc_pred, .8, col="green")
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
    
    segments(trunc_pred,0, trunc_pred, .8, col="green")
    axis(2, at=seq(.1,.8, by=.1), las=1)
    
    segments(0,.1, 0, .8, lty=3)
    for (i in 1:length(quantilecuts)) {
      segments(quantilecuts[i],.1, quantilecuts[i], .8, lty=3)
    }
    text(quantilemean, 0.83, paste("Q", 1:q, sep=""), cex=.8)
    
    lines(seq(0.02468, .55, by=.005), yy.er.idealicu[,1], col="#00a1d5")
    lines(seq(0.02468, .55, by=.005), yy2.er.idealicu[,1], col="#b24745")
    
    legend(.36, .35, inset=.05, legend=c("Early strategy", "Delayed strategy"),
           lty=c(1,1), pch=c(18, 4), col=c("#00a1d5", "#b24745"), lwd=1, box.lty=0, cex=.7)
    
    
    ### Hazard Ratio Plots ###
    ## HR plot : akiki & ideal-icu
    plot(1, xlim=c(0, .6), ylim = c(.3,3.7), yaxt="n", bty="n", type="n", log="y", xlab="", ylab="Hazard Ratio")
    
    abline(v=pred, col="green")
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
    
    abline(v=pred, col="green")
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
    
    abline(v=pred, col="green")
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
    plot(quantilemean, xlab="", ylab="Absolute Risk Reduction", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
    title(xlab="Predicted Probability of RRT Initiation Within 48Hours", cex.lab=.8)
    
    abline(v=pred, col="green")
    segments(0,ARD_overall[1], .6, ARD_overall[1], lty=2, lwd=2)
    lines(c(0,.6), y=c(0,0))
    plotCI(quantilemean, ARDd60_ci[1,], ui=ARDd60_ci[3,], li=ARDd60_ci[2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n", add=TRUE)
    
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
    plotCI(quantilemean, ARDd60_ci_akiki[1,], ui=ARDd60_ci_akiki[3,], li=ARDd60_ci_akiki[2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="", ylab="", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
    title(xlab="Predicted Probability of RRT Initiation Within 48Hours", cex.lab=.8)
    
    abline(v=pred, col="green")
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
    plotCI(quantilemean, ARDd60_ci_idealicu[1,], ui=ARDd60_ci_idealicu[3,], li=ARDd60_ci_idealicu[2,], pch=15, gap=0, sfrac=0.003, col="#79af97", barcol="black", xlab="", ylab="", main="", xlim=c(0, .6), ylim=c(-.35,.3), las=1, bty="n")
    title(xlab="Predicted Probability of RRT Initiation Within 48Hours", cex.lab=.8)
    
    abline(v=pred, col="green")
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

###    
    dev.off()
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 600,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
})