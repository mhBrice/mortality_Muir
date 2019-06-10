

# Function to compute risk ratio + CI + pvalue stars for Cox PH ####

risk_ph <- function(mod, var) {
  
  summ <- summary(mod)
 
  l_var <- match(var, names(mod$coefficients))
  risk_ci <- summ$conf.int[l_var,-2]
  colnames(risk_ci) <- c("risk", "ci.low", "ci.high")
  
  pval <- summ$coefficients[l_var,5]
  star <- stars.pval(pval)

  res <- cbind.data.frame(risk_ci, pval, star)
  
  row.names(res) <- var
  return(res)
}


# function to plot risk ratio ####

plot_risk <- function(mod, var, sp = NULL, lab = NULL, ylab = T) {
 
  # Get estimates, CI and p-value
  summ <- summary(mod)
  
  l_var <- match(var, names(mod$coefficients))
  n_var <- length(l_var)
  risk <- summ$conf.int[l_var,1]
  ci.low <- summ$conf.int[l_var,3]
  ci.high <- summ$conf.int[l_var,4]

  pval <- summ$coefficients[l_var,5]
  star <- stars.pval(pval)
  
  # Model fit
  rsq <- round(summ$rsq[1]*100, 2)
  ploglik <- summ$logtest[3]
  ploglik <- ifelse(ploglik < 0.001, "< 0.001", signif(ploglik, 1))
  
  # Color
  col_pt <- pval
  
  col_pt[which(pval<.05 & risk>1)] <- "dodgerblue4"
  col_pt[which(pval>.05 & risk>1)] <- "#C3D3E2"
  col_pt[which(pval<.05 & risk<1)] <- "#B41414"
  col_pt[which(pval>.05 & risk<1)] <- "#ECC4C4"
  
  # Plot
  plot(risk, ylim = c(.03,22), log = "y", col = "transparent", 
       ann=F, xaxt="n", yaxt="n",  bty = "l", frame.plot = T)
  
  abline(h=1, col = "grey65")
  
  # bars for confidence interval
  arrows(x0 = 1:n_var, 
         y0 = ci.low,
         y1 = ci.high, 
         angle = 90, code = 1, 
         length = 0, col = "grey65", lwd = 1.5, xpd = NA)
  #points
  points(risk, pch = 21, bg = col_pt, col = col_pt)
  
  # text
  text(x=1:n_var, y=rep(0.019,n_var), labels = lab, font=2,
       srt = 90, xpd = NA, adj = 1)

  axis(2, at=seq(0,20,.5), tcl= -0.2, labels=F, col = "grey35")
  axis(2, at=c(.1,1,5,10,20), labels = c(.1,1,5,10,20), las=1, cex.axis = .8)
  
  mtext(sp, adj = 0.95, cex=.85, font = 3, line = -.8)
  
  mtext(bquote("R"[2]~.(rsq)*"%"~italic("  p-value")~.(ploglik)), 
        adj = 0.95, cex=.6, font = 3, line = -2)
  
  if(ylab) mtext("Hazard ratio", 2, adj = 0.05, line = 2.2, cex=.85, font = 2, las=0)
}

