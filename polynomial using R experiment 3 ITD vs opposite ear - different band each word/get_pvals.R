get_pvals <- function(model){
  # only argument is the model that you're evaluating
  # it returns the summary, but with p values
  
  coefs <- as.data.frame(coef(summary(model)))
  coefs$p <- format.pval(2*(1-pnorm(abs(coefs[,"t value"]))), digits=2, eps=0.0001) #also make the p-values a bit more readable
  
  coefs$star <- ifelse(coefs$p < 0.001, "***", 
                       ifelse(coefs$p < 0.01, "**", 
                              ifelse(coefs$p < 0.05, "*", 
                                     ifelse(coefs$p < 0.1, ".", ""))))
  coefs$Estimate <- round(coefs$Estimate,4)
  coefs[,"Std. Error"] <- round(coefs[,"Std. Error"],4)
  coefs[,"t value"] <- round(coefs[,"t value"],2)
  
  return(coefs)
}