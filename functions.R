#==================================================================#
#===========================Functions==============================#
#Optimal Cutpoint Function
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
} 

#ROC Plot function
rocplot=function(pred, truth, ...){
  predob = prediction (pred, truth)
  perf = performance (predob , "tpr", "fpr") 
  plot(perf ,...)}

#One Fourth Root For Examine Plots
one.fourth.root=function(x){
  x^0.25
} #For extermal model diagnostics sourced below
source("~/Desktop/examine.logistic.reg.R") #Model diagnostics

#Plot predicted probabilities and thresholds separated by capsul penetration
plot.aic <- function(model.prob, c,string){
  dx <- density(model.prob[capsule ==0])$x[density(model.prob[capsule ==0])$x >=0 &
                                             density(model.prob[capsule ==0])$x <=1]
  dy <- density(model.prob[capsule ==0])$y[density(model.prob[capsule ==0])$x >=0 &
                                             density(model.prob[capsule ==0])$x <=1]
  zx <- density(model.prob[capsule ==1])$x[density(model.prob[capsule ==1])$x >=0 &
                                             density(model.prob[capsule ==1])$x <=1]
  zy <- density(model.prob[capsule ==1])$y[density(model.prob[capsule ==1])$x >=0 &
                                             density(model.prob[capsule ==1])$x <=1]
  plot(dx,dy , col = "blue", xlim = c(0,1), type = "l", main = string,
       xlab = "Predicted probabilities",
       ylab = "Density") 
  lines(zx, zy, col = "red") 
  abline(v = c, lty =2, col = "green")
  rug(model.prob[capsule ==0], col = "blue")
  rug(model.prob[capsule ==1], col = "red")
  legend("topright", bty = "n", legend = c("No Penetration", "Penetration", "Threshold"), 
         fill = c("blue", "red", "green"))
}
#Odds ratio SE
get.or.se <- function(model) {
  broom::tidy(model) %>% 
    mutate(or = exp(estimate),
           var.diag = diag(vcov(model)),
           or.se = sqrt(or^2 * var.diag)) %>%
    select(or.se) %>% unlist %>% unname
}
