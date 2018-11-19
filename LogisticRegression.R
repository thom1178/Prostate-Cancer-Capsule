# Jamel M. Thomas

#=================================================================#
#===========================Packages==============================#
library(ISLR)
library(pROC)
library(ROCR)
library(corrplot)
library(MASS) # for stepAIC
library(boot) # for cross-validation
library(vioplot) #Violin plots
library(VIM) #missing data plot
library(dplyr)

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
#===================================================================#
#===========================Data Input==============================#
#Prostateate Cancer Data
Prostate <- read.table("~/Desktop/prostate.txt", header = T)

# Prostateate cancer, along with skin cancer, is the most common cancer among men. 
#It can often be treated if detected early enough. 
#The Ohio State University Comprehensive Cancer Center lead a study to determine
# if baseline exam measurements can predict whether a tumor has 
# penetrated the Prostateatic capsule.

#Response Variable is capsule:

head(Prostate) 
dim(Prostate) 
n = dim(Prostate)[1]  # sample size
names(Prostate)
sum(is.na(Prostate))# missing values
which(is.na(Prostate))
set.seed(1)

aggr(Prostate, prop = F, numbers = T, bars = F, col = c("skyblue", "red", "orange")) 
#Plot missing data
Prostate[which(is.na(Prostate)),] #We will need to remove all of these values
Prostate<-na.omit(Prostate); dim(Prostate) # 376  9

str(Prostate)
Prostate$race <- as.factor(Prostate$race)
Prostate$dpros <- as.factor(Prostate$dpros)
Prostate$dcaps <- as.factor(Prostate$dcaps)

str(Prostate)
summary(Prostate)
attach(Prostate)

#=============================================================#
#========================== EDA ==============================#
# PSA is a measure of a protein produced by prostate gland cells. 
# Elevated levels may suggest prostate
# cancer and is thus used as a screening test. 
# The Gleason score is a scale measuring the abnormality of cells.
# Larger values suggest higher risk of cancer

plot(capsule) #Response is 0 or 1. Consider a logistic regression. 
#
table(capsule);
barplot(table(Prostate$capsule), col = c("cyan", "red"), names = c("No", "Yes"), xlab ="Capsule Penetration",
        ylim = c(0,225), cex.axis = 1.5) 

par(mfrow = c(1,2))
# Box-plots of covariates over good
plot(factor(capsule), age, col=c("cyan", "red"), varwidth=T,names = c("No", "Yes"),
     ylab="Age",xlab="Capsule Penetration", cex.lab=1.4) #Does not indicate difference in plot

#Indicate difference in plot
plot(factor(capsule), psa, col=c("cyan", "red"), varwidth=T,names = c("No", "Yes"),
     ylab="PSA",xlab="Capsule Penetration", 
     cex.lab=1.4) 

plot(factor(capsule), gleason, col=c("cyan", "red"), varwidth=T,names = c("No", "Yes"),
     ylab="Gleason",xlab="Capsule Penetration", cex.lab=1.4) #Possible differences

plot(factor(capsule), vol, col=c("cyan", "red"), varwidth=T, names = c("No", "Yes"),
     ylab="Volume",xlab="Capsule Penetration", cex.lab=1.4) #Does not indicate difference in plot but possible
par(mfrow = c(1,1))
barplot(table(capsule, dcaps), beside = T, col = c("cyan", "red"), 
        xlab = "Detection of Capsular Involvement", ylab = "Frequency", 
        names = c("No", "Yes"), las = 1.2, cex.lab = 1.4)
legend("topright", fill = c("cyan", "red"), legend = c("No Penetration", "Penetration"),
       bty = "n", cex = 1.5)

barplot(table(capsule, dpros), beside = T, col = c("cyan", "red"), 
        xlab = "Results of Digital Rectal Exam", ylab = "Frequency", 
        names = c("No Nodule", "Nodule (left)", "Nodule (right)", "Bilobar"), las = 1.2, cex.lab = 1.4)
legend("topright", fill = c("cyan", "red"), legend = c("No Penetration", "Penetration"),
       bty = "n", cex = 1.5)
barplot(table(capsule, race), beside = T, col = c("cyan", "red"), 
        xlab = "Race", ylab = "Frequency", 
        names = c("White Males", "Black Males"), las = 1.2, cex.lab = 1.4)
legend("topright", fill = c("cyan", "red"), legend = c("No Penetration", "Penetration"),
       bty = "n", cex = 1.5)

hist(gleason[capsule==0], main = "", xlab = "Gleason of those without Capsule Penetration", 
     xlim = c(0,10),
     ylim = c(0,100), breaks = 40)
hist(gleason[capsule==1], main = "", xlab = "Gleason of those with Capsule Penetration", 
     xlim = c(0,10),
     ylim = c(0,100), breaks = 40, cex.lab=1.5)

hist(psa[capsule==0], main = "", xlab = "PSA without Capsule Penetration", xlim = c(0,140), 
     ylim = c(0,80),
     breaks = 20, col = "cyan")
hist(psa[capsule==1], main = "", xlab = "PSA with Capsule Penetration", ylim = c(0,80),
     breaks = 20, col = "red")

hist(vol[capsule==0], main = "", xlab = "Volume without Capsule Penetration", xlim = c(0,100),
     breaks = 20, col = "cyan")
hist(vol[capsule==1], main = "", xlab = "Volume with Capsule Penetration", xlim = c(0,100),
     breaks = 20, col = "red")

#Correlation plot based only on numerical variables
Prostate$capsule <- as.numeric(Prostate$capsule)
vol2 <- Prostate$vol[Prostate$vol > 0]
cap2 <- Prostate$capsule[Prostate$vol > 0]
cor(cap2, vol2)

Prostcor = cor(Prostate[,c(2,3,7,8,9)])
corrplot(Prostcor) 
#PSA and Gleason possible mostly correlated. May need to remove one

t.test(log(psa[capsule==0]),log(psa[capsule==1])) 
#Test for differences in psa by capsule
# Enough evidence to suggest differences
t.test(age[capsule==0], age[capsule==1]) #Test for differences in psa by capsule
# No evidence to suggest differences

plot(psa, capsule)


## Chi-square test for categorical variables
chisq.test(x=dcaps, y=capsule) #There is an association between capsule and dcaps

#========================================================================#
#========================== Model building ==============================#

# Consider Logistic Regression on all main effect & interactions
fit.all=glm(capsule~. + .*., family=binomial(link=logit), data=Prostate[,-1]) 
summary(fit.all)

#Backwards Model selection
aic.models <- stepAIC(fit.all)

#AIC Model includes:
  # age
  # race
  # dpros
  # dcaps
  # psa
  # vol
  # gleason
  # age:gleason
  # race:volume
  # dpros:vol
  # dcaps:psa
  # dcaps:gleason

# AIC Model
fit.aic <- glm(aic.models$formula,
                family = binomial(link = logit), data = Prostate[,-1])
summary(fit.aic)
fit.aic
# Due to correlations with PSA, we will remove Gleason and it's interactions

fit.aic2 <- glm(formula = capsule ~ age + race + dpros + dcaps + psa + vol + 
                       + race:vol + dpros:vol + dcaps:psa,  
                         family = binomial(link = logit),
                              data = Prostate[,-1])
summary(fit.aic2)
# age is the least significant, consider a model without it.
fit.aic3 <- glm(formula = capsule ~ race + dpros + dcaps + psa + vol + 
                   race:vol + dpros:vol ,  
                family = binomial(link = logit),
                data = Prostate[,-1])
summary(fit.aic3)
# dpros:volume is the least significant
#       consider removing it.
fit.aic4 <- glm(formula = capsule ~  race + dpros + psa + vol + 
                   race:vol,
                family = binomial(link = logit),
                data = Prostate[,-1])
summary(fit.aic4)
anova(fit.aic3, fit.aic4, test = "Chisq") #Sufficient evidence to remove this interaction
# race:volume is significant at the alpha = 0.15 level.
# We consider this a first final model:

#=================================================================================#
#========================== Model AIC 4 Diagnostics ==============================#
probsAIC4 = predict.glm(fit.aic4, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC4, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC4)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c4 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint by min sq error 

table((probsAIC4 > c4), capsule) #Table with optimal cutpoint
plot.aic(probsAIC4, c4, string = "")
#We would like to lower Predicted False but actually true: False Negatives

#=================================================================================#
#========================== Model AIC 5 Diagnostics ==============================#

#We Consider the model without the interactions of Race:Volume
fit.aic5 <- glm(formula = capsule ~  race + dpros + psa + vol,  
                family = binomial(link = logit), data = Prostate[,-1])
summary(fit.aic5)

#Same Analysis as above:
probsAIC5 = predict.glm(fit.aic5, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC5, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC5)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c5 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC5 > c5), capsule) #Table with optimal cutpoint
plot.aic(probsAIC5, c5, string = "")
anova(fit.aic4, fit.aic5, test = "Chisq")
#The table prediction Yield the same results as before,
#Chi Sq test indicates no statistical differences between the two
#We continue with this as the new final model. 

#We then consider removing Race:
#=================================================================================#
#========================== Model AIC 6 Diagnostics ==============================#

#We Consider the model without the interactions of Race
fit.aic6 <- glm(formula = capsule ~ dpros + psa + vol,
                family = binomial(link = logit), data = Prostate[,-1])
summary(fit.aic6)

#Same Analysis as above:
probsAIC6 = predict.glm(fit.aic6, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC6, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1, by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC6)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c6 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC6 > c6), capsule) #Table with optimal cutpoint
plot.aic(probsAIC6, c6, "Model 1")
anova(fit.aic5, fit.aic6, test = "Chisq")

#The model differences are not significant, so we will remove Race and 
#   and proceed with model 6

#This is our first final model:

#=================================================================================#
#========================== Model AIC 7 Diagnostics ==============================#
# Consider if we took out psa instead of gleason:
fit.aic7 <- glm(formula = capsule ~ age + race + dpros + dcaps  + vol + gleason +
                  + race:vol + dpros:vol,  
                family = binomial(link = logit),data = Prostate[,-1])
summary(fit.aic7)

#Same Analysis as above:
probsAIC7 = predict.glm(fit.aic7, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC7, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC7)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c7 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC7 > c7), capsule) #Table with optimal cutpoint

anova(fit.aic6, fit.aic7, test = "Chisq")


interaction.plot(dpros, vol, capsule)


#Consider removing: Age
#=================================================================================#
#========================== Model AIC 8 Diagnostics ==============================#

fit.aic8 <- glm(formula = capsule ~ race + dpros + dcaps  + vol + gleason +
                  + race:vol + dpros:vol,  
                family = binomial(link = logit),data = Prostate[,-1])
summary(fit.aic8)

#Same Analysis as above:
probsAIC8 = predict.glm(fit.aic8, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC8, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC8)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c8 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC8 > c), capsule) #Table with optimal cutpoint

anova(fit.aic7, fit.aic8, test = "Chisq")

plot.aic(probsAIC8, c8, string = "")
table((probsAIC8 > c8), capsule) #Table with optimal cutpoint


#The chi-sq test indicates that the two models are not statisticlly different at 
# the alpha = 0.05 level.

#=================================================================================#
#========================== Model AIC 9 Diagnostics ==============================#

#We Consider the model without the interactions of race:volume
fit.aic9 <- glm(formula = capsule ~ race + dpros + dcaps  + vol + gleason +
                 + dpros:vol,  
                family = binomial(link = logit),data = Prostate[,-1])
summary(fit.aic9)

#Same Analysis as above:
probsAIC9 = predict.glm(fit.aic9, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC9, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC9)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c9 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC9 > c9), capsule) #Table with optimal cutpoint

anova(fit.aic7, fit.aic9, test = "Chisq")

#It follows that we now consider taking out : race
#=================================================================================#
#========================== Model AIC 10 Diagnostics ==============================#

#We Consider the model without the interactions of race
fit.aic10 <- glm(formula = capsule ~ dpros + gleason + vol + dpros:vol,  
                family = binomial(link = logit),data = Prostate[,-1])
summary(fit.aic10)

#Same Analysis as above:
probsAIC10 = predict.glm(fit.aic10, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC10, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC10)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c10 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC10 > c10), capsule) #Table with optimal cutpoint

anova(fit.aic9, fit.aic10, test = "Chisq")

#We keep this model and consider taking out: dpros:volume
#=================================================================================#
#========================== Model AIC 11 Diagnostics ==============================#

#We Consider the model without the interactions of Dpros:volume
fit.aic11 <- glm(formula = capsule ~ dpros + gleason + vol ,  
                 family = binomial(link = logit),data = Prostate[,-1])
summary(fit.aic11)

#Same Analysis as above:
probsAIC11 = predict.glm(fit.aic11, data = Prostate, type="response")

# ROC curve
#Prediction function for Optimal cutpoint and ROC curves
pred = prediction(probsAIC11, Prostate[,"capsule"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")

plot(roc.perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7)) 
abline(a=0, b= 1, col="gray")

auc(roc(Prostate$capsule,probsAIC11)) #Area under the ROC Curve

# Compute optimal cutoff
print(opt.cut(roc.perf, pred))
c11 <- opt.cut(roc.perf, pred)[3] #Optimal Cutpoint

table((probsAIC11 > c11), capsule) #Table with optimal cutpoint

anova(fit.aic10, fit.aic11, test = "Chisq")
#The models are statistically different, we must keep dpros:volume 

#Therefore our final model is Model 10

#Final Model1 includes: model6
# dpros
# psa
# vol

#Final Model2 includes: model 10
# dpros
# gleason
# vol
# dpros:vol

fit.final1 <- fit.aic6
plot.aic(probsAIC6, c6, string = "Model 1")

fit.final2 <- fit.aic10
plot.aic(probsAIC10, c10, string = "Model 2")
table((probsAIC10 > c10), capsule) #Table with optimal cutpoint

#Consider threshold of .15
plot.aic(probsAIC6,.2, string = "Model 1") #Threshold 
table((probsAIC6 > .2), capsule) #Table with optimal cutpoint

plot.aic(probsAIC10,.25, string = "Model 2") #Threshold 
table((probsAIC10 > .25), capsule) #Table with optimal cutpoint

#The final model is the one with dpros, psa, and volume.
#=================================================================================#
#========================== Final Model Diagnostics ==============================#

# Take a look at diagnostics for final PK data model.
#Create EVPs by bining continuous covariates. 
g  <- 3 # Number of categories
psa_interval <- cut(psa, quantile(psa, 0:g/g), include.lowest = T)
levels(psa_interval)
quantile(vol, 0:g/g)
#Volume has a lot of zeros, need to bind manually
vol_interval <- cut(vol, c(-1, 0, 20, 30, 100)) #for g=5
levels(vol_interval)

w <- aggregate( formula = capsule ~ dpros + gleason + vol_interval + dpros:vol_interval, data = Prostate, FUN = sum)
n <- aggregate( formula = capsule ~ dpros + gleason + vol_interval + dpros:vol_interval, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule, digits = 2))
mod.prelim <- glm(formula = capsule/trials ~dpros + gleason + vol_interval + dpros:vol_interval, 
                  family = binomial(link = logit), data = w.n, weights = trials)

save1=examine.logistic.reg(mod.prelim, identify.points=T, 
                           scale.n=one.fourth.root, scale.cookd=sqrt)
#The Residuals vs estimated probabilites are centered around zero. 
# They indicate mostly homoscedastic behaviour and 
# points 62, and 73 may be infuential.

#The Cooks distance vs. leverage indicates a large cooks distance for point 62, at a t

names(save1)

# Let's examine the outliers/influential points idenitfied by the plots
# Store prediction, residual, Cook's D, and leverage values
w.n.diag1=data.frame(w.n, pi.hat=round(save1$pi.hat, 2), std.res=round(save1$stand.resid, 2), 
                     cookd=round(save1$cookd, 2), h=round(save1$h, 2))
p=length(mod.prelim$coef) # number of parameters in model (# coefficients)
# Set thresholds for residuals, Cooks' distance, and leverage 
# Recall now working with EVPs, not complete data set.  
# So denominators here are the #EVPs=124, not n=1438
ck.out=abs(w.n.diag1$std.res)>2 | w.n.diag1$cookd>4/nrow(w.n) | w.n.diag1$h > 3*p/nrow(w.n)
extract.EVPs=w.n.diag1[ck.out, ]
extract.EVPs[order(extract.EVPs$vol),]  # order by vol
w.n.diag1[c(24, 36, 38, 58), ]
#=====================================================================================#
#========================== Coefficient Interpretations ==============================#
round(exp(fit.final2$coefficients), digits = 4)
get.or.se(fit.aic10)#standard errors of OR

# Intercept
# DPros2 : The odds of capusle penetration of unilobar nodule on left increases by 3.61x 
#         compared to those  without a nodule (baseline), assuming all other variables are fixed.
# DPros3 : The odds of capusle penetration of unilobar nodule on right increases by 4.1x  
#         compared to those without a nodule, assuming all other variables are fixed.
# DPros4 : The odds of capusle penetration of bilobar nodule increases by 12.15x compared to those 
#         Without a nodule (baseline), assuming all other variables are fixed.
# Gleason:   For every 1 unit increase in age (1-10), the odds of capusle penetration increases by 3.34x
#         Assuming all other variables are fixed.
# Vol:   For every 1 unit increase in vol (cm^3), the odds of capsule penetration 
#          increases by 0.85% Assuming all other variables are fixed.
# Dpros2*Vol:  The difference between the odds ratio corresponding to a change in volume by 
#           1 cm^3 amongst the baseline and the the odds ratio corresponding to an increase in
#             volume by 1 year amongst those with a left nodule. 


#detach(Prostate)
#rm(list = ls())

