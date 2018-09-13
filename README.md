# Prostate-Cancer-Capsule
Prostate cancer is one type of cancer found in males that can be treated if detected in the early stages of
development. In this paper, we obtain prostate cancer screening data from Ohio State University and build
a logistic regression model to classify the capsule penetration of a tumor. We begin with an exploratory
analysis of the data, then build a model with all covariates and interactions considered. From this model,
we utilize backward stepwise model selection with AIC. From this model, we obtain covariates that are
correlated, so we build two models and proceed with model selecting individually. We compare the two
models, and discuss the issues with both, and selecting our final model based on the lowest misclassification
rate. We then discuss the diagnostic plots of the model and interpret the odds ratios. This final model
contains the result of a digital rectal exam, with four level of nodules, the volume of a tumor, Gleason
Score, and an interaction term of the volume of a tumor and the results of the digital rectal exam. This
model performs only slightly better than its competitor in the paper, and we discuss some other models to
improve prediction accuracy.
