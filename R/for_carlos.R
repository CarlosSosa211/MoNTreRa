#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
## Path setting:
rm(list=ls())
options(java.parameters = "-Xmx4g")
## Required functions and packages
library(xlsx)
library(DMwR)
library(cvAUC)
library(pROC)
library(boot)
library(tidyverse)
library(caret)
#Path to store results (optional)
#File to read: Column_First= Patients ID, Column_Last=event (0 or 1), Other columns= input data to the model
#FILE=read.xlsx("../../Carlos/VPIPRO/Khemara/Tumor_76_ADC_T2_simp.xlsx",header = T, sheetIndex = 1)

data = read.table("/home/calcul/Carlos/Results/Recurrence/simp/rec_summary_8wTumVol.csv", sep = ',', header = TRUE)
for(i in 1:76){
  if(data$bio_rec[i] == 0){
    data$bio_rec[i] = "N"
  }
  else{
    data$bio_rec[i] = "Y"
  }
}
data$bio_rec = factor(data$bio_rec) 

reg = glm(bio_rec ~ max_tum_area + ADC_ave + T2w_ave + total_dose, data = data, family = "binomial")
summary(reg)
#Returns estimates and p-values
result.d = summary(reg)$coefficients[2,]
#Computes Odds Ratios with CIs
exp(result.d[1] - 1.96 * result.d[2])
exp(result.d[1])
exp(result.d[1] + 1.96 * result.d[2])
prob.pred = predict(reg, type = "response")
#Compute AUC
roc.AUC = roc(data$bio_rec ~ prob.pred, ci = TRUE)
roc.AUC
#Simple AUC plot
plot.roc(roc.AUC, col = "blue", lty = 1, main = "ROC AUC")

reg = glm(bio_rec ~ noHypNec_vascDensNoPref_038_ADCT2w, data = data,
          family = "binomial")
summary(reg)
result.d = summary(reg)$coefficients[2,]
exp(result.d[1] - 1.96 * result.d[2])
exp(result.d[1])
exp(result.d[1] + 1.96 * result.d[2])
prob.pred = predict(reg, type = "response")
cl <- ifelse(prob.pred > 0.5, "M", "R")
roc.AUC = roc(data$bio_rec ~ prob.pred, ci = TRUE)
roc.AUC
plot.roc(roc.AUC, add = TRUE, col = "green", lty = 1, main = "ROC AUC")

legend("bottomright", legend = c("Pre-treatment features and total dose",
                                 "In silico tumor area at t = 8 weeks"),
       col = c("blue", "green"), lty = 1:1, cex = 0.8)




reg = glm(FILE$bio_rec ~ FILE$max_tum_area + FILE$ADC_ave + FILE$T2w_ave + FILE$total_dose, data = FILE,
          family = "binomial")
summary(reg)
#Returns estimates and p-values
result.d = summary(reg)$coefficients[2,]
#Computes Odds Ratios with CIs
exp(result.d[1] - 1.96 * result.d[2])
exp(result.d[1])
exp(result.d[1] + 1.96 * result.d[2])
prob.pred = predict(reg, type = "response")
#Compute AUC
roc.AUC = roc(FILE$bio_rec ~ prob.pred, ci = TRUE)
roc.AUC
#Simple AUC plot
plot.roc(roc.AUC, col = "blue", lty = 1, main = "ROC AUC")
k = 3
kfCV = cv.glm(data = FILE, glmfit = reg, K=k)
kfCV$delta

myControl = trainControl(method = "repeatedcv", number = 3, repeats = 1000,
                         summaryFunction = twoClassSummary, classProbs = TRUE,
                         savePredictions = TRUE, sampling = "down")

model = train(bio_rec ~ tum_vol + T2w_ave + ADC_ave,
              data, method = "glm", trControl = myControl)
roc = (roc(predictor = model$pred$Y, response = model$pred$obs))
auc(predictor = model$pred$Y, response = model$pred$obs)
plot.roc(roc, col = "blue", lty = 1, main = "ROC AUC", xaxs = "i", yaxs = "i",
         print.auc = TRUE)

# model = train(bio_rec ~ init_tum_area_ADC,
#               data, method = "glm", trControl = myControl)
# roc = (roc(predictor = model$pred$Y, response = model$pred$obs))
# auc(predictor = model$pred$Y, response = model$pred$obs)
# plot.roc(roc, add = TRUE, col = "violet", lty = 1, main = "ROC AUC")

model = train(bio_rec ~ tum_area_from_vol + dens_ADCT2w,
              data, method = "glm", trControl = myControl)
roc = (roc(predictor = model$pred$Y, response = model$pred$obs))
auc(predictor = model$pred$Y, response = model$pred$obs)
plot.roc(roc, add = TRUE, col = "orange", lty = 1, main = "ROC AUC",
         print.auc = TRUE)

model = train(bio_rec ~ init_tum_area_tumVolADCT2w,
              data, method = "glm", trControl = myControl)
roc = (roc(predictor = model$pred$Y, response = model$pred$obs))
auc(predictor = model$pred$Y, response = model$pred$obs)
plot.roc(roc, add = TRUE, col = "green", lty = 1, main = "ROC AUC",
         print.auc = TRUE)

model = train(bio_rec ~ TTum330_alphaG1120_ADCT2w,
              data, method = "glm", trControl = myControl)
roc = (roc(predictor = model$pred$Y, response = model$pred$obs))
auc(predictor = model$pred$Y, response = model$pred$obs)
plot.roc(roc, add = TRUE, col = "violet", lty = 1, main = "ROC AUC",
         print.auc = TRUE)




