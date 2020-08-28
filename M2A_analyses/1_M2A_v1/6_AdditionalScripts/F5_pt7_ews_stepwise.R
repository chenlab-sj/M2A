# module load R/3.6.1

# backward stepwise cox model selection for relevant alt promoters and DE genes

library('generics')
library('broom')
library('ggsci')
library('ggsignif')
library('polynom')
library('ggpubr')
library('exactRankTests')
library('mvtnorm')
library('maxstat')
library('KMsurv')
library('km.ci')
library('survMisc')
library('survminer')
library(My.stepwise)
library(survival)
library(dplyr)

Variables = c('TP53_mutation', 'STAG2_mutation',
			'TEX40_ENST00000328404.6', 'TNS1_ENST00000446903.1',
			'RET_ENST00000479913.1', 'CASZ1_ENST00000496432.2',
			'SLC27A6_ENST00000508645.1', 'CALCB')

DF = read.table("StepwiseModel_Input.txt", header=TRUE, sep="\t")
SurvObj <- Surv(time=DF$OS..months., event=DF$Patient.Status)
Formula = as.formula(paste("SurvObj ~", paste(Variables[1:8],collapse=" + ")))
FullCoxphModel <- coxph(Formula, data=DF)
# 
StepResults <- step(FullCoxphModel,direction="backward")
summary(StepResults)
StepCoxForestPlot <- ggforest(StepResults, main = "Ewing's sarcoma hazard ratio analysis", data=DF, cpositions = c(0.005, 0.25, 0.38), fontsize=.7)
ggsave("05_Final_Coxph_Results_Step.png", plot=StepCoxForestPlot, device="png", dpi = 300)

