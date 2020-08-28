#!/usr/bin/python

'''
Purpose:
Initial test of EWS survivorship with filtered promoters OR DE genes
'''

import os
import glob
import lifelines
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------
def loadAltProm(FileName):
	DF = pd.read_csv(FileName,header="infer",sep="\t")
	return DF.Transcript.values

#------------------------------------------------------------------------------
def loadClinical(FileName):
	# load clinical
	ClinicalDF = pd.read_csv(FileName,sep="\t", header="infer")
	ClinicalDF = ClinicalDF[["EWS ID", "TP53_mutation", "STAG2_mutation", "OS (months)", "Patient Status"]]
	ClinicalDF.set_index(["EWS ID"], inplace=True, drop=True)
	return ClinicalDF

#------------------------------------------------------------------------------
def loadPredData(FileName, Transcripts, ClinicalDF):
	# load and filter transcripts
	DF = pd.read_csv(FileName,header="infer",sep="\t")
	DF = DF[DF.Transcript.isin(Transcripts)]
	# index by promoter + gene name
	DF["Promoter"] = DF.GeneName.astype(str)+"_"+DF.Transcript.astype(str)
	DF.set_index("Promoter",inplace=True,drop=True)
	# keep only sample columns
	DF = DF[[c for c in list(DF) if "EWS_T" in c and "Median_Pred_H3K27ac" in c]]
	RenameCols = ["_".join(c.split("_")[:2]) for c in list(DF)]
	DF.columns= RenameCols
	# combine with Clinical data
	DF = DF.T
	DF.index.name = ClinicalDF.index.name
	CDF = pd.concat([ClinicalDF,DF],axis=1,join="inner")
	# replace text with numerical binary
	CDF = CDF.replace("no",0).replace("yes",1)
	for Term, Num in zip(["Alive","DOD","DOC","DCD"],[0,1,1,1]):
		CDF = CDF.replace(Term,Num)
	print(CDF)
	return CDF

#------------------------------------------------------------------------------
def testCoxModel(FitDF, Promoter):
	# build/fit model
	CPH = CoxPHFitter()
	CPH.fit(FitDF, duration_col="OS (months)", event_col="Patient Status")
	return CPH.summary

#------------------------------------------------------------------------------
def runExperimentCox(PredDF):
	StaticCols = ["TP53_mutation", "STAG2_mutation", "OS (months)", "Patient Status"]
	PromoterList = [c for c in list(PredDF) if c not in StaticCols]
	# first, iterate over each promoter
	ResultsList = []
	for Promoter in PromoterList:
		TotalCols = StaticCols+[Promoter]
		TempPredDF = PredDF[TotalCols]
		PartialResultsDF = testCoxModel(TempPredDF, Promoter)
		PartialResultsDF["Test"] = Promoter
		PartialResultsDF["Variables"] = ["TP53_mutation", "STAG2_mutation", Promoter]
		ResultsList.append(PartialResultsDF)
	ResultsDF = pd.concat(ResultsList,axis=0)
	FDR_Results = multipletests(ResultsDF["p"].values,alpha=.05,method="fdr_bh")
	ResultsDF["reject"] = FDR_Results[0]
	ResultsDF["p.adj"] = FDR_Results[1]
	ResultsDF.to_csv("CoxModel_All_Results_Statistics.txt",
						header=True,sep="\t",index=False)
	Parsed_Results = ResultsDF[(ResultsDF["reject"]==True)&(~ResultsDF.Variables.isin(["TP53_mutation","STAG2_mutation"]))]
	Parsed_Results["Variables"] = Parsed_Results["Variables"].str.split("_").str[1]
	Parsed_Results.to_csv("CoxModel_Parsed_Results_Statistics.txt",
						header=True,sep="\t",index=False)
	print(Parsed_Results)
	print("Finished Running Cox Model Experiment")

#------------------------------------------------------------------------------
def main():
	Transcripts = loadAltProm("FilteredList.txt")
	ClinicalDF = loadClinical("StJudeClinical_wMatchup.txt")
	PredDF = loadPredData("EWS_Ensemble_Combined_Predictions.txt",Transcripts, ClinicalDF)
	runExperimentCox(PredDF)

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()