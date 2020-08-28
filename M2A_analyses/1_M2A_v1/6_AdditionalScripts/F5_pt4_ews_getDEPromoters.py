#!/usr/bin/python

'''
Usage:

'''

import os
import time
import numpy as np
import pandas as pd
import glob
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from scipy.stats import spearmanr
import statsmodels
from statsmodels.stats.multitest import fdrcorrection
# supressing copy in-place warnings
pd.options.mode.chained_assignment = None

#------------------------------------------------------------------------------
def loadData():
	# load clinical data
	ClinicalDF = pd.read_csv("StJudeClinical_wMatchup.txt",
						header="infer",sep="\t")
	ClinicalDF = ClinicalDF[["EWS ID", "TP53_mutation"]]
	#ClinicalDF = ClinicalDF[ClinicalDF["Relapse"]!="NA"]
	ClinicalDF["TP53"] = np.where(ClinicalDF["TP53_mutation"]=="yes",1,0)
	ClinicalDF.set_index(["EWS ID"],inplace=True,drop=True)
	ClinicalDF = ClinicalDF.T
	ClinicalDF.index.name = "Promoter"
	print(ClinicalDF)
	# load prediction data and clinical data
	PredDF = pd.read_csv("EWS_Combined_Predictions.txt",
						header="infer",sep="\t")
	PredDF.drop(columns=["Chr","Start","End","Strand"],inplace=True)
	# get multipromoter genes only
	PredDF = PredDF[PredDF.duplicated(subset=["GeneName"],keep=False)]
	# set index
	PredDF["Promoter"] = PredDF["GeneName"] +"_"+ PredDF["Transcript"]+"_"+ PredDF["Gene"]
	PredDF.set_index(["Promoter"],inplace=True,drop=True)	
	# filter Predictions for protein coding genes only
	PredDF = PredDF[PredDF["GeneType"]=="protein_coding"]
	# remove all non-pred cols
	PredCols = 	[c for c in list(PredDF) if "EWS_T" in c]
	PredOnlyDF = PredDF[PredCols]
	Cols = ["_".join(c.split("_")[:2]) for c in list(PredOnlyDF)]
	PredOnlyDF.columns = Cols 
	#combine with clinical DF
	CombinedDF = pd.concat([PredOnlyDF, ClinicalDF],axis=0,join="inner")
	# rename samples with class (survival, death)
	CombinedDF.columns = [Sample+"_"+str(int(i)) for Sample,i in zip(list(CombinedDF),CombinedDF.loc["TP53"])]
	CombinedDF = CombinedDF.iloc[:-2]
	print(CombinedDF)
	return CombinedDF

#------------------------------------------------------------------------------
def getColList(DFCols, ID):
	GroupColList = [c for c in DFCols if str(c.split("_")[-1])==str(ID)]
	print(str(ID)+" ColList ("+str(len(GroupColList))+"):", GroupColList)
	return GroupColList

#------------------------------------------------------------------------------
def getTtest(DF, Group1ID, Group2ID, FDR):
	#--------------------------------------------------------------------------
	def getPval(DF, Group1ID, Group2ID):
		#----------------------------------------------------------------------
		def calcPval(Row, Group1Cols, Group2Cols):
			# divide into two groups
			Group1Values = Row[Group1Cols].astype(np.float16).values
			Group2Values = Row[Group2Cols].astype(np.float16).values
			# return test statistics and pval (independent)
			Pval = ttest_ind(Group1Values, 
								Group2Values,
								equal_var=False, 
								nan_policy="propagate").pvalue
			# note: omit drops nans, and continues with ttest
			# propagate returns nan for any nans found
			return Pval

		#----------------------------------------------------------------------
		print("Caluclating Pval")
		PredCols = list(DF)
		# group 1 Column list
		Group1ColList = getColList(PredCols, Group1ID)
		# group 2 Column list
		Group2ColList = getColList(PredCols, Group2ID)
		# Paired Students t-test for pvalues
		DF["Pval"] = DF.apply(calcPval, args=(Group1ColList, Group2ColList), axis=1)
		TtestDF = DF[~pd.isnull(DF["Pval"])]
		# report removed genes due to ttest ommission
		NanOmmited = DF[pd.isnull(DF["Pval"])]
		NanOmmited["Pval"] = 1
		TtestDF = pd.concat([TtestDF, NanOmmited],axis=0)
		print("NaNs removed from ttest: "+str(NanOmmited.shape[0]))
		return TtestDF

	#--------------------------------------------------------------------------		
	def getAdjustedP(TtestDF, FDR):
		print("FDR adjusting Pval")
		# multiple test correction # i or p == "fdr_bh"
		Rejected, AdjustedPval = fdrcorrection(TtestDF["Pval"], 
												alpha=FDR,
												method="indep",
												is_sorted=False)
		# add additional columns from analysis
		TtestDF["Rejected"] = Rejected
		TtestDF["p.adj"] = AdjustedPval
		return TtestDF

	#--------------------------------------------------------------------------
	def formatDF(PadjDF, Group1ID, Group2ID):
		print("Formatting Pval DF")
		# get H3K27ac Cols
		PredCols = list(PadjDF)
		# group 1 Column list
		Group1ColList = getColList(PredCols, Group1ID)
		# group 2 Column list
		Group2ColList = getColList(PredCols, Group2ID)
		# Add mean value columns by group (group1 vs group2)
		Group1Name = "Mean_"+str(Group1ID)
		Group2Name = "Mean_"+str(Group2ID)
		PadjDF[Group1Name] = PadjDF[Group1ColList].mean(axis=1)
		PadjDF[Group2Name] = PadjDF[Group2ColList].mean(axis=1)
		# Reorder Dataframe
		ColOrder = [Group1Name, Group2Name, "Pval", "p.adj", "Rejected"]
		RawDataCols = Group1ColList + Group2ColList
		RawDataDF = PadjDF[RawDataCols]
		PadjDF = PadjDF[ColOrder]
		# sort pvals descending
		PadjDF.sort_values(by=["p.adj"], ascending=True, inplace=True)
		return PadjDF, RawDataDF
	
	#------------------------------------------------------------------------------
	def getUpDownReg(DF, Group1ID, Group2ID):
		Group1Name = "Mean_"+str(Group1ID)
		Group2Name = "Mean_"+str(Group2ID)
		# only passing FDR tests
		SigDF = DF[DF["Rejected"] == True]
		ExpSigDF = SigDF[np.logical_or(SigDF[Group1Name] > -1, SigDF[Group2Name] > -1)]
		# Look at absolute difference greater than 1
		UpDF = ExpSigDF[ExpSigDF[Group1Name] >= ExpSigDF[Group2Name]+1]
		UpDF["DE_Up"] = Group1ID
		DownDF = ExpSigDF[ExpSigDF[Group2Name] >= ExpSigDF[Group1Name]+1]
		DownDF["DE_Up"] = Group2ID
		# combine these two
		CombinedDF = pd.concat([UpDF, DownDF], axis=0)
		return CombinedDF

	#--------------------------------------------------------------------------
	# remove nans
	DF.dropna(inplace=True)
	# get Pvalues of t-test between two desired groups
	TtestDF = getPval(DF, Group1ID, Group2ID)
	# adjust pvalues for FDR
	PadjDF = getAdjustedP(TtestDF, FDR)
	# caluclate group means, sort cols and rows by p.adj
	FormatPadjDF, RawDataDF = formatDF(PadjDF, Group1ID, Group2ID)
	# return DE transcripts by thresholding, at least one expressed/active
	DE_PadjDF = getUpDownReg(FormatPadjDF, Group1ID, Group2ID)
	return FormatPadjDF, DE_PadjDF, RawDataDF

#------------------------------------------------------------------------------
def main():
	FDR = .05
	# start timer
	T0 = time.time()
	# load raw H3K27ac preds
	PredH3K27ac_DF = loadData()
	# there are multiple groups to test, 
	# TEST1: DE in groups: -------------------------------------------
	Group1ID = "1" # Event happend
	Group2ID = "0" # Event didn't happen
	ResultsName = "Pred_Mutation_ttest.txt"
	ALL_Survival, DE_Survival, RAW_Survival = getTtest(PredH3K27ac_DF.copy(), Group1ID, Group2ID, FDR)
	ALL_Survival.to_csv("Mean_ALL_"+ResultsName, header=True,index=True,sep="\t")
	print("Finished: "+ResultsName)
	#
	T1 = time.time()
	print("Total Time: ", str(int(T1-T0))+"s")

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()