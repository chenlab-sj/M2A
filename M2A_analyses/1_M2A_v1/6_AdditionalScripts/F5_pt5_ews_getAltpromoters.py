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
def getRankOnePromoters():
	# load actual
	DF = pd.read_csv("ALL_Pred_Mutation_ttest.txt",header="infer",sep="\t")
	# Edit columns
	DF[["GeneName","Transcript","Gene"]] =  pd.DataFrame(DF.Promoter.str.split("_",expand=True,n=2))
	# add transcript orders
	DF["NoMutation_Pred_Order"] = DF.groupby(["Gene"])["Mean_0"].rank(ascending=False)
	DF["Mutation_Pred_Order"] = DF.groupby(["Gene"])["Mean_1"].rank(ascending=False)
	# Promoters to Keep: RanK = 1
	DF = DF[(DF["NoMutation_Pred_Order"]==1)|(DF["Mutation_Pred_Order"]==1)]
	DF.set_index(["Gene","Transcript","GeneName"],inplace=True,drop=False)
	DF = DF[["NoMutation_Pred_Order", "Mutation_Pred_Order"]]
	print("RankOne DF:")
	print(DF)
	return DF

#------------------------------------------------------------------------------
def loadPredData():
	# load prediction data and clinical data
	PredDF = pd.read_csv("EWS_Combined_Predictions.txt",
						header="infer",sep="\t")
	# set index to remove promoter from the operable df
	PredDF["Promoter"] = PredDF["GeneName"] +"_"+ PredDF["Transcript"]+"_"+ PredDF["Gene"]
	PredDF.set_index(["Promoter"],inplace=True,drop=True)
	# remove GMM cols
	PredCols = 	[c for c in list(PredDF) if "EWS_T" in c]
	PredOnlyDF = PredDF[PredCols]
	Cols = ["_".join(c.split("_")[:2]) for c in list(PredOnlyDF)]
	PredOnlyDF.columns = Cols 
	#######
	# load clinical data
	ClinicalDF = pd.read_csv("StJudeClinical_wMatchup.txt",
						header="infer",sep="\t")
	ClinicalDF = ClinicalDF[["EWS ID", "TP53_mutation"]]
	ClinicalDF["TP53"] = np.where(ClinicalDF["TP53_mutation"]=="yes",1,0)
	ClinicalDF.set_index(["EWS ID"],inplace=True,drop=True)
	ClinicalDF = ClinicalDF.T
	ClinicalDF.index.name = "Promoter"
	print(ClinicalDF)
	#######
	#combine with clinical DF
	CombinedDF = pd.concat([PredOnlyDF, ClinicalDF],axis=0,join="inner")	
	# set columns based on patient status
	CombinedDF.columns = [Sample+"_"+str(int(i)) for Sample,i in zip(list(CombinedDF),CombinedDF.loc["TP53"])]
	# remove survival data
	CombinedDF = CombinedDF.iloc[:-2]
	# reset index and expand data
	CombinedDF.reset_index(inplace=True,drop=False)
	CombinedDF[["GeneName","Transcript","Gene"]] =  pd.DataFrame(CombinedDF.Promoter.str.split("_",expand=True,n=2))
	# drop old index
	CombinedDF.drop(columns=["Promoter"],inplace=True)
	# and adopt new expanded index
	CombinedDF.set_index(["Gene","Transcript","GeneName"],inplace=True, drop=True)
	# note we don't need to filter all genes out in this data frame, because the unique index join (inner)
	# will remove extraneous genes/transcripts/etc.
	print("Preds data shape")
	print(CombinedDF)
	return CombinedDF

#------------------------------------------------------------------------------
def calcPromoterDiff(RankOneDF, PredDF):
	# build dict of all data
	UseCols = ["Gene", "Transcript", "GeneName", "Regressor_Pred"]
	DFList = []
	for Sample in list(PredDF):
		# split on filename for identifiers
		GroupID = Sample.split("_")[-1]
		SampleName = "_".join(Sample.split("_")[:-1])
		ColName_Ac = Sample+"_Pred_H3K27ac"
		print("Processing: ", ColName_Ac)
		# split DF into single column and merge with RANK one
		DF = PredDF[[Sample]]
		DF.columns = [ColName_Ac]
		DF = DF.merge(RankOneDF, left_index=True,right_index=True,how="inner")
		# iterate over gene groups
		GroupList = []
		for Group in DF.groupby(["Gene"]):
			TempDF = Group[1]
			Survived_Promoter = TempDF[TempDF["NoMutation_Pred_Order"]==1][ColName_Ac].values
			Died_Promoter = TempDF[TempDF["Mutation_Pred_Order"]==1][ColName_Ac].values
			Diff = Survived_Promoter - Died_Promoter
			if Diff:
				ColName_Diff = Sample+"_Diff"
				TempDF[ColName_Diff] = np.repeat(Diff,TempDF.shape[0])
				GroupList.append(TempDF)
			else:
				ColName_Diff = Sample+"_Diff"
				TempDF[ColName_Diff] = np.repeat(0.0,TempDF.shape[0])
				GroupList.append(TempDF)
		DF = pd.concat(GroupList,axis=0)		
		DF.reset_index(drop=False,inplace=True)
		DF.set_index(["Gene","Transcript","GeneName", 
					"NoMutation_Pred_Order", "Mutation_Pred_Order"],
					inplace=True,drop=True)
		# add to list
		DFList.append(DF)
	CombinedDF = pd.concat(DFList, axis=1)
	return CombinedDF

#------------------------------------------------------------------------------
def getColList(ColsList, ID):
	GroupColList = [c for c in ColsList if ID in c.split("_")[-2]]
	return GroupColList

#------------------------------------------------------------------------------
def getTtest(DF, Type, Group1ID, Group2ID, FDR):
	#--------------------------------------------------------------------------
	def getPval(Type, DF, Group1ID, Group2ID):
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
		# group 1 Column list
		Group1ColList = getColList(list(DF), Group1ID)
		# group 2 Column list
		Group2ColList = getColList(list(DF), Group2ID)
		# Paired Students t-test for pvalues
		DF["Pval_"+Type] = DF.apply(calcPval, args=(Group1ColList, Group2ColList), axis=1)
		TtestDF = DF[~pd.isnull(DF["Pval_"+Type])]
		# report removed genes due to ttest ommission
		NanOmmited = DF[pd.isnull(DF["Pval_"+Type])]
		NanOmmited["Pval_"+Type] = 1
		TtestDF = pd.concat([TtestDF, NanOmmited],axis=0)
		print("NaNs removed from ttest: "+str(NanOmmited.shape[0]))
		return TtestDF

	#--------------------------------------------------------------------------		
	def getAdjustedP(Type, TtestDF, FDR):
		print("FDR adjusting Pval")
		# multiple test correction # i or p == "fdr_bh"
		Rejected, AdjustedPval = fdrcorrection(TtestDF["Pval_"+Type], 
												alpha=FDR,
												method="indep",
												is_sorted=False)
		# add additional columns from analysis
		TtestDF["Rejected_"+Type] = Rejected
		TtestDF["p.adj_"+Type] = AdjustedPval
		return TtestDF

	#--------------------------------------------------------------------------
	def formatDF(Type, PadjDF, Group1ID, Group2ID):
		print("Formatting Pval DF")
		# get H3K27ac Cols
		TypeCols = [c for c in list(PadjDF) if Type in c]
		# group 1 Column list
		Group1ColList = getColList(TypeCols, Group1ID)
		# group 2 Column list
		Group2ColList = getColList(TypeCols, Group2ID)
		# Add mean value columns by group (group1 vs group2)
		Group1Name = "Mean_"+str(Group1ID)+"_"+Type
		Group2Name = "Mean_"+str(Group2ID)+"_"+Type
		PadjDF[Group1Name] = PadjDF[Group1ColList].mean(axis=1)
		PadjDF[Group2Name] = PadjDF[Group2ColList].mean(axis=1)
		# Reorder Dataframe
		ColOrder = [Group1Name, Group2Name, "Pval_"+Type, "p.adj_"+Type, "Rejected_"+Type]
		RawDataCols = Group1ColList + Group2ColList
		RawDataDF = PadjDF[RawDataCols]
		PadjDF = PadjDF[ColOrder]
		# sort pvals descending
		PadjDF.sort_values(by=["p.adj_"+Type], ascending=True, inplace=True)
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
	# get DF Type:
	TypeCols = [c for c in list(DF) if Type in c]
	DF = DF[TypeCols]
	# remove duplicates for FPKM
	if Type == "Diff":
		DF.reset_index(inplace=True,drop=False)
		Duplicated = DF.duplicated(subset=["Gene","GeneName"]+TypeCols, 
									keep="first")
		DF = DF[~Duplicated]
		DF.drop(["Transcript", "NoMutation_Pred_Order",
					"Mutation_Pred_Order"], axis=1, inplace=True)
		DF.set_index(["Gene", "GeneName"], drop=True,inplace=True)
	# get Pvalues of t-test between two desired groups
	TtestDF = getPval(Type, DF, Group1ID, Group2ID)
	# adjust pvalues for FDR
	PadjDF = getAdjustedP(Type, TtestDF, FDR)
	# caluclate group means, sort cols and rows by p.adj
	FormatPadjDF, RawDataDF = formatDF(Type, PadjDF, Group1ID, Group2ID)
	return FormatPadjDF

#------------------------------------------------------------------------------
def main():
	Group1ID = "1"
	Group2ID = "0"
	FDR = .05
	# rank promoter activity
	RankOneDF = getRankOnePromoters()
	# load preds dataframe
	PredDF = loadPredData()
	# for each sample, find difference in these promoters
	CombinedDF = calcPromoterDiff(RankOneDF, PredDF)
	CombinedDF.to_csv("1_Predicted_Promoter_Differences.txt",header=True,sep="\t",index=True)
	# get Ttest pval, adjusted
	DIFF_K27ac_ALL_ERMSvsARMS = \
				getTtest(CombinedDF.copy(), "Diff", Group1ID, Group2ID, FDR)
	DIFF_K27ac_ALL_ERMSvsARMS.to_csv("2_Predicted_TTEST_Promoter_Differences.txt",header=True,sep="\t",index=True)		
	print("Pipeline Completed")


#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    