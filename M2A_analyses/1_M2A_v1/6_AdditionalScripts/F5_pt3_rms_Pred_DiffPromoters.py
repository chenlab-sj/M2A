#!/usr/bin/python

'''
get matching predictions for ground truth "actual" determined values
'''

import os
import time
import numpy as np
import pandas as pd
import glob
from scipy.stats import ranksums
import statsmodels
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import mannwhitneyu
# supressing copy in-place warnings
pd.options.mode.chained_assignment = None


#------------------------------------------------------------------------------
def loadActualData(FileName):
	AltPromDF = pd.read_csv(FileName, header="infer", sep="\t")
	#AltPromDF = AltPromDF[AltPromDF["Rejected_Act-Diff"] == True]
	AltPromDF = AltPromDF[["Gene", "Transcript", "GeneName", "ARMS_ActMean","ERMS_ActMean",
							"ARMS_ACTUAL_ORDER", "ERMS_ACTUAL_ORDER",
							"Pval_Act-Diff", "p.adj_Act-Diff",
							"Mean_ARMS_Act-Diff","Mean_ERMS_Act-Diff"]]
	AltPromDF.set_index(["Gene", "Transcript", "GeneName"], inplace=True, drop=True)
	return AltPromDF

#------------------------------------------------------------------------------
def getTypeLookUp(FileID):
	TypeLookUpTable = {'Mast111':"ERMS", 
						'Mast139':"ERMS", 
						'Mast71':"ERMS", 
						'Mast74':"ERMS", 
						'Mast39':"ERMS",
						'Mast85':"ERMS",
						"Mast100":"ERMS",
						'Rh73':"ERMS", 
						'Rh70':"ERMS",
						'Rh74':"ERMS",
						'Rh75':"ERMS",
						'Mast118':"ARMS", 
						'Mast60':"ARMS",
						'Mast95':"ARMS",
						'Mast35':"ARMS"}
	if FileID not in list(TypeLookUpTable.keys()):
		Value = "NA" # skip MAST161, MAST100, Rh73, etc?
	else:
		Value = TypeLookUpTable[FileID]
	return Value

#------------------------------------------------------------------------------
def loadPredData(Dir):
	# load gene target list
	#TargetDF = pd.read_csv("pax3_fox01_targets.txt",header="infer",sep="\t")
	# load prediction files
	RegexPath = Dir+"/*_predictions.txt"
	PredFilePaths = glob.glob(RegexPath)
	# build dict of all data
	UseCols = ["Gene", "Transcript", "GeneName" ,"Regressor_Pred"]
	DFList = []
	for FP in PredFilePaths:
		# split on filename for identifiers
		FileIDList = FP.split("/")[-1].split("_")
		FileID = FileIDList[1]
		RHBType = getTypeLookUp(FileID)
		if RHBType != "NA":
			ColName_PredAc = "_".join([RHBType,FileID,"Pred-H3K27ac"])
			# load DF
			DF = pd.read_csv(FP, header="infer", usecols=UseCols, sep="\t")
			print("loading: ",ColName_PredAc,DF.shape)
			DF = DF[UseCols]
			DF.columns = UseCols[:-1] + [ColName_PredAc]
			# Set Index
			DF.set_index(["Gene", "Transcript", "GeneName"], inplace=True, drop=True)
			# add to list
			DFList.append(DF)
	H3K27ac_DF = pd.concat(DFList, axis=1)
	return H3K27ac_DF

#------------------------------------------------------------------------------
def getSubtypeMean(PredDF):
	ERMSPredCols = [c for c in list(PredDF) if "ERMS" in c and "Pred-H3K27ac" in c]
	ARMSPredCols = [c for c in list(PredDF) if "ARMS" in c and "Pred-H3K27ac" in c]
	# get subtype promoter mean
	PredDF["ARMS_PredMean"] = PredDF[ERMSPredCols].mean(axis=1)
	PredDF["ERMS_PredMean"] = PredDF[ARMSPredCols].mean(axis=1)
	return PredDF

#------------------------------------------------------------------------------
def filterPredData(ActualDF,SubTypeMeanDF):
	FilterDF = pd.concat([ActualDF,SubTypeMeanDF],axis=1,join="inner")
	return FilterDF

#------------------------------------------------------------------------------
def calcPromoterDiff(RankOneDF):
	# get actual H3K27ac samples
	SampleList = [c for c in list(RankOneDF) if "Pred-H3K27ac" in c]
	DFList = []
	for Sample in SampleList:
			DF = RankOneDF[[Sample, "ARMS_ACTUAL_ORDER", "ERMS_ACTUAL_ORDER",
						"ARMS_PredMean","ERMS_PredMean",
						"ARMS_ActMean","ERMS_ActMean",
							"Pval_Act-Diff", "p.adj_Act-Diff",
							"Mean_ARMS_Act-Diff","Mean_ERMS_Act-Diff"]]
			GroupList = []
			for Group in DF.groupby(["Gene"]):
				TempDF = Group[1]
				ARMS_Promoter = TempDF[TempDF["ARMS_ACTUAL_ORDER"]==1][Sample].values
				ARMS_Promoter = ARMS_Promoter.mean()
				ERMS_Promoter = TempDF[TempDF["ERMS_ACTUAL_ORDER"]==1][Sample].values
				ERMS_Promoter = ERMS_Promoter.mean()
				Diff = ARMS_Promoter - ERMS_Promoter
				if Diff:
					ColName_Diff = Sample.replace("H3K27ac","Diff")
					TempDF[ColName_Diff] = np.repeat(Diff,TempDF.shape[0])
					GroupList.append(TempDF)
				else:
					ColName_Diff = Sample.replace("H3K27ac","Diff")
					TempDF[ColName_Diff] = np.repeat(0.0,TempDF.shape[0])
					GroupList.append(TempDF)
			#------				
			DF = pd.concat(GroupList,axis=0)		
			DF.reset_index(drop=False,inplace=True)
			DF.set_index(["Gene", "Transcript", "GeneName", 
				"ARMS_ACTUAL_ORDER", "ERMS_ACTUAL_ORDER",
						"ARMS_PredMean","ERMS_PredMean",
						"ARMS_ActMean","ERMS_ActMean",
							"Pval_Act-Diff", "p.adj_Act-Diff",
							"Mean_ARMS_Act-Diff","Mean_ERMS_Act-Diff"], 
				inplace=True,drop=True)
			# add to list
			DFList.append(DF)
	CombinedDF = pd.concat(DFList, axis=1)
	return CombinedDF

#------------------------------------------------------------------------------
def getColList(ColsList, ID):
	GroupColList = [c for c in ColsList if ID in c.split("_")[0]]
	return GroupColList

#------------------------------------------------------------------------------
def getWilcRS(DF, Type, Group1ID, Group2ID, FDR):
	#--------------------------------------------------------------------------
	def getPval(Type, DF, Group1ID, Group2ID):
		#----------------------------------------------------------------------
		def calcPval(Row, Group1Cols, Group2Cols):
			# divide into two groups
			print(Row[Group1Cols].values)
			print(Row[Group2Cols].values)
			Group1Values = np.round(Row[Group1Cols].astype(np.float16).values,3)#.astype(np.float16).values
			Group2Values = np.round(Row[Group2Cols].astype(np.float16).values,3)#.astype(np.float16).values
			Stats = mannwhitneyu(Group1Values, 
								Group2Values,
								alternative="two-sided")
			return Stats[1] # pval

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
		PadjDF["Abs_"+Group1Name] = PadjDF[Group1ColList].abs().mean(axis=1)
		PadjDF["Abs_"+Group2Name] = PadjDF[Group2ColList].abs().mean(axis=1)
		# Reorder Dataframe
		ColOrder = ["Abs_"+Group1Name, "Abs_"+Group2Name, Group1Name, Group2Name, 
					"Pval_"+Type, "p.adj_"+Type, "Rejected_"+Type]
		RawDataCols = Group1ColList + Group2ColList
		RawDataDF = PadjDF[RawDataCols]
		PadjDF = PadjDF[ColOrder]
		# sort pvals descending
		PadjDF.sort_values(by=["p.adj_"+Type], ascending=True, inplace=True)
		return PadjDF, RawDataDF

	#--------------------------------------------------------------------------
	# remove nans
	#DF.dropna(inplace=True)
	# get DF Type:
	TypeCols = [c for c in list(DF) if Type in c]
	DF = DF[TypeCols]
	# remove duplicates
	DF.reset_index(inplace=True,drop=False)
	####
	IndexDF = DF[["Gene", "Transcript", "GeneName", 
				"ARMS_ACTUAL_ORDER", "ERMS_ACTUAL_ORDER",
				"ARMS_PredMean","ERMS_PredMean",
				"ARMS_ActMean","ERMS_ActMean",
				"Pval_Act-Diff", "p.adj_Act-Diff",
				"Mean_ARMS_Act-Diff","Mean_ERMS_Act-Diff"]]
	IndexDF.set_index(["Gene"], drop=True,inplace=True)
	###
	Duplicated = DF.duplicated(subset=["Gene","GeneName"]+TypeCols, keep="first")
	DF = DF[~Duplicated]
	DF.drop(["Transcript", "ARMS_ACTUAL_ORDER", "ERMS_ACTUAL_ORDER",
				"ARMS_PredMean","ERMS_PredMean",
				"ARMS_ActMean","ERMS_ActMean",
				"Pval_Act-Diff", "p.adj_Act-Diff",
				"Mean_ARMS_Act-Diff","Mean_ERMS_Act-Diff"], axis=1, inplace=True)
	DF.set_index(["Gene"], drop=True,inplace=True)
	####
	print("Pval pipeline dimensions:",DF.shape)
	# get Pvalues of t-test between two desired groups
	TtestDF = getPval(Type, DF, Group1ID, Group2ID)
	# adjust pvalues for FDR
	PadjDF = getAdjustedP(Type, TtestDF, FDR)
	# caluclate group means, sort cols and rows by p.adj
	FormatPadjDF, RawDataDF = formatDF(Type, PadjDF, Group1ID, Group2ID)
	# recombine data
	CDF = pd.concat([IndexDF,FormatPadjDF],axis=1,join="inner")
	return CDF

#------------------------------------------------------------------------------
def main():
	FDR = .1
	# load actual
	ActualDF = loadActualData("2b_RMS_AltPromoter_RESULTS_DFw_PVAL.txt")
	# load  preds
	PredDF = loadPredData("2_NBL_transfer_predictions")
	# get subtype mean preds
	SubTypeMeanDF = getSubtypeMean(PredDF)
	# filter pred data by transcripts in actual
	FilteredPred = filterPredData(ActualDF,SubTypeMeanDF)
	# get pred diff using actual definitions 
	AltPromDF = calcPromoterDiff(FilteredPred)
	AltPromDF.to_csv("3extra_RMS_PRED_DIFF_RESULTS_DF.txt",header=True,index=True,sep="\t")
	# get pval for these pred diffs
	PvalDF = getWilcRS(AltPromDF, "Pred-Diff", "ARMS","ERMS", FDR)
	PvalDF.to_csv("3b_RMS_ACTUAL_vs_PRED_RESULTS_DF.txt",header=True,index=True,sep="\t")

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()                                                                                                                                                                                                                                                                                                                                                                                                     