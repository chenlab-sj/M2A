#!/usr/bin/python

'''
determine significance of actual values
'''

import os
import time
import numpy as np
import pandas as pd
import glob
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import mannwhitneyu

#------------------------------------------------------------------------------
def loadData(FileName):
	AltPromDF = pd.read_csv(FileName, header="infer", sep="\t")
	AltPromDF.set_index(["Gene", "Transcript", "GeneName", "Mapability", 
					"Chr", "Start", "End"], inplace=True, drop=True)
	return AltPromDF

#------------------------------------------------------------------------------
def getActive(AltPromDF):
	ActiveDF = AltPromDF[(AltPromDF.ARMS_ActMean>1)&(AltPromDF.ERMS_ActMean>1)]
	Genes = ActiveDF.index.get_level_values("Gene")
	print("Active promoter gene promoters", ActiveDF.shape)
	print("Active promoter genes", len(Genes.unique()))
	return ActiveDF

#------------------------------------------------------------------------------
#ARMS_ACTUAL_ORDER	ERMS_ACTUAL_ORDER
def getAltProm(ActiveDF):
	PrimAltPromDF = ActiveDF[ActiveDF.ARMS_ACTUAL_ORDER!=ActiveDF.ERMS_ACTUAL_ORDER]
	Genes = PrimAltPromDF.index.get_level_values("Gene")
	print("Primary alternate promoter gene promoters", PrimAltPromDF.shape)
	print("Primary alternate promoter genes", len(Genes.unique()))
	return PrimAltPromDF

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
			Group1Values = np.round(Row[Group1Cols].values,3)#.astype(np.float16).values
			Group2Values = np.round(Row[Group2Cols].values,3)#.astype(np.float16).values
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
	Duplicated = DF.duplicated(subset=["Gene","GeneName"]+TypeCols, keep="first")
	DF = DF[~Duplicated]
	DF.drop(["Transcript", "Mapability", "Chr", "Start", "End"], axis=1, inplace=True)
	DF.set_index(["Gene","GeneName"], drop=True,inplace=True)
	####
	print("Pval pipeline dimensions:",DF.shape)
	# get Pvalues of t-test between two desired groups
	TtestDF = getPval(Type, DF, Group1ID, Group2ID)
	# adjust pvalues for FDR
	PadjDF = getAdjustedP(Type, TtestDF, FDR)
	# caluclate group means, sort cols and rows by p.adj
	FormatPadjDF, RawDataDF = formatDF(Type, PadjDF, Group1ID, Group2ID)
	return FormatPadjDF

#------------------------------------------------------------------------------
def main():
	FDR = .1
	# load raw actual and preds
	AltPromDF = loadData("1b_RMS_AltPromoter_RESULTS_DF.txt")
	# different primary promoters
	PrimAltPromDF = getAltProm(AltPromDF)
	# get pvalue from wilcoxon rank sum test
	PvalDF = getWilcRS(PrimAltPromDF, "Act-Diff", "ARMS","ERMS",FDR)
	#PvalDF.to_csv("2_RMS_AltPromoter_RESULTS_DFw_PVAL.txt",header="infer",sep="\t",index=True)
	PvalDF.reset_index(drop=False,inplace=True)
	print(PvalDF[PvalDF["Rejected_Act-Diff"]==True]["Gene"].unique().shape)
	PvalDF.set_index(["Gene","GeneName"], drop=True,inplace=True)
	# combine with promoter level details
	PrimAltPromDF.reset_index(inplace=True,drop=False)
	PrimAltPromDF.set_index(["Gene","GeneName"], drop=True,inplace=True)
	CDF = PvalDF.merge(PrimAltPromDF,left_index=True,right_index=True,how="outer")
	CDF.to_csv("2b_RMS_AltPromoter_RESULTS_DFw_PVAL.txt",header="infer",sep="\t",index=True)
#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()