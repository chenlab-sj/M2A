#!/usr/bin/python

'''
In observed data, find diff primary promtoers
'''

import os
import time
import numpy as np
import pandas as pd
import glob

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
def loadData(Dir):
	# load prediction files
	RegexPath = Dir+"/*_predictions.txt"
	PredFilePaths = glob.glob(RegexPath)
	# build dict of all data
	UseCols = ["Gene", "Transcript", "GeneName", "Mapability", "Chr", "Start", "End", "H3K27Ac"]
	DFList = []
	for FP in PredFilePaths:
		# split on filename for identifiers
		FileIDList = FP.split("/")[-1].split("_")
		FileID = FileIDList[1]
		RHBType = getTypeLookUp(FileID)
		if RHBType != "NA":
			ColName_Ac = "_".join([RHBType,FileID,"Act-H3K27ac"])
			# load DF
			DF = pd.read_csv(FP, header="infer", usecols=UseCols, sep="\t")
			print("loading: ", ColName_Ac,DF.shape)
			DF = DF[UseCols]
			DF.columns = UseCols[:-1] + [ColName_Ac]
			# Set Index
			DF.set_index(["Gene", "Transcript", "GeneName", "Mapability", 
							"Chr", "Start", "End"], inplace=True, drop=True)
			# add to list
			DFList.append(DF)
	H3K27ac_DF = pd.concat(DFList, axis=1)
	return H3K27ac_DF

#------------------------------------------------------------------------------
def getMultiPromoterGenes(H3K27ac_DF):
	Genes = H3K27ac_DF.index.get_level_values("Gene")
	MultiPromDF = H3K27ac_DF[Genes.duplicated(keep=False)]
	Genes = MultiPromDF.index.get_level_values("Gene")
	print("Multi promoter gene promoters", MultiPromDF.shape)
	print("Multi promoter genes", len(Genes.unique()))
	return MultiPromDF

#------------------------------------------------------------------------------
def getMapabilityFilterGenes(MultiPromDF):
	Map = MultiPromDF.index.get_level_values("Mapability")
	MapDF = MultiPromDF[Map>0.75]
	Genes = MapDF.index.get_level_values("Gene")
	print("Map>.75 promoter gene promoters", MapDF.shape)
	print("Multi promoter genes", len(Genes.unique()))
	return MapDF

#------------------------------------------------------------------------------
def getSubtypeMean(MapDF):
	# get subtype columns
	ERMSActCols = [c for c in list(MapDF) if "ERMS" in c and "Act-H3K27ac" in c]
	ARMSActCols = [c for c in list(MapDF) if "ARMS" in c and "Act-H3K27ac" in c]
	ERMSPredCols = [c for c in list(MapDF) if "ERMS" in c and "Pred-H3K27ac" in c]
	ARMSPredCols = [c for c in list(MapDF) if "ARMS" in c and "Pred-H3K27ac" in c]
	# get subtype promoter mean
	MapDF["ARMS_ActMean"] = MapDF[ARMSActCols].mean(axis=1)
	MapDF["ERMS_ActMean"] = MapDF[ERMSActCols].mean(axis=1)
	# get all multi promoters expressed in at least 1 subtype:
	ActiveDF = MapDF[(MapDF.ARMS_ActMean>1)|(MapDF.ERMS_ActMean>1)]
	# get unique genes with active promoters in at least 1 subtype
	Genes = ActiveDF.index.get_level_values("Gene")
	# sanity check the number of active genes (by gene ID, not GeneName)
	print("Actual active genes in at least 1 subtype (subtype mean > 1 )", len(Genes.unique()))
	return ActiveDF

#------------------------------------------------------------------------------
def getAlternatePromoters(SubTypeMeanDF, Type, FDR):
	#--------------------------------------------------------------------------
	def getRankOnePromoters(DF,Type):
		if Type == "PRED":
			SampleType = "Pred"
		else:
			SampleType = "Act"
		DF.reset_index(inplace=True,drop=False)
		# add transcript orders
		DF["ARMS_"+Type+"_ORDER"] = DF.groupby(["Gene"])["ARMS_"+SampleType+"Mean"].rank(ascending=False, method="min")
		DF["ERMS_"+Type+"_ORDER"] = DF.groupby(["Gene"])["ERMS_"+SampleType+"Mean"].rank(ascending=False, method="min")
		# Promoters to Keep: RanK = 1
		DF = DF[(DF["ARMS_"+Type+"_ORDER"]==1)|(DF["ERMS_"+Type+"_ORDER"]==1)]
		DF.set_index(["Gene", "Transcript", "GeneName", "Mapability", 
			"Chr", "Start", "End"], inplace=True, drop=True)
		DF = DF[["ARMS_"+Type+"_ORDER", "ERMS_"+Type+"_ORDER"]]
		return DF

	#------------------------------------------------------------------------------
	def calcPromoterDiff(RankOneDF, SubTypeMeanDF, Type):
		if Type == "PRED":
			SampleType = "Pred-H3K27ac"
			ARMS_Mean = "ARMS_PredMean"
			ERMS_Mean = "ERMS_PredMean"
		else:
			SampleType = "Act-H3K27ac"
			ARMS_Mean = "ARMS_ActMean"
			ERMS_Mean = "ERMS_ActMean"
		# combine DFs 
		CDF = pd.concat([RankOneDF, SubTypeMeanDF],axis=1,join="inner")
		# get actual H3K27ac samples
		SampleList = [c for c in list(SubTypeMeanDF) if SampleType in c]
		DFList = []
		for Sample in SampleList:
				DF = CDF[[Sample, "ARMS_"+Type+"_ORDER", "ERMS_"+Type+"_ORDER",
							ARMS_Mean,ERMS_Mean]]
				GroupList = []
				for Group in DF.groupby(["Gene"]):
					TempDF = Group[1]
					ARMS_Promoter = TempDF[TempDF["ARMS_"+Type+"_ORDER"]==1][Sample].values
					ARMS_Promoter = ARMS_Promoter.mean()
					ERMS_Promoter = TempDF[TempDF["ERMS_"+Type+"_ORDER"]==1][Sample].values
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
				DF.set_index(["Gene", "Transcript", "GeneName", "Mapability", 
					"Chr", "Start", "End", "ARMS_"+Type+"_ORDER", "ERMS_"+Type+"_ORDER",
					ARMS_Mean,ERMS_Mean], 
					inplace=True,drop=True)
				# add to list
				DFList.append(DF)
		CombinedDF = pd.concat(DFList, axis=1)
		return CombinedDF

	#------------------------------------------------------------------------------
	def getSubtypeDiffAve(CombinedDF, Type):
		if Type == "PRED":
			SampleType = "Pred-Diff"
		else:
			SampleType = "Act-Diff"
		ERMSCols = [c for c in list(CombinedDF) if "ERMS" in c and SampleType in c]
		ARMSCols = [c for c in list(CombinedDF) if "ARMS" in c and SampleType in c]
		# get ave:
		CombinedDF["ARMS_MEAN_"+SampleType] = CombinedDF[ARMSCols].mean(axis=1)
		CombinedDF["ERMS_MEAN_"+SampleType] = CombinedDF[ERMSCols].mean(axis=1)
			# sanity check gene count actual with alternate promoter
		print("Samples: "+Type)
		DiffDF = CombinedDF[CombinedDF["ARMS_MEAN_"+SampleType]>0]
		Genes = DiffDF.index.get_level_values("Gene")
		# sanity check the number of active genes (by gene ID, not GeneName)
		print("genes with different primary promoters", len(Genes.unique()))
		print(CombinedDF)
		return CombinedDF

	#--------------------------------------------------------------------------
	RankOneDF = getRankOnePromoters(SubTypeMeanDF.copy(),Type)
	# for each sample, find difference in these promoters
	CombinedDF = calcPromoterDiff(RankOneDF, SubTypeMeanDF, Type)
	# get subtype promtoer diff mean
	DiffMeanDF = getSubtypeDiffAve(CombinedDF, Type)
	return DiffMeanDF

#------------------------------------------------------------------------------
def main():
	FDR = .1
	# load raw actual and preds
	H3K27ac_DF = loadData("2_NBL_transfer_predictions")
	# get multi promoters
	MultiPromDF = getMultiPromoterGenes(H3K27ac_DF)
	# filter by mapability
	MapDF = getMapabilityFilterGenes(MultiPromDF)
	# determine subtype mean H3K27ac for each promoter
	SubTypeMeanDF = getSubtypeMean(MapDF)
	# get alternate promoters 
	AltPromDF = getAlternatePromoters(SubTypeMeanDF,"ACTUAL",FDR)
	AltPromDF.to_csv("1b_RMS_AltPromoter_RESULTS_DF.txt",header=True,index=True,sep="\t")
	print(PDZRN3DF)
#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()                                                                                                                                                                                                                                                                                                                                                                                                     