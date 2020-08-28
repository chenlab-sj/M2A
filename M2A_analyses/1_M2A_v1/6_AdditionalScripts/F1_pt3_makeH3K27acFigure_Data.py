#!/usr/bin/python
'''
Generate Figure plot input, including promoter class occupancy
'''
import numpy as np
np.random.seed(9)

import os
import sys
import glob
import pandas as pd
from itertools import combinations
from sklearn.mixture import GaussianMixture
#------------------------------------------------------------------------------
def loadData(FileRegex):
	# first load cancer concensus gene list
	CCDF = pd.read_csv("Census_allThu Feb 20 17_12_20 2020.tsv",header="infer",sep="\t")
	CCDF = CCDF[["Gene Symbol", "Tier"]]
	CCDF.set_index("Gene Symbol",inplace=True,drop=True)
	#
	IndexCols = ["EnsmblID_T","EnsmblID_G","Gene"]
	IndexCols2 = ["EnsmblID_T","EnsmblID_G","Gene", "Tier"]
	DropCols = ["Strand","Chr","Start","End","RStart","REnd"]
	# second load 2000bp window promoter region data
	PromoterDF = pd.read_csv("2000bpWindowRespVal_ALL/CombinedResults_2000bp.txt",
							header="infer",sep="\t")
	PromoterDF.set_index("Gene",inplace=True,drop=False)
	# combine CC with Promoter values
	TempDF = PromoterDF.merge(CCDF,left_index=True,right_index=True,how="left")
	TempDF = TempDF[~TempDF.duplicated(["EnsmblID_G"],keep=False)]
	Index = list(set(TempDF.index)&set(CCDF.index))
	CCDF = CCDF.ix[Index]
	# now combine final CCDF with promoter df
	PromoterDF = PromoterDF.merge(CCDF,left_index=True,right_index=True,how="left")
	PromoterDF["Tier"].replace(np.nan,0,inplace=True)
	PromoterDF["Tier"].replace(2,1,inplace=True)
	# reset index
	PromoterDF.set_index(IndexCols,inplace=True,drop=True)
	PromoterDF.drop(columns=DropCols,inplace=True)
	###### 
	FileList = glob.glob(FileRegex)
	# iterate eover files and combine
	DFList = []
	for Fle in FileList:
		Sample = Fle.split("_")[-4]
		DF = pd.read_csv(Fle,header="infer",sep="\t")
		DF.drop(columns=DropCols,inplace=True)
		DF.set_index(IndexCols,inplace=True,drop=True)
		DF = DF.merge(PromoterDF[[Sample, "Tier"]], left_index=True,right_index=True,how="inner")
		DF.reset_index(inplace=True,drop=False)
		DF.set_index(IndexCols2, inplace=True,drop=True)
		GMM = GaussianMixture(n_components=2, covariance_type='full').fit(DF[Sample].values.reshape(-1,1))
		if GMM.means_[0] > GMM.means_[1]:
			ClassPreds = GMM.predict(DF[Sample].values.reshape(-1,1)).flatten()
			DF["GMM-Class"] = np.where(ClassPreds==1,0,1)
		else:
			DF["GMM-Class"] = GMM.predict(DF[Sample].values.reshape(-1,1)).flatten()
		DF.columns = [Sample+"_"+c for c in list(DF)]
		DFList.append(DF)
	CDF = pd.concat(DFList,axis=1)
	return CDF

#------------------------------------------------------------------------------
def main():
	DF1 = loadData("H3K27ac_CellLine_Data/*FinalSignal.txt")
	DF1.to_csv("H3K27ac_CellLine_Data/H3K27ac_CellLine_Processed_FigureInput.txt",header=True,index=True,sep="\t")
	DF2 = loadData("H3K27ac_PDX_Data/*FinalSignal.txt")
	DF2.to_csv("H3K27ac_PDX_Data/H3K27ac_PDX_Processed_FigureInput.txt",header=True,index=True,sep="\t")

#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()