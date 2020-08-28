#!/usr/bin/python
'''
Plot data and make venn diagram
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
def sortDF(DF, DataName):
	# parse the samples, activatoin classes, and PromoterRegion column names
	SampleList = np.unique([c.split("_")[0] for c in list(DF)]).tolist()
	ClassGMMList = np.array([c for c in list(DF) if "_GMM-Class" in c])
	PromoterSampleList = np.array([n+"_"+n for n in SampleList])
	# print a report of the samples being prcessed
	print(SampleList)
	# hard code the anotations and the use-index "ParseList"
	AnotList = ["(+,+,+)", "(+,+,-)", "(+,-,+)", "(+,-,-)", "(-,+,+)", "(-,+,-)", "(-,-,+)"]
	ParseList = np.array([[1,1,1],[1,1,0],[1,0,1],[1,0,0],[0,1,1],[0,1,0],[0,0,1]])
	#
	#
	# Iterate over ParseList, slice dataframe based on class activations from parse list
	DFList = []
	HLineList = []
	HSum = 0
	RowY = []
	CC_AllActive = []
	CC_DiffActive = []
	for Parse in ParseList:
		# get active GMM columns
		ActClassGMMCols = ClassGMMList[Parse==1]
		# get inactive GMM columns
		NonActClassGMMCols = ClassGMMList[Parse==0]
		# return all active promoter columns
		SampleCols = PromoterSampleList[Parse==1]
		# slice Activity DF by active and non active columns from Parse list
		nDF = DF[(DF[ActClassGMMCols].mean(axis=1)==1)&(DF[NonActClassGMMCols].sum(axis=1)==0)]
		# determine the mean value of all active promoter columns
		nDF["Mean"] = nDF[SampleCols].mean(axis=1)
		# sort by the mean of active promoter columns
		nDF.sort_values("Mean", ascending=False,inplace=True)
		# remove non-activity windowed columns
		nDF.drop(columns=PromoterSampleList.tolist()+ClassGMMList.tolist(),inplace=True)
		# print report
		print(nDF.shape)
		#
		# retain all positional information for annotation
		# get the horizontal group division line
		HSum = HSum + nDF.shape[0]
		RowY.append(HSum-(nDF.shape[0]/2))
		HLineList.append(HSum)
		# capture CC genes
		if sum(Parse) == 3:
			nDF["All_vs_DA"] = 0
			CC_AllActive.append(nDF.index.get_level_values(3).values.sum())
		else:
			nDF["All_vs_DA"] = 1
			CC_DiffActive.append(nDF.index.get_level_values(3).values.sum())
		# keep sliced dataframe
		DFList.append(nDF)
	#print CC gene report
	print("DataName:",DataName)
	print("CC_AllActive", CC_AllActive)
	print("CC_DiffActive", CC_DiffActive)
	print("Sum CC_AllActive:",sum(CC_AllActive))
	print("Sum CC_DiffActive:",sum(CC_DiffActive))
	Percentage = sum(CC_DiffActive)/sum(CC_AllActive + CC_DiffActive)
	print("Percentage:", Percentage)
	#----------------
	# merge all dataframes in order
	SortDF = pd.concat(DFList,axis=0)
	SortDF.to_csv(DataName+".txt",header=True,sep="\t",index=True)
	# remove the index. no longer needed, seaborn will try to print it
	DropCols = ["EnsmblID_T","EnsmblID_G","Gene", "Tier", "Mean", "All_vs_DA"]
	SortDF.reset_index(inplace=True,drop=False)	
	TempSortDF = SortDF.sort_values(["Gene","Mean"],ascending=[True,False])
	TempSortDF = TempSortDF[(TempSortDF["Tier"]==1)&(~TempSortDF["Gene"].duplicated())]
	CC_RowY = list(TempSortDF.index)
	SortDF.drop(columns=DropCols,inplace=True)
	SortDF.index.name=""
	# calculate all positional information for annotation
	# iterate over samples for vertical line and Sample name placement
	VLineList = []
	ColX = []
	VSum = 0
	for Sample in SampleList:
		VSum = VSum + len([c for c in list(SortDF) if Sample in c])
		ColX.append(VSum-(len([c for c in list(SortDF) if Sample in c])/2))
		VLineList.append(VSum)
	######
	return SortDF, DFList, [HLineList, VLineList, SampleList, ColX, AnotList, RowY, CC_RowY]

#------------------------------------------------------------------------------
def combineDataTypes(DF1,DF2):
	# Retain all columns for the first sample from two groups: CellLine and PDX
	DFList = []
	for DF in [DF1,DF2]:
		SampleList = np.unique([c.split("_")[0] for c in list(DF)]).tolist()
		DFCols = [c for c in list(DF) if SampleList[0] in c]
		DF = DF[DFCols]
		DFList.append(DF)
	# Combine this list into as single DF
	DF = pd.concat(DFList,axis=1,join="inner")
	# parse the samples, activation classes, and PromoterRegion column names
	SampleList = np.unique([c.split("_")[0] for c in list(DF)]).tolist()
	print("Combined SampleList", SampleList)
	ClassGMMList = np.array([c for c in list(DF) if "_GMM-Class" in c])
	PromoterSampleList = np.array([n+"_"+n for n in SampleList])
	# hard code the anotations and the use-index "ParseList"
	AnotList = ["(+,+)", "(+,-)", "(-,+)"]
	ParseList = np.array([[1,1],[1,0],[0,1]])
	#
	#
	# Iterate over ParseList, slice dataframe based on class activations from parse list
	DFList = []
	HLineList = []
	HSum = 0
	RowY = []
	for Parse in ParseList:
		# get active GMM columns
		ActClassGMMCols = ClassGMMList[Parse==1]
		# get inactive GMM columns
		NonActClassGMMCols = ClassGMMList[Parse==0]
		# return all active promoter columns
		SampleCols = PromoterSampleList[Parse==1]
		# slice Activity DF by active and non active columns from Parse list
		nDF = DF[(DF[ActClassGMMCols].mean(axis=1)==1)&(DF[NonActClassGMMCols].sum(axis=1)==0)]
		# determine the mean value of all active promoter columns
		nDF["Mean"] = nDF[SampleCols].mean(axis=1)
		# sort by the mean of active promoter columns
		nDF.sort_values("Mean", ascending=False,inplace=True)
		# remove non-activity windowed columns
		nDF.drop(columns=["Mean"]+PromoterSampleList.tolist()+ClassGMMList.tolist(),inplace=True)
		# print report
		print(nDF.shape)
		#
		# retain all positional information for annotation
		# get the horizontal group division line
		HSum = HSum + nDF.shape[0]
		RowY.append(HSum-(nDF.shape[0]/2))
		HLineList.append(HSum)
		# keep sliced dataframe
		DFList.append(nDF)
	#----------------
	# merge all dataframes in order
	SortDF = pd.concat(DFList,axis=0)
	# remove the index. no longer needed, seaborn will try to print it
	DropCols = ["EnsmblID_T","EnsmblID_G","Gene", "Tier"]
	SortDF.reset_index(inplace=True,drop=False)
	CC_RowY = SortDF.index[SortDF["Tier"]==1].tolist()
	SortDF.drop(columns=DropCols,inplace=True)
	SortDF.index.name=""
	#
	#
	# calculate all positional information for annotation
	# iterate over samples for vertical line and Sample name placement
	VLineList = []
	ColX = []
	VSum = 0
	for Sample in SampleList:
		VSum = VSum + len([c for c in list(SortDF) if Sample in c])
		ColX.append(VSum-(len([c for c in list(SortDF) if Sample in c])/2))
		VLineList.append(VSum)
	######
	return SortDF, DFList, [HLineList, VLineList, SampleList, ColX, AnotList, RowY, CC_RowY]

#------------------------------------------------------------------------------
def makeCombinedHeatmap(DFList, AnnotList_List, RawDFList, PlotName):
	def makeBoxplot(CellLineDF, PDX_DF, Ax):
		#--------------------------------------------------------------------------
		def getData(CellLineDF, PDX_DF):
			# combine both dataframes
			CellLine_Cols = [c for c in list(CellLineDF) if "GMM" in c]
			PDX_Cols = [c for c in list(PDX_DF) if "GMM" in c]
			####-------------
			CDF = CellLineDF.merge(PDX_DF,left_index=True,right_index=True,how="inner")
			Cols = [c for c in list(CDF) if "GMM" in c]
			CDF = CDF[Cols]
			#-------------
			ResultsList = []
			for Combo in combinations(list(CDF),2):
				print(Combo)
				Combo = list(Combo)
				if Combo[0] in PDX_Cols and Combo[1] in PDX_Cols:
					Comparison = "PDX"
				elif Combo[0] in CellLine_Cols and Combo[1] in CellLine_Cols:
					Comparison = "Cell Line"
				else:
					Comparison = "PDX-CellLine"
				#---------------
				DF = CDF[Combo]
				Denominator = DF[DF[Combo].sum(axis=1)>0].shape[0]
				Numerator = DF[DF[Combo].sum(axis=1)==1].shape[0]
				Value = Numerator/Denominator
				ResultsList.append([Combo[0], Combo[1], Comparison, Value])
			ResultsDF = pd.DataFrame(ResultsList,columns=["Sample1", "Sample2", 
									"Group Comparison", "Fraction of Sample-Specific Active Promoters"])
			ResultsDF.to_csv("For_group_comparison.txt",header=True,sep="\t")
			return ResultsDF
		#--------------------------------------------------------------------------
		def plotBoxplot(DF, Ax):
			# Bring in matplotlib
			import matplotlib
			matplotlib.use('Agg')
			import matplotlib.pyplot as plt
			from matplotlib import colors
			matplotlib.style.use('ggplot')
			import seaborn as sns
			from matplotlib.gridspec import GridSpec
			sns.set(font_scale=5)
			Ax.set_facecolor('white')
			DF.replace("PDX-CellLine", "PDX\nvs\nCellLine", inplace=True)
			MyPal = {"Cell Line": "Blue", "PDX": "Red", "PDX\nvs\nCellLine":"Purple"}
			sns.boxplot(x="Group Comparison", y="Fraction of Sample-Specific Active Promoters", 
						ax=Ax, data=DF, palette=MyPal, 
						order=["Cell Line", "PDX", "PDX\nvs\nCellLine"])

			for patch in Ax.artists:
				r, g, b, a = patch.get_facecolor()
				patch.set_facecolor((r, g, b, .3))

			Ax.set_xlabel("", fontsize=50)
			Ax.set_ylabel("Fraction of Sample-\nSpecific Active Promoters",fontsize=50)
			Ax.tick_params(axis="y",labelsize=45)
			Ax.tick_params(axis="x",labelsize=45)#labelrotation=45,

		#------
		ResultsDF = getData(CellLineDF, PDX_DF)
		plotBoxplot(ResultsDF, Ax)
	#--------------------------------------------------------------------------	
	def makeVenn(DF, Ax):
		#--------------------------------------------------------------------------
		def getData(DF):
			# format a temp dataframe
			SampleCols = [c for c in list(DF) if "GMM" in c]
			Samples = [c.split("_")[0] for c in SampleCols]
			TempDF = DF[SampleCols]
			TempDF.columns = Samples
			TempDF.reset_index(inplace=True)
			# iterate over samples
			SampleList = []
			ActivePromoterList = []
			print("processing VENN samples:")
			for Sample in list(Samples):
				ActivePromoters = set(TempDF[TempDF[Sample]==1]["EnsmblID_T"].values.tolist())
				print("\t", Sample)
				SampleList.append(Sample)
				ActivePromoterList.append(ActivePromoters)
			return SampleList, ActivePromoterList
		#--------------------------------------------------------------------------
		def plotVenn(SampleList, ActivePromoterList, Ax):
			# Bring in matplotlib
			import matplotlib
			matplotlib.use('Agg')
			import matplotlib.pyplot as plt
			from matplotlib import colors
			matplotlib.style.use('ggplot')
			import seaborn as sns
			from matplotlib.gridspec import GridSpec
			from matplotlib_venn import venn2, venn3, venn3_circles
			sns.set(font_scale=3.5)
			Ax.set_facecolor('white')
			#--------
			venn3(ActivePromoterList,set_labels=SampleList,ax=Ax)
		#----------------------------------------------------------------------
		SampleList, ActivePromoterList = getData(DF)
		plotVenn(SampleList, ActivePromoterList, Ax)

	#--------------------------------------------------------------------------	
	# Bring in matplotlib
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib import colors
	matplotlib.style.use('ggplot')
	import seaborn as sns
	from matplotlib.gridspec import GridSpec
	sns.set(font_scale=2)
	#
	ColorList = ["Blues","Reds", "Purples"]
	Fig = plt.figure(num=None, figsize=(45, 45), dpi=80, facecolor='w', edgecolor='k')
	#Fig, Ax = plt.subplots(ncols=3,figsize=(35, 35), dpi=80, facecolor='w', edgecolor='k')
	# create axes space
	ncols = 3
	nrows = 2
	GS = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[12, 12, 8], height_ratios=[3,1])
	HeatMapAxList = []
	Row2AxList = []
	for i in range(0,(ncols*nrows),1):
		NewAx = Fig.add_subplot(GS[i])		
		if i in [0,1,2]:
			HeatMapAxList.append(NewAx)
		else:
			Row2AxList.append(NewAx)
	# adjust subplots before laying down annotations or figures (diagrams)
	Fig.subplots_adjust(left=.1, right=.9, top=.9, bottom=.1, wspace=.3, hspace=.1)
	#plt.tight_layout()
	###########
	for i, Ax in enumerate(HeatMapAxList):		
		DF = DFList[i]
		AnnotList = AnnotList_List[i]
		# parse annotations placement
		HLineList, VLineList, = AnnotList[:2]
		ColText, ColX = AnnotList[2:4]
		RowText, RowY = AnnotList[4:6]
		CC_RowY = AnnotList[-1]
		print("Plotting Combined Heatmap GMM version 1")
		print(HLineList)
		print(VLineList)
		# plot everything
		sns.heatmap(DF,
					ax=Ax,
					cmap=ColorList[i],
					annot=False,
					vmin=0,
					vmax=3,
					yticklabels=False,
					xticklabels=False,
					cbar_kws={#"orientation": "horizontal", 
							"location":"top",
							"use_gridspec":False,
							"shrink": .75}, 
					cbar=True, 
					square=False, 
					robust=False)
		Ax.hlines(HLineList[:-1], *Ax.get_xlim(),linestyles="dashed", linewidth=4)
		Ax.vlines(VLineList[:-1], *Ax.get_ylim(),linestyles="solid", linewidth=5)	
		# annotate columns: samples
		for n, txt in enumerate(ColText):
			TrCO1 = Ax.transData.transform((ColX[n], Ax.get_ylim()[1]))
			Inv = Fig.transFigure.inverted()
			TrCO2 = Inv.transform(TrCO1)
			TrCO2 = (TrCO2[0], TrCO2[1]+.003)
			plt.annotate(txt, TrCO2,  xycoords="figure fraction",
						va="bottom", ha="center", fontweight="bold",
						fontsize=35)
		# annotate rows: groups
		for n, txt in enumerate(RowText):
			TrCO1 = Ax.transData.transform((Ax.get_xlim()[0], RowY[n]))
			Inv = Fig.transFigure.inverted()
			TrCO2 = Inv.transform(TrCO1)
			TrCO2 = (TrCO2[0]-.05, TrCO2[1])
			plt.annotate(txt, TrCO2, xycoords="figure fraction",#transform=plt.gcf().transFigure,
						va="center", ha="left", fontweight="bold",
						fontsize=30)
		# annotate rows: cancer concensus genes
		ColPosList = list(range(40,len(ColText)*40,40))
		for ColPos in ColPosList:
			for RowPos in CC_RowY:
				StartVs1 = Ax.transData.transform((ColPos-2, RowPos))
				EndVs1 = Ax.transData.transform((ColPos+2, RowPos))
				Inv = Fig.transFigure.inverted()
				StartVs2 = Inv.transform(StartVs1)
				EndVs2 = Inv.transform(EndVs1)
				#plt.annotate("", StartVs2, EndVs2, xycoords="figure fraction",
				#			arrowprops=dict(arrowstyle='<->',color="black"))  
	# plot boxplot
	makeVenn(RawDFList[0], Row2AxList[0])
	makeVenn(RawDFList[1], Row2AxList[1])
	makeBoxplot(RawDFList[0], RawDFList[1], Row2AxList[2])
	#
	Fig.savefig(PlotName) 
#------------------------------------------------------------------------------
def loadData(FilePath):
	IndexCols = ["EnsmblID_T","EnsmblID_G","Gene", "Tier"]
	DF = pd.read_csv(FilePath,header="infer",sep="\t")
	DF.set_index(IndexCols, inplace=True,drop=True)
	return DF

#------------------------------------------------------------------------------
def main():
	CL_DF = loadData("H3K27ac_CellLine_Data/H3K27ac_CellLine_Processed_FigureInput.txt")
	PDX_DF = loadData("H3K27ac_PDX_Data/H3K27ac_PDX_Processed_FigureInput.txt")
	SortCL_DF, CL_DFList1, AnnotCL_DFList = sortDF(CL_DF, "H3K27ac_CellLine_SortedData")
	SortPDX_DF, PDX_DFList2, AnnotPDX_DFList = sortDF(PDX_DF, "H3K27ac_PDX_SortedData")
	CDF, CDFList, Annot_CDF_List = combineDataTypes(CL_DF, PDX_DF)
	makeCombinedHeatmap([SortCL_DF, SortPDX_DF, CDF,], [AnnotCL_DFList, AnnotPDX_DFList, Annot_CDF_List], [CL_DF, PDX_DF], "Fig1_Combined.png")

#
if __name__ == '__main__':
	main()