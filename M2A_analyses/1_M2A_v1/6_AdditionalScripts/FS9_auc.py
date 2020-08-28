#!/usr/bin/python
'''
'''


import glob
import numpy as np
import pandas as pd
from math import ceil,floor,sqrt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import MinMaxScaler
import warnings
warnings.filterwarnings('ignore')
#------------------------------------------------------------------------------
def loadData(FilePath, Type):
	DF = pd.read_csv(FilePath,header="infer",sep="\t")
	DF = DF[DF["GeneType"]=="protein_coding"]
	IndexCols = ["Transcript", "Gene", "GeneType"]
	IndexCols = list(set(IndexCols)&set(list(DF)))
	DF.set_index(IndexCols,inplace=True,drop=True)
	DF.dropna(inplace=True)
	return DF

#------------------------------------------------------------------------------
def makeAUCROCFig(DF, PlotName, Type):
	# Bring in matplotlib
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib import colors
	matplotlib.style.use('ggplot')
	import seaborn as sns
	from matplotlib.gridspec import GridSpec
	sns.set(font_scale=2)
	from string import ascii_lowercase
	#
	#
	#
	# unique list of samples
	Samples = np.unique([c.split("_")[0] for c in list(DF)])
	# determine dimensions of subplots
	SR = sqrt(len(Samples))
	if ((SR - floor(SR)) == 0):
		ncols = nrows = int(SR)
	else:
		nrows = int(ceil(SR))
		ncols = int(ceil(len(Samples)/nrows))
	# initialize figure and subplots
	Fig = plt.figure(num=None, figsize=(125, 125), dpi=80, facecolor='w', edgecolor='k')
	GS = GridSpec(ncols=ncols, nrows=nrows)#, width_ratios=1, height_ratios=1)
	# iterate over subplots, samples in DF
	ResultsList = []
	for i, (Sample, Letter) in enumerate(zip(Samples,ascii_lowercase)):
		# add subplot
		Ax = Fig.add_subplot(GS[i])		
		# get relevant columns from DF
		FPKMCol = [c for c in list(DF) if "FPKM" in c and Sample in c][0]
		PredCol = [c for c in list(DF) if "Pred" in c and Sample in c and Type in c][0]
		ActCol = [c for c in list(DF) if "Actual" in c and Sample in c and Type in c][0]
		#####		
		print("#####\nFPKM: "+FPKMCol)
		print("PredCol: " + PredCol)
		print("ActCol: " + ActCol)
		# iterate over predicted, actual values for an AUC ROC comparison
		for ColType, ColName in zip(["Predicted", "Actual"], [PredCol, ActCol]):
			# reset index to get access to gene
			TempDF = DF.reset_index(inplace=False,drop=False)
			# sort on genes and colname values
			TempDF.sort_values(["Gene",ColName],ascending=[True,False],inplace=True)
			# keep first duplicate (largest value)
			Duplicates = TempDF.duplicated("Gene", keep="first")
			# remove duplicates
			TempDF = TempDF[~Duplicates]
			# caluculate ROC, AUC
			AUC = roc_auc_score(TempDF[FPKMCol].astype(int).values, TempDF[ColName].values)
			FPR,TPR,_ = roc_curve(TempDF[FPKMCol].astype(int).values, TempDF[ColName].values)
			Label = ColType+"\nAUC: "+str(round(AUC,3))
			plt.plot(FPR, TPR, label=Label, linewidth=25)
			ResultsList.append([Sample, Type, ColType, AUC])
		plt.grid()
		plt.xlim(0.0, 1.0)
		plt.ylim(0.0, 1.0)
		plt.plot([0, 1], [0, 1], 'k--',linewidth=15)
		plt.xlabel('False positive rate',fontsize=200)
		plt.ylabel('True positive rate',fontsize=200)
		plt.title("("+Letter+")  "+Sample,fontsize=250, loc="left")
		plt.legend(loc="lower right", prop={"size":175})#'best')		
	# save table to file
	ResultsDF = pd.DataFrame(ResultsList,columns=["Sample", "Type", "Pred/Actual", "AUC"])
	ResultsDF.to_csv("Table_"+PlotName+".txt",header=True,index=False,sep="\t")
	# adjust subplots before laying down annotations or figures (diagrams)
	Fig.subplots_adjust(left=.03, right=.975, top=.975, bottom=.025, wspace=.25, hspace=.25)
	#Fig.tight_layout()
	Fig.savefig(PlotName+".png")
	plt.clf()

#------------------------------------------------------------------------------
def main():
	Type = "H3K27ac"
	ENCODE = loadData("Final_Curated_ENCODE_Datasets.txt", Type)
	makeAUCROCFig(ENCODE, Type+"_ENCODE_AUC_ROC_Combined", Type)

if __name__ == '__main__':
	main()