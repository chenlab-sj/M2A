#!/usr/bin/python

'''
Create and plot RMS tSNE analysis for observed, predicted(vanilla) and predicted(transfer)

'''

import numpy as np
np.random.seed(123456)
import pandas as pd
from sklearn.manifold import TSNE as sklearn_TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
# supressing copy in-place warnings
pd.options.mode.chained_assignment = None

#------------------------------------------------------------------------------
def loadData(FilePath):
	DF = pd.read_csv(FilePath,header="infer",sep="\t")
	DF.set_index("Transcript",inplace=True,drop=True)
	return DF

#------------------------------------------------------------------------------
def reduceDimensions(DF,Perplexity):
	PredDF = DF.copy()
	MMS = StandardScaler()
	Scaled = MMS.fit_transform(PredDF.values)
	X_embedded = sklearn_TSNE(n_components=2,method="barnes_hut",perplexity=Perplexity).fit_transform(Scaled.T)
	TSNE = pd.DataFrame(X_embedded, columns=["tSNE-1","tSNE-2"])
	TSNE["SampleName"] = PredDF.T.index
	return TSNE

#------------------------------------------------------------------------------
def formatTSNE(DF_tSNE):
	DF_tSNE[['Sample','Experiment',"Sample Type"]] = DF_tSNE.SampleName.str.split("_",expand=True,)
	DF_tSNE.drop(columns=["SampleName"],inplace=True)
	return DF_tSNE

#------------------------------------------------------------------------------
def plotTSNE(DF, Perplexity):
	# Bring in matplotlib
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib import colors
	matplotlib.style.use('default')
	from matplotlib.gridspec import GridSpec
	plt.rcParams.update({'font.size': 22})
	#---
	PlotName = "P"+str(Perplexity)+"_PLOT_TSNE.png"
	Fig = plt.figure(num=None, figsize=(24, 8), dpi=80, facecolor='w', edgecolor='k')
	GS = GridSpec(ncols=3, nrows=1)#, width_ratios=1, height_ratios=1)
	#--
	Marks = ["D","o","s"]
	Experiments = ["Actual","Baseline","Transfer"]
	for Idx in [0,1,2]:
		Mark = Marks[Idx]
		Experiment = Experiments[Idx]
		Ax = Fig.add_subplot(GS[Idx])	
		for Subtype, Color in zip(["ARMS","ERMS","Spindle"],["orange","blue","gray"]):
			TempDF = DF[(DF["Sample Type"]==Subtype)&(DF["Experiment"]==Experiment)]
			plt.scatter(TempDF["tSNE-1"], TempDF["tSNE-2"], label=Subtype, s=175, marker=Mark, facecolors=Color, edgecolors='black',linewidth=1)
			Ax.set_ylabel("tSNE-2")
			Ax.set_xlabel("tSNE-1")
			if Experiment == "Actual":
				Ax.legend()
				Ax.set_title("Observed promoter activity landscape\nof RMS tumors",fontweight="bold",fontsize=20)
			elif Experiment =="Baseline":
				Ax.set_title("Promoter activity landscape of RMS tumors\ninferred from the baseline model",fontweight="bold",fontsize=20)
			elif Experiment =="Transfer":
				Ax.set_title("Promoter activity landscape of RMS tumors\ninferred from the transfer model",fontweight="bold",fontsize=20)
	#----
	print("Plotting:", PlotName)
	plt.tight_layout()
	Fig.subplots_adjust(top=0.9)
	Fig.savefig(PlotName)
#------------------------------------------------------------------------------
def main():
	Perplexity = 6
	DF = loadData("RMS_tSNE_dataset.txt")
	DF_tSNE = reduceDimensions(DF,Perplexity)
	DF_tSNE = formatTSNE(DF_tSNE)
	plotTSNE(DF_tSNE, Perplexity)
#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()