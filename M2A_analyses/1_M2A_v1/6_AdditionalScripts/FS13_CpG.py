#!/usr/bin/python

'''

'''
import os
import time
import math
import glob
import argparse
import numpy as np
import pandas as pd

def loadData(WindowSize):
	Root = "CpG/Path/"
	Path = Root+"Sliding_W_TESTING_bp_S"+str(WindowSize)+"bp"
	FilePath = Path+"CpG_Sliding_W_TESTING_bp_S"+str(WindowSize)+"bp.txt"
	DF = pd.read_csv(FilePath,header="infer",sep="\t")
	Icols = ["EnsmblID_T","EnsmblID_G","Gene","Strand","Chr","Start","End","BinStart","BinEnd"]
	DF.set_index(Icols,inplace=True,drop=True)
	DF.replace(np.nan,0,inplace=True)
	return DF

def plotHist(LowerDFList,HigherDFList):
	# Bring in matplotlib
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib import colors
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	import seaborn as sns
	from matplotlib.gridspec import GridSpec
	plt.rcParams.update({'font.size': 15})
	sns.set(font_scale=1.5)
	sns.set_style("ticks")

	PlotName = "CpG_Distribution_V6.png"
	MaxSmallerWinCpG = max([df.values.max() for df in LowerDFList])
	MaxLargeWinCpG = max([df.values.max() for df in HigherDFList])
	SmallWinBins = list(range(0,11,1))+[MaxSmallerWinCpG]
	LargeWinBins = list(range(0,11,1))+[MaxLargeWinCpG]
	Fig = plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
	GS = GridSpec(ncols=3, nrows=2)#, width_ratios=1, height_ratios=1)
	# size 1
	for PlotIdx in [0,1,2]:
		if PlotIdx > 0:
			Ax = Fig.add_subplot(GS[0,PlotIdx],sharey=Ax)
		else:
			Ax = Fig.add_subplot(GS[0,PlotIdx])
		##
		Color = ["red","blue","black"][PlotIdx]
		Label = ["100bp", "250bp", "500bp"][PlotIdx]
		##
		Hist, Edges = np.histogram(LowerDFList[PlotIdx].values, SmallWinBins)
		Freq = Hist/float(Hist.sum())
		plt.bar(SmallWinBins[:-1], Freq, width=1, align="edge", ec="white",color=Color,label=Label)
		##
		plt.xlim(0, 11)
		Ax.legend()
		plt.title(Label+" windowed CpG distribution",fontsize=16)
		plt.xticks(np.arange(0,11,1),[str(i) for i in np.arange(0,10,1)]+[" 10+"])
		plt.xlabel("CpGs captured per window",fontsize=15)
		plt.ylabel("Fraction of windows",fontsize=15)
	# size 2
	for PlotIdx in [0,1,2]:
		if PlotIdx > 0:
			Ax = Fig.add_subplot(GS[1,PlotIdx],sharey=Ax)
		else:
			Ax = Fig.add_subplot(GS[1,PlotIdx])
		###
		Color = ["grey", "green", "purple"][PlotIdx]
		Label = ["1000bp", "2500bp", "5000bp"][PlotIdx]
		##
		Hist, Edges = np.histogram(HigherDFList[PlotIdx].values, LargeWinBins)
		Freq = Hist/float(Hist.sum())
		plt.bar(LargeWinBins[:-1], Freq, width=1, align="edge", ec="white",color=Color,label=Label)
		##
		plt.xlim(0, 11)
		Ax.legend()
		plt.title(Label+" windowed CpG distribution",fontsize=16)
		plt.xticks(np.arange(0,11,1),[str(i) for i in np.arange(0,10,1)]+[" 10+"])
		plt.xlabel("CpGs captured per window",fontsize=15)
		plt.ylabel("Fraction of windows",fontsize=15)

	print("Plotting:", PlotName)
	# adjust subplots before laying down annotations or figures (diagrams)
	Fig.tight_layout()
	Fig.savefig(PlotName)
	plt.clf()

def main():
	DF_100bp = loadData("100")
	DF_250bp = loadData("250")
	DF_500bp = loadData("500")
	LowerDFList = [DF_100bp, DF_250bp, DF_500bp]
	DF_1000bp = loadData("1000")
	DF_2500bp = loadData("2500")
	DF_5000bp = loadData("5000")
	HigherDFList = [DF_1000bp, DF_2500bp, DF_5000bp]
	plotHist(LowerDFList, HigherDFList)

if __name__ == '__main__':
	# Parse the arguments
	main()