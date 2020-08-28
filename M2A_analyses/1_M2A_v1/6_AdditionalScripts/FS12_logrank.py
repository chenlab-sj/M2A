#!/usr/bin/python

'''
Purpose:
KaplanMeier Curves, logrank for mutational status

'''
import os
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import numpy as np
import pandas as pd


#-------------------------------------------------------------------------
def plotKaplanMeier(DF, GroupName):
	# Bring in matplotlib
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib import colors
	matplotlib.style.use('seaborn-talk')
	import seaborn as sns
	# fit data
	sns.set(font_scale=.95)
	# plot baseline survival
	KMF = KaplanMeierFitter()
	# get log rank results
	Unique = sorted(DF[GroupName].unique())[::-1]
	G0 = Unique[0]
	G1 = Unique[1]
	Group0 = DF[DF[GroupName]==G0]
	Group1 = DF[DF[GroupName]==G1]
	LR_Results = logrank_test(Group0["OS (months)"].astype(float), 
							Group1["OS (months)"].astype(float), 
							Group0["Patient Status"].astype(int), 
							Group1["Patient Status"].astype(int))
	# build figure
	fig, ax = plt.subplots(1,1,figsize=(7, 5))
	# add logrank to plot
	# iterate over unique groups and plot KM curves
	Median_Survival_List = []
	for Group in Unique:
		TempDF = DF[DF[GroupName]==Group]
		KMF = KaplanMeierFitter()
		KMF.fit(TempDF["OS (months)"].astype(float), 
				TempDF["Patient Status"].astype(int),
				label=Group).plot(ax=ax)
		print("KMF.median_survival_time_", KMF.median_survival_time_)
		Median_Survival_List.append(float(KMF.median_survival_time_))
	Median_Survival_List = sorted(Median_Survival_List)
	Annotation = "logrank: "+str(round(LR_Results.test_statistic,3))+"\n"+\
					"pvalue: "+str(round(LR_Results.p_value,5))+"\n"+\
					"median survival times: "+str(Median_Survival_List)
	ax.text(.05, .05, Annotation, 
			horizontalalignment='left',
			verticalalignment='bottom',
			transform=ax.transAxes)			
	# set labels
	ax.set_xlabel("Time (months)")
	ax.set_ylabel("Overall survival probability")
	plt.tight_layout()
	fig.savefig("KaplanMeierSurvival.png")
	plt.clf()
	plt.close()
	return(LR_Results)

#------------------------------------------------------------------------------
def main():
	InputDF = pd.read_csv("InputData.txt",header="infer",sep="\t")
	plotKaplanMeier(DF, "MutationStatus")

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse the arguments
	main()