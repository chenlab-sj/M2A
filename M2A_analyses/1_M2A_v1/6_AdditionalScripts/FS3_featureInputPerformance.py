#!/usr/bin/python


'''
1) Compare Train feature performance (all features) stratified by resolution and feature type (input vs CNN extracted), across entire train set

2) Compare best features determined by step 1) performance in validation set by sample.
'''

# Set random seed
import numpy as np
np.random.seed(9)
# Standard modules
import sys
import math
import glob
import argparse
import pandas as pd
from os.path import isfile, join
# statistics
from scipy.stats import pearsonr
import scipy.stats
# Feature File Processing
import h5py
from keras.models import Sequential
from keras.models import load_model
from keras.models import Model


#------------------------------------------------------------------------------
def getInputFiles():
	RootPath = "Dir/path"
	InputDirArray = np.array(["4_FeatData_NBL_121_File_0.h5",
								"4_FeatData_NBL_97_File_1.h5",
								"4_FeatData_NBL_109_File_2.h5",
								"4_FeatData_NBL_127_File_3.h5",
								"4_FeatData_NBL_3_File_4.h5",
								"4_FeatData_NBL_78_File_5.h5",
								"4_FeatData_NBL_CHLA90_File_6.h5",
								"4_FeatData_NBL_IMR32_File_7.h5",
								"4_FeatData_NBL_LAN6_File_8.h5",
								"4_FeatData_NBL_NB5_File_9.h5",
								"4_FeatData_NBL_NBLS_File_10.h5",
								"4_FeatData_NBL_SKNBE2_File_11.h5",
								"4_FeatData_NBL_SKNFI_File_12.h5",
								"4_FeatData_NBL_SKNMM_File_13.h5",
								"4_FeatData_NBL_126_File_14.h5",
								"4_FeatData_NBL_268_File_15.h5"])

	PathList = [join(RootPath, SampleFile) for SampleFile in InputDirArray]
	return PathList

#------------------------------------------------------------------------------
def getSampleGroup(Sample):
	SampleGroupDict = {"NBL_121":"Train",
						"NBL_97":"Train",
						"NBL_109":"Train",
						"NBL_127":"Train",
						"NBL_3":"Train",
						"NBL_78":"Train",
						"NBL_CHLA90":"Test",
						"NBL_IMR32":"Test",
						"NBL_LAN6":"Test",
						"NBL_NB5":"Test",
						"NBL_NBLS":"Test",
						"NBL_SKNBE2":"Test",
						"NBL_SKNFI":"Test",
						"NBL_SKNMM":"Test",
						"NBL_126":"Test",
						"NBL_268":"Test"}
	return SampleGroupDict[Sample]

#------------------------------------------------------------------------------
def loadData(InFileList):
	SampleDict = dict()
	for InFile in InFileList:
		Sample = "_".join(InFile.split("/")[-1].split("_")[2:4]).replace(".h5","")
		# Read h5 images and resp vars
		SampleDataHDF = h5py.File(InFile, mode="r")
		IndexKeys = list(SampleDataHDF.keys())
		# Transcript and GeneNames
		GeneNameIdx = [k for k in IndexKeys if "Gene" in k][0]
		GeneNameArray = np.array(SampleDataHDF[GeneNameIdx])
		##
		TranscriptIdx = [k for k in IndexKeys if "EnsmblID_T" in k][0]
		TranscriptArray = np.array(SampleDataHDF[TranscriptIdx])
		# images will always be last in list
		FeatDataIdx = [k for k in IndexKeys if "images" in k][0]
		FeatDataArray = np.array(SampleDataHDF[FeatDataIdx])
		# h2k27ac response var idx
		K27AcIdx = [k for k in IndexKeys if "Ac_div_Input" in k][0]
		H3K27acArray = np.array(SampleDataHDF[K27AcIdx])
		# h23k4me3 response var idx
		K3Me3Idx = [k for k in IndexKeys if "K4Me3_div_Input" in k][0]
		H3K4me3Array = np.array(SampleDataHDF[K3Me3Idx])
		###
		#build dataframe for meta and RV data
		RespVarDF_List = [H3K27acArray, H3K4me3Array, 
						GeneNameArray, TranscriptArray]
		Columns = ["H3K27ac", "H3K4me3", 
					"GeneName", "Transcript"]
		RespVarDF = pd.DataFrame(RespVarDF_List).T
		RespVarDF.columns = Columns
		# replace b'' from text
		for Col in ["GeneName", "Transcript"]:
			RespVarDF[Col] = RespVarDF[Col].astype(str).str.replace("b'","",regex=True)
			RespVarDF[Col] = RespVarDF[Col].astype(str).str.replace("'","",regex=True)
		# append to sample dict
		SampleDict[Sample] = [FeatDataArray, RespVarDF]
		print(Sample, RespVarDF.shape, FeatDataArray.shape)
	return SampleDict 

#------------------------------------------------------------------------------
def getCNNFeats(SampleDict, M2A_Model):
	# parse base M2A model to get access to CNN extracted features:
	FeatExtractModel = Model(input=M2A_Model.input, output=M2A_Model.get_layer('leaky_re_lu_3').output)
	CNNSampleDict = dict()
	SampleList = list(SampleDict.keys())
	for Sample in SampleList:
		Features = SampleDict[Sample][0]
		Features = Features[:,:,:,[0,1,3]]
		CNNFeatures = FeatExtractModel.predict(Features)
		CNNSampleDict[Sample] = [CNNFeatures,SampleDict[Sample][1]]
	return CNNSampleDict

#------------------------------------------------------------------------------
def getTotalPerformance(SampleDict, Type, TestTrain):
	if Type == "Input":
		FIdx = [0,1,3]
		FList = ["Mean", "FracSSD", "Var"]
	elif Type == "CNN":
		FIdx = list(range(0,8,1))
		FList = ["filter_"+str(i) for i in range(1,9,1)]
	# determine type of analysis
	if TestTrain == "Train":
		ResultsList = []
		SampleList = list(SampleDict.keys())
		TrainSamples = [s for s in SampleList if getSampleGroup(s)=="Train"]
		# iterate over each feature across all train samples
		for ResIdx, Res in zip([0,1],["250bp","2500bp"]):
			for WindowIdx in range(0,20,1):
				for FeatIdx, Feat in zip(FIdx,FList):
					TrainSampleFeatureList = []
					TrainSampleRespVarList = []
					for Sample in TrainSamples:
						# get features:
						FeaturesArray = SampleDict[Sample][0] 
						Features = FeaturesArray[:,ResIdx,WindowIdx,FeatIdx].flatten()
						TrainSampleFeatureList.append(Features)
						# get response:
						RespVarDF = SampleDict[Sample][1]
						TrainSampleRespVarList.append(RespVarDF.H3K27ac.values.flatten())
					#-------------
					H3K27acAllSample = np.concatenate(TrainSampleRespVarList,axis=0).flatten()
					FeatureAllSample = np.concatenate(TrainSampleFeatureList,axis=0).flatten()
					H3K27ac_R2 = pearsonr(H3K27acAllSample, FeatureAllSample)[0]**2
					Window = WindowIdx+1
					ResultsList.append([Res, Window, Feat, H3K27ac_R2])
		Columns=["window size","window","feature","performance"]
		PerformanceDF = pd.DataFrame(ResultsList,columns=Columns)
		print(Type, "TRAIN TOTAL PERFORMERS:")
		print(PerformanceDF)
	# determine type of analysis
	elif TestTrain == "Test":
		ResultsList = []
		SampleList = list(SampleDict.keys())
		TestSampleList = [s for s in SampleList if getSampleGroup(s)=="Test"]
		# iterate over every test sample for each feature
		for Sample in TestSampleList:
			for ResIdx, Res in zip([0,1],["250bp","2500bp"]):
				for WindowIdx in range(0,20,1):
					for FeatIdx, Feat in zip(FIdx,FList):
						# get features:
						FeaturesArray = SampleDict[Sample][0] 
						Features = FeaturesArray[:,ResIdx,WindowIdx,FeatIdx].flatten()
						# get response:
						RespVarDF = SampleDict[Sample][1]
						H3K27ac = RespVarDF.H3K27ac.values.flatten()
						H3K27ac_R2 = pearsonr(H3K27ac, Features.flatten())[0]**2
						Window = WindowIdx+1
						ResultsList.append([Sample, Res, Window, Feat, H3K27ac_R2])
						#print(Window, H3K27ac.shape, Features.shape, H3K27ac_R2)
		Columns=["sample","window size","window","feature","performance"]
		PerformanceDF = pd.DataFrame(ResultsList,columns=Columns)
		print(Type, "TEST TOTAL PERFORMERS:")
		print(PerformanceDF)
	return PerformanceDF

#------------------------------------------------------------------------------
def getBestPerformers(DFList, FeatType, TrainTest):
	# train
	if TrainTest=="Train":
		TotalPerformanceDF = DFList[0]
		BestPerformerList = []
		for Res in ["250bp","2500bp"]:
			ResDF = TotalPerformanceDF[TotalPerformanceDF["window size"]==Res]
			MaxPerf = ResDF["performance"].max()#.values[0]
			BestPerformerList.append(ResDF[ResDF["performance"]==MaxPerf])
		BestPerformanceDF = pd.concat(BestPerformerList,axis=0)
		BestPerformanceDF["feature type"] = FeatType
		print(FeatType, "TRAIN BEST PERFORMERS:")
		print(BestPerformanceDF)
	#-----------
	# test		
	elif TrainTest=="Test":
		TotalPerformanceDF = DFList[0]
		BestPerfDF = DFList[1]
		BestPerformerList = []
		for Res in ["250bp","2500bp"]:
			#-----
			# Test performance
			ResDF = TotalPerformanceDF[TotalPerformanceDF["window size"]==Res]
			ResDF["feature performance"] = "other features"
			#-----
			# get best perf from train data
			Res_BestPerfDF = BestPerfDF[BestPerfDF["window size"]==Res]
			BestWindow = Res_BestPerfDF["window"].values[0]
			BestFeat = Res_BestPerfDF["feature"].values[0]
			#-----
			ResDF.loc[(ResDF.window==BestWindow)&(ResDF.feature==BestFeat), 'feature performance'] = "best feature"
			BestPerformerList.append(ResDF)
		BestPerformanceDF = pd.concat(BestPerformerList,axis=0)
		BestPerformanceDF["feature type"] = FeatType
		print(FeatType, "TEST BEST PERFORMERS:")
		print(BestPerformanceDF)
	#------------
	return BestPerformanceDF

#------------------------------------------------------------------------------
def plotPerformance(TestDF, TrainDF):
	# Bring in matplotlib
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib import colors
	matplotlib.style.use('default')
	import seaborn as sns
	from matplotlib.gridspec import GridSpec
	plt.rcParams.update({'font.size': 22})
	sns.set(font_scale=1.25,style="ticks")

	PerfDF = pd.concat([TrainDF,TestDF],axis=0)
	PlotName = "NBL_Input_Feature_Performance_Dist.png"
	Boxplot = sns.catplot(x="window size", y="performance", hue="feature type",
			data=PerfDF, col="sample set", kind="box", palette="colorblind", legend_out=True,
			height=5, aspect=1)
	Boxplot.set_ylabels('Feature performance ($\mathregular{R^2}$)')
	Boxplot.fig.suptitle("Feature performance comparison, Input features vs CNN mapped features\n\n")
	Boxplot.fig.subplots_adjust(top=0.85)
	plt.savefig(PlotName)
	plt.clf()


#------------------------------------------------------------------------------
def main():
	# load model
	ModelPath = "Model.h5"
	M2A_Model = load_model(ModelPath)
	# load data
	InFileList = getInputFiles()
	SampleDict = loadData(InFileList)
	# process cnn features
	SampleCNNDict = getCNNFeats(SampleDict, M2A_Model)
	# determine individual feature predictive power
	Train_Input_TotalPerformanceDF = getTotalPerformance(SampleDict,"Input", "Train")
	Train_CNN_TotalPerformanceDF = getTotalPerformance(SampleCNNDict,"CNN", "Train")
	Test_Input_TotalPerformanceDF = getTotalPerformance(SampleDict,"Input", "Test")
	Test_CNN_TotalPerformanceDF = getTotalPerformance(SampleCNNDict,"CNN", "Test")
	#### to get the best performaers from train in Test
	# 1) Determine the best performers in Train
	Train_Input_BestPerformanceDF = getBestPerformers([Train_Input_TotalPerformanceDF], "input features", "Train")
	Train_CNN_BestPerformanceDF = getBestPerformers([Train_CNN_TotalPerformanceDF], "CNN features", "Train")
	# 2) get the Train performaers in the test set
	Input_BestPerformanceDF = getBestPerformers([Test_Input_TotalPerformanceDF,Train_Input_BestPerformanceDF], "input features", "Test")
	CNN_BestPerformanceDF = getBestPerformers([Test_CNN_TotalPerformanceDF, Train_CNN_BestPerformanceDF], "CNN features", "Test")
	# 3) combine Input and CNN best Test performers
	Test_BestPerformanceDF = pd.concat([Input_BestPerformanceDF,CNN_BestPerformanceDF],axis=0)
	Test_BestPerformanceDF['sample set'] = "Test"
	Test_BestPerformanceDF = Test_BestPerformanceDF[Test_BestPerformanceDF['feature performance']=="best feature"]
	#### to get the total train features:
	Train_CNN_TotalPerformanceDF["feature type"] = "CNN features" 
	Train_Input_TotalPerformanceDF["feature type"] = "input features"
	Train_BestPerformanceDF = pd.concat([Train_CNN_TotalPerformanceDF,Train_Input_TotalPerformanceDF],axis=0)
	Train_BestPerformanceDF['feature performance'] = "all features"
	Train_BestPerformanceDF['sample set'] = "Train"
	# plot performance by resolution
	plotPerformance(Test_BestPerformanceDF, Train_BestPerformanceDF)

# Parse line arguments
if __name__ == '__main__':
	main()