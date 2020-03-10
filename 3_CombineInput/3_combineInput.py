#!/usr/bin/env python
'''
Requires:
python/3.6.5+
sklearn version '0.20.2+'
time
os
h5py
argparse
numpy
pandas


Purpose:
Combine/interleave all features into pseudo image arrays for input to a CNN model

General Notes:
In the event transfer learning is required, this .h5 file output will also include
the response variable.

#------------------------------------------------------------------------------
Input Required:
1) Tab delimited list of methylation features for 250, 2500bp window sizes
	1A) Features (from Windowed M-value): 
		Average, Variation, SSD, Fraction of Region SSD
	1B) Includes positional and gene data:
		"EnsmblID_T", TranscriptID
		"EnsmblID_G", GeneID
		"Gene", GeneName
		"Strand", 
		"Chr", 
		"Start", 
		"End", 
		"RStart", response variable region start
		"REnd", response variable region end

#------------------------------------------------------------------------------
Output:
1)This script produces one HDF5 file for CNN predictions/ training.

The main output dataset is termed "FeatureInput" in the HDF5, and is sorted 
such that:

An m by n array, 
	where m = Windowsizes (resolutions),
	where n = Window position relative to the TSS (W-10,W-9...W-1,{TSS}W1,W2...W10)

"pixel" == Image(i,j) == WindowedFeatures(Windowsize,WindowPosition)
"channels" == features (e.g., Ave,Var,SSD,FracSSD)

[[[Features250bp at W-10],[Features250bp at W-9]...[Features250bp at W9][Features250bp at W10]],
[[Features2500bp at W-10][Features2500bp at W-9]...[Features2500bp at W9][Features2500bp at W10]]]

To do this, both columns and rows are interleaved, then the matrix is reshaped to:
(N,Resolutions,Windows,Features)
N = 96,757 number of promoters from all chromosomes 
	(or however many in the file provided)
Resolutions = 2, or number of windowsizes (250, 2500)
Windows = 20, windows
Features = 4, (Ave, Var, SSD, FracSSD)

Output Notes:
	All meta data/ positional values are also recorded in hdf5, astype(bytes).
'''

import time
import os
import h5py
import argparse
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

#------------------------------------------------------------------------------
def parseArguments(): 
	# Create argument parser
	parser = argparse.ArgumentParser()   
	# Positional mandatory arguments
	parser.add_argument("MethylationFilePath", 
		help="Full path to methylation features file.", type=str)
	parser.add_argument("--ResponseVariablePath", help="Full path to response variable file.", type=str)
	parser.add_argument("--outFileName", help="File name to use for output", type=str)
	parser.add_argument("--outDirectory", help="Directory to write output to", default="output", type=str)

	# Parse arguments
	args = parser.parse_args()
	return args

#------------------------------------------------------------------------------
# 1)
def prepMethylationData(MethylationFilePath, WindowSizeList):
	#--------------------------------------------------------------------------
	# 1.1)
	def loadData(MethylationFilePath):
		ExtraCols = ["EnsmblID_T", "EnsmblID_G", "Gene",
					"Strand", "Chr", "Start", "End", "RStart", "REnd"]
		MethylationDF = pd.read_csv(MethylationFilePath,header="infer",sep="\t")
		MethylationDF.set_index("EnsmblID_T", inplace=True,drop=False)
		MetaDF = MethylationDF[ExtraCols]
		MethylationDF.drop(ExtraCols,inplace=True,axis=1)
		return MethylationDF, MetaDF

	#--------------------------------------------------------------------------
	# 1.2)
	def scaleData(MethylationFeatDF):
		# hard coded parameters, most likely will not change:
		NewMval = 0
		OldMval = np.nan
		MinMval = .1
		MaxMval = 1
		#
		# Note:
		# The full range of a given feature for a given window size is required
		# for an accurate scaler, so the feature values are reshaped/flattened.
		# This requires iteration over unique window size/ feature type.
		#
		# Also: not all versions of MinMaxScaler (sklearn '0.20.2') 
		# ignores np.nans. This functionality is essential for this chunk of 
		# code,(np.nan masking)
		#
		#
		WindowSizeList = np.unique([c.split("_")[0] \
											for c in list(MethylationFeatDF)])
		FeatureTypeList = np.unique(["_".join(c.split("_")[-2:]) \
											for c in list(MethylationFeatDF)])
		ScaledFeatDFList = []
		for WindowSize in WindowSizeList:
			for FeatureType in FeatureTypeList:
				FeatureCols = [c for c in list(MethylationFeatDF) if \
								c.split("_")[0]==str(WindowSize) and \
								"_".join(c.split("_")[-2:])==FeatureType]
				TempMethylDF = MethylationFeatDF[FeatureCols]
				# Scaling
				MMS = MinMaxScaler(feature_range=(MinMval,MaxMval))
				NewShape = (np.prod(TempMethylDF.shape),1)
				ReshapedMethylArray = TempMethylDF.values.reshape(NewShape)
				MMS.fit(ReshapedMethylArray)
				ScaledFeatDF = pd.DataFrame(MMS.transform(TempMethylDF),
												columns=TempMethylDF.columns)
				ScaledFeatDF.index = MethylationFeatDF.index
				ScaledFeatDFList.append(ScaledFeatDF)
		#--------------------------------------------
		ScaledMethylDF = pd.concat(ScaledFeatDFList,axis=1)
		ScaledMethylDF.replace(OldMval,NewMval,inplace=True)
		return ScaledMethylDF

	#-------------------------------------------------------------------------#
	# load data and drop extraneous fields, set index
	MethylationFeatDF, MetaDF = loadData(MethylationFilePath)
	# scale DFs from 0-1, replace np.nan with 0
	ScaledMethylDF = scaleData(MethylationFeatDF)
	return ScaledMethylDF, MetaDF
	
#------------------------------------------------------------------------------
# 2)
def runCombineDataPipeline(ScaledMethylDF, NumOfWin, WindowSizeList):
	#-------------------------------------------------------------------------#
	# 2.1)
	def splitOnWindowSize(ScaledMethylDF, WindowSizeList):
		WindowSizeMethylDFList = []
		for WindowSize in WindowSizeList:
			FeatureCols = [c for c in list(ScaledMethylDF) \
										if c.split("_")[0]==str(WindowSize)]
			WindowSizeMethylDF = ScaledMethylDF[FeatureCols]
			WindowSizeMethylDF.columns = ["_".join(c.split("_")[1:]) \
										for c in list(WindowSizeMethylDF)]
			WindowSizeMethylDFList.append(WindowSizeMethylDF)
		return WindowSizeMethylDFList

	#-------------------------------------------------------------------------#
	# 2.2)
	def interleaveCols(WindowSizeMethylDFList, NumOfWin, WindowSizeList):
		#---------------------------------------------------------------------#
		# 2.2.1)
		def getFeatOrder(Window250bpDF, NumOfWin):
			DF_Feat_List = list(Window250bpDF)
			# sort on window number:
			SortedFeatList = sorted(DF_Feat_List, 
							key=lambda x: int(x.split("_")[0].replace("W","")))
			TotalNumFeatures = len(DF_Feat_List)
			FeatPerWindow = int(TotalNumFeatures / NumOfWin)
			WindowFeatList = [SortedFeatList[i:i+FeatPerWindow] \
							for i in range(0,TotalNumFeatures,FeatPerWindow)]			
			FeatOrder = [sorted(i,key=lambda x:"_".join(x.split("_")[1:]))\
								for i in WindowFeatList]
			FeatOrder = np.array(FeatOrder).flatten().tolist()
			return FeatOrder

		#---------------------------------------------------------------------#
		# FeatureCols are the same between WindowSizeDFs, feature order should
		# be the same for both:
		Window250bpDF = WindowSizeMethylDFList[0]
		FeatureOrder = getFeatOrder(Window250bpDF, NumOfWin)
		# iterate over each WindowSize DF and reorder features
		InterleaveDFList = []
		for WinIdx, WindowSizeMethylDF in enumerate(WindowSizeMethylDFList):
			# reorder these columns
			InterleaveDF = WindowSizeMethylDF[FeatureOrder]
			# Adjust index for row interleaving/sorting later
			# NOTE: only required for multiple resolutions
			InterleaveDF.index = InterleaveDF.index+"_"+str(WinIdx)
			print("Col Interleave "+str(WindowSizeList[WinIdx])+\
					" Resolution shape:", InterleaveDF.shape)
			InterleaveDFList.append(InterleaveDF)
		return InterleaveDFList

	#-------------------------------------------------------------------------#
	# 2.3)
	def interleaveRows(InterleavedColDFList):
		InterleavedRowDF = pd.concat(InterleavedColDFList, axis=0).sort_index()
		print("Example DF Row Interleaved,\n", InterleavedRowDF.iloc[:9])
		return InterleavedRowDF

	#-------------------------------------------------------------------------#
	# 2.4)
	def reshapeDF(InterleavedDF, NumOfWin, WindowSizeList):
		# USING: (n,2,20,4), (N,R,W,1,Wcols)
		R = len(WindowSizeList) # num of window sizes, resolutions
		N = int(InterleavedDF.shape[0] / R) # number of "images", promoters
		W = NumOfWin # number of windows
		Wcols = int(InterleavedDF.shape[1] / W) # pixel channels
		Reshape = (N,R,W,Wcols)
		print("Using Reshape Values:", str(Reshape))
		# Reshape and return...
		ReshapedImageArray = InterleavedDF.values.reshape(Reshape)
		print("ReshapedImageArray", ReshapedImageArray.shape)
		return ReshapedImageArray

	#-------------------------------------------------------------------------#
	# split methylationDF into separate window sizes/ resolutions for 
	# interleaving
	WindowSizeMethylDFList = splitOnWindowSize(ScaledMethylDF, WindowSizeList)
	# interleave features, or columns		
	InterleavedColDFList = interleaveCols(WindowSizeMethylDFList, 
											NumOfWin, 
											WindowSizeList)
	InterleavedDF = interleaveRows(InterleavedColDFList)
	ReshapedImageArray = reshapeDF(InterleavedDF, NumOfWin, WindowSizeList)
	return ReshapedImageArray, InterleavedDF.index

#------------------------------------------------------------------------------
# 4)
def sortMetaDF(SortedIndex, MetaDF):
	#--------------------------------------------------------------------------
	# 4.1)
	def processIdx(SortedIndex):
		# retrieving unique index, removing resolution level data
		IDX = pd.DataFrame(SortedIndex.str.split("_").tolist(), 
														columns = ["1","2"])
		return IDX["1"].unique()
	#--------------------------------------------------------------------------
	FinalIndex = processIdx(SortedIndex)
	NewMetaDF = MetaDF.loc[FinalIndex]
	return NewMetaDF

#------------------------------------------------------------------------------
# 5)
def createHDF5(ReshapedImageArray, MetaDF, OutputFilePath):
	# initial hdf5 file 
	HDF5_File = h5py.File(OutputFilePath, mode='w')
	# create "image" feature dataset
	HDF5_File.create_dataset("FeatureInput", data=ReshapedImageArray)
	# create response variable dataset
	if "log2_ChipDivInput" in list(MetaDF):
		HDF5_File.create_dataset("log2_ChipDivInput", data=MetaDF["log2_ChipDivInput"].values)
		MetaDF.drop(columns=["log2_ChipDivInput"],inplace=True)
	# load in meta data	
	for Meta in list(MetaDF):
		HDF5_File.create_dataset(Meta,data=MetaDF[Meta].values.astype(bytes))
	HDF5_File.close()
	print("Created HDF5 Dataset:", OutputFilePath)

#------------------------------------------------------------------------------
def main(LineArgs):
	T0 = time.time()
	# Variable input
	MethylationFilePath = LineArgs.MethylationFilePath
	Sample = ".".join(MethylationFilePath.split("/")[-1].split(".")[:-1])
	OutputFilePath = LineArgs.outDirectory
	if not os.path.exists(LineArgs.outDirectory):
		os.mkdir(LineArgs.outDirectory)

	if not LineArgs.outFileName:
		OutputFilePath += "/" + Sample+".h5"
	else:
		OutputFilePath += "/" + LineArgs.outFileName
	# hard coded parameters
	NumOfWin = 20
	WindowSizeList = [250, 2500]
	# 1) prep methylation for interleaving feature arrays:
	# 1.1) load data, remove extraneous columns
	# 1.1) slice extraneous columns into MetaData DF:
	#		positional, gene/transcript IDs, etc.
	# 1.2) Scale Features, minmax .1-1, replace np.nan with 0 after scaling
	ScaledMethylDF, MetaDF = prepMethylationData(MethylationFilePath, 
																WindowSizeList)
	# 2) Combines data into reshaped arrays for CNN model input
	# 2.1) Split features by corresponding window size
	# 2.2) Interleave columns, 2.2.1) Get column order by sorting
	# 2.3) Interleave rows (promoters/transcripts) on index (resolutions)
	# 2.4) reshape DF to arrray
	ReshapedImageArray, SortedIndex = runCombineDataPipeline(ScaledMethylDF, 
														NumOfWin, WindowSizeList)
	# 3) load Response variable data if available
	if LineArgs.ResponseVariablePath:
		ResponseVarDF = pd.read_csv(LineArgs.ResponseVariablePath, header="infer",sep="\t")
		ResponseVarDF.set_index("EnsmblID_T",inplace=True,drop=True)
		# keeping only the relevant column (response variable)
		ResponseVarDF = ResponseVarDF[["log2_ChipDivInput"]]
		# concat with MetaDF for easy sorting and addition to H5 file
		MetaDF = pd.concat([MetaDF,ResponseVarDF],axis=1)
	# 4) sort meta df on index to match corresponding promoter data
	SortedMetaDF = sortMetaDF(SortedIndex, MetaDF)
	# 5) convert arrays  and metaDF to hdf5 datasets:
	createHDF5(ReshapedImageArray, SortedMetaDF, OutputFilePath)
	print("Total time for Image Feature Processing,", str(time.time()-T0)+"s")

# Capture command line args, with or without defaults
if __name__ == '__main__':
	# Parse the arguments
	LineArgs = parseArguments()
	main(LineArgs)
