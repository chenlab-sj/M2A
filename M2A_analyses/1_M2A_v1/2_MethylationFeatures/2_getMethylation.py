#!/usr/bin/env python
'''
Requires:
python/3.6.5+
time
argparse
numpy as np
pandas as pd
pandarallel
os
subprocess
warnings


Purpose:
Extract methylation features on a per window, promoter (region) basis.

General Notes:
20 Windows/ promoter/ windowsize
2 window sizes: 250bp, and 2500bp
4 features/ window
total windows/promoter = 40
total features/promoter = 160

Diagram (strand sensitive):
RegionStart|----------|TSS|----------|RegionEnd
          W-10,W-9...W-1,W1,W2 ... W10

#------------------------------------------------------------------------------
Input Required:
1) Methylation bed-like file.
	1A) Tab delimited
	2A) Required headers(order does not matter): 
		chrom, (chromosome ID) 
		pos, (position of the C in the CpG)
		mval, (caluclated mvalue of a given CpG) 

2) Chromosome ID, e.g. 1,2,3...X,Y

3) By default, uses the promoter definitions previously generated:
	2_Promoter_Definitions_hg19.txt

Input Notes:
	Not all WGBS output is the same, thus it is possible that
	np.nan occur in the mvalue fields; CpG sites with associated
	mvalues are not constant, so a given window between two samples
	may have a different number of Mvalues.

#------------------------------------------------------------------------------
Output:
1)This script produces four main features:
	1: Window Average Mval
	2: Window Sum of Squared Differences (SSD) (of Mval)
	3: Window Fraction of SSD (FracSSD)
		defined as WindowSSD/RegionSSD
	4: Window Mvalue variability
		defined as WindowSSD / Number of Mvals

Two window sizes are used to calculate features: 250bp, and 2500bp. Capturing
different "resolutions" actually enables us to capture different trends, and 
also extend the features further away from the TSS. Mulitple resolutions can be
included in the model using a CNN which will learn which features are important
for prediction.

#------------------------------------------------------------------------------
What if...?
	1) Q: the region extends passed chromosome boundaries, especially with larger
		window sizes?
		A: these windows will return "np.nan", and later be converted to 0 value.
			In fact, any window without mvalues will be returned as np.nan,
			including centromeric or gapped regions (low map, NNNNs)
'''

import time
import argparse
import numpy as np
import pandas as pd
from pandarallel import pandarallel
import os
#from contextlib import contextmanager
#import subprocess
import warnings
warnings.filterwarnings("ignore")

#------------------------------------------------------------------------------
def parseArguments(): 
	# Create argument parser
	parser = argparse.ArgumentParser()   
	# Positional mandatory arguments
	parser.add_argument("MethylFilePath", 
		help="Full path to methylation bed-like file.", type=str)
	parser.add_argument("PromoterDefinitions", type=str, help="Path fo Promoter Definitions file")
	parser.add_argument("--nbWorkers", help="No. of threads to use", default=2, type=int)	
	parser.add_argument("--outFileName", help="File name to use for output", type=str)
	parser.add_argument("--outDirectory", help="Directory to write output to", default="output", type=str)

	# Parse arguments
	args = parser.parse_args()
	return args

#------------------------------------------------------------------------------
# 1a)
def getMethylData(MethylFilePath):
	MethylDF = pd.read_csv(MethylFilePath, sep="\t", header="infer")
	MethylDF["chrom"] = MethylDF["chrom"].astype(str)
	MethylDF["mval"] = MethylDF["mval"].astype(np.float64)
	MethylDF["pos"] = MethylDF["pos"].astype(int)
	print("Methylation data loaded, total pos, "+str(MethylDF.shape[0]))
	return MethylDF

#------------------------------------------------------------------------------
# 1b)
def getMethylDataByChr(MethylDF, Chr):
	Chr = Chr.replace("chr","")
	MethylDF = MethylDF[MethylDF["chrom"].astype(str)==str(Chr)]
	MethylDF.sort_values(by=["pos"], inplace=True)
	MethylDF.set_index(keys=["pos"], inplace=True)
	print("chr"+Chr+" MethylDF Mvalues:", MethylDF.shape)
	return MethylDF

#------------------------------------------------------------------------------
# 2)
def getFeatureRegionDelimiters(PromoterDF, NumWindows, WindowSize):
	#--------------------------------------------------------------------------
	# 2.1)
	def getRegion(Row, NumWindows, WindowSize):
		if Row["Strand"] == "+":
			TSS = Row["Start"]
		else:
			TSS = Row["End"]
		# feature region is centered on the TSS
		FeatureRegionStart = TSS - (WindowSize * (NumWindows/2))
		FeatureRegionEnd = TSS + (WindowSize * (NumWindows/2))
		return FeatureRegionStart, FeatureRegionEnd
	#--------------------------------------------------------------------------
	# define feature region on strandedness, centered on TSS
	PromoterDF["FeatureRegion"] = PromoterDF.parallel_apply(getRegion,
											args=(NumWindows,WindowSize),
											axis=1)
	# split these into separate columns
	PromoterDF[["FeatureRegionStart","FeatureRegionEnd"]] = pd.DataFrame(
											PromoterDF.FeatureRegion.values.tolist(), 
											index=PromoterDF.index)
	# remove unsplit column
	PromoterDF.drop(["FeatureRegion"],inplace=True, axis=1)
	return PromoterDF

#------------------------------------------------------------------------------
# 3)
def computeWindowRanges(NumWindows, WindowSize):
	RegionBinSize = NumWindows * WindowSize
	WindowRanges = [(i, i + WindowSize) for i in range(0, RegionBinSize, 
														WindowSize)]
	return WindowRanges

#------------------------------------------------------------------------------
# 4)
def getFeatures(Row, MethylChrDF, WindowRanges):
	#--------------------------------------------------------------------------
	# 4.1)
	def getMvalueFeatures(WinStart, WinStop, MethylChrDF, RegionSSD):
		# determine cpg sites with signal in this window
		#MethylationDF = MethylChrDF[(MethylChrDF["pos"]>=int(WinStart))&\
		#							(MethylChrDF["pos"]<int(WinStop))]
		MethylationDF = MethylChrDF.loc[int(WinStart):int(WinStop)-1]

		# flatten to array of floats
		MvalueArray = MethylationDF["mval"].values
		# remove all potential nans
		MvalueArray = MvalueArray[~np.isnan(MvalueArray)]
		#
		# Check results, if none, report np.nan
		NumOfMvals = int(MvalueArray.size)
		if NumOfMvals > 0:
			# average mvalue
			Ave = np.mean(MvalueArray)
			# note: can't calculate window SSD, etc. if entire region SSD is 0
			if RegionSSD == 0.0:
				SSD = 0
				FracSSD = 0
				Var = 0
			else:
				SSD = np.square(MvalueArray - MvalueArray.mean()).sum()
				FracSSD = SSD / RegionSSD
				Var = SSD / NumOfMvals
		else:
			Ave = np.nan
			SSD = np.nan
			FracSSD = np.nan
			Var = np.nan

		# we really can't justify extreme precision (not useful)
		StatsList = np.round([Ave, FracSSD, Var, SSD],decimals=3).tolist()
		return StatsList

	#--------------------------------------------------------------------------
	# Define limits of feature region
	RegionStart = int(Row["FeatureRegionStart"])
	RegionEnd = int(Row["FeatureRegionEnd"])
	Strand = str(Row["Strand"])
	# calculate directionality of gene 
	# (differentiate upstream vs downstream)
	if Strand == "+":
		WindowRanges = WindowRanges
	else:
		WindowRanges = WindowRanges[::-1]
	# calculate sum of squared differences (SSD) over entire feature region
	# "1" is passed to act as the pseudo region SSD
	RegionResults = getMvalueFeatures(RegionStart, RegionEnd, MethylChrDF, 1)
	RegionSSD = RegionResults[-1]
	# Calculate each features for each window,
	#DF = MethylChrDF[(MethylChrDF["pos"]>=int(RegionStart))&\
	#			(MethylChrDF["pos"]<int(RegionEnd))]
	DF = MethylChrDF.loc[int(RegionStart):int(RegionEnd)-1]

	WindowResultsList = []
	for WindowRange in WindowRanges:
		FeatureResults = getMvalueFeatures(WindowRange[0]+RegionStart, 
								WindowRange[1]+RegionStart,
								DF, RegionSSD) 
		WindowResultsList = WindowResultsList + FeatureResults
	return WindowResultsList

#------------------------------------------------------------------------------
# 5)
def buildHeader(NumWindows, WindowSize):
	# N+1 to account for upper bounds
	L=int(-(NumWindows/2))
	R=int((NumWindows/2) + 1)
	RangeList = list(range(L,0,1)) + list(range(1,R,1))
	Header = []
	for x in RangeList:
		Header = Header + [str(WindowSize)+"_W"+str(x)+"_M_Ave", 
							str(WindowSize)+"_W"+str(x)+"_M_FracSSD",
							str(WindowSize)+"_W"+str(x)+"_M_Var",
							str(WindowSize)+"_W"+str(x)+"_M_SSD"]
	return Header

#------------------------------------------------------------------------------
def main(LineArgs):
	T0 = time.time()
	# variable input
	MethylFilePath = LineArgs.MethylFilePath
	# hard coded input, subject to change:
	NumWindows = 20
	WindowSizesList = [250, 2500]
	PromoterDefinitionPath = LineArgs.PromoterDefinitions #"../0_PromoterDefinitions/output/2_Promoter_Definitions_hg19.txt"
	# 1) load methylation data
	MethylDF = getMethylData(MethylFilePath)
	# load promoter definitions from file
	PromoterDF = pd.read_csv(PromoterDefinitionPath, header="infer",sep="\t")
	# Iterate over Chrs
	ChrResultsDFList = []
	for Chr in PromoterDF["Chr"].unique():
		print("Processing Features for: "+str(Chr))
		T1 = time.time()
		ChrPromoterDF = PromoterDF[PromoterDF["Chr"]==Chr]
		# iterate over window sizes to generate corresponding features
		ChrMethylDF = getMethylDataByChr(MethylDF, Chr)
		for WindowSize in WindowSizesList:
			T2 = time.time()
			print("\tWindow: "+str(WindowSize))
			# 2) define feature regions
			PromoterRegionDF = getFeatureRegionDelimiters(ChrPromoterDF, NumWindows, 
															WindowSize)
			# 3) generate window ranges for feature calculations
			WindowRanges = computeWindowRanges(NumWindows, WindowSize)
			# 4) generate features by promoter, for a given window size:
			ChrPromoterDF["FeatureOutput"] = PromoterRegionDF.parallel_apply(getFeatures, 
											args=(ChrMethylDF, WindowRanges), axis=1)
			T3 = time.time()
			print("\tTime to compute Window: ", T3-T2)
			# 5) build corresponding feature names (WindowSize, Window, Feature)
			Header = buildHeader(NumWindows, WindowSize)
			# split feature signal into separate signal types
			ChrPromoterDF[Header] = pd.DataFrame(ChrPromoterDF.FeatureOutput.values.tolist(), 
										index=ChrPromoterDF.index)
			# remove unparsed col "FeatureOutput"
			ChrPromoterDF.drop(["FeatureOutput", "FeatureRegionStart","FeatureRegionEnd"], 
							inplace=True, axis=1)
		ChrResultsDFList.append(ChrPromoterDF)
		# report time taken for one chr
		T4 = time.time()
		Time1 = T4 - T1
		print("\tTime to complete "+str(Chr)+", "+str(Time1)+"s")
	# completed features printed to file
	OutputFilePath = LineArgs.outDirectory
	if not os.path.exists(LineArgs.outDirectory):
		os.mkdir(LineArgs.outDirectory)

	if not LineArgs.outFileName:
		OutputFilePath += "/" + "_".join(["Features", MethylFilePath.split("/")[-1]])
	else:
		OutputFilePath += "/" + LineArgs.outFileName

	print ("Printing to file", OutputFilePath)
	# concat results DFs
	PromoterDF = pd.concat(ChrResultsDFList, axis=0)
	PromoterDF.to_csv(OutputFilePath, header=True, sep="\t",index=False)
	# report total time to generate feature for a given chromosome
	T3 = time.time()
	TotalTime = T3 - T0
	print ("Total time to complete,",str(TotalTime)+"s")

# Capture command line args, with or without defaults
if __name__ == '__main__':
	LineArgs = parseArguments()
	par = pandarallel.initialize(nb_workers=LineArgs.nbWorkers)
	# Parse the arguments
	main(LineArgs)
