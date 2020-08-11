#!/usr/bin/env python
'''
Purpose:
Calculate the response variable (histone enrichment) for each promoter region, 
for transfer learning.

General Notes:
1) Get the sum of the signal in promoter regions for:
	A) Chip-seq (H3K27ac or H3K4me3) termed ChIPSum
	B) Input (Chip-seq control) termed InputSum
2) calculate the chip-seq signal enrichment (response variable) as:
	log2_ChipDivInput = log2((ChipSumSignal+Alpha)/(InputSumSignal+Alpha))
Where 
	Alpha = ["Input_SumSig"].quantile(.25)
	log2_ChipDivInput = the response variable for transfer learning

#------------------------------------------------------------------------------
Input required:
1) ChipSeq.bw and matching Input.bw
2) Parsed Transcript/promoter definition file in this format:
	EnsmblID_T	EnsmblID_G	Gene	Strand	Chr	Start	End	RStart	REnd
	(e.g. 2_Promoter_Definitions_hg19.txt)
Where:
	EnsmblID_T = Ensemble transcript ID (unique)
	EnsmblID_G = Ensemble gene ID (not unique)
	Gene = human readable gene name (abbrev, not unique)
	Chr = "chr1, chr2, chrX, etc."
	RStart = TSS - 1000bp
	REnd = TSS + 1000bp
	
#------------------------------------------------------------------------------
Output:
A flat file, tab delimited containing fields:
EnsmblID_T	EnsmblID_G	Gene	Strand	Chr	Start	End	RStart	REnd log2_ChipDivInput

#------------------------------------------------------------------------------
What if...?
	1) Q: the region contains no signal?
		A: All promoter regions without signal are considered 0 values (even after 
		log2 transformation). Therefore, we use a small value calculated from 
		the input quantile(.25) as an adjustment to both the numerator and 
		denominator (see General Notes). This means regions without signal = 
		log2(ChipSum+Alpha/Input+Alpha), or log2(Alpha/Alpha) = log2(1) = 0
'''

import os
import time
import argparse
import numpy as np
import pandas as pd
import pyBigWig

#------------------------------------------------------------------------------
def parseArguments(): 
	# Create argument parser
	parser = argparse.ArgumentParser()   
	# Positional mandatory arguments
	parser.add_argument("ChIP_Path", 
		help="FullPath to ChIP File (bw).", type=str)
	parser.add_argument("Input_Path", 
		help="FullPath to Input File (bw).", type=str)
	parser.add_argument("PromoterDefinitions", type=str, help="Path fo Promoter Definitions file")
	parser.add_argument("--outFileName", help="File name to use for output", type=str)
	parser.add_argument("--outDirectory", help="Directory to write output to", default="output", type=str)
	# Parse arguments
	args = parser.parse_args()
	return args

#------------------------------------------------------------------------------
# 3) Caluclate signal from promoter region
def getSignal(Row, ChipBHW, Inp_BHW):
	Start = int(Row["RStart"])
	End = int(Row["REnd"])
	Chr = str(Row["Chr"])
	#
	ChIP = np.array(ChipBHW.values(Chr, Start, End))
	Input = np.array(Inp_BHW.values(Chr, Start, End))
	#
	ChIP_SumSig = ChIP[~np.isnan(ChIP)].sum()
	Input_SumSig = Input[~np.isnan(Input)].sum()
	return [ChIP_SumSig, Input_SumSig]

#------------------------------------------------------------------------------
# 6) Calculate the response variable from the windowed ChIP_SumSig and Input_SumSig values
def scaleRespVarData(PromoterDF):
	Alpha = PromoterDF["Input_SumSig"].quantile(.25)
	# Add alpha
	PromoterDF["Chip_SumSigAlpha"] = PromoterDF["ChIP_SumSig"].astype(float) + Alpha
	PromoterDF["Input_SumSigAlpha"] = PromoterDF["Input_SumSig"].astype(float) + Alpha
	# calculate response variables
	PromoterDF["log2_ChipDivInput"] = \
					np.log2((PromoterDF["Chip_SumSigAlpha"].astype(float)/\
					PromoterDF["Input_SumSigAlpha"].astype(float)))
	# ensure the only non float value is np.nan
	PromoterDF.replace([np.inf, -np.inf], np.nan, inplace=True)
	# all np.nan are considered 0.0, i.e. np.log2(1)
	# this would be equivalent to log2(0+Alpha/0+Alpha)
	PromoterDF.replace(np.nan, 0.0, inplace=True)
	return PromoterDF

#------------------------------------------------------------------------------
def main(LineArgs):
	T0 = time.time()
	# Return Sample input and output file paths
	Sample = os.path.splitext(os.path.basename(LineArgs.ChIP_Path))[0]
	# determine output directory
	OutputFilePath = LineArgs.outDirectory
	if not os.path.exists(LineArgs.outDirectory):
		os.mkdir(LineArgs.outDirectory)
	# determine fullpath to output file
	if not LineArgs.outFileName:
		OutputFilePath += "/" + Sample + ".txt"
	else:
		OutputFilePath += "/" + LineArgs.outFileName	
	# 1) Load transcript definitions:
	PromoterDefFile = LineArgs.PromoterDefinitions 
	PromoterDF = pd.read_csv(PromoterDefFile, header="infer", sep="\t")
	# 2) Load ChIP-seq data
	ChIP_BHW = pyBigWig.open(LineArgs.ChIP_Path)
	Inp_BHW = pyBigWig.open(LineArgs.Input_Path)
	# 3) Calculate signal from promoter region
	PromoterDF["SignalOutput"] = PromoterDF.apply(getSignal, 
									args=(ChIP_BHW,Inp_BHW), 
									axis=1)
	# 4) Split results into separate signal types
	SigOutCols = ["ChIP_SumSig", "Input_SumSig"]
	PromoterDF[SigOutCols] = pd.DataFrame(PromoterDF.SignalOutput.values.tolist(), 
							index= PromoterDF.index)
	# 5) Calculate the response variable from the windowed ChIP_SumSig and Input_SumSig values
	ScaledRespVarDF = scaleRespVarData(PromoterDF)
	# 6) Remove unparse col of signal types
	ScaledRespVarDF.drop(["SignalOutput", "ChIP_SumSig", "Input_SumSig",
							"Chip_SumSigAlpha", "Input_SumSigAlpha"], 
							inplace=True, axis=1)
	# 7) Save to file
	print("Printing to file", OutputFilePath)
	ScaledRespVarDF.to_csv(OutputFilePath, header=True,sep="\t",index=False)
	#
	T1 = time.time()
	Time = T1 - T0
	print("Total time to complete,", str(Time)+"s")

# Capture command line args, with or without defaults
if __name__ == '__main__':
	# Parse the arguments
	LineArgs = parseArguments()
	main(LineArgs)
