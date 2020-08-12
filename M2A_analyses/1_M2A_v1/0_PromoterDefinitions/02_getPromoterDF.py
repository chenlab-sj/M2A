#!/usr/bin/env python
'''
Requires:
module load python/3.7.0

Purpose:
Format GFF file for feature generation.
'''

import os
import time
import argparse
import pandas as pd

#------------------------------------------------------------------------------
def parseArguments(): 
	# Create argument parser
	parser = argparse.ArgumentParser()   
	# Positional mandatory arguments
	parser.add_argument("GFFFilePath", 
		help="Full path to GFF file.", type=str)
	parser.add_argument("--outFileName", help="File name to use for output", default="Promoter_Definitions.txt", type=str)
	parser.add_argument("--outDirectory", help="Directory to write output to", default="output", type=str)

	# Parse arguments
	args = parser.parse_args()
	return args


#------------------------------------------------------------------------------
def getParsedTranscripts(GffFile, RespVarWin):
	#-------------------------------------------------------------------------#
	def loadGFFData(GffFile):
		# cols in gff:
		ColNames = ['Chr', 
					"CallType", 
					"FeatType", 
					"Start", 
					"End", 
					"Blank1", 
					"Strand", 
					"Blank2", 
					"Attributes"]
		# read file, remove unneeded cols
		GFF_DF = pd.read_csv(GffFile, header=None, sep="\t", names=ColNames)
		print("Raw Transcript Data Shape", GFF_DF.shape)
		GFF_DF.drop(["CallType", 
					"FeatType",
					"Blank1",
					"Blank2"], 
					inplace=True, axis=1)
		return GFF_DF

	#-------------------------------------------------------------------------#
	def parseGFFAttributes(DF):
		AttrCols = ["EnsmblID_T", "EnsmblID_G", "Gene"]
		ParseStrList = ["ID=", "Parent=", "gene_name="]
		DF[AttrCols] = pd.DataFrame(DF["Attributes"].str.split(";",2,expand=True))
		for RepString, Col in zip(ParseStrList, AttrCols):
			DF[Col].replace(RepString, "", regex=True, inplace=True)
		# remove pre-parsed col
		DF.drop(["Attributes"], inplace=True, axis=1)
		# reorder columns
		DF = DF[["EnsmblID_T", 
				"EnsmblID_G", 
				"Gene", 
				"Strand", 
				"Chr", 
				"Start", 
				"End"]]
		return DF

	#-------------------------------------------------------------------------#
	def getResponseRegion(Row, RespVarWin):
		if Row["Strand"] == "+":
			Start = Row["Start"]
		else:
			Start = Row["End"]
		RStart = int(Start - (RespVarWin/2))
		REnd = int(Start + (RespVarWin/2))
		return RStart, REnd

	#-------------------------------------------------------------------------#
	# simple read csv, remove unneeded cols
	GFF_DF = loadGFFData(GffFile)
	# get gene name and ensemble ids
	TranscriptDF = parseGFFAttributes(GFF_DF)
	# define feature region on strandedness
	TranscriptDF["Bins"] = TranscriptDF.apply(getResponseRegion,
											args=(RespVarWin,),
											axis=1)
	Cols = ["RStart", "REnd"]
	TranscriptDF[Cols] = pd.DataFrame(TranscriptDF.Bins.values.tolist(), 
									index=TranscriptDF.index)
	TranscriptDF.drop(["Bins"], inplace=True, axis=1)
	return TranscriptDF

#------------------------------------------------------------------------------
def main(LineArgs):
	T0 = time.time()
	RespVarWin = 2000
	GffFile = LineArgs.GFFFilePath #"output/1_Gencode_parsedGFF.txt"
	TranscriptDF = getParsedTranscripts(GffFile, RespVarWin)
	if not os.path.exists(LineArgs.outDirectory):
		os.mkdir(LineArgs.outDirectory)
	TranscriptDF.to_csv(os.path.join(LineArgs.outDirectory, LineArgs.outFileName), header=True, index=False, sep="\t")
	T1 = time.time()
	Time = T1 - T0
	print("Total time to complete,", str(Time)+"s")

# Capture command line args, with or without defaults
if __name__ == '__main__':
	# Parse the arguments
	LineArgs = parseArguments()
	main(LineArgs)
