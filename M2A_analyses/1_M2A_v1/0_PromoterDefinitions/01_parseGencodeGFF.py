#!/usr/bin/env python
'''
Purpose:
Return Desired transcripts with < 50% Response Variable region overlap 

Definition of Overlap:
	1) Response variable is cumulative normalized H3K27ac signal 
		1A) Response variable region is 2000bp centered on the TSS.
	2) To avoid exact duplicate records, allowing only 50% overlap.
	3) This means at most 1000bp overlap in response variable region
	4) Because promoters are directional, Negative and Positive Strand
	   transcripts were considered seperately.

Notes:
1) Gencode hg19 was used for transcript isoform definitions
'''

import os
import argparse
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------
def parseArguments(): 
	# Create argument parser
	parser = argparse.ArgumentParser()   
	# Positional mandatory arguments
	parser.add_argument("GFFFilePath", 
		help="Full path to GFF file.", type=str)
	parser.add_argument("--outFileName", help="File name to use for output", default="1_Gencode_parsedGFF.txt", type=str)
	parser.add_argument("--outDirectory", help="Directory to write output to", default="output", type=str)

	# Parse arguments
	args = parser.parse_args()
	return args

def getParsed(ChrDF, StartPos, Col):
	# Grab <=50% overlapping promoter regions
	IndexList = []
	for i, Row in ChrDF.iterrows():
		TSS = Row[Col]
		if TSS > StartPos:
			StartPos = TSS + 1000
			# place non-overlap promoters in list
			IndexList.append(i)	
	return ChrDF.loc[IndexList]


def parseGFFPromoter(DF, Strand, Overlap):
	KeepChrs = ["chr"+str(x) for x in range(1,23)]
	DF = DF[DF["strand"]==Strand]
	# iterate by chromosome to avoid positional overlap
	ChrDFList = []
	for Chr in KeepChrs:
		ChrDF = DF[DF["chr"]==Chr]
		# Determine strandedness
		if Strand == "-":
			StartPos = -1e12
			ChrDF["end"] = ChrDF["end"] * -1
			ChrDF.sort_values(["end"], ascending=True, inplace=True)
			ParsedDF = getParsed(ChrDF, StartPos, "end")
			ParsedDF["end"] = ParsedDF["end"] * -1
		elif Strand == "+":
			StartPos = 0
			ChrDF.sort_values(["start"], ascending=True, inplace=True)
			ParsedDF = getParsed(ChrDF, StartPos, "start")
		ChrDFList.append(ParsedDF)
	# Combine results
	ParsedGff = pd.concat(ChrDFList, axis=0)
	print("Unique "+Strand+" Strand; <1000bp Overlap Data shape", ParsedGff.shape)
	return ParsedGff


def getFormatedGFF(ParsedDF):
	# split info col
	ID = ParsedDF["ID"].str.split(";",expand=True)
	# get ID, Parent gene, and gene name
	ID = ID[[0,1,6]].apply(lambda x: ";".join(x), axis=1)
	# drop original
	ParsedDF.drop(["ID"], axis=1, inplace=True)
	# replace with new
	ParsedDF["ID"] = ID
	return ParsedDF


def main(LineArgs):
	# Populate vars with args
	GffFile = LineArgs.GFFFilePath

	# load gff3 file
	MyCols = ["chr", "def", "type", "start", "end", "n1", "strand", "n2", "ID"]
	GffData = pd.read_csv(GffFile, sep="\t",names=MyCols, skiprows=range(0,7))
	print("Raw GFF Shape", GffData.shape)
	# Filter correct chrs
	KeepChrs = ["chr"+str(x) for x in range(1,23)]
	GffData = GffData[GffData["chr"].isin(KeepChrs)]
	# Get transcripts
	GffData = GffData[GffData["type"]=="transcript"]
	# Parse with no more than 50% promoter overlap
	Overlap = 1000 # response var == 2000
	PosGFF = parseGFFPromoter(GffData, "+", Overlap)
	NegGFF = parseGFFPromoter(GffData, "-", Overlap)
	ParsedGff = pd.concat([PosGFF, NegGFF],axis=0)
	#ParsedGff.index = range(0,ParsedGff.shape[0])
	FormatedGFF = getFormatedGFF(ParsedGff)
	print("Total Unique Promotor; <1000bp Overlap Data shape", FormatedGFF.shape)
	OutFile = os.path.join(LineArgs.outDirectory, LineArgs.outFileName)
	if not os.path.exists(LineArgs.outDirectory):
		os.mkdir(LineArgs.outDirectory)

	FormatedGFF.to_csv(OutFile, header=None, index=None, sep="\t")

# Capture command line args, with or without defaults
if __name__ == '__main__':
	# Parse the arguments
	LineArgs = parseArguments()
	main(LineArgs)
