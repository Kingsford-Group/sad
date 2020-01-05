#!@PYTHON3@

import sys
import numpy as np
import statsmodels.stats.multitest

def ReadInputFile(filename):
	PvalueInfo = []
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		PvalueInfo.append( [strs[0], float(strs[1]), float(strs[2]), float(strs[3]), float(strs[4]), float(strs[5]), int(strs[6])] )
	fp.close()
	return PvalueInfo


def ReadInputFile_long(filename):
	PvalueInfo = []
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		# Name   Coverage        AnomalyScorePos RegionStartPos  RegionEndPos    PValue_Pos      AdjPValue_Pos   AnomalyScoreNeg RegionStartNeg  RegionEndNeg    PValue_Neg      AdjPValue_Neg   MinAdjPValue    Choice
		PvalueInfo.append( [strs[0], float(strs[1]), float(strs[2]), int(strs[3]), int(strs[4]), float(strs[5]), float(strs[6]), float(strs[7]), int(strs[8]), int(strs[9]), float(strs[10]), float(strs[11]), float(strs[12]), int(strs[13])] )
	fp.close()
	return PvalueInfo


def WriteOutputFile(filename, PvalueInfo):
	fp = open(filename, 'w')
	fp.write("# Name\tCoverage\tDeletionScorePos\tDeletionScoreNeg\tRawPvalue\tAdjustedPvalue\tChoice\n")
	for t in PvalueInfo:
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(t[0], t[1], t[2], t[3], t[4], t[5], t[6]))
	fp.close()


def WriteOutputFile_long(filename, PvalueInfo):
	PvalueInfo.sort(key = lambda x:x[12])
	fp = open(filename, 'w')
	fp.write("# Name\tCoverage\tAnomalyScorePos\tRegionStartPos\tRegionEndPos\tPValue_Pos\tAdjPValue_Pos\tAnomalyScoreNeg\tRegionStartNeg\tRegionEndNeg\tPValue_Neg\tAdjPValue_Neg\tMinAdjPValue\tChoice\n")
	for t in PvalueInfo:
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11], t[12], t[13]))
	fp.close()


def Adjust(PvalueInfo):
	oldP = np.array([t[4] for t in (PvalueInfo)])
	_,adjP,_,_ = statsmodels.stats.multitest.multipletests(oldP, method="fdr_bh")
	for i in range(len(PvalueInfo)):
		PvalueInfo[i][5] = adjP[i]
	return PvalueInfo


def Adjust_long(PvalueInfo):
	oldP = np.array([min(t[5], t[10]) for t in (PvalueInfo)])
	_,adjP,_,_ = statsmodels.stats.multitest.multipletests(oldP, method="fdr_bh")
	for i in range(len(PvalueInfo)):
		PvalueInfo[i][12] = adjP[i]
	return PvalueInfo


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python3 AdjustPValue.py <mode(0 for short, 1 for long)> <test_pvalue_overall> <test_pvalue_overall_sort>")
	else:
		mode = sys.argv[1]
		InputFile = sys.argv[2]
		OutputFile = sys.argv[3]

		if mode == '0':
			PvalueInfo = ReadInputFile(InputFile)
			PvalueInfo = Adjust(PvalueInfo)
			WriteOutputFile(OutputFile, PvalueInfo)
		elif mode == '1':
			PvalueInfo = ReadInputFile_long(InputFile)
			#PvalueInfo = Adjust_long(PvalueInfo)
			WriteOutputFile_long(OutputFile, PvalueInfo)
