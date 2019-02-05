#!/bin/python
# The script answers the following question:
# for under-expressed region predicted by SAD, does scallop/stringtie also predicts a transcripts corresponding to the deletion and contains all other splicing sites?

import sys
import copy
import pickle
from pathlib import Path
from TranscriptClass import *


class SADpred(object):
	def __init__(self, _TransID, _GeneID, _Score, _Region, _IsOver, _AdjPvalue):
		self.TransID = _TransID
		self.GeneID = _GeneID
		self.Score = _Score
		self.Region = _Region
		self.IsOver = _IsOver
		self.AdjPvalue = _AdjPvalue
		self.Correspondence = None


def ReadSADprediction(filename, pthresh=0.01):
	pred = []
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if float(strs[9]) < pthresh:
			tmp = SADpred(strs[0], strs[11], float(strs[2]), [int(strs[3]), int(strs[4])], strs[10]=="1", float(strs[9]))
			if tmp.IsOver:
				tmp.Score = float(strs[5])
				tmp.Region = [int(strs[6]), int(strs[7])]
			pred.append(tmp)
	fp.close()
	return pred


def ReadSADprediction_long(filename, pthresh=0.01):
	pred = []
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if float(strs[12]) < pthresh:
			tmp = SADpred(strs[0], strs[14], float(strs[2]), [int(strs[3]), int(strs[4])], strs[13]=="1", float(strs[12]))
			if tmp.IsOver:
				tmp.Score = float(strs[7])
				tmp.Region = [int(strs[8]), int(strs[9])]
			assert(tmp.Region[0] >= 0 and tmp.Region[1] >= tmp.Region[0])
			pred.append(tmp)
	fp.close()
	return pred


def ConvertRegion2GenomeCoord(t, region):
	region_convert = []
	coveredlen = 0
	if t.Strand:
		for e in t.Exons:
			if max(coveredlen, region[0]) < min(coveredlen+e[1]-e[0], region[1]):
				overlap_start = max(coveredlen, region[0])
				overlap_end = min(coveredlen+e[1]-e[0], region[1])
				region_convert.append( [e[0]+overlap_start-coveredlen, e[0]+overlap_end-coveredlen] )
			coveredlen += e[1] - e[0]
	else:
		for e in t.Exons:
			if max(coveredlen, region[0]) < min(coveredlen+e[1]-e[0], region[1]):
				overlap_start = max(coveredlen, region[0])
				overlap_end = min(coveredlen+e[1]-e[0], region[1])
				region_convert.append( [e[1]-(overlap_end-coveredlen), e[1]-(overlap_start-coveredlen)] )
			coveredlen += e[1] - e[0]
	region_convert = [(e[0],e[1]) for e in region_convert]
	assert(np.sum([e[1]-e[0] for e in region_convert]) == region[1] - region[0])
	return region_convert


def CheckAssembledTranscript(t_asm, t_ref, region_convert, isoverexpress = False, window = 5):
	if t_asm.Chr != t_ref.Chr or max(t_asm.StartPos, t_ref.StartPos) > min(t_asm.EndPos, t_ref.EndPos):
		return False
	oldStrand = t_asm.Strand
	# considering that transcriptome assembly may have wrong inferred strand, make it the same as reference
	if t_ref.Strand:
		t_asm.Exons.sort(key = lambda x:x[0])
	else:
		t_asm.Exons.sort(key = lambda x:x[0])
		t_asm.Exons = t_asm.Exons[::-1]
	# test whether the assembled transcripts correspond to SAD predicted over and under expression
	if not isoverexpress: # deletion
		# calculate length of overlapping with deleted region_convert
		overlaplength = 0
		for e in t_asm.Exons:
			overlaplength += np.sum( [max(0, min(e[1],e2[1]) - max(e[0],e2[0])) for e2 in region_convert] )
		overlapratio = 1.0*overlaplength / np.sum([e[1]-e[0] for e in region_convert])
		# label whether each splicing site in t_ref is covered by t_asm
		CoverSplicing = []
		RelatedDeletion = []
		if t_ref.Strand:
			for i in range(len(t_ref.Exons)-1):
				exist = np.sum( [abs(t_asm.Exons[j][1] - t_ref.Exons[i][1]) < window and abs(t_asm.Exons[j+1][0] - t_ref.Exons[i+1][0]) < window for j in range(len(t_asm.Exons)-1)] )
				CoverSplicing.append( (exist > 0) )
				related = np.sum([max(t_ref.Exons[i][0],region_convert[j][0]) < min(t_ref.Exons[i][1],region_convert[j][1]) for j in range(len(region_convert))]) \
					+ np.sum([max(t_ref.Exons[i+1][0],region_convert[j][0]) < min(t_ref.Exons[i+1][1],region_convert[j][1]) for j in range(len(region_convert))])
				RelatedDeletion.append( (related > 0) )
		else:
			for i in range(len(t_ref.Exons)-1):
				exist = np.sum( [abs(t_asm.Exons[j+1][1] - t_ref.Exons[i+1][1]) < window and abs(t_asm.Exons[j][0] - t_ref.Exons[i][0]) < window for j in range(len(t_asm.Exons)-1)] )
				CoverSplicing.append( (exist > 0) )
				related = np.sum([max(t_ref.Exons[i][0],region_convert[j][0]) < min(t_ref.Exons[i][1],region_convert[j][1]) for j in range(len(region_convert))]) \
					+ np.sum([max(t_ref.Exons[i+1][0],region_convert[j][0]) < min(t_ref.Exons[i+1][1],region_convert[j][1]) for j in range(len(region_convert))])
				RelatedDeletion.append( (related > 0) )
		# check if assembled transcripts contains all splicing junctions non-related to deletion
		CoverSplicing = np.array(CoverSplicing)
		RelatedDeletion = np.array(RelatedDeletion)
		splicing_match = np.all(CoverSplicing[np.where(np.logical_not(RelatedDeletion))[0]])
		return overlapratio < 0.5 and splicing_match
	else: # over-expression
		# calculate length of overlapping with the non-over-expressed region_convert
		overlaplength_over = 0
		for e in t_asm.Exons:
			overlaplength_over += np.sum( [max(0, min(e[1],e2[1]) - max(e[0],e2[0])) for e2 in region_convert] )
		overlaplength_all = 0
		for e in t_asm.Exons:
			overlaplength_all += np.sum( [max(0, min(e[1],e2[1]) - max(e[0],e2[0])) for e2 in t_ref.Exons] )
		# ratio of overlapping with non-over-expressed region
		overlapratio = (overlaplength_all - overlaplength_over) / (np.sum([e[1]-e[0] for e in t_ref.Exons]) - np.sum([e[1]-e[0] for e in region_convert]))
		# label whether each splicing site in t_ref is covered by t_asm
		CoverSplicing = []
		RelatedOverExp = []
		if t_ref.Strand:
			for i in range(len(t_ref.Exons)-1):
				exist = np.sum( [abs(t_asm.Exons[j][1] - t_ref.Exons[i][1]) < window and abs(t_asm.Exons[j+1][0] - t_ref.Exons[i+1][0]) < window for j in range(len(t_asm.Exons)-1)] )
				CoverSplicing.append( (exist > 0) )
				related1 = np.sum([max(t_ref.Exons[i][0],region_convert[j][0]) < min(t_ref.Exons[i][1],region_convert[j][1]) for j in range(len(region_convert))])
				related2 = np.sum([max(t_ref.Exons[i+1][0],region_convert[j][0]) < min(t_ref.Exons[i+1][1],region_convert[j][1]) for j in range(len(region_convert))])
				RelatedOverExp.append( (related1 > 0 and related2 > 0) )
		else:
			for i in range(len(t_ref.Exons)-1):
				exist = np.sum( [abs(t_asm.Exons[j+1][1] - t_ref.Exons[i+1][1]) < window and abs(t_asm.Exons[j][0] - t_ref.Exons[i][0]) < window for j in range(len(t_asm.Exons)-1)] )
				CoverSplicing.append( (exist > 0) )
				related1 = np.sum([max(t_ref.Exons[i][0],region_convert[j][0]) < min(t_ref.Exons[i][1],region_convert[j][1]) for j in range(len(region_convert))])
				related2 = np.sum([max(t_ref.Exons[i+1][0],region_convert[j][0]) < min(t_ref.Exons[i+1][1],region_convert[j][1]) for j in range(len(region_convert))])
				RelatedOverExp.append( (related1 > 0 and related2 > 0) )
		# check if the assembled transcripts contains all splicing junctions within the over-expression region
		CoverSplicing = np.array(CoverSplicing)
		RelatedOverExp = np.array(RelatedOverExp)
		splicing_match = np.all( CoverSplicing[np.where(RelatedOverExp)[0]] )
		return overlapratio < 0.5 and splicing_match


def FindPotentialCorrespondence(RefTranscripts, AsmTranscripts, SADpredictions):
	NewSADpredictions = []
	for pred in SADpredictions:
		t_ref = RefTranscripts[pred.TransID]
		region_convert = ConvertRegion2GenomeCoord(t_ref, pred.Region)
		for tname, t_asm in AsmTranscripts.items():
			if CheckAssembledTranscript(t_asm, t_ref, region_convert, pred.IsOver):
				pred.Correspondence = tname
				break
		NewSADpredictions.append(pred)
	return NewSADpredictions


def WriteCorrespondence(RefTranscripts, AsmTranscripts, SADpredictions, outputfile):
	fp = open(outputfile, 'w')
	fp.write("# TransID\tGeneID\tScore\tRegionStart\tRegionEnd\tIsOverExpression\tCorrespondingAssembly\tEitherEnds\n")
	for pred in SADpredictions:
		t_ref = RefTranscripts[pred.TransID]
		translength = np.sum([e[1]-e[0] for e in t_ref.Exons])
		isend = (pred.Region[0] < translength*0.05) or (pred.Region[1] > translength*0.95)
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pred.TransID, pred.GeneID, pred.Score, pred.Region[0], pred.Region[1], int(pred.IsOver), pred.Correspondence, int(isend)))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python3 FindSADpredinAssembly.py <mode(0 for short, 1 for long)> <RefGTFfile> <AssembleGTFfile> <SADpredFile> <OutputFile>")
	else:
		mode = sys.argv[1]
		RefGTFfile = sys.argv[2]
		AssembleGTFfile = sys.argv[3]
		SADpredFile = sys.argv[4]
		OutputFile = sys.argv[5]

		RefTranscripts = ReadGTF(RefGTFfile)
		AsmTranscripts = ReadGTF(AssembleGTFfile)
		if mode == "0":
			SADpredictions = ReadSADprediction(SADpredFile)
		elif mode == "1":
			SADpredictions = ReadSADprediction_long(SADpredFile)

		print("length of SAD prediction = {}".format(len(SADpredictions)))
		SADpredictions = FindPotentialCorrespondence(RefTranscripts, AsmTranscripts, SADpredictions)
		WriteCorrespondence(RefTranscripts, AsmTranscripts, SADpredictions, OutputFile)