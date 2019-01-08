#!/bin/python

import sys


def ReadIsoformFamily(gtffile):
	GeneTransMap={}
	TransGeneMap={}
	fp=open(gtffile, 'r')
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[2]=="transcript":
			s=line.index("gene_id")
			t=line.index(";", s+1)
			geneid=line[(s+9):(t-1)]
			s=line.index("transcript_id")
			t=line.index(";", s+1)
			transid=line[(s+15):(t-1)]
			TransGeneMap[transid]=geneid
			if geneid in GeneTransMap:
				GeneTransMap[geneid].append(transid)
			else:
				GeneTransMap[geneid]=[transid]
	fp.close()
	return [GeneTransMap, TransGeneMap]


def CollapsePrediction(inpredfile, outpredfile, TransGeneMap):
	fpin = open(inpredfile, 'r')
	fpout = open(outpredfile, 'w')
	PredictedGenes = {'0'}
	for line in fpin:
		if line[0] == '#':
			fpout.write(line.strip()+"\tGeneName\n")
		else:
			transname = line.strip().split("\t")[0]
			if not (TransGeneMap[transname] in PredictedGenes):
				PredictedGenes.add(TransGeneMap[transname])
				fpout.write(line.strip()+"\t"+TransGeneMap[transname]+"\n")
	fpin.close()
	fpout.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python3 Pred_Trans2Gene.py <GTFfile> <InPredFile> <OutPredFile>")
	else:
		GTFfile = sys.argv[1]
		InPredFile = sys.argv[2]
		OutPredFile = sys.argv[3]

		[GeneTransMap, TransGeneMap] = ReadIsoformFamily(GTFfile)
		CollapsePrediction(InPredFile, OutPredFile, TransGeneMap)