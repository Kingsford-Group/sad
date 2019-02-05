#!/bin/python

import sys
import gzip
import numpy as np
import pickle

class Transcript_t(object):
	def __init__(self, _TransID, _GeneID, _Chr, _Strand, _StartPos, _EndPos):
		self.TransID=_TransID
		self.GeneID=_GeneID
		self.Chr = _Chr
		self.Strand = _Strand
		self.StartPos = _StartPos
		self.EndPos = _EndPos
		self.Exons = [] # each exon is a tuple of two integers
	def __eq__(self, other):
		if isinstance(other, Transcript_t):
			return (self.Chr==other.Chr and self.Strand==other.Strand and len(self.Exons)==len(other.Exons) and \
				min([self.Exons[i]==other.Exons[i] for i in range(len(self.Exons))])!=0)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, Transcript_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos<other.StartPos
			else:
				return self.EndPos<other.EndPos
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, Transcript_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos>other.StartPos
			else:
				return self.EndPos<other.EndPos
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result
		

def GetFeature(line, key):
	s=line.index(key)
	t=line.index(";", s+1)
	return line[(s+len(key)+2):(t-1)]


def ReadGTF(gtffile):
	Transcripts={}
	strand=""
	fp=open(gtffile, 'r')
	tmptransname=""
	tmptranscript=None
	extraExons = []
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[2]=="transcript":
			if tmptransname!="" and not (tmptranscript is None):
				Transcripts[tmptransname]=tmptranscript
			tmptransname=GetFeature(line, "transcript_id")
			tmpgeneid=GetFeature(line, "gene_id")
			tmptranscript=Transcript_t(tmptransname, tmpgeneid, strs[0], (strs[6]=="+"), int(strs[3])-1, int(strs[4]))
		elif strs[2]=="exon":
			thistransid=GetFeature(line, "transcript_id")
			if thistransid == tmptransname and not (tmptranscript is None):
				tmptranscript.Exons.append((int(strs[3])-1, int(strs[4])))
			else:
				extraExons.append([thistransid, int(strs[3])-1, int(strs[4])])
	if tmptransname!="" and not (tmptranscript is None):
		Transcripts[tmptransname]=tmptranscript
	for e in extraExons:
		assert(e[0] in Transcripts)
		Transcripts[e[0]].Exons.append((e[1],e[2]))
	for t in Transcripts:
		Transcripts[t].Exons.sort(key=lambda x:x[0])
		if not Transcripts[t].Strand:
			Transcripts[t].Exons = Transcripts[t].Exons[::-1]
	fp.close()
	return Transcripts


def Map_Gene_Trans(Transcripts):
	GeneTransMap={}
	TransGeneMap={}
	for v in Transcripts.values():
		TransGeneMap[v.TransID]=v.GeneID
		if v.GeneID in GeneTransMap:
			GeneTransMap[v.GeneID].append(v.TransID)
		else:
			GeneTransMap[v.GeneID]=[v.TransID]
	return [GeneTransMap, TransGeneMap]


def ReadSimulatedStartPos(targetseqfile, filename):
	# read transcripts length in the target fasta file
	TransLength={}
	Count={}
	fp=open(targetseqfile, 'r')
	prevname=""
	seqlen=0
	for line in fp:
		if line[0]=='>':
			strs=line.strip().split(" ")
			if strs[0][1:] != prevname and prevname!="":
				TransLength[prevname]=seqlen
			prevname=strs[0][1:]
			seqlen=0
		else:
			seqlen+=len(line)-1
	TransLength[prevname]=seqlen
	fp.close()
	# read fastq.gz file and store the read positions
	TrueRaw={}
	fp=gzip.open(filename, 'rb')
	file_content=fp.read()
	fp.close()
	ReadsInfo=file_content.decode('utf-8')
	ReadsInfo=ReadsInfo.strip().split("\n")
	for i in range(len(ReadsInfo)):
		if i%2!=0:
			continue
		strs=ReadsInfo[i].split(";")
		transname=strs[0].split("/")[1].split(" ")[0]
		strs[1]=strs[1][(strs[1].index(":")+1):]
		pos=int(strs[1].split("-")[0])
		if transname in TrueRaw:
			TrueRaw[transname].append(pos)
		else:
			TrueRaw[transname]=[pos]
	# convert the read positions into counts at each position
	newTrueRaw = {}
	for t,v in TrueRaw.items():
		assert(t in TransLength)
		tmp = np.zeros(TransLength[t])
		for pos in v:
			if pos<=0:
				print([t, np.where(v==pos)[0][0]])
			assert(pos>0)
			tmp[pos-1] += 1
		newTrueRaw[t] = tmp
	return newTrueRaw


def SimulationCount(targetseqfile, simufile, binsize=50):
	TransLength={}
	Count={}
	fp=open(targetseqfile, 'r')
	prevname=""
	seqlen=0
	for line in fp:
		if line[0]=='>':
			strs=line.strip().split(" ")
			if strs[0][1:] != prevname and prevname!="":
				TransLength[prevname]=seqlen
			prevname=strs[0][1:]
			seqlen=0
		else:
			seqlen+=len(line)-1
	TransLength[prevname]=seqlen
	fp.close()
	# read simulated fasta.gz file
	fp=gzip.open(simufile, 'rb')
	file_content=fp.read()
	fp.close()
	ReadsInfo=file_content.decode('utf-8')
	ReadsInfo=ReadsInfo.strip().split("\n")
	for i in range(len(ReadsInfo)):
		if i%2!=0:
			continue
		strs=ReadsInfo[i].split(";")
		transname=strs[0].split("/")[1].split(" ")[0]
		strs[1]=strs[1][(strs[1].index(":")+1):]
		pos=int(strs[1].split("-")[0])
		ind=int(pos/binsize)
		if transname in Count:
			Count[transname][ind]+=1
		else:
			tmp=np.zeros((int(TransLength[transname]/binsize)+1), dtype=np.int)
			tmp[ind]+=1
			Count[transname]=tmp
	return Count


def WriteSimulationCount(Count, outputfile):
	fp=open(outputfile, 'w')
	fp.write("Name\tbinstart\tbinend\tnsimulated\n")
	for k,v in Count.items():
		for i in range(len(v)):
			fp.write("{}\t{}\t{}\t{}\n".format(k, i*50, i*50+50, v[i]))
	fp.close()


def ReadSalmonCount(salmonquant):
	SalmonCount={}
	fp=open(salmonquant, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		SalmonCount[strs[0]]=[int(strs[1]), float(strs[2]), float(strs[3]), float(strs[4])] # Length, EffectiveLength, TPM, NumReads
	fp.close()
	return SalmonCount


def WriteSalmonComparison(SalmonCount, SimuCount, OUTfile):
	overlap=set(SalmonCount.keys())&set(SimuCount.keys())
	diffTPM=[]
	for k in overlap:
		[translen, efflen, tpm, numreads]=SalmonCount[k]
		truecount=SimuCount[k]
		diffTPM.append([k, abs(numreads-truecount)/efflen])
	diffTPM.sort(key=lambda x:x[1], reverse=True)
	fp=open(OUTfile, 'w')
	fp.write("# Name\tLength\tEffectiveLength\tTPM\tNumReads\tnsimulated\tdiffTPM\n")
	for e in diffTPM:
		[translen, efflen, tpm, numreads]=SalmonCount[e[0]]
		truecount=SimuCount[e[0]]
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(e[0], translen, efflen, tpm, numreads, truecount, e[1]))
	fp.close()


def ReadKallistoCount(kallistoquant):
	KallistoCount={}
	fp=open(kallistoquant, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		KallistoCount[strs[0]] = [int(strs[1]), float(strs[2]), float(strs[4]), float(strs[3])]
	fp.close()
	return KallistoCount


def WriteKallistoComparison(KallistoCount, SimuCount, OUTfile):
	overlap=set(KallistoCount.keys())&set(SimuCount.keys())
	diffTPM=[]
	for k in overlap:
		[translen, efflen, tpm, numreads]=KallistoCount[k]
		truecount=SimuCount[k]
		diffTPM.append([k, abs(numreads-truecount)/efflen])
	diffTPM.sort(key=lambda x:x[1], reverse=True)
	fp=open(OUTfile, 'w')
	fp.write("# Name\tLength\tEffectiveLength\tTPM\tNumReads\tnsimulated\tdiffTPM\n")
	for e in diffTPM:
		[translen, efflen, tpm, numreads]=KallistoCount[e[0]]
		truecount=SimuCount[e[0]]
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(e[0], translen, efflen, tpm, numreads, truecount, e[1]))
	fp.close()


def MergedGeneDiffTPM(SalmonCount, SimuCount, GeneTransMap):
	DiffGeneTPM={}
	for g,v in GeneTransMap.items():
		diff=0
		for t in v:
			if (t in SalmonCount) and (t in SimuCount):
				[translen, efflen, tpm, numreads]=SalmonCount[t]
				truecount=SimuCount[t]
				if abs(tpm - truecount/efflen) > diff:
					diff = abs(tpm - truecount/efflen)
		DiffGeneTPM[g]=diff
	return DiffGeneTPM


def WriteDiffGeneTPM(DiffGeneTPM, OUTfile):
	fp=open(OUTfile, 'w')
	fp.write("# GeneID\tdiffTPM\n")
	for g,v in DiffGeneTPM.items():
		fp.write("{}\t{}\n".format(g, v))
	fp.close()


def WriteTrueExpression(SimuCount, OUTfile):
	fp=open(OUTfile, 'w')
	fp.write("Name\tnsimulated\n")
	for t,v in SimuCount.items():
		truecount = v
		fp.write("{}\t{}\n".format(t, truecount))
	fp.close()


if __name__=="__main__":
	if len(sys.argv)==1:
		print("python3 GroupTruthExp.py <target_sequence> <sample_01_1.fasta.gz> <salmonquant/kallistoquant> <gtffile> <outputprefix>")
	else:
		TargetseqFile=sys.argv[1]
		FastaGZFile=sys.argv[2]
		QuantFile=sys.argv[3]
		GTFfile=sys.argv[4]
		OUTprefix=sys.argv[5]

		Transcripts=ReadGTF(GTFfile)
		[GeneTransMap, TransGeneMap]=Map_Gene_Trans(Transcripts)

		Simulated = ReadSimulatedStartPos(TargetseqFile, FastaGZFile)
		pickle.dump(Simulated, open(OUTprefix+"simupos.pickle", 'wb'))
		SimuCount = {t:np.sum(v) for t,v in Simulated.items()}
		for tid,t in Transcripts.items():
			if not (tid in SimuCount):
				SimuCount[t] = 0
		WriteTrueExpression(SimuCount, OUTprefix+"expression_truth.txt")

		fp=open(QuantFile, 'r')
		line=fp.readline()
		fp.close()
		if line.strip().split("\t")[0]=="Name":
			SalmonCount=ReadSalmonCount(QuantFile)
			WriteSalmonComparison(SalmonCount, SimuCount, OUTprefix+"expression_difference.txt")
		elif line.strip().split("\t")[0]=="target_id":
			KallistoCount=ReadKallistoCount(QuantFile)
			WriteKallistoComparison(KallistoCount, SimuCount, OUTprefix+"expression_difference.txt")

		# DiffGeneTPM=MergedGeneDiffTPM(SalmonCount, SimuCount, GeneTransMap)
		# WriteDiffGeneTPM(DiffGeneTPM, OUTprefix+"gene_expression_difference.txt")
