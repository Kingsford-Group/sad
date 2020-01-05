#!@PYTHON3@

import numpy as np

class Transcript_t(object):
	def __init__(self, _TransID, _GeneID, _GeneName, _Chr, _Strand, _StartPos, _EndPos):
		self.TransID=_TransID
		self.GeneID=_GeneID
		self.GeneName = _GeneName
		self.Chr = _Chr
		self.Strand = _Strand
		self.StartPos = _StartPos
		self.EndPos = _EndPos
		self.Exons = [] # each exon is a tuple of two integers
		self.CDS = []
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
	extraCDS = []
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[2]=="transcript":
			if tmptransname!="" and not (tmptranscript is None):
				Transcripts[tmptransname]=tmptranscript
			tmptransname=GetFeature(line, "transcript_id")
			tmpgeneid=GetFeature(line, "gene_id")
			tmpgenename = ""
			if "gene_name" in line:
				tmpgenename = GetFeature(line, "gene_name")
			tmptranscript=Transcript_t(tmptransname, tmpgeneid, tmpgenename, strs[0], (strs[6]=="+"), int(strs[3])-1, int(strs[4]))
		elif strs[2]=="exon":
			thistransid=GetFeature(line, "transcript_id")
			if thistransid == tmptransname and not (tmptranscript is None):
				tmptranscript.Exons.append((int(strs[3])-1, int(strs[4])))
			else:
				extraExons.append([thistransid, int(strs[3])-1, int(strs[4])])
		elif strs[2] == "CDS":
			thistransid=GetFeature(line, "transcript_id")
			if thistransid == tmptransname and not (tmptranscript is None):
				tmptranscript.CDS.append((int(strs[3])-1, int(strs[4])))
			else:
				extraCDS.append([thistransid, int(strs[3])-1, int(strs[4])])
	if tmptransname!="" and not (tmptranscript is None):
		Transcripts[tmptransname]=tmptranscript
	for e in extraExons:
		assert(e[0] in Transcripts)
		Transcripts[e[0]].Exons.append((e[1],e[2]))
	for e in extraCDS:
		assert(e[0] in Transcripts)
		Transcripts[e[0]].CDS.append((e[1],e[2]))
	for t in Transcripts:
		Transcripts[t].Exons.sort(key=lambda x:x[0])
		Transcripts[t].CDS.sort(key = lambda x:x[0])
		if not Transcripts[t].Strand:
			Transcripts[t].Exons = Transcripts[t].Exons[::-1]
			Transcripts[t].CDS = Transcripts[t].CDS[::-1]
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
	for g,v in GeneTransMap.items():
		sortedv = sorted(v)
		GeneTransMap[g] = sortedv
	return [GeneTransMap, TransGeneMap]


def GetTransLength(Transcripts):
	TransLength = {t:np.sum([e[1]-e[0] for e in v.Exons]) for t,v in Transcripts.items()}
	return TransLength


def ReadTranscriptFasta(filename, namesplitter = " "):
	TransSequence = {}
	fp = open(filename, 'r')
	name = ""
	seq = ""
	for line in fp:
		if line[0] == '>':
			if len(name) != 0:
				TransSequence[name] = seq
			name = line.strip().split(namesplitter)[0][1:]
			seq = ""
		else:
			seq += line.strip()
	if len(name) != "":
		TransSequence[name] = seq
	fp.close()
	return TransSequence


def GetSharedExon(Transcripts, tnames):
	assert(len(tnames) != 0)
	Chr = Transcripts[tnames[0]].Chr
	Strand = Transcripts[tnames[0]].Strand
	for t in tnames:
		assert(Chr == Transcripts[t].Chr)
		assert(Strand == Transcripts[t].Strand)
	allexons = []
	for t in tnames:
		allexons += Transcripts[t].Exons
	allexons.sort(key=lambda x:x[0])
	uniqexons = []
	for e in allexons:
		if len(uniqexons) == 0 or uniqexons[-1][1] < e[0]:
			uniqexons.append([e[0], e[1]])
		else:
			uniqexons[-1][1] = max(uniqexons[-1][1], e[1])
	if not Strand:
		uniqexons = uniqexons[::-1]
	return uniqexons


def JaccardDistance(Transcripts, t1, t2):
	jacard = 0
	if Transcripts[t1].Chr != Transcripts[t2].Chr or Transcripts[t1].Strand != Transcripts[t2].Strand:
		return jacard
	t1_exons = Transcripts[t1].Exons
	t2_exons = Transcripts[t2].Exons
	# calculate union exons
	tmpunion_exons = [[e[0], e[1]] for e in t1_exons] + [[e[0], e[1]] for e in t2_exons]
	tmpunion_exons.sort(key = lambda x:x[0])
	union_exons = []
	for e in tmpunion_exons:
		if len(union_exons) == 0 or union_exons[-1][1] < e[0]:
			union_exons.append(e)
		else:
			union_exons[-1][1] = max(union_exons[-1][1], e[1])
	union_exons = [(e[0], e[1]) for e in union_exons]
	union_length = np.sum([e[1]-e[0] for e in union_exons])
	t1_length = np.sum([e[1]-e[0] for e in t1_exons])
	t2_length = np.sum([e[1]-e[0] for e in t2_exons])
	jacard = (union_length - t1_length + union_length - t2_length) / union_length
	return jacard
