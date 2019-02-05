#!/bin/python

import sys
import numpy as np
import re

# The functions below are simulating deleted transcripts and fusion transcripts
def ReadGenome(fafile):
	genome={}
	fp=open(fafile,'r')
	line=fp.readline().strip()
	tmpseq=''
	tmpname=''
	while line!='':
		if line[0]=='>':
			if len(tmpseq)!=0:
				genome[tmpname]=tmpseq
			tmpseq=''
			tmpname=re.split(" |\|", line)[0][1:]
		else:
			tmpseq+=line
		line=fp.readline().strip()
	genome[tmpname]=tmpseq
	fp.close()
	return genome


def ReadQuant(quantfile):
	fp=open(quantfile, 'r')
	linecount=0
	TPM={}
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		TPM[strs[0]]=float(strs[4])/float(strs[1])
	fp.close()
	return TPM


def SimuFusion(genome, TPM, n=50, mincutthresh=20):
	TransNames=np.random.permutation(list(TPM.keys()))
	fusion={}
	count=0
	for i in range(int(len(TransNames)/2)):
		t1=TransNames[(2*i)]
		t2=TransNames[(2*i+1)]
		seq1=genome[t1]
		seq2=genome[t2]
		if len(seq1)<2*mincutthresh or len(seq2)<2*mincutthresh:
			continue
		cut1=np.random.randint(mincutthresh, len(seq1)-mincutthresh)
		cut2=np.random.randint(mincutthresh, len(seq2)-mincutthresh)
		newseq=seq1[:cut1]+seq2[cut2:]
		newname="Fusion"+str(count).zfill(4)+" "+t1+":0:"+str(cut1)+" "+t2+":"+str(cut2)+":"+str(len(seq2))
		fusion[newname]=newseq
		# record number of fusion event created
		count+=1
		if count==n:
			break
	return fusion


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


def SimuDeletion_MoreRandom(genome, GeneTransMap, TPM, n=100):
	Removelist = []
	GeneList = list(GeneTransMap.keys())
	permute = np.random.permutation(len(GeneList))
	for i in permute:
		g = GeneList[i]
		if not (g in GeneTransMap):
			print(g)
		v = [t for t in GeneTransMap[g] if (t in TPM) and (TPM[t]>0)]
		if len(v) == 0:
			continue
		ind = np.random.randint(len(v))
		Removelist.append(v[ind])
		if len(Removelist) >= n:
			break
	RemoveSet = set(Removelist)
	newgenome = {}
	for t,v in genome.items():
		if not (t in RemoveSet):
			newgenome[t] = v
	return [Removelist, newgenome]


def WriteNewGenome(genome, outfile):
	fp=open(outfile, 'w')
	for k,v in genome.items():
		fp.write(">"+k+"\n")
		i=0
		while i<len(v):
			fp.write(v[i:min(i+80, len(v))]+"\n")
			i+=80
	fp.close()


def WriteList(List, outfile):
	fp=open(outfile, 'w')
	for i in range(len(List)):
		fp.write(List[i]+"\n")
	fp.close()


def WriteTheoreticalExpression(TPM, genome, fusion, outfile):
	# print original NumReads
	originalNumReads = np.sum([v*len(genome[k]) for k,v in TPM.items()])
	print("original NumReads = "+str(originalNumReads))
	# write new theoretical NumReads
	mean_logTPM = np.mean([np.log(x) for x in TPM.values() if x>0])
	std_logTPM = np.std([np.log(x) for x in TPM.values() if x>0])
	fp=open(outfile, 'w')
	fp.write("# Name\tNumReads\n")
	sumnumreads = 0
	for k,v in TPM.items():
		numreads=int(round(len(genome[k])*v))
		sumnumreads += numreads
		fp.write("{}\t{}\n".format(k, numreads))
	print("Without fusion = "+str(sumnumreads))
	for k,v in fusion.items():
		logtpm = np.random.normal(loc=mean_logTPM, scale=std_logTPM)
		tpm = np.exp(logtpm)
		numreads=int(round(len(v)*tpm))
		sumnumreads += numreads
		fp.write("{}\t{}\n".format(k, numreads))
	print("With fusion = "+str(sumnumreads))
	fp.close()


# functions below are evaluating the EMD of affected transcripts
def ReadDeletion(filename):
	fp=open(filename, 'r')
	removelist=[]
	for line in fp:
		strs=line.strip().split("\t")
		removelist.append(strs[0])
	fp.close()
	return removelist


def ReadFusion(filename):
	fusionlist=[]
	fp=open(filename, 'r')
	for line in fp:
		if line[:7]=='>Fusion':
			strs=line.strip().split(" ")
			mergedname=strs[0][1:]
			tname1=strs[1].split(":")[0]
			tname2=strs[2].split(":")[0]
			fusionlist.append([mergedname, tname1, tname2])
	fp.close()
	return fusionlist


def ReadQuantError(filename):
	QuantError={}
	fp=open(filename, 'r')
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		QuantError[strs[0]]=[float(strs[4]), float(strs[5])] # pair<salmon quantification, simulated reads>
	fp.close()
	return QuantError


def WriteEvaluation(RANKfile, OUTfile, removelist, fusionlist, TransGeneMap):
	fpin=open(RANKfile, 'r')
	fpout=open(OUTfile, 'w')
	fpout.write("Name\tEMD\tRank\tHit\tEventName\n")
	Hit_remove=[0]*len(removelist)
	Hit_fusion=[0]*len(fusionlist)
	linecount=0
	for line in fpin:
		if line[0]=='#':
			continue
		linecount+=1
		strs=line.strip().split("\t")
		events=[]
		for i in range(len(removelist)):
			# if Hit_remove[i]==0 and TransGeneMap[strs[0]]==TransGeneMap[removelist[i]]:
			if TransGeneMap[strs[0]]==TransGeneMap[removelist[i]]:
				Hit_remove[i]=1
				events.append(removelist[i])
		for i in range(len(fusionlist)):
			# if Hit_fusion[i]==0 and TransGeneMap[strs[0]]==TransGeneMap[fusionlist[i][1]]:
			if TransGeneMap[strs[0]]==TransGeneMap[fusionlist[i][1]]:
				Hit_fusion[i]=1
				events.append(fusionlist[i][0])
			# elif Hit_fusion[i]==0 and TransGeneMap[strs[0]]==TransGeneMap[fusionlist[i][2]]:
			elif TransGeneMap[strs[0]]==TransGeneMap[fusionlist[i][2]]:
				Hit_fusion[i]=1
				events.append(fusionlist[i][0])
		# quanterror=QuantError[strs[0]]
		# if (quanterror[0]-quanterror[1])/max(1,quanterror[1]) < 1.0/40:
		# 	events=[]
		if len(events)==0:
			fpout.write("{}\t{}\t{}\t{}\t.\n".format(strs[0], strs[1], linecount, 0))
		else:
			fpout.write("{}\t{}\t{}\t{}\t{}\n".format(strs[0], strs[1], linecount, 1, ",".join(events)))
	fpin.close()
	fpout.close()


def WriteGeneEvaluation(RANKfile, OUTfile, removelist, fusionlist, TransGeneMap):
	fpin=open(RANKfile, 'r')
	fpout=open(OUTfile, 'w')
	fpout.write("# GeneName\tEMD\tRank\tHit\tEventName\n")
	linecount=0
	for line in fpin:
		if line[0]=='#':
			continue
		linecount+=1
		strs=line.strip().split("\t")
		events=[]
		for i in range(len(removelist)):
			if strs[0]==TransGeneMap[removelist[i]]:
				events.append(removelist[i])
		for i in range(len(fusionlist)):
			if strs[0]==TransGeneMap[fusionlist[i][1]]:
				events.append(fusionlist[i][0])
			elif strs[0]==TransGeneMap[fusionlist[i][2]]:
				events.append(fusionlist[i][0])
		# quanterror=QuantError[strs[0]]
		# if (quanterror[0]-quanterror[1])/max(1,quanterror[1]) < 1.0/40:
		# 	events=[]
		if len(events)==0:
			fpout.write("{}\t{}\t{}\t{}\t.\n".format(strs[0], strs[1], linecount, 0))
		else:
			fpout.write("{}\t{}\t{}\t{}\t{}\n".format(strs[0], strs[1], linecount, 1, ",".join(events)))
	fpin.close()
	fpout.close()


if __name__=="__main__":
	if len(sys.argv)==1:
		print("python3 Simulation.py <FAfile> <GTFfile> <QUANTfile> <OUTprefix> [<nfusion> <ndeletion>]")
		print("python3 Simulation.py <evatrans or evagenes> <SimuRemove> <SimuFasta> <RANKfile> <GTFfile> <OUTfile>")
	elif sys.argv[1][:3]!="eva":
		FAfile=sys.argv[1]
		GTFfile=sys.argv[2]
		QUANTfile=sys.argv[3]
		OUTprefix=sys.argv[4]
		nfusion=50
		if len(sys.argv)>5:
			nfusion=int(sys.argv[5])
		ndeletion=150
		if len(sys.argv)>6:
			ndeletion=int(sys.argv[6])

		genome=ReadGenome(FAfile)
		TPM=ReadQuant(QUANTfile)
		# removing ERCC spike in
		TPM = {k:v for k,v in TPM.items() if not ("ERCC" in k)}
		# adding transcripts missing TPM
		for k in genome.keys():
			if not (k in TPM):
				TPM[k] = 0
		[GeneTransMap, TransGeneMap]=ReadIsoformFamily(GTFfile)
		fusion=SimuFusion(genome, TPM, nfusion)
		[Removelist, genome2]=SimuDeletion_MoreRandom(genome, GeneTransMap, TPM, ndeletion)
		genome.update(fusion)
		WriteNewGenome(genome, OUTprefix+"_target.fa")
		WriteNewGenome(genome2, OUTprefix+"_reference.fa")
		WriteList(Removelist, OUTprefix+"_removelist.txt")
		WriteTheoreticalExpression(TPM, genome, fusion, OUTprefix+"_theoexp.txt")
	else:
		SimuRemove=sys.argv[2]
		SimuFasta=sys.argv[3]
		RANKfile=sys.argv[4]
		GTFfile=sys.argv[5]
		OUTfile=sys.argv[6]

		[GeneTransMap, TransGeneMap]=ReadIsoformFamily(GTFfile)
		removelist=ReadDeletion(SimuRemove)
		fusionlist=ReadFusion(SimuFasta)
		if sys.argv[1]=="evatrans":
			WriteEvaluation(RANKfile, OUTfile, removelist, fusionlist, TransGeneMap)
		elif sys.argv[1]=="evagenes":
			WriteGeneEvaluation(RANKfile, OUTfile, removelist, fusionlist, TransGeneMap)
