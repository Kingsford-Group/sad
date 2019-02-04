#!/bin/python

import sys


def GetFeature(line, key):
	s=line.index(key)
	t=line.index(";", s+1)
	return line[(s+len(key)+2):(t-1)]


def ReadTransIDs(infasta):
	transids = []
	fp = open(ingtf, 'r')
	for line in fp:
		if line[0] == '>':
			transids.append(line.strip().split("\t"))
	fp.close()
	transids = set(transids)
	return transids


def ExtractGTF(transids, ingtf, outgtf):
	fpin = open(ingtf, 'r')
	fpout = open(outgtf,'w')
	count = 0
	for line in fpin:
		if line[0] == '#':
			continue
		if "transcript_id" in line:
			if GetFeature(line, "transcript_id") in transids:
				count += 1
				fpout.write(line)
	fpin.close()
	fpout.close()
	print("Finish extracting on GTF file. {} line written.".format(count))


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python ExtractPCAnnotation.py <infasta> <ingtf> <outgtf>")
	else:
		infasta = sys.argv[1]
		ingtf = sys.argv[2]
		outgtf = sys.argv[3]

		transids = ReadTransIDs(infasta)
		ExtractFasta(transids, ingtf, outgtf)
