#!/bin/python

import sys
import struct
import numpy as np
import pandas as pd
from plotnine import *
from scripts.TranscriptClass import *
from matplotlib import pyplot as plt
from pathlib import Path

plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams.update({'font.size': 13})
# plt.rcParams.update({'font.size': 15})

def ReadRawCorrection(filename):
	ExpectedProb={}
	fp=open(filename, 'rb')
	numtrans=struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen=struct.unpack('i', fp.read(4))[0]
		seqlen=struct.unpack('i', fp.read(4))[0]
		name=""
		correction=np.zeros(seqlen)
		for j in range(namelen):
			name+=struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			correction[j]=struct.unpack('d', fp.read(8))[0]
		ExpectedProb[name]=correction
	fp.close()
	print("Finish reading theoretical distributions for {} transcripts.".format(len(ExpectedProb)))
	return ExpectedProb


def ReadAjustedTheoDistribution(filename, ExpectedProb):
	AdjExpected = {}
	fp = open(filename, 'rb')
	numtrans=struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen=struct.unpack('i', fp.read(4))[0]
		seqlen=struct.unpack('i', fp.read(4))[0]
		name=""
		correction=np.zeros(seqlen)
		for j in range(namelen):
			name+=struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			correction[j]=struct.unpack('d', fp.read(8))[0]
		if name in ExpectedProb:
			# expanding to the full length
			oriexp = ExpectedProb[name]
			adjexp = np.zeros(len(oriexp))
			for j in range(len(correction)):
				ind_start = int(j*len(oriexp)/len(correction))
				ind_end = int((j+1)*len(oriexp)/len(correction))
				if np.sum(oriexp[ind_start:ind_end]) != 0:
					adjexp[ind_start:ind_end] = correction[j] * oriexp[ind_start:ind_end] / np.sum(oriexp[ind_start:ind_end])
				else:
					adjexp[ind_start:ind_end] = correction[j] / (ind_end-ind_start)
			AdjExpected[name] = adjexp
	fp.close()
	print("Finish reading adjusted theoretical distribution, {} out of {} transcripts are expanded to the original length.".format(len(AdjExpected), numtrans))
	return AdjExpected


def ReadCovariance(filename):
	Covariances = []
	fp = open(filename, 'rb')
	numclass = struct.unpack('i', fp.read(4))[0]
	for i in range(numclass):
		matrixsize = struct.unpack('i', fp.read(4))[0]
		cov = np.zeros((matrixsize, matrixsize))
		for j in range(matrixsize):
			for k in range(matrixsize):
				cov[j,k] = struct.unpack('d', fp.read(8))[0]
		Covariances.append(cov)
	LenClass = {}
	numtrans = struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen = struct.unpack('i', fp.read(4))[0]
		name = ""
		for j in range(namelen):
			name += struct.unpack('c', fp.read(1))[0].decode('utf-8')
		classid = struct.unpack('i', fp.read(4))[0]
		LenClass[name] = classid
	fp.close()
	print("Finish reading {} covariance matrix, and LenClass of {} transcripts.".format(len(Covariances), len(LenClass)))
	return [Covariances, LenClass]


def ReadRawStartPos(filename):
	TrueRaw={}
	fp=open(filename, 'rb')
	numtrans=struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen=struct.unpack('i', fp.read(4))[0]
		seqlen=struct.unpack('i', fp.read(4))[0]
		name=""
		poses=np.zeros(seqlen, dtype=np.int)
		counts=np.zeros(seqlen)
		for j in range(namelen):
			name+=struct.unpack('c', fp.read(1))[0].decode('utf-8')
		for j in range(seqlen):
			poses[j]=struct.unpack('i', fp.read(4))[0]
		for j in range(seqlen):
			counts[j]=struct.unpack('d', fp.read(8))[0]
		tmp=np.zeros(poses[-1]+1)
		for j in range(len(poses)):
			tmp[poses[j]] = counts[j]
		TrueRaw[name]=tmp
	fp.close()
	print("Finish reading actual distribution for {} transcripts.".format(len(TrueRaw)))
	return TrueRaw


def ReadLPAdjustment(prefix):
	# read the names of adjusted transcripts in order
	AdjustmentList = []
	fp = open(prefix, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		AdjustmentList.append(strs[0])
	fp.close()
	# read the binary file of adjusted observed
	AdjObserved = {}
	fp = open(prefix+"_dist", 'rb')
	numtrans=struct.unpack('i', fp.read(4))[0]
	for i in range(numtrans):
		namelen=struct.unpack('i', fp.read(4))[0]
		seqlen=struct.unpack('i', fp.read(4))[0]
		name=""
		for j in range(namelen):
			name+=struct.unpack('c', fp.read(1))[0].decode('utf-8')
		counts = np.zeros(seqlen)
		for j in range(seqlen):
			counts[j]=struct.unpack('d', fp.read(8))[0]
		AdjObserved[AdjustmentList[i]] = counts
	fp.close()
	print("Finish reading LP adjusted distribution for {} transcripts.".format(len(AdjObserved)))
	return AdjObserved


def ConvertCoordinate(exp, Transcripts, sharedexons, t):
	if not Transcripts[t].Strand:
		for i  in range(1, len(sharedexons)):
			assert(sharedexons[i][0] < sharedexons[i-1][0])
	exons = Transcripts[t].Exons
	newexp = np.zeros(np.sum([e[1]-e[0] for e in sharedexons]))
	pos_ori = 0
	pos_new = 0
	idx_ori = 0
	idx_new = 0
	if Transcripts[t].Strand:
		for idx_new in range(len(sharedexons)):
			while idx_ori < len(exons) and exons[idx_ori][1] < sharedexons[idx_new][0]:
				pos_ori += exons[idx_ori][1] - exons[idx_ori][0]
				idx_ori += 1
			while idx_ori < len(exons) and max(exons[idx_ori][0], sharedexons[idx_new][0]) <= min(exons[idx_ori][1], sharedexons[idx_new][1]):
				overlap_start = max(exons[idx_ori][0], sharedexons[idx_new][0])
				overlap_end = min(exons[idx_ori][1], sharedexons[idx_new][1])
				newexp[(pos_new+overlap_start-sharedexons[idx_new][0]):(pos_new+overlap_end-sharedexons[idx_new][0])] = exp[(pos_ori+overlap_start-exons[idx_ori][0]):(pos_ori+overlap_end-exons[idx_ori][0])]
				pos_ori += exons[idx_ori][1] - exons[idx_ori][0]
				idx_ori += 1
			pos_new += sharedexons[idx_new][1] - sharedexons[idx_new][0]
	else:
		for idx_new in range(len(sharedexons)):
			while idx_ori < len(exons) and exons[idx_ori][0] > sharedexons[idx_new][1]:
				pos_ori += exons[idx_ori][1] - exons[idx_ori][0]
				idx_ori += 1
			while idx_ori < len(exons) and max(exons[idx_ori][0], sharedexons[idx_new][0]) <= min(exons[idx_ori][1], sharedexons[idx_new][1]):
				overlap_start = max(exons[idx_ori][0], sharedexons[idx_new][0])
				overlap_end = min(exons[idx_ori][1], sharedexons[idx_new][1])
				newexp[(pos_new+sharedexons[idx_new][1]-overlap_end):(pos_new+sharedexons[idx_new][1]-overlap_start)] = exp[(pos_ori+exons[idx_ori][1]-overlap_end):(pos_ori+exons[idx_ori][1]-overlap_start)]
				pos_ori += exons[idx_ori][1] - exons[idx_ori][0]
				idx_ori += 1
			pos_new += sharedexons[idx_new][1] - sharedexons[idx_new][0]
	return newexp


def PlotSingleDist(thisProb, tname="", _color=None, binsize=50, _ratio=1/3, _contain_x_title=True, _contain_y_title=True):
	P = np.zeros(int(len(thisProb)/binsize) + 1)
	for i in range(len(thisProb)):
		ind = int(i/binsize)
		P[ind] += thisProb[i]
	P /= np.sum(P)
	xaxis=np.zeros(len(P))
	for i in range(len(xaxis)):
		xaxis[i]=50*i
	d = {"pos":xaxis, "prob":P}
	df = pd.DataFrame(data=d)
	if _color is None:
		p = ggplot(aes(x='pos', y='prob'), data=df) + geom_point() + theme_light() + theme(aspect_ratio=_ratio) + ggtitle(tname)
	else:
		p = ggplot(aes(x='pos', y='prob'), data=df) + geom_point(color=_color) + theme_light() + theme(aspect_ratio=_ratio) + ggtitle(tname)
	if _contain_x_title:
		p += xlab("position in transcript")
	if _contain_y_title:
		p += ylab("coverage probability")
	return p


def PlotDist(thisExpectedProb, thisTrueRaw, tname="", binsize=50):
	assert(len(thisExpectedProb)==len(thisTrueRaw))
	P=np.zeros(int(len(thisExpectedProb)/binsize)+1)
	Q=np.zeros(int(len(thisExpectedProb)/binsize)+1)
	for i in range(len(thisExpectedProb)):
		ind=int(i/binsize)
		P[ind]+=thisExpectedProb[i]
		Q[ind]+=thisTrueRaw[i]
	P /= np.sum(P)
	Q /= np.sum(Q)
	xaxis=np.zeros(len(P))
	for i in range(len(xaxis)):
		xaxis[i]=binsize*i
	d = {"pos":np.hstack((xaxis,xaxis)), "prob":np.hstack((Q, P)), "Label":(["observed"]*len(P)+["theoretical"]*len(Q))}
	df = pd.DataFrame(data=d)
	p = ggplot(aes(x='pos', y='prob', color='Label'), data=df) + geom_point(alpha=0.5) + theme_light() + theme(aspect_ratio=1/3) + xlab("position in transcript") + ylab("coverage probability") + ggtitle(tname) + scale_color_manual(values=["#F8766D", "#00BFC4"])
	return p


def PlotSingle(dist, tname="", binsize=50, is_exp=True, legend_bottom=False, coordinate="transcript", removetick=True):
	length = len(dist)
	P = np.zeros(int(len(dist)/binsize)+1)
	for i in range(len(P)-1):
		a = i*binsize
		b = max((i+1)*binsize, len(P))
		P[i] = np.sum(dist[a:b]) + 1e-10
	P /= np.sum(P)
	xaxis=np.zeros(len(P))
	for i in range(len(xaxis)):
		xaxis[i]=binsize*i
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['mathtext.fontset'] = 'dejavuserif'
	if not (legend_bottom is None) or legend_bottom:
		fig, ax = plt.subplots(1, figsize=(8,3))
	else:
		fig, ax = plt.subplots(1, figsize=(10,3))
	if is_exp:
		ax.scatter(xaxis, P, color="#00BFC4", alpha=0.7, label="expected")
	else:
		ax.scatter(xaxis, P, color="#F8766D", alpha=0.7, label="observed")
	ax.grid(alpha=0.1)
	if legend_bottom:
		ax.legend(bbox_to_anchor=(0.15, -0.5, 0.7, 0.1), loc=3, ncol=1, mode="expand", borderaxespad=0., frameon=False)
	else:
		ax.legend(bbox_to_anchor=(1, 0.7), loc=2, mode="expand", borderaxespad=0., frameon=False)
	ax.set_xlabel("position in "+coordinate)
	ax.set_ylabel("coverage")
	ax.set_title(tname)
	ax.set_ylim(bottom=0)
	if not (legend_bottom is None) or legend_bottom:
		fig.subplots_adjust(left=0.1, right=0.975, bottom = 0.3)
	else:
		fig.subplots_adjust(left=0.1, right=0.8, bottom = 0.175)
	if removetick:
		ax.set_xticks([])
		ax.set_yticks([])
	return fig


def PlotDist_ribbon(exp, obs, cov, tname="", binsize=50, legend_bottom=False, coordinate="transcript", removetick=False):
	assert(len(exp) == len(obs))
	length = len(exp)
	P = np.zeros(int(len(exp)/binsize)+1)
	Q = np.zeros(int(len(exp)/binsize)+1)
	Std = np.zeros(int(len(exp)/binsize)+1)
	for i in range(len(P)-1):
		a = i*binsize
		b = max((i+1)*binsize, len(P))
		P[i] = np.sum(exp[a:b]) + 1e-10
		Q[i] = np.sum(obs[a:b]) + 1e-10
		a_effective = len(np.where(exp[:a] != 0)[0])
		b_effective = len(np.where(exp[:b] != 0)[0])
		var = 0
		if a_effective != b_effective:
			(frac_a, int_a) = np.modf(1.0 * a_effective / length * cov.shape[0])
			(frac_b, int_b) = np.modf(1.0 * b_effective / length * cov.shape[1])
			int_a = int(int_a)
			int_b = int(int_b)
			var = (1.0 - frac_a) * cov[int_a:(int_a+1), int_a:(int_a+1)] + (frac_b) * cov[int_b:(int_b+1), int_b:(int_b+1)]
			var += np.sum(cov[int_a:int_b, int_a:int_b])
		Std[i] = np.sqrt(var)
	P /= np.sum(P)
	Q /= np.sum(Q)
	xaxis=np.zeros(len(P))
	for i in range(len(xaxis)):
		xaxis[i]=binsize*i
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['mathtext.fontset'] = 'dejavuserif'
	if not (legend_bottom is None) or legend_bottom:
		fig, ax = plt.subplots(1, figsize=(8,3))
	else:
		fig, ax = plt.subplots(1, figsize=(10,3))
	ax.scatter(xaxis, P, color="#00BFC4", alpha=0.7, label="expected")
	ax.scatter(xaxis, Q, color="#F8766D", alpha=0.7, label="observed")
	ax.fill_between( xaxis, np.maximum(P-Std, 0), np.minimum(P+Std, 1), color="#00BFC4", alpha=0.1, label="std of Gaussian error" )
	ax.grid(alpha=0.1)
	if not (legend_bottom is None):
		if legend_bottom:
			ax.legend(bbox_to_anchor=(0.1, -0.5, 0.8, 0.1), loc=3, ncol=3, mode="expand", borderaxespad=0., frameon=False)
		else:
			ax.legend(bbox_to_anchor=(1, 0.7), loc=2, mode="expand", borderaxespad=0., frameon=False)
	ax.set_xlabel("position in "+coordinate)
	ax.set_ylabel("coverage")
	ax.set_title(tname)
	ax.set_ylim(bottom=0)
	if not (legend_bottom is None) or legend_bottom:
		fig.subplots_adjust(left=0.1, right=0.975, bottom = 0.3)
	else:
		fig.subplots_adjust(left=0.1, right=0.8, bottom = 0.15)
	if removetick:
		ax.set_xticks([])
		ax.set_yticks([])
	return fig


def PlotDist_ribbon_v2(ax, exp, obs, cov, tname="", binsize=50, legend_bottom=False, coordinate="transcript", removetick=False, ytop_lim=None):
	assert(len(exp) == len(obs))
	length = len(exp)
	P = np.zeros(int(len(exp)/binsize)+1)
	Q = np.zeros(int(len(exp)/binsize)+1)
	Std = np.zeros(int(len(exp)/binsize)+1)
	for i in range(len(P)-1):
		a = i*binsize
		b = max((i+1)*binsize, len(P))
		P[i] = np.sum(exp[a:b]) + 1e-10
		Q[i] = np.sum(obs[a:b]) + 1e-10
		a_effective = len(np.where(exp[:a] != 0)[0])
		b_effective = len(np.where(exp[:b] != 0)[0])
		var = 0
		if a_effective != b_effective:
			(frac_a, int_a) = np.modf(1.0 * a_effective / length * cov.shape[0])
			(frac_b, int_b) = np.modf(1.0 * b_effective / length * cov.shape[1])
			int_a = int(int_a)
			int_b = int(int_b)
			var = (1.0 - frac_a) * cov[int_a:(int_a+1), int_a:(int_a+1)] + (frac_b) * cov[int_b:(int_b+1), int_b:(int_b+1)]
			var += np.sum(cov[int_a:int_b, int_a:int_b])
		Std[i] = np.sqrt(var)
	P /= np.sum(P)
	Q /= np.sum(Q)
	xaxis=np.zeros(len(P))
	for i in range(len(xaxis)):
		xaxis[i]=binsize*i
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['mathtext.fontset'] = 'dejavuserif'
	ax.scatter(xaxis, P, color="#00BFC4", alpha=0.7, label="expected")
	ax.scatter(xaxis, Q, color="#F8766D", alpha=0.7, label="observed")
	ax.fill_between( xaxis, np.maximum(P-Std, 0), np.minimum(P+Std, 1), color="#00BFC4", alpha=0.1, label="std of Gaussian error" )
	ax.grid(alpha=0.1)
	if not (legend_bottom is None):
		if legend_bottom:
			print("legend bottom")
			ax.legend(bbox_to_anchor=(0.1, -0.4, 0.8, 0.1), loc=3, ncol=3, mode="expand", borderaxespad=0., frameon=False)
		else:
			ax.legend(bbox_to_anchor=(1, 0.7), loc=2, mode="expand", borderaxespad=0., frameon=False)
	ax.set_xlabel("position in "+coordinate)
	ax.set_ylabel("coverage")
	ax.set_title(tname)
	ax.set_ylim(bottom=0)
	if not (ytop_lim is None):
		ax.set_ylim(top=ytop_lim)
	if removetick:
		ax.set_xticks([])
		ax.set_yticks([])
	return ax


def PlotDist2(thisExpectedProb1, thisTrueRaw1, thisExpectedProb2, thisTrueRaw2, tnames, binsize=50):
	assert(len(thisExpectedProb1)==len(thisTrueRaw1) and len(thisExpectedProb2)==len(thisTrueRaw2))
	P1=np.zeros(int(len(thisExpectedProb1)/binsize)+1)
	Q1=np.zeros(int(len(thisExpectedProb1)/binsize)+1)
	for i in range(len(thisExpectedProb1)):
		ind=int(i/binsize)
		P1[ind]+=thisExpectedProb1[i]
		Q1[ind]+=thisTrueRaw1[i]
	P1 /= np.sum(P1)
	Q1 /= np.sum(Q1)
	# transcript 2
	P2=np.zeros(int(len(thisExpectedProb2)/binsize)+1)
	Q2=np.zeros(int(len(thisExpectedProb2)/binsize)+1)
	for i in range(len(thisExpectedProb2)):
		ind=int(i/binsize)
		P2[ind]+=thisExpectedProb2[i]
		Q2[ind]+=thisTrueRaw2[i]
	P2 /= np.sum(P2)
	Q2 /= np.sum(Q2)
	# make axis
	xaxis1=np.zeros(len(P1))
	for i in range(len(xaxis1)):
		xaxis1[i]=binsize*i
	xaxis2=np.zeros(len(P2))
	for i in range(len(xaxis2)):
		xaxis2[i]=binsize*i
	# panda table
	d = {"pos":np.hstack((xaxis1,xaxis1, xaxis2, xaxis2)), "prob":np.hstack((Q1, P1, Q2, P2)), \
		"Label":(["observed"]*len(P1)+["theoretical"]*len(Q1)+["observed"]*len(P2)+["theoretical"]*len(Q2)), \
		"Transcript":([tnames[0]]*2*len(P1)+[tnames[1]]*2*len(P1))}
	df = pd.DataFrame(data=d)
	p = ggplot(aes(x='pos', y='prob', color='Label'), data=df) + geom_point(alpha=0.5) + theme_light() + theme(aspect_ratio=1/3) + xlab("position in gene") + ylab("coverage probability") + facet_wrap("Transcript", nrow=2) + scale_color_manual(values=["#F8766D", "#00BFC4"])
	return p