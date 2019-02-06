/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __Transcript_H__
#define __Transcript_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
#include <mutex>
#include "boost/algorithm/string.hpp"

using namespace std;


class Exon_t
{
public:
	int32_t StartPos;
	int32_t EndPos;

public:
	Exon_t(){};
	Exon_t(int32_t StartPos, int32_t EndPos): StartPos(StartPos), EndPos(EndPos) {};

	bool operator < (const Exon_t& rhs) const {
		if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else
			return EndPos < rhs.EndPos;
	};

	bool operator == (const Exon_t& rhs) const {
		return (StartPos == rhs.StartPos) && (EndPos == rhs.EndPos);
	};
};

class Transcript_t
{
public:
	string TransID;
	string GeneID;
	string Chr;
	int32_t StartPos;
	int32_t EndPos;
	bool Strand;
	vector<Exon_t> Exons;

public:
	Transcript_t(){};
	Transcript_t(string TransID, string GeneID, string Chr, int32_t StartPos, int32_t EndPos, bool Strand):
		TransID(TransID), GeneID(GeneID), Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand)
	{
		Exons.clear();
	};

	bool operator < (const Transcript_t& rhs) const {
		if (Chr != rhs.Chr)
			return Chr < rhs.Chr;
		else if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else if (EndPos != rhs.EndPos)
			return EndPos < rhs.EndPos;
		else
			return Strand < rhs.Strand;
	};

	static bool CompTransID(const Transcript_t& lhs, const Transcript_t& rhs) {
		return lhs.TransID < rhs.TransID;
	};

	void AddExon (int32_t StartPos, int32_t EndPos) {
		Exon_t tmp(StartPos, EndPos);
		Exons.push_back(tmp);
	};

	void AddExon (Exon_t exon){
		Exons.push_back(exon);
	};

	void SortExons() {
		assert(Exons.size() > 0);
		sort(Exons.begin(), Exons.end());
		if (!Strand)
			reverse(Exons.begin(), Exons.end());
	};
};

// small function to extract tag information in GTF file
string GetFeature(string line, string key);

// read transcript info from GTF file
void ReadGTF (string GTFfile, vector<Transcript_t>& Transcripts);

// trim transcripts: salmon mark some transcripts as duplicated and remove them from quantification
// We also remove them from any later steps according to what salmon removed
void TrimTranscripts (const map<string,double>& SalmonExp, vector<Transcript_t>& Transcripts);

// get GeneID and TransID mapping
void Map_Gene_Trans (const vector<Transcript_t>& Transcripts, map<string,int32_t>& TransIndex, 
	map<string,string>& TransGeneMap, map< string,vector<string> >& GeneTransMap);

// collect non-redundant unique exonic positions for each genes
void CombineGeneLevelExons (const vector<Transcript_t>& Transcripts, const map<string,int32_t>& TransIndex, 
	const map< string,vector<string> >& GeneTransMap, map<string, vector<Exon_t> >& GeneLevelExons);

#endif
