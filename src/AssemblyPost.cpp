#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <map>
#include "boost/algorithm/string.hpp"

using namespace std;

struct Exon_t{
	string Chr;
	int32_t StartPos;
	int32_t EndPos;

	Exon_t(){};
	Exon_t(string Chr, int32_t StartPos, int32_t EndPos): Chr(Chr), StartPos(StartPos), EndPos(EndPos) {};

	bool operator < (const Exon_t& rhs) const {
		if (Chr != rhs.Chr)
			return Chr < rhs.Chr;
		else if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else
			return EndPos < rhs.EndPos;
	};
};


struct Junction_t{
	string Chr;
	int32_t Pos1;
	int32_t Pos2;

	Junction_t(){};
	Junction_t(string Chr, int32_t Pos1, int32_t Pos2): Chr(Chr), Pos1(Pos1), Pos2(Pos2) {};

	bool operator < (const Junction_t& rhs) const {
		if (Chr != rhs.Chr)
			return Chr < rhs.Chr;
		else if (Pos1 != rhs.Pos1)
			return Pos1 < rhs.Pos1;
		else
			return Pos2 < rhs.Pos2;
	};

	bool operator == (const Junction_t& rhs) const {
		return (Chr==rhs.Chr && Pos1==rhs.Pos1 && Pos2==rhs.Pos2);
	};
};


class Transcript_t{
public:
	string GeneID;
	string TransID;
	string Chr;
	int32_t StartPos;
	int32_t EndPos;
	bool Strand;
	vector<Exon_t> Exons;

public:
	Transcript_t(){};
	Transcript_t(string GeneID, string TransID, string Chr, int32_t StartPos, int32_t EndPos, bool Strand):
		GeneID(GeneID), TransID(TransID), Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand) {};

	bool operator < (const Transcript_t& rhs) const {
		if (Chr != rhs.Chr)
			return Chr < rhs.Chr;
		else if (StartPos != rhs.StartPos)
			return StartPos < rhs.StartPos;
		else
			return EndPos < rhs.EndPos;
	};

	static bool CompGeneID (const Transcript_t& lhs, const Transcript_t& rhs) {
		return lhs.GeneID < rhs.GeneID;
	};

	static bool CompTransID (const Transcript_t& lhs, const Transcript_t& rhs) {
		return lhs.TransID < rhs.TransID;
	};

	void SortExons(){
		sort(Exons.begin(), Exons.end());
		if (!Strand)
			reverse(Exons.begin(), Exons.end());
	};
};


string GetFeature(string line, string key)
{
	size_t s = line.find(key);
	s += key.size() + 2;
	if (s == string::npos)
		cout << "Error: key \"" << key << "\" not find in string \"" << line << "\"" << endl;
	assert(s != string::npos);
	size_t t = line.find_first_of(";", s+1);

	return line.substr(s, t-s-1);
};


void ReadGTF(string gtffile, vector<Transcript_t>& Transcripts)
{
	Transcripts.clear();

	vector< pair<string, Exon_t> > extraexons;

	ifstream input(gtffile);
	string line;
	string prevtransid;
	while (getline(input, line)){
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if (strs[2] == "transcript"){
			string transid = GetFeature(line, "transcript_id");
			string geneid = GetFeature(line, "gene_id");
			Transcript_t tmp(geneid, transid, strs[0], stoi(strs[3])-1, stoi(strs[4]), (strs[6]=="+"));
			Transcripts.push_back(tmp);
			// add to previous info
			prevtransid = transid;
		}
		else if (strs[2] == "exon"){
			string transid = GetFeature(line, "transcript_id");
			Exon_t tmpexon(strs[0], stoi(strs[3])-1, stoi(strs[4]));
			if (transid == prevtransid)
				Transcripts.back().Exons.push_back(tmpexon);
			else
				extraexons.push_back( make_pair(transid, tmpexon) );
		}
	}
	input.close();

	map<string,int32_t> AllTransID;
	for (int32_t i=0; i<Transcripts.size(); i++)
		AllTransID[Transcripts[i].TransID] = i;

	for (vector< pair<string, Exon_t> >::iterator itextra = extraexons.begin(); itextra != extraexons.end(); itextra++) {
		map<string,int32_t>::iterator itmap = AllTransID.find(itextra->first);
		assert(itmap != AllTransID.end());
		int32_t ind = itmap->second;
		assert(Transcripts[ind].TransID == itextra->first);
		Transcripts[ind].Exons.push_back(itextra->second);
	}

	for (vector<Transcript_t>::iterator it = Transcripts.begin(); it != Transcripts.end(); it++)
		it->SortExons();
};


void CalculateTransLength(const vector<Transcript_t>& Transcripts, map<string,int32_t>& TransLength)
{
	TransLength.clear();
	for(vector<Transcript_t>::const_iterator it = Transcripts.cbegin(); it != Transcripts.cend(); it++){
		int32_t sumexonlen = 0;
		for(vector<Exon_t>::const_iterator itexon = it->Exons.cbegin(); itexon != it->Exons.cend(); itexon++)
			sumexonlen += itexon->EndPos - itexon->StartPos;
		TransLength[it->TransID] = sumexonlen;
	}
};


void CollectJunction(const vector<Transcript_t>& Transcripts, vector<Junction_t>& Junctions)
{
	Junctions.clear();
	for (vector<Transcript_t>::const_iterator it = Transcripts.cbegin(); it != Transcripts.cend(); it++){
		if (it->Strand){
			for (int32_t i=1; i < it->Exons.size(); i++){
				Junction_t tmp(it->Chr, (it->Exons[i-1]).EndPos, (it->Exons[i]).StartPos);
				Junctions.push_back(tmp);
			}
		}
		else{
			for (int32_t i=1; i < it->Exons.size(); i++){
				Junction_t tmp(it->Chr, (it->Exons[i]).EndPos, (it->Exons[i-1]).StartPos);
				Junctions.push_back(tmp);
			}
		}
	}

	sort(Junctions.begin(), Junctions.end());
	vector<Junction_t>::iterator itend = unique(Junctions.begin(), Junctions.end());
	Junctions.resize(distance(Junctions.begin(), itend));
	Junctions.reserve(Junctions.size());
};


void ProcessTmapFile(string tmapfile, vector<Transcript_t>& AsmTranscripts, const map<string,int32_t>& TransLength, map<string,string>& Correspondence, map<string,double>& Cov)
{
	Correspondence.clear();
	Cov.clear();

	ifstream input(tmapfile);
	string line;
	int32_t linecount = 0;
	while (getline(input, line)) {
		linecount++;
		if (linecount == 1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if (strs[2]!="u" && strs[2]!="r" && strs[2]!="=" && strs[2]!="s" && strs[2]!="i" && strs[1]!="-"){
			Correspondence[strs[4]] = strs[1];
			Cov[strs[4]] = stod(strs[7]);
		}
		else if(strs[2]=="="){
			map<string,int32_t>::const_iterator itmap = TransLength.find(strs[1]);
			assert(itmap != TransLength.cend());
			int32_t len_original = itmap->second;
			int32_t len_assemble = stoi(strs[9]);
			if (abs(len_original - len_assemble) > 200) {
			//if(!(len_original > 0.5*len_assemble && len_assemble > 0.5*len_original)){
				Correspondence[strs[4]] = strs[1];
				Cov[strs[4]] = stod(strs[8]);
			}
		}
	}
	input.close();

	vector<Transcript_t> newAsmTranscripts;
	for (vector<Transcript_t>::iterator it=AsmTranscripts.begin(); it != AsmTranscripts.end(); it++){
		map<string,string>::iterator itmap = Correspondence.find(it->TransID);
		if (itmap != Correspondence.end())
			newAsmTranscripts.push_back((*it));
	}
	newAsmTranscripts.reserve(newAsmTranscripts.size());
	AsmTranscripts = newAsmTranscripts;
};


void UpdateAssemblyStrandInfo(vector<Transcript_t>& AsmTranscripts, const vector<Transcript_t>& RefTranscripts, const map<string,string>& Correspondence)
{
	map<string,int32_t> RefAllTransID;
	for (int32_t i=0; i<RefTranscripts.size(); i++)
		RefAllTransID[RefTranscripts[i].TransID] = i;

	for (vector<Transcript_t>::iterator it = AsmTranscripts.begin(); it != AsmTranscripts.end(); it++){
		map<string,string>::const_iterator itmap1 = Correspondence.find(it->TransID);
		assert(itmap1 != Correspondence.cend());
		map<string,int32_t>::const_iterator itmap2 = RefAllTransID.find(itmap1->second);
		assert(itmap2 != RefAllTransID.cend());

		int32_t ind = itmap2->second;
		it->Strand = RefTranscripts[ind].Strand;
	}
};


void LabelingMultiExon(const vector<Transcript_t>& AsmTranscripts, map<string,bool>& MultiExonLabel)
{
	MultiExonLabel.clear();
	for (vector<Transcript_t>::const_iterator it = AsmTranscripts.begin(); it != AsmTranscripts.end(); it++){
		assert( it->Exons.size() > 0);
		MultiExonLabel[it->TransID] = (it->Exons.size()>1);
	}
};


void LabelingNovelJunction(const vector<Transcript_t>& AsmTranscripts, const vector<Junction_t>& Junctions, map<string,bool>& NovelJunctionLabel)
{
	int32_t thresh = 5;
	NovelJunctionLabel.clear();
	for (vector<Transcript_t>::const_iterator it = AsmTranscripts.begin(); it != AsmTranscripts.end(); it++){
		if (it->Exons.size() <= 1)
			NovelJunctionLabel[it->TransID] = false;
		else{
			vector<Junction_t> tmpjunctions;
			if (it->Strand){
				for (int32_t i=1; i < it->Exons.size(); i++){
					Junction_t tmp(it->Chr, (it->Exons[i-1]).EndPos, (it->Exons[i]).StartPos);
					tmpjunctions.push_back(tmp);
				}
			}
			else{
				for (int32_t i=1; i < it->Exons.size(); i++){
					Junction_t tmp(it->Chr, (it->Exons[i]).EndPos, (it->Exons[i-1]).StartPos);
					tmpjunctions.push_back(tmp);
				}
			}
			// check each junction in this assembled transcript, to see whether the junction is novel compared to all reference junctions
			bool hasnovel = false;
			for (int32_t i=0; i<tmpjunctions.size(); i++){
				vector<Junction_t>::const_iterator lb = lower_bound(Junctions.cbegin(), Junctions.cend(), tmpjunctions[i]);
				bool ishit = false;
				for (vector<Junction_t>::const_iterator itjunc = lb; itjunc != Junctions.cend() && itjunc->Chr == tmpjunctions[i].Chr && 
					itjunc->Pos1 < tmpjunctions[i].Pos1+thresh; itjunc++)
				{
					if (itjunc->Chr == tmpjunctions[i].Chr && abs(itjunc->Pos1-tmpjunctions[i].Pos1) + abs(itjunc->Pos2-tmpjunctions[i].Pos2) < thresh){
						ishit = true;
						break;
					}
				}
				for (vector<Junction_t>::const_iterator itjunc = lb; itjunc != Junctions.cbegin()-1 && itjunc->Chr == tmpjunctions[i].Chr &&
					itjunc->Pos1 > tmpjunctions[i].Pos1-thresh; itjunc--)
				{
					if (itjunc->Chr == tmpjunctions[i].Chr && abs(itjunc->Pos1-tmpjunctions[i].Pos1) + abs(itjunc->Pos2-tmpjunctions[i].Pos2) < thresh){
						ishit = true;
						break;
					}
				}
				if (!ishit){
					hasnovel = true;
					break;
				}
			}

			NovelJunctionLabel[it->TransID] = hasnovel;
		}
	}
};


void WriteRanking(const map<string,string>& Correspondence, const map<string,double>& Cov, const map<string,bool>& MultiExonLabel, 
	const map<string,bool>& NovelJunctionLabel, string outputfile, bool onlymultiexon=false, bool onlynoveljunc=false)
{
	vector< pair<string,double> > VectorCov;
	for (map<string,double>::const_iterator itmap = Cov.cbegin(); itmap != Cov.cend(); itmap++)
		VectorCov.push_back( make_pair(itmap->first, itmap->second) );
	sort(VectorCov.begin(), VectorCov.end(), [](pair<string,double> a, pair<string,double> b){return a.second > b.second;} );

	ofstream output(outputfile);
	output << "# RefTransID\tAssemblyCov\tAssemblyTransID\n";
	for (vector< pair<string,double> >::iterator it = VectorCov.begin(); it != VectorCov.end(); it++){
		map<string,string>::const_iterator itcorres = Correspondence.find(it->first);
		map<string,bool>::const_iterator itmulti = MultiExonLabel.find(it->first);
		map<string,bool>::const_iterator itnovel = NovelJunctionLabel.find(it->first);
		if (onlymultiexon && (!itmulti->second || itnovel->second))
			continue;
		if (onlynoveljunc && !itnovel->second)
			continue;
		output << (itcorres->second) << "\t" << (it->second) <<"\t" << (it->first) << endl;
	}
};


int32_t main(int32_t argc, char* argv[])
{
	if (argc == 1)
		printf("assemblypost <RefGTFfile> <AsmGTFfile> <TmapFile> <OutputPrefix>\n");
	else{
		string RefGTFfile = string(argv[1]);
		string AsmGTFfile = string(argv[2]);
		string TmapFile = string(argv[3]);
		string OutputPrefix = string(argv[4]);

		vector<Transcript_t> RefTranscripts;
		vector<Transcript_t> AsmTranscripts;
		ReadGTF(RefGTFfile, RefTranscripts);
		ReadGTF(AsmGTFfile, AsmTranscripts);

		map<string,int32_t> RefTransLength;
		CalculateTransLength(RefTranscripts, RefTransLength);

		vector<Junction_t> RefJunctions;
		map<string,string> Correspondence;
		map<string,double> Cov;
		map<string,bool> MultiExonLabel;
		map<string,bool> NovelJunctionLabel;

		CollectJunction(RefTranscripts, RefJunctions);
		ProcessTmapFile(TmapFile, AsmTranscripts, RefTransLength, Correspondence, Cov);
		UpdateAssemblyStrandInfo(AsmTranscripts, RefTranscripts, Correspondence);
		LabelingMultiExon(AsmTranscripts, MultiExonLabel);
		LabelingNovelJunction(AsmTranscripts, RefJunctions, NovelJunctionLabel);

		WriteRanking(Correspondence, Cov, MultiExonLabel, NovelJunctionLabel, OutputPrefix+"_overall_novelisoforms.txt");
		WriteRanking(Correspondence, Cov, MultiExonLabel, NovelJunctionLabel, OutputPrefix+"_multiexon_novelisoforms.txt", true, false);
		WriteRanking(Correspondence, Cov, MultiExonLabel, NovelJunctionLabel, OutputPrefix+"_noveljunc_novelisoforms.txt", false, true);
	}
};
