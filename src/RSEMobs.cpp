/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
#include <tuple>
#include "boost/algorithm/string.hpp"
#include "htslib/sam.h"

using namespace std;

// about transcript
class Transcript_t
{
public:
	string TransName;
	string Chr;
	int32_t StartPos;
	int32_t EndPos;
	bool Strand;
	vector< pair<int32_t,int32_t> > Exons;
public:
	Transcript_t(){};
	Transcript_t(string TransName, string Chr, int32_t StartPos, int32_t EndPos, bool Strand):
		TransName(TransName), Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand) {};
	void SortExons(){
		if(Strand)
			sort(Exons.begin(), Exons.end(), [](pair<int32_t,int32_t> a, pair<int32_t,int32_t> b){return a.first<b.first;} );
		else
			sort(Exons.begin(), Exons.end(), [](pair<int32_t,int32_t> a, pair<int32_t,int32_t> b){return a.first>b.first;} );
	};
};


string GetFeature(string line, string key)
{
	size_t s=line.find(key);
	if(s==string::npos)
		return "";
	size_t t=line.find_first_of(";", s+1);
	if(t==string::npos)
		return "";
	return line.substr(s+key.size()+2, t-s-key.size()-3);
};


void ReadGTF(string GTFfile, map<string, Transcript_t>& Transcripts)
{
	Transcripts.clear();

	// first pass collect all transcripts
	ifstream input(GTFfile);
	string line;
	while(getline(input, line)){
		if(line[0]=='#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if (strs[2] == "transcript") {
			string transname = GetFeature(line, "transcript_id");
			Transcript_t tmp(transname, strs[0], stoi(strs[3])-1, stoi(strs[4]), (strs[6]=="+"));
			Transcripts[transname] = tmp;
		}
	}
	input.close();

	// second pass collect exon of each transcript
	input.open(GTFfile, std::ifstream::in);
	map<string, Transcript_t>::iterator it = Transcripts.begin();
	while(getline(input, line)) {
		if(line[0]=='#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if (strs[2]=="exon"){
			string transname=GetFeature(line, "transcript_id");
			if ((it->second).TransName != transname)
				it = Transcripts.find(transname);
			if (it == Transcripts.end()) {
				cout << "Warning: exons contain unrecognized transcript id " << transname << ", ignoring this exon.\n";
				continue;
			}
			assert(strs[0]==(it->second).Chr);
			assert((strs[6]=="+") == (it->second).Strand);
			assert((it->second).TransName==transname);
			(it->second).Exons.push_back(make_pair(stoi(strs[3])-1, stoi(strs[4])));
		}
	}
	input.close();

	// sort exons
	for (map<string, Transcript_t>::iterator it = Transcripts.begin(); it != Transcripts.end(); it++)
		(it->second).SortExons();
};


struct Brdy{
	int32_t TID;
	int32_t Pos;
	bool IsEnd;
	double Weight;

	Brdy(){};
	Brdy(int32_t TID, int32_t Pos, bool IsEnd, double Weight): TID(TID), Pos(Pos), IsEnd(IsEnd), Weight(Weight){};

	bool operator < (const Brdy& rhs) const {
		if(TID!=rhs.TID)
			return TID<rhs.TID;
		else if(Pos!=rhs.Pos)
			return Pos<rhs.Pos;
		else
			return IsEnd<rhs.IsEnd;
	};
};


bool SpanSplicing(int32_t StartPos, int32_t EndPos, const Transcript_t& t, int32_t threshold = 5)
{
	// loop over the exons, and check whether read/mate align blocks hit the junction
	int32_t coveredlen = 0;
	for (uint32_t i = 0; i < t.Exons.size() - 1; i++) {
		// junctions are considered as the end of the exon
		coveredlen += t.Exons[i].second - t.Exons[i].first;
		if (StartPos < coveredlen - threshold && EndPos > coveredlen + threshold)
			return true;
	}
	return false;
};


int32_t ReadBAMStartPos(string bamfile, const map<string, Transcript_t>& transcripts, vector<string>& TransNames, vector<int32_t>& TransLength, vector<Brdy>& ReadBrdy, vector<Brdy>& JunctionBrdy)
{
	// initialize read alignment start vector
	ReadBrdy.clear();
	TransNames.clear();
	TransLength.clear();

	samFile * bamreader=sam_open(bamfile.c_str(), "r");
	bam_hdr_t *header=NULL;
	bam1_t *b=bam_init1();

	// open bam file
	if(bamreader==NULL){
		cout<<"cannot open bam file\n";
		return -1;
	}
	if((header=sam_hdr_read(bamreader))==0){
		cout<<"cannot open header\n";
		return -1;
	}

	// get header info
	for(int32_t i = 0; i < header->n_targets; i++) {
		int32_t len = (header->target_len)[i];
		string name = (header->target_name)[i];
		TransNames.push_back(name);
		TransLength.push_back(len);
	}

	// map of (readname, tid, pos) to record whether the alignment has been processed
	map< tuple<string,int32_t,int32_t,int32_t>, bool > alignname;
	map< tuple<string,int32_t,int32_t,int32_t>, bool > alignname_junction;
	// read bam record
	while(sam_read1(bamreader, header, b)>0){
		if((b)->core.flag&BAM_FUNMAP)
			continue;
		// only retrieve read alignment if the (readname, tid) pair has not been processed before
		string readname = bam_get_qname(b);
		int32_t tid = b->core.tid;
		int32_t fragstart = min(b->core.pos, b->core.mpos);
		int32_t matestart = max(b->core.pos, b->core.mpos);
		// try to find the alignment info in alignname
		map< tuple<string,int32_t,int32_t,int32_t>, bool >::iterator itmap = alignname.find( make_tuple(readname, tid, fragstart, matestart) );
		if (itmap == alignname.cend()) {
			// add to Brdy
			uint8_t* s = bam_aux_get(b, "ZW");
			double w = bam_aux2f(s);
			Brdy tmpstart(tid, fragstart, false, w);
			assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
			ReadBrdy.push_back(tmpstart);
			// add to alignname map
			alignname[make_tuple(readname, tid, fragstart, matestart)] = true;
		}

		// try to see whether the alignment spanning junction
		string transname( header->target_name[tid] );
		map<string, Transcript_t>::const_iterator ittrans = transcripts.find(transname);
		assert(ittrans != transcripts.cend());
		if (SpanSplicing(fragstart, bam_endpos(b), ittrans->second)) {
			// check whether it has been recorded in the JunctionBrdy
			itmap = alignname_junction.find( make_tuple(readname, tid, fragstart, matestart) );
			if (itmap == alignname_junction.cend()) {
				// add to Brdy
				uint8_t* s = bam_aux_get(b, "ZW");
				double w = bam_aux2f(s);
				Brdy tmpstart(tid, fragstart, false, w);
				JunctionBrdy.push_back(tmpstart);
				alignname_junction[make_tuple(readname, tid, fragstart, matestart)] = true;
			}
		}

	}

	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(bamreader);

	// calculate sum of weight
	double sum_weight = 0;
	for (vector<Brdy>::iterator it = ReadBrdy.begin(); it != ReadBrdy.end(); it++)
		sum_weight += it->Weight;
	cout << "sum of weight = " << sum_weight << endl;


	ReadBrdy.shrink_to_fit();
	sort(ReadBrdy.begin(), ReadBrdy.end());

	JunctionBrdy.shrink_to_fit();
	sort(JunctionBrdy.begin(), JunctionBrdy.end());

	return 0;
};


int32_t WriteReadBrdy(const vector<string>& TransNames, const vector<int32_t>& TransLength, const vector<Brdy>& ReadBrdy, ofstream& ss)
{
	assert(is_sorted(ReadBrdy.cbegin(), ReadBrdy.cend()));
	
	int32_t numtrans=1;
	for(vector<Brdy>::const_iterator it=ReadBrdy.cbegin(); it!=ReadBrdy.cend(); it++)
		if(it!=ReadBrdy.cbegin() && it->TID!=(it-1)->TID)
			numtrans++;
	ss.write((char*)(&numtrans), sizeof(int32_t));

	vector<int32_t> poses;
	vector<double> counts;
	for(vector<Brdy>::const_iterator it=ReadBrdy.cbegin(); it!=ReadBrdy.cend(); it++){
		if(it!=ReadBrdy.cbegin() && it->TID!=(it-1)->TID){
			if(poses.back()!=TransLength[(it-1)->TID]-1){
				poses.push_back(TransLength[(it-1)->TID]-1);
				counts.push_back(0);
			}

			string name=TransNames[(it-1)->TID];
			int32_t namelen=name.size();
			assert(counts.size()==poses.size());
			int32_t vectorlen=poses.size();
			ss.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
			ss.write(reinterpret_cast<char*>(&vectorlen), sizeof(int32_t));
			ss.write(name.c_str(), namelen*sizeof(char));
			ss.write(reinterpret_cast<char*>(poses.data()), vectorlen*sizeof(int32_t));
			ss.write(reinterpret_cast<char*>(counts.data()), vectorlen*sizeof(double));

			poses.clear();
			counts.clear();
		}
		if(poses.size()==0 || it->Pos!=poses.back()){
			poses.push_back(it->Pos);
			counts.push_back(it->Weight);
		}
		else
			counts.back()+=it->Weight;
	}
	vector<Brdy>::const_iterator it=ReadBrdy.cend(); it--;
	if(poses.back()!=TransLength[it->TID]-1){
		poses.push_back(TransLength[it->TID]-1);
		counts.push_back(0);
	}
	string name=TransNames[it->TID];
	int32_t namelen=name.size();
	assert(counts.size()==poses.size());
	int32_t vectorlen=poses.size();
	ss.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
	ss.write(reinterpret_cast<char*>(&vectorlen), sizeof(int32_t));
	ss.write(name.c_str(), namelen*sizeof(char));
	ss.write(reinterpret_cast<char*>(poses.data()), vectorlen*sizeof(int32_t));
	ss.write(reinterpret_cast<char*>(counts.data()), vectorlen*sizeof(double));

	return 0;
};


int32_t main(int32_t argc, char* argv[]){
	if(argc==1){
		printf("rsemobs <gtffile> <bamfile> <outfile>\n");
	}
	else{
		string GTFfile(argv[1]);
		string BamFile(argv[2]);
		string OutFile(argv[3]);

		string OutJunctionFile;
		size_t pos_name = OutFile.find_last_of("/");
		size_t suffixpos = OutFile.find_last_of(".");
		if(suffixpos != string::npos && suffixpos > pos_name)
			OutJunctionFile = OutFile.substr(0, suffixpos) + "_junction"+OutFile.substr(suffixpos);
		else
			OutJunctionFile = OutFile+"_junction";

		map<string, Transcript_t> transcripts;
		ReadGTF(GTFfile, transcripts);

		vector<string> TransNames;
		vector<int32_t> TransLength;
		vector<Brdy> ReadBrdy;
		vector<Brdy> JunctionBrdy;

		ReadBAMStartPos(BamFile, transcripts, TransNames, TransLength, ReadBrdy, JunctionBrdy);
		ofstream ss(OutFile, ios::out | ios::binary);
		WriteReadBrdy(TransNames, TransLength, ReadBrdy, ss);

		ofstream ss_junction(OutJunctionFile, ios::out | ios::binary);
		WriteReadBrdy(TransNames, TransLength, JunctionBrdy, ss_junction);
	}
}

