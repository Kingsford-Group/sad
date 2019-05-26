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


int32_t ReadBAMStartPos(string bamfile, vector<string>& TransNames, vector<int32_t>& TransLength, vector<Brdy>& ReadBrdy)
{
	// initialize read alignment start vector
	ReadBrdy.clear();
	ReadBrdy.reserve(65536);
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
			if(ReadBrdy.capacity()==ReadBrdy.size())
				ReadBrdy.reserve(ReadBrdy.size()*2);
			// add to alignname map
			alignname[make_tuple(readname, tid, fragstart, matestart)] = true;
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


	ReadBrdy.reserve(ReadBrdy.size());
	sort(ReadBrdy.begin(), ReadBrdy.end());

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
		printf("rsemobs <bamfile> <outfile>\n");
	}
	else{
		string BamFile(argv[1]);
		string OutFile(argv[2]);

		vector<string> TransNames;
		vector<int32_t> TransLength;
		vector<Brdy> ReadBrdy;

		ReadBAMStartPos(BamFile, TransNames, TransLength, ReadBrdy);
		ofstream ss(OutFile, ios::out | ios::binary);
		WriteReadBrdy(TransNames, TransLength, ReadBrdy, ss);
	}
}

