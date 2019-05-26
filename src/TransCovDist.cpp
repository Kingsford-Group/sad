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
#include "boost/algorithm/string.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/iostreams/filter/gzip.hpp"
#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/copy.hpp"
#include "boost/math/distributions/binomial.hpp"
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


class Alignment_t{
public:
	int32_t TID;
	int32_t StartPos;
	int32_t	EndPos;
	int32_t mate_StartPos;
	int32_t mate_EndPos;

public:
	Alignment_t(){};
	Alignment_t(int32_t TID, int32_t StartPos, int32_t EndPos, int32_t mate_StartPos):
		TID(TID), StartPos(StartPos), EndPos(EndPos), mate_StartPos(mate_StartPos) {mate_EndPos = -1;};
};


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


vector<Alignment_t> MergeAlignments(const vector<Alignment_t>& alignments)
{
	vector<Alignment_t> merged;
	for (vector<Alignment_t>::const_iterator it = alignments.cbegin(); it != alignments.cend(); it++) {
		bool has_mate_inserted = false;
		for (vector<Alignment_t>::iterator itmate = merged.begin(); itmate != merged.end(); itmate++) {
			if (itmate->TID == it->TID && itmate->mate_StartPos == it->StartPos) {
				itmate->mate_EndPos = it->EndPos;
				has_mate_inserted = true;
				break;
			}
		}
		if (!has_mate_inserted)
			merged.push_back(*it);
	}
	return merged;
};


bool SpanSplicing(const Alignment_t& ali, const Transcript_t& t, int32_t threshold = 5)
{
	// loop over the exons, and check whether read/mate align blocks hit the junction
	int32_t coveredlen = 0;
	for (int32_t i = 0; i < t.Exons.size(); i++) {
		// junctions are considered as the end of the exon
		coveredlen += t.Exons[i].second - t.Exons[i].first;
		if (ali.StartPos < coveredlen - threshold && ali.EndPos > coveredlen + threshold)
			return true;
		else if (ali.mate_StartPos < coveredlen - threshold && ali.mate_EndPos > coveredlen + threshold)
			return true;
	}
	return false;
};


void ReadLowQualReadNames(string file1name, string file2name, vector<string>& LowQualReadNames, bool Phred33=true, int32_t LowThresh=3, double LenPropThresh=0.5)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Start reading fastq.gz file."<<endl;

	LowQualReadNames.clear();

	std::ifstream file1(file1name, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf1;
    inbuf1.push(boost::iostreams::gzip_decompressor());
    inbuf1.push(file1);
    std::istream input1(&inbuf1); //Convert streambuf to istream

    std::ifstream file2(file2name, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf2;
    inbuf2.push(boost::iostreams::gzip_decompressor());
    inbuf2.push(file2);
    std::istream input2(&inbuf2); //Convert streambuf to istream

	string line1, line2;
	string readname = "";
	while(getline(input1, line1) && getline(input2, line2)){
		assert(line1[0] == '@' && line2[0] == '@');
		vector<string> strs;
		boost::split(strs, line1, boost::is_any_of("\t "));
		readname = strs[0].substr(1);
		for(int32_t i=0; i<3; i++){
			getline(input1, line1);
			getline(input2, line2);
		}
		// low quality defined as Phred score < 4 for > 50% length
		int32_t lowcount1 = 0, lowcount2 = 0;
		for(int32_t i=0; i<line1.size(); i++){
			if(Phred33 && line1[i] < '!'+LowThresh)
				lowcount1++;
			else if(!Phred33 && line1[i] < '@'+LowThresh)
				lowcount1++;
		}
		for(int32_t i=0; i<line2.size(); i++){
			if(Phred33 && line2[i] < '!'+LowThresh)
				lowcount2++;
			else if(!Phred33 && line2[i] < '@'+LowThresh)
				lowcount2++;
		}
		if(1.0*lowcount1/line1.size() > LenPropThresh || 1.0*lowcount2/line2.size() > LenPropThresh)
			LowQualReadNames.push_back(readname);
	}
	file1.close();
	file2.close();

	sort(LowQualReadNames.begin(), LowQualReadNames.end());

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finding "<<(LowQualReadNames.size())<<" low quality reads."<<endl;
};

void ReadSalmonQuant(string quantfile, const map<string,int32_t>& Trans, vector<double>& SalmonQuant, vector<int32_t>& TransLength)
{
	SalmonQuant.clear();
	SalmonQuant.assign(Trans.size(), 0);
	TransLength.clear();
	TransLength.assign(Trans.size(), 0);

	ifstream input(quantfile);
	string line;
	int32_t linecount=0;
	while(getline(input, line)){
		linecount++;
		if(linecount==1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		map<string,int32_t>::const_iterator itmap=Trans.find(strs[0]);
		if(itmap==Trans.cend())
			cout<<(Trans.size())<<"\t"<<strs[0]<<endl;
		assert(itmap!=Trans.cend());
		SalmonQuant[itmap->second] = stod(strs[4]);
		TransLength[itmap->second] = stoi(strs[1]);
	}
	for (vector<int32_t>::iterator it = TransLength.begin(); it != TransLength.end(); it++)
		assert((*it) > 0);
	input.close();
};


void GetEqTrans(string eqclassfile, map<string,int32_t>& Trans, vector<string>& TransNames, map< vector<int32_t>,int32_t >& EqTransID, map< pair<int32_t,int32_t>,double >& Aux){
	EqTransID.clear();
	Aux.clear();

	ifstream input(eqclassfile);
	string line;

	int32_t linecount=0;
	int32_t numtrans=0;
	int32_t numeqclass=0;
	while(getline(input,line)){
		linecount++;
		if(linecount==1)
			numtrans=stoi(line);
		else if(linecount==2)
			numeqclass=stoi(line);
		else if(linecount<=2+numtrans){
			Trans[line]=linecount-3;
			TransNames.push_back(line);
		}
		else if(linecount>2+numtrans){
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of("\t"));

			int32_t numrelatedtid=stoi(strs[0]);
			assert(strs.size() == 2*numrelatedtid+2);
			vector<int32_t> TIDs(numrelatedtid);
			for(int32_t i=0; i<numrelatedtid; i++)
				TIDs[i]=stoi(strs[1+i]);
			EqTransID[TIDs]=linecount-3-numtrans;
			for(int32_t i=0; i<numrelatedtid; i++)
				Aux[make_pair(linecount-3-numtrans, TIDs[i])] = stod(strs[1+numrelatedtid+i]);
		}
	}

	input.close();
};

void GetSalmonWeightAssign(map< pair<int32_t,int32_t>,double >& WeightAssign, const map< pair<int32_t,int32_t>,double >& Aux, 
	const map< vector<int32_t>,int32_t >& EqTransID, const vector<double>& SalmonQuant)
{
	WeightAssign.clear();
	for(map< vector<int32_t>,int32_t >::const_iterator it=EqTransID.cbegin(); it!=EqTransID.cend(); it++){
		const vector<int32_t>& TIDs=it->first;
		int32_t eqID=it->second;
		vector<double> auxs;

		for(int32_t i=0; i<TIDs.size(); i++){
			map< pair<int32_t,int32_t>,double >::const_iterator itaux = Aux.find(make_pair(eqID, TIDs[i]));
			if(itaux == Aux.cend())
				cout<<"Error\n";
			auxs.push_back(itaux->second);
		}

		double denom=0.0;
		for(int32_t i=0; i<TIDs.size(); i++)
			denom += SalmonQuant[TIDs[i]] * auxs[i];
		if(denom>0){
			for(int32_t i=0; i<TIDs.size(); i++)
				WeightAssign[make_pair(eqID, TIDs[i])] = SalmonQuant[TIDs[i]] * auxs[i] / denom;
		}
		else{
			for(int32_t i=0; i<TIDs.size(); i++)
				WeightAssign[make_pair(eqID, TIDs[i])] = 0;
		}
	}
};

string GetFeature(string line, string key)
{
	size_t s=line.find(key);
	if(s==string::npos)
		return "";
	size_t t=line.find_first_of(";", s+1);
	if(t==string::npos)
		return "";
	// XXX check whether this is correct
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

void FLDKDE(vector<int32_t>& RawFLD, vector<double>& FLD, int32_t kernel_n=10, double kernel_p=0.5){
	// initialize fragment length distribution FLD
	FLD.assign(RawFLD.size(), 0);
	// calculate binomial kernel
	boost::math::binomial bino(kernel_n, kernel_p);
	vector<double> kernel(kernel_n+1, 0);
	for(int32_t i=0; i<kernel_n+1; i++)
		kernel[i]=pdf(bino, i);
	// calculate FLD based on kernel
	for(int32_t i=0; i<RawFLD.size(); i++){
		if(RawFLD[i]==0)
			continue;
		int32_t offset=max(0, i-kernel_n/2);
		while(offset<=i+kernel_n/2 && offset<RawFLD.size()){
			FLD[offset]+=RawFLD[i]*kernel[offset-i+kernel_n/2];
			offset++;
		}
	}
	double sum=0;
	for(int32_t i=0; i<RawFLD.size(); i++){
		FLD[i]+=1e-8;
		sum+=FLD[i];
	}
	for(int32_t i=0; i<RawFLD.size(); i++)
		FLD[i]/=sum;
};

void ReadProcessFLD(string filename, vector<double>& FLD){
	boost::iostreams::filtering_istream fpin;
	fpin.push(boost::iostreams::file_source(filename, std::ios_base::in | std::ios_base::binary));
	vector<int32_t> RawFLD(1001, 0);
	fpin.read((char*)&RawFLD[0], 1001*sizeof(int32_t));
	FLDKDE(RawFLD, FLD);
};

template <typename T>
T getMax(vector<T> x){
	assert(x.size()>0);
	T maxvalue=x[0];
	for(const T& v : x){
		if(v>maxvalue)
			maxvalue=v;
	}
	return maxvalue;
};

template <typename T>
T getMin(vector<T> x){
	assert(x.size()>0);
	T minvalue=x[0];
	for(const T& v : x){
		if(v<minvalue)
			minvalue=v;
	}
	return minvalue;
};

pair<string,int32_t> GetGenomicPosition(const vector<Transcript_t>& Transcripts, int32_t ind, int32_t pos)
{
	const Transcript_t& t=Transcripts[ind];
	int32_t genomepos=-1;
	if(t.Strand){
		int32_t coveredlen=0;
		for(vector< pair<int32_t,int32_t> >::const_iterator it=t.Exons.cbegin(); it!=t.Exons.cend(); it++){
			if(coveredlen+it->second-it->first < pos)
				coveredlen+=it->second-it->first;
			else if(coveredlen+it->second-it->first >= pos){
				genomepos = it->first+(pos-coveredlen);
				break;
			}
		}
	}
	else{
		int32_t coveredlen=0;
		for(vector< pair<int32_t,int32_t> >::const_iterator it=t.Exons.cbegin(); it!=t.Exons.cend(); it++){
			if(coveredlen+it->second-it->first < pos)
				coveredlen+=it->second-it->first;
			else if(coveredlen+it->second-it->first >= pos){
				genomepos = it->second-(pos-coveredlen);
				break;
			}
		}
	}
	assert(genomepos!=-1);
	return make_pair(t.Chr, genomepos);
};


int32_t WriteReadBrdy(map<string,int32_t>& Trans, vector<int32_t>& TransLength, vector<Brdy>& ReadBrdy, ofstream& ss)
{
	// vector of Trans names
	vector<string> TransName(Trans.size());
	for(map<string,int32_t>::iterator it=Trans.begin(); it!=Trans.end(); it++)
		TransName[it->second]=it->first;

	// sort read boundaries
	ReadBrdy.reserve(ReadBrdy.size());
	sort(ReadBrdy.begin(), ReadBrdy.end());

	int32_t numtrans=1;
	for(vector<Brdy>::iterator it=ReadBrdy.begin(); it!=ReadBrdy.end(); it++)
		if(it!=ReadBrdy.begin() && it->TID!=(it-1)->TID)
			numtrans++;
	ss.write((char*)(&numtrans), sizeof(int32_t));

	vector<int32_t> poses;
	vector<double> counts;
	for(vector<Brdy>::iterator it=ReadBrdy.begin(); it!=ReadBrdy.end(); it++){
		if(it!=ReadBrdy.begin() && it->TID!=(it-1)->TID){
			if(poses.back()!=TransLength[(it-1)->TID]-1){
				poses.push_back(TransLength[(it-1)->TID]-1);
				counts.push_back(0);
			}

			string name=TransName[(it-1)->TID];
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
	vector<Brdy>::iterator it=ReadBrdy.end(); it--;
	if(poses.back()!=TransLength[it->TID]-1){
		poses.push_back(TransLength[it->TID]-1);
		counts.push_back(0);
	}
	string name=TransName[it->TID];
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


int32_t ReadBAMStartPos(string bamfile, const vector<string>& LowQualReadNames, map<string,int32_t>& Trans, map< vector<int32_t>,int32_t >& EqTransID, 
	vector<int32_t>& TransLength, map< pair<int32_t,int32_t>,double >& WeightAssign, const map<string, Transcript_t>& transcripts, 
	ofstream& ss, ofstream& ss_junction){
	vector<Brdy> ReadBrdy;
	vector<Brdy> JunctionBrdy;

	// vector of Trans names
	vector<string> TransName(Trans.size());
	for(map<string,int32_t>::iterator it=Trans.begin(); it!=Trans.end(); it++)
		TransName[it->second]=it->first;

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
	// read bam records
	string prevreadname="";
	// vector<int32_t> readtids;
	// vector< pair<int32_t,int32_t> > alignments;
	vector<Alignment_t> alignments;
	while(sam_read1(bamreader, header, b)>0){
		// read current info
		string readname=bam_get_qname(b);
		map<string,int32_t>::iterator itmap = Trans.find((string)header->target_name[b->core.tid]);
		if (itmap == Trans.end())
			continue;
		int32_t tid=itmap->second;
		// pair<int32_t,int32_t> curalignment(b->core.pos, bam_endpos(b));
		// if((b)->core.flag&BAM_FUNMAP)
		// 	curalignment.second=curalignment.first;
		Alignment_t curalignment(tid, b->core.pos, bam_endpos(b), b->core.mpos);
		// add to group or process coverage
		if(readname==prevreadname)
			alignments.push_back(curalignment);
		else{
			if(prevreadname.size()!=0 && !binary_search(LowQualReadNames.cbegin(), LowQualReadNames.cend(), prevreadname)){
				if(alignments.size()>0){
					vector<Alignment_t> merged = MergeAlignments(alignments);
					// find the read in eq class
					vector<int32_t> tidgroup;
					for (vector<Alignment_t>::const_iterator it = merged.cbegin(); it != merged.cend(); it++)
						tidgroup.push_back(it->TID);
					sort(tidgroup.begin(), tidgroup.end());
					vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
					tidgroup.resize(distance(tidgroup.begin(), itend));

					int32_t eqID=EqTransID[tidgroup];
					for(int32_t i=0; i<merged.size(); i++){
						double w=WeightAssign[make_pair(eqID, merged[i].TID)];
						
						// coverage is counted as fragment, not each end in read pair
						int32_t fragstart;
						fragstart=min(merged[i].StartPos, merged[i].mate_StartPos);
						Brdy tmpstart(merged[i].TID, fragstart, false, w);
						assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
						ReadBrdy.push_back(tmpstart);

						// check whether hit splicing junction
						string transname(header->target_name[merged[i].TID]);
						map<string, Transcript_t>::const_iterator ittrans = transcripts.find(transname);
						assert(ittrans != transcripts.cend());
						if (SpanSplicing(merged[i], ittrans->second))
							JunctionBrdy.push_back(tmpstart);
					}
				}
			}
			prevreadname=readname;
			alignments.clear();
			alignments.push_back(curalignment);
		}
	}
	// store the last read alignments
	if(prevreadname.size()!=0 && !binary_search(LowQualReadNames.cbegin(), LowQualReadNames.cend(), prevreadname)){
		if(alignments.size()>0){
			vector<Alignment_t> merged = MergeAlignments(alignments);
			// find the read in the eq class
			vector<int32_t> tidgroup;
			for (vector<Alignment_t>::const_iterator it = merged.cbegin(); it != merged.cend(); it++)
				tidgroup.push_back(it->TID);
			sort(tidgroup.begin(), tidgroup.end());
			vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
			tidgroup.resize(distance(tidgroup.begin(), itend));

			int32_t eqID=EqTransID[tidgroup];
			for(int32_t i=0; i<merged.size(); i+=2){
				double w=WeightAssign[make_pair(eqID,merged[i].TID)];
				
				int32_t fragstart;
				fragstart=min(merged[i].StartPos, merged[i].mate_StartPos);
				Brdy tmpstart(merged[i].TID, fragstart, false, w);
				assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
				// Brdy tmpstart(readtids[i], fragmiddle, false, w);
				ReadBrdy.push_back(tmpstart);

				// check whether hit splicing junction
				string transname(header->target_name[merged[i].TID]);
				map<string, Transcript_t>::const_iterator ittrans = transcripts.find(transname);
				assert(ittrans != transcripts.cend());
				if (SpanSplicing(merged[i], ittrans->second))
					JunctionBrdy.push_back(tmpstart);
			}
		}
	}

	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(bamreader);

	int32_t flag1 = WriteReadBrdy(Trans, TransLength, ReadBrdy, ss);
	int32_t flag2 = WriteReadBrdy(Trans, TransLength, JunctionBrdy, ss_junction);
};


int32_t ReadBAMSingleMapStartPos2(string bamfile, const vector<string>& LowQualReadNames, map<string,int32_t>& Trans, map< vector<int32_t>,int32_t >& EqTransID, 
	vector<int32_t>& TransLength, map< pair<int32_t,int32_t>,double >& WeightAssign, ofstream& ss, ofstream& ss_bad){

	vector<Brdy> ReadBrdy;
	vector<Brdy> BadReadBrdy;

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
	// read bam records
	string prevreadname="";
	vector<int32_t> allreadtids;
	vector<int32_t> readtids;
	vector< pair<int32_t,int32_t> > alignments;
	while(sam_read1(bamreader, header, b)>0){
		// read current info
		string readname=bam_get_qname(b);
		map<string,int32_t>::iterator itmap = Trans.find((string)header->target_name[b->core.tid]);
		if (itmap == Trans.end())
			continue;
		int32_t tid=itmap->second;
		pair<int32_t,int32_t> curalignment(b->core.pos, bam_endpos(b));
		// add to group or process coverage
		if(readname==prevreadname){
			allreadtids.push_back(tid);
			if((b)->core.flag&BAM_FUNMAP){
				curalignment.second=curalignment.first;
				readtids.push_back(tid);
				alignments.push_back(curalignment);
			}
		}
		else{
			if(prevreadname.size()!=0 && !binary_search(LowQualReadNames.cbegin(), LowQualReadNames.cend(), prevreadname)){
				if(alignments.size()>0){
					// find the read in eq class
					vector<int32_t> tidgroup=allreadtids;
					sort(tidgroup.begin(), tidgroup.end());
					vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
					tidgroup.resize(distance(tidgroup.begin(), itend));

					int32_t eqID=EqTransID[tidgroup];
					for(int32_t i=0; i<alignments.size(); i++){
						double w=WeightAssign[make_pair(eqID,readtids[i])];
						
						// coverage is counted as fragment, not each end in read pair
						int32_t fragstart=alignments[i].first;
						Brdy tmpstart(readtids[i], fragstart, false, w);
						assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
						ReadBrdy.push_back(tmpstart);
					}
				}
			}
			else if(prevreadname.size()!=0){
				if(alignments.size()>0){
					// find the read in eq class
					vector<int32_t> tidgroup=allreadtids;
					sort(tidgroup.begin(), tidgroup.end());
					vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
					tidgroup.resize(distance(tidgroup.begin(), itend));

					int32_t eqID=EqTransID[tidgroup];
					for(int32_t i=0; i<alignments.size(); i++){
						double w=WeightAssign[make_pair(eqID,readtids[i])];
						
						// coverage is counted as fragment, not each end in read pair
						int32_t fragstart=alignments[i].first;
						Brdy tmpstart(readtids[i], fragstart, false, w);
						assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
						BadReadBrdy.push_back(tmpstart);
					}
				}
			}
			prevreadname=readname;
			allreadtids.clear();
			readtids.clear();
			alignments.clear();
			allreadtids.push_back(tid);
			if((b)->core.flag&BAM_FUNMAP){
				curalignment.second=curalignment.first;
				readtids.push_back(tid);
				alignments.push_back(curalignment);
			}
		}
	}
	// store the last read alignments
	if(prevreadname.size()!=0 && !binary_search(LowQualReadNames.cbegin(), LowQualReadNames.cend(), prevreadname)){
		if(alignments.size()>0){
			vector<int32_t> tidgroup=allreadtids;
			sort(tidgroup.begin(), tidgroup.end());
			vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
			tidgroup.resize(distance(tidgroup.begin(), itend));

			int32_t eqID=EqTransID[tidgroup];
			for(int32_t i=0; i<alignments.size(); i++){
				double w=WeightAssign[make_pair(eqID,readtids[i])];
				
				int32_t fragstart=alignments[i].first;
				Brdy tmpstart(readtids[i], fragstart, false, w);
				assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
				ReadBrdy.push_back(tmpstart);
			}
		}
	}
	else if(prevreadname.size()!=0){
		if(alignments.size()>0){
			vector<int32_t> tidgroup=allreadtids;
			sort(tidgroup.begin(), tidgroup.end());
			vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
			tidgroup.resize(distance(tidgroup.begin(), itend));

			int32_t eqID=EqTransID[tidgroup];
			for(int32_t i=0; i<alignments.size(); i++){
				double w=WeightAssign[make_pair(eqID,readtids[i])];
				
				int32_t fragstart=alignments[i].first;
				Brdy tmpstart(readtids[i], fragstart, false, w);
				assert(!std::isnan(tmpstart.Weight) && !std::isinf(tmpstart.Weight));
				BadReadBrdy.push_back(tmpstart);
			}
		}
	}

	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(bamreader);

	int32_t flag1 = WriteReadBrdy(Trans, TransLength, ReadBrdy, ss);
	int32_t flag2 = WriteReadBrdy(Trans, TransLength, BadReadBrdy, ss_bad);
	return min(flag1, flag2);
};

int32_t ReadSalmonFragLen(string bamfile, const vector<string>& LowQualReadNames, vector<double>& FLD, map<string,int32_t>& Trans, map< vector<int32_t>,int32_t >& EqTransID, 
	vector<int32_t>& TransLength, map< pair<int32_t,int32_t>,double >& WeightAssign, ofstream& ss, double cutoff=0.005)
{
	// vector of Trans names
	vector<string> TransName(Trans.size());
	for(map<string,int32_t>::iterator it=Trans.begin(); it!=Trans.end(); it++)
		TransName[it->second]=it->first;

	vector<Brdy> ReadBrdy;
	ReadBrdy.reserve(65536);

	// process lower FLD length and upper FLD length corresponding to distribution cutoff
	int32_t fldLow=0, fldHigh=FLD.size();
	double fldcummulative=0;
	bool lb=false, ub=false;
	for(int32_t i=0; i<FLD.size(); i++){
		if(lb && ub)
			break;
		fldcummulative += FLD[i];
		if(fldcummulative >= cutoff && !lb){
			fldLow = i;
			lb = true;
		}
		if(fldcummulative >= 1-cutoff && !ub){
			fldHigh = i-1;
			ub = true;
		}
	}

	// bam file bandler
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
	// read bam records
	string prevreadname="";
	vector<int32_t> readtids;
	vector< pair<int32_t,int32_t> > alignments;
	while(sam_read1(bamreader, header, b)>0){
		// read current info
		string readname=bam_get_qname(b);
		map<string,int32_t>::iterator itmap = Trans.find((string)header->target_name[b->core.tid]);
		if (itmap == Trans.end())
			continue;
		int32_t tid=itmap->second;
		pair<int32_t,int32_t> curalignment(b->core.pos, bam_endpos(b));

		if(readname==prevreadname){
			if(readtids.size()!=0 && tid==readtids.back()){
				alignments.back().first = min(alignments.back().first, curalignment.first);
				alignments.back().second = max(alignments.back().second, curalignment.second);
			}
			else{
				readtids.push_back(tid);
				alignments.push_back(curalignment);
			}
		}
		else{
			if(prevreadname.size()!=0 && !binary_search(LowQualReadNames.cbegin(), LowQualReadNames.cend(), prevreadname)){
				if(alignments.size()>0){
					// find the read in eq class
					vector<int32_t> tidgroup=readtids;
					sort(tidgroup.begin(), tidgroup.end());
					vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
					tidgroup.resize(distance(tidgroup.begin(), itend));

					int32_t eqID=EqTransID[tidgroup];

					assert(readtids.size() == alignments.size());
					for(int32_t i=0; i<readtids.size(); i++){
						double w=WeightAssign[make_pair(eqID,readtids[i])];
						if(alignments[i].second-alignments[i].first < fldLow || alignments[i].second-alignments[i].first > fldHigh){
							Brdy tmp(readtids[i], alignments[i].first, false, w);
							assert(!std::isnan(tmp.Weight) && !std::isinf(tmp.Weight));
							ReadBrdy.push_back(tmp);
						}
					}
				}
			}
			prevreadname=readname;
			readtids.clear();
			alignments.clear();

			readtids.push_back(tid);
			alignments.push_back(curalignment);
		}
	}
	// store the last alignment group
	if(prevreadname.size()!=0 && !binary_search(LowQualReadNames.cbegin(), LowQualReadNames.cend(), prevreadname)){
		if(alignments.size()>0){
			vector<int32_t> tidgroup=readtids;
			sort(tidgroup.begin(), tidgroup.end());
			vector<int32_t>::iterator itend=unique(tidgroup.begin(), tidgroup.end());
			tidgroup.resize(distance(tidgroup.begin(), itend));

			int32_t eqID=EqTransID[tidgroup];
			
			assert(readtids.size() == alignments.size());
			for(int32_t i=0; i<readtids.size(); i++){
				double w=WeightAssign[make_pair(eqID,readtids[i])];
				if(alignments[i].second-alignments[i].first < fldLow || alignments[i].second-alignments[i].first > fldHigh){
					Brdy tmp(readtids[i], alignments[i].first, false, w);
					assert(!std::isnan(tmp.Weight) && !std::isinf(tmp.Weight));
					ReadBrdy.push_back(tmp);
				}
			}
		}
	}

	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(bamreader);

	// compute fragment length mean and std for each transcript
	ReadBrdy.reserve(ReadBrdy.size());
	sort(ReadBrdy.begin(), ReadBrdy.end());

	int32_t numtrans=1;
	for(vector<Brdy>::iterator it=ReadBrdy.begin(); it!=ReadBrdy.end(); it++)
		if(it!=ReadBrdy.begin() && it->TID!=(it-1)->TID)
			numtrans++;
	ss.write((char*)(&numtrans), sizeof(int32_t));

	vector<int32_t> poses;
	vector<double> counts;
	for(vector<Brdy>::iterator it=ReadBrdy.begin(); it!=ReadBrdy.end(); it++){
		if(it!=ReadBrdy.begin() && it->TID!=(it-1)->TID){
			if(poses.back()!=TransLength[(it-1)->TID]-1){
				poses.push_back(TransLength[(it-1)->TID]-1);
				counts.push_back(0);
			}

			string name=TransName[(it-1)->TID];
			int32_t namelen=name.size();
			assert(counts.size()==poses.size());
			int32_t vectorlen=poses.size();
			ss.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
			ss.write(reinterpret_cast<char*>(&vectorlen), sizeof(int32_t));
			ss.write(name.c_str(), namelen*sizeof(char));
			ss.write(reinterpret_cast<char*>(poses.data()), vectorlen*sizeof(int32_t));
			ss.write(reinterpret_cast<char*>(counts.data()), vectorlen*sizeof(double));

			counts.clear();
			poses.clear();
		}
		if(poses.size()==0 || it->Pos!=poses.back()){
			poses.push_back(it->Pos);
			counts.push_back(it->Weight);
		}
		else
			counts.back()+=it->Weight;
	}
	vector<Brdy>::iterator it=ReadBrdy.end(); it--;
	if(poses.back()!=TransLength[it->TID]-1){
		poses.push_back(TransLength[it->TID]-1);
		counts.push_back(0);
	}
	string name=TransName[it->TID];
	int32_t namelen=name.size();
	assert(counts.size()==poses.size());
	int32_t vectorlen=poses.size();
	ss.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
	ss.write(reinterpret_cast<char*>(&vectorlen), sizeof(int32_t));
	ss.write(name.c_str(), namelen*sizeof(char));
	ss.write(reinterpret_cast<char*>(poses.data()), vectorlen*sizeof(int32_t));
	ss.write(reinterpret_cast<char*>(counts.data()), vectorlen*sizeof(double));
};

int32_t main(int32_t argc, char* argv[]){
	if(argc==1){
		printf("transcovdist <mode> <gtffile> <quantfile> <eqfile> <bamfile> <outfile> (--fld <FLDfile if mode 2>\n");
		printf("\tmode 0: all salmon mapped reads\n");
		printf("\tmode 1: single-end mapped redas by salmon\n");
		printf("\tmode 2: fragment length mean and std\n");
	}
	else{
		string GTFfile(argv[2]);
		string QuantFile(argv[3]);
		string EqFile(argv[4]);
		string BamFile(argv[5]);
		string OutFile(argv[6]);

		string FLDfile="";
		for(int32_t i=7; i<argc; i+=2){
			if(string(argv[i])=="--fld")
				FLDfile=string(argv[i+1]);
		}

		map<string, Transcript_t> transcripts;
		ReadGTF(GTFfile, transcripts);

		vector<string> LowQualReadNames;

		map<string,int32_t> Trans;
		vector<string> TransNames;
		map< vector<int32_t>,int32_t > EqTransID;
		map< pair<int32_t,int32_t>,double > Aux;
		map< pair<int32_t,int32_t>,double > WeightAssign;
		vector<double> SalmonQuant;
		vector<int32_t> TransLength;
		ofstream ss(OutFile, ios::out | ios::binary);

		GetEqTrans(EqFile, Trans, TransNames, EqTransID, Aux);
		ReadSalmonQuant(QuantFile, Trans, SalmonQuant, TransLength);
		GetSalmonWeightAssign(WeightAssign, Aux, EqTransID, SalmonQuant);

		/*GetEqTrans(EqFile, Trans, TransNames, EqTransID);
		GetSalmonWeightAssign(AssignFile, WeightAssign);*/

		if(atoi(argv[1])==0) {
			string OutJunctionFile;
			size_t suffixpos = OutFile.find_last_of(".");
			if(suffixpos != string::npos)
				OutJunctionFile = OutFile.substr(0, suffixpos) + "_junction"+OutFile.substr(suffixpos);
			else
				OutJunctionFile = OutFile+"_junction";
			ofstream ss_junction(OutJunctionFile, ios::out | ios::binary);
			ReadBAMStartPos(BamFile, LowQualReadNames, Trans, EqTransID, TransLength, WeightAssign, transcripts, ss, ss_junction);
		}
		else if(atoi(argv[1])==1){
			string OutBadFile;
			size_t suffixpos = OutFile.find_last_of(".");
			if(suffixpos != string::npos)
				OutBadFile = OutFile.substr(0, suffixpos) + "_bad"+OutFile.substr(suffixpos);
			else
				OutBadFile = OutFile+"_bad";
			ofstream ss_bad(OutBadFile, ios::out | ios::binary);
			ReadBAMSingleMapStartPos2(BamFile, LowQualReadNames, Trans, EqTransID, TransLength, WeightAssign, ss, ss_bad);
		}
		else if(atoi(argv[1])==2){
			assert(FLDfile!="");
			vector<double> FLD;
			ReadProcessFLD(FLDfile, FLD);
			ReadSalmonFragLen(BamFile, LowQualReadNames, FLD, Trans, EqTransID, TransLength, WeightAssign, ss);
		}
		ss.close();
	}
}
