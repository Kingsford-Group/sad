/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
#include <mutex>
#include <numeric>
#include "boost/algorithm/string.hpp"

using namespace std;


void ReadTransSequence(string filename, vector<string>& TransSequences, vector<string>& TransNames){
	ifstream input(filename);
	string line;
	string curname="";
	string curseq="";
	while(getline(input, line)){
		if(line[0]=='>'){
			if(curname!=""){
				TransSequences.push_back(curseq);
				TransNames.push_back(curname);
			}
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" |"));
			curname=strs[0].substr(1);
			curseq="";
		}
		else
			curseq+=line;
	}
	TransSequences.push_back(curseq);
	TransNames.push_back(curname);
};


vector<double> ReadRSEMrspd(string filename)
{
	// initialze result variable
	vector<double> rspd;
	// read file
	ifstream input(filename);
	string line;
	string prevline;
	bool is_rspd_line = false;

	while(getline(input, line)) {
		if (is_rspd_line) {
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of("\t "));
			assert( strs.size() == stoi(prevline) );

			for (uint32_t i = 0; i < strs.size(); i++)
				rspd.push_back( stod(strs[i]) );
			rspd.reserve( rspd.size() );
			is_rspd_line = false;
		}

		// check this line and previous line to see whether it is the rspd section.
		if (prevline != "1") {
			prevline = line;
			continue;
		}
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t "));
		if (strs.size() == 1 && strs[0] != "") 
			is_rspd_line = true;
		prevline = line;
	}
	input.close();

	if (rspd.size() == 0)
		cerr << "Cannot find the rspd information in file " << filename << endl;

#ifndef NDEBUG
	double sum_rspd = std::accumulate(rspd.begin(), rspd.end(), 0.0);
#endif
	assert( fabs(sum_rspd - 1) < 1e-4);

	return rspd;
};


vector<double> BiasCorrectTrans(string seq, const vector<double>& rspd)
{
	assert( rspd.size() > 0 );
	assert( seq.size() >= rspd.size() );

	// initialze result vector
	vector<double> correction(seq.size(), 0);
	int32_t idx = 0;
	int32_t start = (int32_t)round(1.0 * idx * seq.size() / rspd.size());
	int32_t end = (int32_t)round(1.0 * (idx + 1) * seq.size() / rspd.size());
	for (int32_t i = 0; i < (int32_t)seq.size(); i++) {
		if (i >= end) {
			idx += 1;
			start = (int32_t)round(1.0 * idx * seq.size() / rspd.size());
			end = (int32_t)round(1.0 * (idx + 1) * seq.size() / rspd.size());
			assert( end <= seq.size() );
			double tmpsum = std::accumulate(correction.begin(), correction.begin() + i, 0.0);
			double tmptruth = std::accumulate(rspd.begin(), rspd.begin() + idx, 0.0);
			if (fabs(tmpsum - tmptruth) > 1e-8)
				cout << "watch here\n";
		}
		correction[i] = rspd[idx] / (end - start);
	}
	// check the sum of correction vector
	double sum = std::accumulate(correction.begin(), correction.end(), 0.0);
	if (fabs(sum - 1) >= 1e-4)
		cout << "watch here\n";
	assert( fabs(sum - 1) < 1e-4 );
	return correction;
};


void WriteRawCorrection(string outfile, map< string,vector<double> >& Corrections){
	ofstream fpout(outfile, ios::out | ios::binary);
	int32_t numtrans=Corrections.size();
	fpout.write(reinterpret_cast<char*>(&numtrans), sizeof(int32_t));
	for(map< string,vector<double> >::iterator it=Corrections.begin(); it!=Corrections.end(); it++){
		int32_t namelen=it->first.size();
		int32_t seqlen=it->second.size();
		string name=it->first;
		fpout.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
		fpout.write(reinterpret_cast<char*>(&seqlen), sizeof(int32_t));
		fpout.write(name.c_str(), namelen*sizeof(char));
		fpout.write(reinterpret_cast<char*>(it->second.data()), seqlen*sizeof(double));
	}
	fpout.close();
};


int32_t main(int32_t argc, char* argv[]){
	if(argc==1){
		printf("readrsembias correction <ModelFile> <TransFasta> <OutputFile>\n");
	}
	else{
		string ModelFile(argv[2]);
		string TransFasta(argv[3]);
		string OutputFile(argv[4]);

		// read read start position distribution (rspd)
		vector<double> rspd = ReadRSEMrspd(ModelFile);

		// read transcript sequences
		vector<string> TransSequences;
		vector<string> TransNames;
		ReadTransSequence(TransFasta, TransSequences, TransNames);
		cout<<"num transcripts = "<<(TransSequences.size())<<endl;

		// correction
		map<string, vector<double> > Corrections;
		for(uint32_t i=0; i<TransSequences.size(); i++){
			if(TransSequences[i].size() < rspd.size()){
				cout<<TransNames[i]<<"\t"<<(TransSequences[i].size())<<endl;
				continue;
			}
			vector<double> correction = BiasCorrectTrans(TransSequences[i], rspd);
			Corrections[TransNames[i]] = correction;

		}

		WriteRawCorrection(OutputFile, Corrections);
	}
}
