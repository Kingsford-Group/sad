/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "Transcript.hpp"

// small function to extract tag information in GTF file
string GetFeature(string line, string key)
{
	size_t s = line.find(key);
	if (s == string::npos){
		cout << "Cannot find key \"" << key <<"\" in string \"" << line << "\"\n";
		return "";
	}
	size_t t = line.find_first_of(";", s+1);
	return line.substr(s + key.size() + 2, t - s - key.size() - 3);
};

// read transcript info from GTF file
void ReadGTF (string GTFfile, vector<Transcript_t>& Transcripts)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Start GTF file."<<endl;
	// clear outputs and intermediate variables
	// ExtraExons is used to store exons who is not directly following its corresponding transcript
	Transcripts.clear();
	vector< pair<string,Exon_t> > ExtraExons;
	// reading GTF file
	ifstream input(GTFfile);
	string line;
	string prevtransid = "";
	while(getline(input, line)) {
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if (strs[2] == "transcript"){
			assert(strs[6] == "+" || strs[6] == "-");
			string transid = GetFeature(line, "transcript_id");
			string geneid = GetFeature(line, "gene_id");
			Transcript_t tmp(transid, geneid, strs[0], stoi(strs[3])-1, stoi(strs[4]), (strs[6] == "+"));
			Transcripts.push_back(tmp);
			// aftering add to Transcripts, update previous trans id info
			prevtransid = transid;
		}
		else if (strs[2] == "exon"){
			string transid = GetFeature(line, "transcript_id");
			Exon_t tmpexon(stoi(strs[3])-1, stoi(strs[4]));
			if (transid == prevtransid)
				Transcripts.back().AddExon(tmpexon);
			else
				ExtraExons.push_back( make_pair(transid, tmpexon) );
		}
	}
	input.close();
	// adding extra exon to corresponding transcript
	map<string,int32_t> TransIndex;
	for (int32_t i=0; i<Transcripts.size(); i++)
		TransIndex[Transcripts[i].TransID] = i;
	for (vector< pair<string,Exon_t> >::iterator it = ExtraExons.begin(); it != ExtraExons.end(); it++){
		map<string,int32_t>::iterator itmap = TransIndex.find(it->first);
		assert(itmap != TransIndex.end());
		Transcripts[itmap->second].AddExon(it->second);
	}
	// sort exons for each transcript
	for (vector<Transcript_t>::iterator it = Transcripts.begin(); it != Transcripts.end(); it++)
		it->SortExons();

	Transcripts.reserve(Transcripts.size());

	// time info
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish reading GTF. Number transcripts = "<<(Transcripts.size())<<endl;

	return;
};

// trim transcripts: salmon mark some transcripts as duplicated and remove them from quantification
// We also remove them from any later steps according to what salmon removed
void TrimTranscripts (const map<string,double>& SalmonExp, vector<Transcript_t>& Transcripts)
{
	vector<Transcript_t> newTranscripts;
	newTranscripts.reserve(Transcripts.size());
	for (vector<Transcript_t>::iterator it = Transcripts.begin(); it != Transcripts.end(); it++){
		map<string,double>::const_iterator itmap = SalmonExp.find(it->TransID);
		if (itmap != SalmonExp.cend())
			newTranscripts.push_back((*it));
	}
	newTranscripts.reserve(newTranscripts.size());
	Transcripts = newTranscripts;

	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Trim transcripts: " << (Transcripts.size()) << " transcripts are kept\n";
	return;
};

// get GeneID and TransID mapping
void Map_Gene_Trans (const vector<Transcript_t>& Transcripts, map<string,int32_t>& TransIndex, 
	map<string,string>& TransGeneMap, map< string,vector<string> >& GeneTransMap)
{
	// clear variables
	TransIndex.clear();
	TransGeneMap.clear();
	GeneTransMap.clear();
	// TransIndex
	for (int32_t i=0; i < Transcripts.size(); i++)
		TransIndex[ Transcripts[i].TransID ] = i;
	// TransGeneMap
	vector< pair<string,string> > GeneTransIDs; // pair < GeneID, TransID >
	for (vector<Transcript_t>::const_iterator it = Transcripts.cbegin(); it != Transcripts.cend(); it++){
		GeneTransIDs.push_back( make_pair(it->GeneID, it->TransID) );
		TransGeneMap[it->TransID] = it->GeneID;
	}
	// GeneTransMap
	sort(GeneTransIDs.begin(), GeneTransIDs.end(), [](pair<string,string> a, pair<string,string> b){return a.first < b.first;} );
	string prevgeneid = GeneTransIDs.front().first;
	vector<string> TransIDs;
	for (vector< pair<string,string> >::iterator it = GeneTransIDs.begin(); it != GeneTransIDs.end(); it++){
		if (it->first != prevgeneid){
			assert(TransIDs.size() > 0);
			GeneTransMap[prevgeneid] = TransIDs;
			prevgeneid = it->first;
			TransIDs.clear();
		}
		TransIDs.push_back(it->second);
	}
	if (TransIDs.size() > 0)
		GeneTransMap[prevgeneid] = TransIDs;

	return;
};

// collect non-redundant unique exonic positions for each genes
void CombineGeneLevelExons (const vector<Transcript_t>& Transcripts, const map<string,int32_t>& TransIndex, 
	const map< string,vector<string> >& GeneTransMap, map<string, vector<Exon_t> >& GeneLevelExons)
{
	GeneLevelExons.clear();
	for (map< string,vector<string> >::const_iterator it = GeneTransMap.cbegin(); it != GeneTransMap.cend(); it++){
		const vector<string>& trans = it->second;
		// collect exons, redundant and overlapping
		vector<Exon_t> allexons;
		bool gene_strand = true;
		for (const string& t : trans) {
			map<string,int32_t>::const_iterator itidx = TransIndex.find(t);
			assert(itidx != TransIndex.cend());
			// add to variable
			if (t == trans[0])
				gene_strand = Transcripts[itidx->second].Strand;
			else
				assert(gene_strand == Transcripts[itidx->second].Strand);
			allexons.insert( allexons.end(), Transcripts[itidx->second].Exons.cbegin(), Transcripts[itidx->second].Exons.cend() );
		}
		// sort and select unique exons
		sort(allexons.begin(), allexons.end());
		vector<Exon_t> uniqueexons;
		for (vector<Exon_t>::iterator itexon = allexons.begin(); itexon != allexons.end(); itexon++) {
			assert( uniqueexons.size() == 0 || itexon->StartPos >= uniqueexons.back().StartPos);
			if (uniqueexons.size() == 0 || itexon->StartPos > uniqueexons.back().EndPos)
				uniqueexons.push_back(*itexon);
			else
				uniqueexons.back().EndPos = max(uniqueexons.back().EndPos, itexon->EndPos);
		}
		assert(uniqueexons.size() != 0);
		if (!gene_strand)
			reverse(uniqueexons.begin(), uniqueexons.end());
		// add to gene exon map
		GeneLevelExons[it->first] = uniqueexons;
	}
	assert(GeneLevelExons.size() == GeneTransMap.size());
};

