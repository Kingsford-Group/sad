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
#include "Transcript.hpp"

using namespace std;


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


void ReadRemoveList(string filename, vector<string>& RemovedTrans)
{
	RemovedTrans.clear();
	ifstream input(filename);
	string line;
	while (getline(input, line)){
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of(" \t"));
		RemovedTrans.push_back(strs[0]);
	}
	input.close();
};


void ReadFusionList(string filename, vector<string>& FusedTrans_ID, vector< pair<Junction_t,Junction_t> >& FusedTrans)
{
	FusedTrans.clear();
	FusedTrans_ID.clear();
	ifstream input(filename);
	string line;
	while (getline(input, line)) {
		if (line.substr(0,7) == ">Fusion") {
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));
			vector<string> fuse1;
			boost::split(fuse1, strs[1], boost::is_any_of(":"));
			assert(fuse1.size() == 3);
			vector<string> fuse2;
			boost::split(fuse2, strs[2], boost::is_any_of(":"));
			assert(fuse2.size() == 3);
			Junction_t tmp1(fuse1[0], stoi(fuse1[1]), stoi(fuse1[2]));
			Junction_t tmp2(fuse2[0], stoi(fuse2[1]), stoi(fuse2[2]));
			FusedTrans.push_back( make_pair(tmp1, tmp2) );
			FusedTrans_ID.push_back(strs[0].substr(1));
		}
	}
	input.close();
};


// separate transcripts into two vectors, one within the reference, the other within remove list
void SeparateTranscripts(const vector<Transcript_t>& Transcripts, const map<string,string>& TransGeneMap, 
	vector<Transcript_t>& Transcripts_ref, vector<Transcript_t>& Transcripts_novel, 
	const vector<string>& RemovedTrans, const vector< pair<Junction_t,Junction_t> >& FusedTrans)
{
	Transcripts_ref.clear();
	Transcripts_novel.clear();
	// collected genes related to removal and fusion
	vector<string> RelatedGenes;
	for (int32_t i = 0; i < RemovedTrans.size(); i++){
		map<string,string>::const_iterator itmap = TransGeneMap.find(RemovedTrans[i]);
		assert( itmap != TransGeneMap.cend() );
		RelatedGenes.push_back(itmap->second);
	}
	for (vector< pair<Junction_t,Junction_t> >::const_iterator it = FusedTrans.cbegin(); it != FusedTrans.cend(); it++){
		map<string,string>::const_iterator itmap = TransGeneMap.find((it->first).Chr);
		assert( itmap != TransGeneMap.cend() );
		RelatedGenes.push_back(itmap->second);

		itmap = TransGeneMap.find((it->second).Chr);
		assert( itmap != TransGeneMap.cend() );
		RelatedGenes.push_back(itmap->second);
	}
	sort(RelatedGenes.begin(), RelatedGenes.end());
	vector<string>::iterator itunique = unique(RelatedGenes.begin(), RelatedGenes.end());
	RelatedGenes.resize( distance(RelatedGenes.begin(), itunique) );
	// sort remove list
	vector<string> RemovedTrans_sort(RemovedTrans.cbegin(), RemovedTrans.cend());
	sort(RemovedTrans_sort.begin(), RemovedTrans_sort.end());
	// sort fusion list
	vector<Junction_t> FusedTrans_sort;
	for(vector< pair<Junction_t,Junction_t> >::const_iterator it = FusedTrans.cbegin(); it != FusedTrans.cend(); it++){
		FusedTrans_sort.push_back(it->first);
		FusedTrans_sort.push_back(it->second);
	}
	sort(FusedTrans_sort.begin(), FusedTrans_sort.end());
	vector<string> FusedTransNames;
	for (int32_t i = 0; i < FusedTrans_sort.size(); i++)
		FusedTransNames.push_back( FusedTrans_sort[i].Chr ); // Junction_t.Chr is used as transcript name here
	// separating 
	for (int32_t i = 0; i < Transcripts.size(); i++) {
		if (binary_search(RelatedGenes.begin(), RelatedGenes.end(), Transcripts[i].GeneID)) {
			Transcript_t t = Transcripts[i];
			// check if the transcript is removed
			if (binary_search(RemovedTrans_sort.begin(), RemovedTrans_sort.end(), t.TransID))
				Transcripts_novel.push_back( t );
			else
				Transcripts_ref.push_back( t );
			// check if we also need to add to the fusion part to novel
			if (binary_search(FusedTransNames.begin(), FusedTransNames.end(), t.TransID)) {
				vector<string>::iterator itfuse = lower_bound(FusedTransNames.begin(), FusedTransNames.end(), t.TransID);
				const Junction_t& fuse = FusedTrans_sort[ distance(FusedTransNames.begin(), itfuse) ];
				// trim transcript based on the part of fusion
				// XXX this is probably wrong
				vector<Exon_t> newExons;
				int32_t curlength = 0;
				if (t.Strand) {
					for (int32_t j = 0; j < t.Exons.size(); j++){
						if (curlength <= fuse.Pos1 && curlength + t.Exons[j].EndPos - t.Exons[j].StartPos > fuse.Pos1){
							int32_t newexon_start = t.Exons[j].StartPos + fuse.Pos1 - curlength;
							int32_t newexon_end = min(t.Exons[j].EndPos, t.Exons[j].StartPos + fuse.Pos2 - curlength);
							Exon_t tmpexon(newexon_start, newexon_end);
							newExons.push_back(tmpexon);
						}
						else if (curlength <= fuse.Pos2 && curlength + t.Exons[j].EndPos - t.Exons[j].StartPos > fuse.Pos2){
							int32_t newexon_start = max(t.Exons[j].StartPos, t.Exons[j].StartPos + fuse.Pos1 - curlength);
							int32_t newexon_end = t.Exons[j].StartPos + fuse.Pos2 - curlength;
							Exon_t tmpexon(newexon_start, newexon_end);
							newExons.push_back(tmpexon);
						}
						else if (curlength > fuse.Pos1 && curlength + t.Exons[j].EndPos - t.Exons[j].StartPos <= fuse.Pos2)
							newExons.push_back(t.Exons[j]);
						curlength += t.Exons[j].EndPos - t.Exons[j].StartPos;
					}
				}
				else{
					for (int32_t j = 0; j < t.Exons.size(); j++){
						if (curlength <= fuse.Pos1 && curlength + t.Exons[j].EndPos - t.Exons[j].StartPos > fuse.Pos1){
							int32_t newexon_end = t.Exons[j].EndPos - fuse.Pos1 + curlength;
							int32_t newexon_start = max(t.Exons[j].StartPos, t.Exons[j].EndPos - fuse.Pos2 + curlength);
							Exon_t tmpexon(newexon_start, newexon_end);
							newExons.push_back(tmpexon);
						}
						else if (curlength <= fuse.Pos2 && curlength + t.Exons[j].EndPos - t.Exons[j].StartPos > fuse.Pos2){
							int32_t newexon_end = min(t.Exons[j].EndPos, t.Exons[j].EndPos - fuse.Pos1 + curlength);
							int32_t newexon_start = t.Exons[j].EndPos - fuse.Pos2 + curlength;
							Exon_t tmpexon(newexon_start, newexon_end);
							newExons.push_back(tmpexon);
						}
						else if (curlength > fuse.Pos1 && curlength + t.Exons[j].EndPos - t.Exons[j].StartPos <= fuse.Pos2)
							newExons.push_back(t.Exons[j]);
						curlength += t.Exons[j].EndPos - t.Exons[j].StartPos;
					}
				}
				t.Exons = newExons;
				t.TransID = "Fusion_" + t.TransID;
				// add to novel transcripts
				Transcripts_novel.push_back( t );
			}
		}
	}
};


void LabelNovelJunction(const vector<Transcript_t>& Transcripts_novel, const vector<Junction_t>& Junctions, 
	map<string,bool>& NovelJunctionLabel)
{
	int32_t thresh = 5;
	NovelJunctionLabel.clear();
	for (vector<Transcript_t>::const_iterator it = Transcripts_novel.begin(); it != Transcripts_novel.end(); it++){
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


void WriteCategory(string OutputPrefix, const vector<string>& RemovedTrans, const vector<string>& FusedTrans_ID, 
	const vector< pair<Junction_t,Junction_t> >& FusedTrans, const map<string,bool>& NovelJunctionLabel)
{
	// writing remove list to separate files
	ofstream output_rem_novel(OutputPrefix+"removelist_noveljunc.txt");
	ofstream output_rem_exist(OutputPrefix+"removelist_existjunc.txt");
	for (int32_t i = 0; i < RemovedTrans.size(); i++){
		map<string,bool>::const_iterator itmap = NovelJunctionLabel.find(RemovedTrans[i]);
		assert(itmap != NovelJunctionLabel.cend());
		if (itmap->second)
			output_rem_novel << RemovedTrans[i] << endl;
		else
			output_rem_exist << RemovedTrans[i] << endl;
	}
	output_rem_novel.close();
	output_rem_exist.close();
	// writing fusion to separate files
	assert(FusedTrans_ID.size() == FusedTrans.size());
	ofstream output_fusion_novel(OutputPrefix+"fusion_noveljunc.txt");
	ofstream output_fusion_exist(OutputPrefix+"fusion_existjunc.txt");
	for (int32_t i = 0; i < FusedTrans.size(); i++) {
		map<string,bool>::const_iterator itmap1 = NovelJunctionLabel.find("Fusion_"+(FusedTrans[i].first.Chr));
		assert(itmap1 != NovelJunctionLabel.cend());
		map<string,bool>::const_iterator itmap2 = NovelJunctionLabel.find("Fusion_"+(FusedTrans[i].second.Chr));
		assert(itmap2 != NovelJunctionLabel.cend());
		if (itmap1->second || itmap2->second){
			output_fusion_novel <<">"<< FusedTrans_ID[i] <<" "<< (FusedTrans[i].first.Chr) <<":"<< (FusedTrans[i].first.Pos1) <<":"<< (FusedTrans[i].first.Pos2);
			output_fusion_novel <<" "<< (FusedTrans[i].second.Chr) <<":"<< (FusedTrans[i].second.Pos1) <<":"<< (FusedTrans[i].second.Pos2) << endl;
		}
		if (!itmap1->second || !itmap2->second){
			output_fusion_exist <<">"<< FusedTrans_ID[i] <<" "<< (FusedTrans[i].first.Chr) <<":"<< (FusedTrans[i].first.Pos1) <<":"<< (FusedTrans[i].first.Pos2);
			output_fusion_exist <<" "<< (FusedTrans[i].second.Chr) <<":"<< (FusedTrans[i].second.Pos1) <<":"<< (FusedTrans[i].second.Pos2) << endl;
		}
	}
	output_fusion_novel.close();
	output_fusion_exist.close();
};


int32_t main(int32_t argc, char* argv[])
{
	if (argc == 1)
		printf("CategorizeSimulation <FullGTF> <RemoveFile> <FusionFile> <OutputPrefix>\n");
	else{
		string FullGTF(argv[1]);
		string RemoveFile(argv[2]);
		string FusionFile(argv[3]);
		string OutputPrefix(argv[4]);

		vector<Transcript_t> Transcripts;
		vector<string> RemovedTrans;
		vector<string> FusedTrans_ID;
		vector< pair<Junction_t, Junction_t> > FusedTrans;

		ReadGTF(FullGTF, Transcripts);
		ReadRemoveList(RemoveFile, RemovedTrans);
		ReadFusionList(FusionFile, FusedTrans_ID, FusedTrans);

		map<string,int32_t> TransIndex;
		map<string,string> TransGeneMap;
		map< string,vector<string> > GeneTransMap;
		Map_Gene_Trans(Transcripts, TransIndex, TransGeneMap, GeneTransMap);

		vector<Transcript_t> Transcripts_ref;
		vector<Transcript_t> Transcripts_novel;
		vector<Junction_t> AllJunctions;
		SeparateTranscripts(Transcripts, TransGeneMap, Transcripts_ref, Transcripts_novel, RemovedTrans, FusedTrans);
		CollectJunction(Transcripts_ref, AllJunctions);

		map<string,bool> NovelJunctionLabel;
		LabelNovelJunction(Transcripts_novel, AllJunctions, NovelJunctionLabel);
		WriteCategory(OutputPrefix, RemovedTrans, FusedTrans_ID, FusedTrans, NovelJunctionLabel);
	}
};