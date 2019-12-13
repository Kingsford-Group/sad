/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
#include <stdexcept>

#include "LPReassign.hpp"
#include "Transcript.hpp"
#include "DistTest.hpp"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

using namespace std;

template <class T>
void reorder(vector<T>& arr, vector<int32_t>& index)
{
	assert(arr.size() == index.size());
	vector<T> newarr;
	for (int32_t i = 0; i < index.size(); i++)
		newarr.push_back(arr[index[i]]);
	arr = newarr;
};

LPReassign_t::LPReassign_t(const map<string,int32_t>& _TransIndex, const map<string,string>& _TransGeneMap, 
	const map<string, vector<string> >& _GeneTransMap, const vector<Transcript_t>& Transcripts, 
	const vector< vector<double> >& Expected, const vector< vector<double> >& Observed)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Initializing LP reassignment."<<endl;

	// copy gene and transcripts mapping
	TransIndex = _TransIndex;
	TransGeneMap.assign(_TransGeneMap.size(), "");
	for(map<string,int32_t>::const_iterator itidx = _TransIndex.cbegin(); itidx != _TransIndex.cend(); itidx++){
		map<string,string>::const_iterator itmap = _TransGeneMap.find(itidx->first);
		assert(itmap != _TransGeneMap.cend());
		TransGeneMap[itidx->second] = itmap->second;
	}
	GeneTransMap.clear();
	for(map<string, vector<string> >::const_iterator itgene = _GeneTransMap.cbegin(); itgene != _GeneTransMap.cend(); itgene++) {
		vector<int32_t> tmp;
		for(vector<string>::const_iterator ittrans = (itgene->second).cbegin(); ittrans != (itgene->second).cend(); ittrans++) {
			map<string,int32_t>::const_iterator itidx = _TransIndex.find(*ittrans);
			assert(itidx != _TransIndex.cend());
			tmp.push_back(itidx->second);
		}
		GeneTransMap[itgene->first] = tmp;
	}
	// clear variables
	PositionExistence.clear();
	for (int32_t i = 0; i < TransIndex.size(); i++){
		vector<bool> tmp(Expected[i].size(), 0);
		PositionExistence.push_back(tmp);
	}
	// convert distributions to full length gene-level exons
	// and construct the PositionExistence indicator
	map< string, vector<Exon_t> > GeneLevelExons;
	CombineGeneLevelExons(Transcripts, _TransIndex, _GeneTransMap, GeneLevelExons);
	assert(GeneLevelExons.size() == _GeneTransMap.size());
	for (map< string,vector<Exon_t> >::const_iterator it = GeneLevelExons.cbegin(); it != GeneLevelExons.cend(); it++) {
		const vector<Exon_t>& totalexons = it->second;
		const vector<string>& trans = _GeneTransMap.at(it->first);
		// calculate total length of unique exons
		int32_t totallen = 0;
		for (int32_t i = 0; i < totalexons.size(); i++)
			totallen += totalexons[i].EndPos - totalexons[i].StartPos;
		// mapping the exons of each transcript to this totalexons position
		for (const string& t : trans) {
			map<string,int32_t>::const_iterator itidx = TransIndex.find(t);
			assert(itidx != TransIndex.cend());
			const vector<Exon_t>& thisexons = Transcripts[itidx->second].Exons;
			// set temporary variables
			vector<bool> tmpExistence(totallen, false);
			int32_t ind_total = 0;
			int32_t ind_this = 0;
			int32_t pos_total = 0;
			if (Transcripts[itidx->second].Strand){
				for(ind_total = 0; ind_total < totalexons.size(); ind_total++) {
					assert(ind_this >= thisexons.size() || totalexons[ind_total].EndPos < thisexons[ind_this].StartPos || 
						(totalexons[ind_total].StartPos <= thisexons[ind_this].StartPos && totalexons[ind_total].EndPos >= thisexons[ind_this].EndPos));
					while (ind_this < thisexons.size() && thisexons[ind_this].StartPos <= totalexons[ind_total].EndPos){
						int32_t overlap_start = max(thisexons[ind_this].StartPos, totalexons[ind_total].StartPos);
						int32_t overlap_end = min(thisexons[ind_this].EndPos, totalexons[ind_total].EndPos);
						for (int32_t i = overlap_start; i < overlap_end; i++){
							tmpExistence[pos_total + i - totalexons[ind_total].StartPos] = 1;
						}
						ind_this ++;
					}
					pos_total += totalexons[ind_total].EndPos - totalexons[ind_total].StartPos;
				}
			}
			else{
				for(ind_total = 0; ind_total < totalexons.size(); ind_total++) {
					assert(ind_this >= thisexons.size() || totalexons[ind_total].StartPos > thisexons[ind_this].EndPos || 
						(totalexons[ind_total].StartPos <= thisexons[ind_this].StartPos && totalexons[ind_total].EndPos >= thisexons[ind_this].EndPos));
					while (ind_this < thisexons.size() && thisexons[ind_this].EndPos > totalexons[ind_total].StartPos){
						int32_t overlap_start = max(thisexons[ind_this].StartPos, totalexons[ind_total].StartPos);
						int32_t overlap_end = min(thisexons[ind_this].EndPos, totalexons[ind_total].EndPos);
						for (int32_t i = overlap_start; i < overlap_end; i++){
							tmpExistence[pos_total + totalexons[ind_total].EndPos - i - 1] = 1;
						}
						ind_this ++;
					}
					pos_total += totalexons[ind_total].EndPos - totalexons[ind_total].StartPos;
				}
			}
			// sanity check length
			int32_t len_processed = 0;
			for (int32_t i = 0; i < tmpExistence.size(); i++)
				len_processed += tmpExistence[i];
			int32_t len_theo = 0;
			for (int32_t i = 0; i < thisexons.size(); i++)
				len_theo += thisexons[i].EndPos - thisexons[i].StartPos;
			assert(len_processed == len_theo);

			// add to PositionExistence
			PositionExistence[itidx->second] = tmpExistence;
		}
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Initializing: step 1 calculating position existence matrix. Finished\n";

	// process full length expected and observed distribution
	// expected distribution sum to 1, observed distribution keeps its original value
	ExpectedFullNorm.clear();
	ObservedFull.clear();
	for (int32_t i = 0; i < PositionExistence.size(); i++) {
		// expected
		Eigen::VectorXd tmpexp = Eigen::VectorXd::Zero(PositionExistence[i].size());
		int32_t count = 0;
		for (int32_t j = 0; j < PositionExistence[i].size(); j++) {
			if (PositionExistence[i][j]) {
				tmpexp(j) = Expected[i][count];
				count++;
			}
		}
		assert(count == Expected[i].size());
		tmpexp /= tmpexp.sum();
		ExpectedFullNorm.push_back(tmpexp);
		// observed
		Eigen::VectorXd tmpobs = Eigen::VectorXd::Zero(PositionExistence[i].size());
		count = 0;
		for (int32_t j = 0; j < PositionExistence[i].size(); j++) {
			if (PositionExistence[i][j]) {
				tmpobs(j) = Observed[i][count];
				count++;
			}
		}
		assert(count == Observed[i].size());
		ObservedFull.push_back(tmpobs);
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Initializing: step 2 calculating full exp and obs in new coordinate. Finished\n";

	// binning distributions
	ExpectedBinNorm.clear();
	ObservedBin.clear();
	for (int32_t i = 0; i < PositionExistence.size(); i++) {
		int32_t oldlen = ExpectedFullNorm[i].size();
		int32_t newlen = int32_t(ceil(1.0 * oldlen / Bin_Size));
		Eigen::VectorXd tmpexp = Eigen::VectorXd::Zero(newlen);
		Eigen::VectorXd tmpobs = Eigen::VectorXd::Zero(newlen);
		for (int32_t j = 0; j < newlen; j++){
			tmpexp(j) = ExpectedFullNorm[i].segment(Bin_Size*j, min(Bin_Size,oldlen-Bin_Size*j)).sum();
			tmpobs(j) = ObservedFull[i].segment(Bin_Size*j, min(Bin_Size,oldlen-Bin_Size*j)).sum();
		}
		ExpectedBinNorm.push_back(tmpexp);
		ObservedBin.push_back(tmpobs);
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Initializing: step 3 binning exp and obs. Finished\n";

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finished Initializing LP reassignment.\n";
};


void LPReassign_t::InitializeJunction(const vector<Transcript_t>& Transcripts, const vector< vector<double> >& junc_obs)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Adding junction support."<<endl;

	// clear variables
	Gene_Junctions.clear();
	JunctionExistence.clear();
	for (int32_t i = 0; i < TransIndex.size(); i++) {
		Eigen::VectorXi tmp;
		JunctionExistence.push_back(tmp);
	}
	JunctionObserveFull.clear();
	JunctionObserveBin.clear();
	JunctionRelevance.clear();
	Original_junction_trans_full.clear();
	Original_junction_trans_bin.clear();

	// convert the observed to gene-level coordinate in full
	assert( junc_obs.size() == TransIndex.size() );
	for (int32_t i = 0; i < junc_obs.size(); i++) {
		Eigen::VectorXd tmpobs = Eigen::VectorXd::Zero(ObservedFull[i].size());
		int32_t count = 0;
		for (int32_t j = 0; j < PositionExistence[i].size(); j++) {
			if (PositionExistence[i][j]) {
				tmpobs(j) = junc_obs[i][count];
				count++;
			}
		}
		Original_junction_trans_full.push_back( tmpobs );
	}
	// binning for Original_junction_trans_bin
	for (int32_t i = 0; i < PositionExistence.size(); i++) {
		int32_t oldlen = Original_junction_trans_full[i].size();
		int32_t newlen = int32_t(ceil(1.0 * oldlen / Bin_Size));
		Eigen::VectorXd tmpobs = Eigen::VectorXd::Zero(newlen);
		for (int32_t j = 0; j < newlen; j++)
			tmpobs(j) = Original_junction_trans_full[i].segment(Bin_Size*j, min(Bin_Size,oldlen-Bin_Size*j)).sum();
		Original_junction_trans_bin.push_back(tmpobs);
	}

	// listing out junctions
	for (map<string, vector<int32_t> >::const_iterator it = GeneTransMap.cbegin(); it != GeneTransMap.cend(); it++) {
		const string& g = it->first;
		const vector<int32_t>& tids = it->second;
		vector<Junction_t> junctions;
		// loop over each transcript to add the junctions
		for (const int32_t& tid : tids) {
			const Transcript_t& t = Transcripts[tid];
			if (t.Strand) {
				for (int32_t i = 1; i < t.Exons.size(); i++) {
					Junction_t tmp(t.Chr, t.Exons[i-1].EndPos, t.Exons[i].StartPos);
					junctions.push_back(tmp);
				}
			}
			else {
				for (int32_t i = 1; i < t.Exons.size(); i++) {
					Junction_t tmp(t.Chr, t.Exons[i].EndPos, t.Exons[i-1].StartPos);
					junctions.push_back(tmp);
				}
			}
		}
		// find unique junctions
		sort(junctions.begin(), junctions.end());
		junctions.resize( distance(junctions.begin(), unique(junctions.begin(), junctions.end()) ) );
		junctions.reserve( junctions.size() );
		// update the gene to junction map
		Gene_Junctions[g] = junctions;
		// update the junction existence matrix of each transcript
		Eigen::VectorXi sumtmp = Eigen::VectorXi::Zero(junctions.size());
		for (const int32_t& tid : tids) {
			vector<int32_t> this_junction_existence(junctions.size(), 0);
			const Transcript_t& t = Transcripts[tid];
			if (t.Strand) {
				for (int32_t i = 1; i < t.Exons.size(); i++) {
					Junction_t tmp(t.Chr, t.Exons[i-1].EndPos, t.Exons[i].StartPos);
					for (int32_t j = 0; j < junctions.size(); j++) {
						if (tmp == junctions[j]) {
							this_junction_existence[j] = 1;
							break;
						}
					}
				}
			}
			else {
				for (int32_t i = 1; i < t.Exons.size(); i++) {
					Junction_t tmp(t.Chr, t.Exons[i].EndPos, t.Exons[i-1].StartPos);
					for (int32_t j = 0; j < junctions.size(); j++) {
						if (tmp == junctions[j]) {
							this_junction_existence[j] = 1;
							break;
						}
					}
				}
			}
			Eigen::VectorXi tmp = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(this_junction_existence.data(), this_junction_existence.size());
			JunctionExistence[tid] = tmp;
			sumtmp += tmp;
		}
		for (int32_t i = 0; i < sumtmp.size(); i++) {
			assert(sumtmp[i] > 0);
		}
	}

	// based on the number of junctions per gene, initialize junction observed
	for (map< string, vector<Junction_t> >::const_iterator it = Gene_Junctions.cbegin(); it != Gene_Junctions.cend(); it++) {
		const string& g = it->first;
		vector< Eigen::VectorXd > v_obs;
		int32_t tid = GeneTransMap[g].front();
		// for JunctionObserveFull
		int32_t npos = ExpectedFullNorm[tid].size();
		for (int32_t i = 0; i < (it->second).size(); i++) {
			Eigen::VectorXd tmp = Eigen::VectorXd::Zero(npos);
			v_obs.push_back(tmp);
		}
		JunctionObserveFull[g] = v_obs;
	}
	// loop over obs per transcript, and add to obs per junction
	for (int32_t i = 0; i < junc_obs.size(); i++) {
		// skip the records with no coverage
		const vector<double>& this_obs = junc_obs[i];
		double sum = std::accumulate(this_obs.cbegin(), this_obs.cend(), 0.0);
		if (sum < 1e-8)
			continue;
		// for the ones with coverage, loop over each exon to see whether the succeeding junction has read support
		string g = TransGeneMap[i];
		const Transcript_t& t = Transcripts[i];
		const vector<Junction_t>& junctions = Gene_Junctions[g];
		// variables to update
		vector< Eigen::VectorXd >& obs_full = JunctionObserveFull[g];
		assert(i == TransIndex.at(t.TransID));
		int32_t covered_length = 0;
		for (vector<Exon_t>::const_iterator itexon = t.Exons.cbegin(); itexon != t.Exons.cend(); itexon++) {
			sum = std::accumulate(this_obs.cbegin() + covered_length, this_obs.cbegin() + covered_length + itexon->EndPos - itexon->StartPos, 0.0);
			int32_t end_length = covered_length + itexon->EndPos - itexon->StartPos;
			if (sum > 1e-8) {
				assert( itexon + 1 != t.Exons.cend() );
				// find the index of junction in the junction list
				Junction_t tmp;
				if (t.Strand)
					tmp.Init(t.Chr, itexon->EndPos, (itexon+1)->StartPos);
				else
					tmp.Init(t.Chr, (itexon+1)->EndPos, itexon->StartPos);
				for (int32_t j = 0; j < junctions.size(); j++) {
					if (tmp == junctions[j]) {
						int32_t count = 0;
						for (int32_t pos = 0; pos < PositionExistence[i].size(); pos++) {
							if (PositionExistence[i][pos]) {
								if (count >= covered_length && count < end_length)
									obs_full[j][pos] += this_obs[count];
								count ++;
							}
						}
						break;
					}
				}
			}
			covered_length += itexon->EndPos - itexon->StartPos;
		}
	}

	// bin JunctionObserveFull to get JunctionObserveBin
	for (map< string, vector< Eigen::VectorXd > >::const_iterator it = JunctionObserveFull.cbegin(); it != JunctionObserveFull.cend(); it++) {
		const string& g = it->first;
		const vector< Eigen::VectorXd >& obs_full = it->second;
		int32_t oldlen = 0;
		int32_t newlen = 0;
		if (obs_full.size() > 0) {
			oldlen = obs_full[0].size();
			newlen = int32_t(ceil(1.0 * oldlen / Bin_Size));
		}
		// new result
		vector< Eigen::VectorXd> obs_bin;
		for (int32_t i = 0; i < obs_full.size(); i++) {
			assert(obs_full[i].size() == oldlen);
			Eigen::VectorXd tmp = Eigen::VectorXd::Zero(newlen);
			for (int32_t j = 0; j < newlen; j++){
				tmp(j) = obs_full[i].segment(Bin_Size*j, min(Bin_Size,oldlen-Bin_Size*j)).sum();
			}
			obs_bin.push_back(tmp);
		}
		JunctionObserveBin[g] = obs_bin;
	}

	// indicate relevant positions: the bin that contain the junction plus the bins with nonzero coverage
	for (map< string, vector<Junction_t> >::const_iterator it = Gene_Junctions.cbegin(); it != Gene_Junctions.cend(); it++) {
		const string& g = it->first;
		const vector<Junction_t>& junctions = it->second;
		const vector< Eigen::VectorXd >& obs_bin = JunctionObserveBin[g];
		const vector<int32_t> tids = GeneTransMap[g];
		// result variable
		vector< Eigen::VectorXi > relevance;
		for (int32_t i = 0; i < junctions.size(); i++) {
			// find the tid that contain the junction
			int32_t tid = -1;
			for (int32_t j = 0; j < tids.size(); j++) {
				assert(JunctionExistence[tids[j]].size() == junctions.size());
				if (JunctionExistence[tids[j]][i]) {
					tid = tids[j];
					break;
				}
			}
			assert(tid != -1);
			// retrieve the trancript info of tid
			const Transcript_t& t = Transcripts[tid];
			// find the corresponding junction position in transcript tid
			int32_t pos_in_trans = -1;
			int32_t covered_length = 0;
			for (vector<Exon_t>::const_iterator itexon = t.Exons.cbegin(); itexon != t.Exons.cend(); itexon++) {
				if (itexon + 1 == t.Exons.cend())
					break;
				Junction_t tmp;
				if (t.Strand)
					tmp.Init(t.Chr, itexon->EndPos, (itexon+1)->StartPos);
				else
					tmp.Init(t.Chr, (itexon+1)->EndPos, itexon->StartPos);
				if (tmp == junctions[i])
					pos_in_trans = covered_length + itexon->EndPos - itexon->StartPos;
				covered_length += itexon->EndPos - itexon->StartPos;
			}
			assert(pos_in_trans != -1);
			// convert to gene-level coordinate
			int32_t pos_in_gene = 0;
			int32_t count = 0;
			for (; pos_in_gene < PositionExistence[tid].size(); pos_in_gene++) {
				if (PositionExistence[tid][pos_in_gene])
					count ++;
				if (count == pos_in_trans)
					break;
			}
			pos_in_gene++;
			assert(pos_in_gene < PositionExistence[tid].size());
			// convert to binned position
			int32_t pos_in_bin = pos_in_gene / Bin_Size;
			Eigen::VectorXi tmp = Eigen::VectorXi::Zero(obs_bin[0].size());
			tmp[pos_in_bin] = 1;
			if (pos_in_bin > 0)
				tmp[pos_in_bin - 1] = 1;
			relevance.push_back(tmp);
		}
		JunctionRelevance[g] = relevance;
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finished adding junction support.\n";
};


LPReassign_t::LPReassign_t(const map<string,int32_t>& _TransIndex, const map<string,string>& _TransGeneMap, 
	const map<string, vector<string> >& _GeneTransMap, const vector<Transcript_t>& Transcripts, 
	const vector< vector<double> >& Expected, const vector< vector<double> >& Observed, int32_t _Bin_Size, int32_t _Num_Threads)
{
	Bin_Size = _Bin_Size;
	Num_Threads = _Num_Threads;
	LPReassign_t(_TransIndex, _TransGeneMap, _GeneTransMap, Transcripts, Expected, Observed);
};

//****************** test *********************
#ifdef HAVE_GUROBI
vector<double> Quantify_singlecase_grb(Eigen::MatrixXd& exp, Eigen::VectorXd& obs)
{
	assert(exp.rows() == obs.size());
	// create GUROBI environment
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	// silence std out parameter
	model.set(GRB_IntParam_OutputFlag, 0);
	// add variables
	vector<GRBVar> X;
	vector<GRBVar> C;
	for (int32_t i = 0; i < exp.cols(); i++) {
		// lower bound, upper bound, objective coefficient, type, name
		GRBVar x = model.addVar(0, obs.sum(), 0, GRB_CONTINUOUS, "x"+to_string(i));
		X.push_back(x);
	}
	for (int32_t i = 0; i < exp.rows(); i++) {
		GRBVar c = model.addVar(0, obs.sum(), 1, GRB_CONTINUOUS, "c"+to_string(i));
		C.push_back(c);
	}
	// set optimization direction to minimize
	model.set(GRB_IntAttr_ModelSense, 1);
	model.update();
	// set up constraints
	// |  T_{d-by-n}  -E_{d-by-d} |   |X|   <=   |  Y_{d-by-1}  |
	//                                |C|
	for (int32_t i = 0; i < exp.rows(); i++) {
		GRBLinExpr obj = 0.0;
		for (int32_t j = 0; j < exp.cols(); j++)
			obj += exp(i,j)*X[j];
		obj += -1*C[i];
		// lhs expression, sense, rhs value, name
		model.addConstr(obj, GRB_LESS_EQUAL, obs(i), "const"+to_string(i));
	}
	// | -T_{d-by-n}  -E_{d-by-d} |   |X|   <=   | -Y_{d-by-1}  |
	//                                |C|
	for (int32_t i = 0; i < exp.rows(); i++) {
		GRBLinExpr obj = 0.0;
		for (int32_t j = 0; j < exp.cols(); j++)
			obj += -exp(i,j)*X[j];
		obj += -1*C[i];
		// lhs expression, sense, rhs value, name
		model.addConstr(obj, GRB_LESS_EQUAL, -obs(i), "const"+to_string(i+exp.rows()));
	}
	// |    -E_{(n+d)-by-(n+d)}   |   |X|   <=   |0_{(n+d)-by-1}|
	//                                |C|
	// this can be ignored, since non-negative contraints are specified in addVar part
	// solve the LP model
	model.optimize();
	// collect quantification alpha
	vector<double> alpha(exp.cols(), 0);
	for (int32_t i = 0; i < exp.cols(); i++)
		alpha[i] = X[i].get(GRB_DoubleAttr_X);
	// for 0 alphas, assign to a very small positive value, in order to assign position-wise number of reads
	for (int32_t i = 0; i < exp.cols(); i++) {
		assert( alpha[i] >= -1e-4);
		if (alpha[i] <= 0)
			alpha[i] = 1e-8;
	}
	return alpha;
};
#endif

#ifdef HAVE_CLP
vector<double> Quantify_singlecase_clp(Eigen::MatrixXd& exp, Eigen::VectorXd& obs) {
  throw std::runtime_error("CLP not yet implemented");
}
#endif

vector<double> LPReassign_t::Quantify_singlecase(Eigen::MatrixXd& exp, Eigen::VectorXd& obs) {
  #if defined HAVE_GUROBI
  return Quantify_singlecase_grb(exp, obs);
  #elif defined HAVE_CLP
  return Quantify_singlecase_clp(exp, obs);
  #else
  throw std::runtime_error("No linear solver available");
  #endif
}

#ifdef HAVE_GUROBI
vector<double> Quantify_singlecase_junction_grb(Eigen::MatrixXd& exp, Eigen::VectorXd& obs, const vector< Eigen::VectorXi >& relevance, 
	const vector< Eigen::VectorXd >& obs_bin, const vector< Eigen::VectorXi >& existence)
{
	assert(exp.rows() == obs.size());
	// create GUROBI environment
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	// silence std out parameter
	model.set(GRB_IntParam_OutputFlag, 0);
	// add variables
	vector<GRBVar> X;
	vector<GRBVar> C;
	for (int32_t i = 0; i < exp.cols(); i++) {
		// lower bound, upper bound, objective coefficient, type, name
		GRBVar x = model.addVar(0, obs.sum(), 0, GRB_CONTINUOUS, "x"+to_string(i));
		X.push_back(x);
	}
	for (int32_t i = 0; i < exp.rows(); i++) {
		GRBVar c = model.addVar(0, obs.sum(), 1, GRB_CONTINUOUS, "c"+to_string(i));
		C.push_back(c);
	}
	// variable for junctions
	vector< vector<GRBVar> > C_junction;
	for (int32_t i = 0; i < obs_bin.size(); i++) {
		vector<GRBVar> tmp;
		for (int32_t j = 0; j < exp.rows(); j++) {
			GRBVar c = model.addVar(0, obs.sum(), 1, GRB_CONTINUOUS, "c_junction_"+to_string(i)+"_"+to_string(j));
			tmp.push_back(c);
		}
		C_junction.push_back(tmp);
	}
	// set optimization direction to minimize
	model.set(GRB_IntAttr_ModelSense, 1);
	model.update();
	// set up constraints
	// |  T_{d-by-n}  -E_{d-by-d} |   |X|   <=   |  Y_{d-by-1}  |
	//                                |C|
	for (int32_t i = 0; i < exp.rows(); i++) {
		GRBLinExpr obj = 0.0;
		for (int32_t j = 0; j < exp.cols(); j++)
			obj += exp(i,j)*X[j];
		obj += -1*C[i];
		// lhs expression, sense, rhs value, name
		model.addConstr(obj, GRB_LESS_EQUAL, obs(i), "const"+to_string(i));
	}
	// | -T_{d-by-n}  -E_{d-by-d} |   |X|   <=   | -Y_{d-by-1}  |
	//                                |C|
	for (int32_t i = 0; i < exp.rows(); i++) {
		GRBLinExpr obj = 0.0;
		for (int32_t j = 0; j < exp.cols(); j++)
			obj += -exp(i,j)*X[j];
		obj += -1*C[i];
		// lhs expression, sense, rhs value, name
		model.addConstr(obj, GRB_LESS_EQUAL, -obs(i), "const"+to_string(i+exp.rows()));
	}
	// comparing the linear combination of expected distribution to the support of each junction
	for (int32_t i = 0; i < obs_bin.size(); i++) {
		// | relevance diagonal * T_{d-by-n} * existence diagonal  -E_{d-by-d} |    |       X      |    <=     |obs_bin_{junction i}|
		//                                                                          |C_{junction i}|
		Eigen::MatrixXd diag_relevance = Eigen::MatrixXd::Zero(relevance[i].size(), relevance[i].size());
		for (int32_t j = 0; j < relevance[i].size(); j++)
			diag_relevance(j,j) = relevance[i](j);
		Eigen::MatrixXd diag_existence = Eigen::MatrixXd::Zero(exp.cols(), exp.cols());
		assert( existence.size() == exp.cols() );
		for (int32_t j = 0; j < existence.size(); j++)
			diag_existence(j,j) = existence[j](i);
		assert(diag_relevance.sum() > 0);
		Eigen::MatrixXd m = diag_relevance * exp * diag_existence;
		assert(m.rows() == exp.rows() && m.cols() == exp.cols());
		for (int32_t j = 0; j < m.rows(); j++) {
			GRBLinExpr obj = 0.0;
			for (int32_t k = 0; k < m.cols(); k++)
				obj += m(j, k) * X[k];
			obj += -1 * C_junction[i][j];
			// lhs expression, sense, rhs value, name
			model.addConstr(obj, GRB_LESS_EQUAL, obs_bin[i](j), "const_junction_"+to_string(i)+"_"+to_string(j));
		}
		// | -relevance diagonal * T_{d-by-n} * existence diagonal  -E_{d-by-d} |    |       X      |    <=     |-obs_bin_{junction i}|
		//                                                                           |C_{junction i}|
		for (int32_t j = 0; j < m.rows(); j++) {
			GRBLinExpr obj = 0.0;
			for (int32_t k = 0; k < m.cols(); k++)
				obj += -m(j, k) * X[k];
			obj += -1 * C_junction[i][j];
			// lhs expression, sense, rhs value, name
			model.addConstr(obj, GRB_LESS_EQUAL, -obs_bin[i](j), "const_junction_"+to_string(i)+"_"+to_string(j+exp.rows()));
		}
	}
	// solve the LP model
	model.optimize();
	// collect quantification alpha
	vector<double> alpha(exp.cols(), 0);
	for (int32_t i = 0; i < exp.cols(); i++)
		alpha[i] = X[i].get(GRB_DoubleAttr_X);
	// for 0 alphas, assign to a very small positive value, in order to assign position-wise number of reads
	for (int32_t i = 0; i < exp.cols(); i++) {
		assert( alpha[i] >= -1e-4);
		if (alpha[i] <= 0)
			alpha[i] = 1e-8;
	}
	return alpha;
};
#endif

#ifdef HAVE_CLP
vector<double> Quantify_singlecase_junction_clp(Eigen::MatrixXd& exp, Eigen::VectorXd& obs, const vector< Eigen::VectorXi >& relevance, 
                                                const vector< Eigen::VectorXd >& obs_bin, const vector< Eigen::VectorXi >& existence) {
  throw std::runtime_error("CLP not yet implemented");
}
#endif

vector<double> LPReassign_t::Quantify_singlecase_junction(Eigen::MatrixXd& exp, Eigen::VectorXd& obs, const vector< Eigen::VectorXi >& relevance, 
                                                          const vector< Eigen::VectorXd >& obs_bin, const vector< Eigen::VectorXi >& existence) {
#if defined HAVE_GUROBI
  return Quantify_singlecase_junction_grb(exp, obs, relevance, obs_bin, existence);
#elif defined HAVE_CLP
  return Quantify_singlecase_junction_clp(exp, obs, relevance, obs_bin, existence);
#else
  throw std::runtime_error("No linear solver available");
#endif
}

//****************** end test *********************


/* vector<double> LPReassign_t::Quantify_singlecase(Eigen::MatrixXd& exp, Eigen::VectorXd& obs)
{
	assert(exp.rows() == obs.size());
	// create LP problem in glpk
	glp_prob *lp=glp_create_prob();
	glp_set_prob_name(lp, "Reassign");
	glp_set_obj_dir(lp, GLP_MIN);
	// d positions, n transcripts
	// constraints in matrix  |  T_{d-by-n}  -E_{d-by-d} |   |X|      |  Y_{d-by-1}  |
	//                        | -T_{d-by-n}  -E_{d-by-d} | * |C|  <=  | -Y_{d-by-1}  |
	//                        |    -E_{(n+d)-by-(n+d)}   |            |0_{(n+d)-by-1}|
	glp_add_cols(lp, exp.cols()+exp.rows());
	for(int32_t i = 0; i < exp.cols(); i++) {
		glp_set_col_name(lp, i+1, ("x"+to_string(i)).c_str());
		glp_set_col_bnds(lp, i+1, GLP_LO, 0, 100);
		glp_set_obj_coef(lp, i+1, 0);
		glp_set_col_kind(lp, i+1, GLP_CV);
	}
	for(int32_t i = 0; i < exp.rows(); i++) {
		glp_set_col_name(lp, i+1+exp.cols(), ("c"+to_string(i)).c_str());
		glp_set_col_bnds(lp, i+1+exp.cols(), GLP_LO, 0, 100);
		glp_set_obj_coef(lp, i+1+exp.cols(), 1);
		glp_set_col_kind(lp, i+1+exp.cols(), GLP_CV);
	}
	// adding constraints
	// matrix coefficient of non-zero elements
	int* ia = new int[2*exp.rows()*exp.cols()+3*exp.rows()+exp.cols()+1];
	int* ja = new int[2*exp.rows()*exp.cols()+3*exp.rows()+exp.cols()+1];
	double* ar = new double[2*exp.rows()*exp.cols()+3*exp.rows()+exp.cols()+1];
	// specifying values of non-zero element of matrix coefficient
	// two sub-matrix T and -T
	int32_t count = 0;
	for(int32_t i = 0; i < exp.rows(); i++)
		for(int32_t j = 0; j < exp.cols(); j++) {
			ia[count + 1] = i + 1; ja[count + 1] = j + 1; ar[count + 1] = exp(i,j);
			ia[count + 2] = i + exp.rows() + 1; ja[count + 2] = j + 1; ar[count + 2] = -exp(i,j);
			count += 2;
		}
	// two matrix -E _{d-by-d} and -E_{d-by-d}
	for(int32_t i = 0; i < exp.rows(); i++) {
		ia[count + 1] = i + 1; ja[count + 1] = i + exp.cols() + 1; ar[count + 1] = -1;
		ia[count + 2] = i + exp.rows() + 1; ja[count + 2] = i + exp.cols() + 1; ar[count + 2] = -1;
		count += 2;
	}
	// the bottom matrix -E_{(n+d)-by-(n+d)}
	for(int32_t i = 0; i < exp.rows()+exp.cols(); i++) {
		ia[count + 1] = i + 2*exp.rows() + 1; ja[count + 1] = i + 1; ar[count + 1] = -1;
		count ++;
	}
	// adding constraints
	// bound of each constraints
	glp_add_rows(lp, 3*exp.rows()+exp.cols());
	for(int32_t i = 0; i < exp.rows(); i++){
		glp_set_row_name(lp, i+1, ("c"+to_string(i+1)).c_str());
		glp_set_row_bnds(lp, i+1, GLP_UP, -1, obs(i));
		glp_set_row_name(lp, i+exp.rows()+1, ("c"+to_string(i+exp.rows()+1)).c_str());
		glp_set_row_bnds(lp, i+exp.rows()+1, GLP_UP, -1, -obs(i));
	}
	for(int32_t i = 0; i < exp.rows()+exp.cols(); i++) {
		glp_set_row_name(lp, i+2*exp.rows()+1, ("c"+to_string(i+2*exp.rows()+1)).c_str());
		glp_set_row_bnds(lp, i+2*exp.rows()+1, GLP_UP, -1, 0);
	}
	// load matrix
	glp_load_matrix(lp, count, ia, ja, ar);
	// solve with parameters
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;
	// parm.presolve = GLP_ON;
	int err = glp_simplex(lp, &parm);
	// prepare result
	if (err != 0){
		cout << "GLP_EBADB\t"<< (err==GLP_EBADB) << endl;
		cout << "GLP_ESING\t"<< (err==GLP_ESING) << endl;
		cout << "GLP_ECOND\t"<< (err==GLP_ECOND) << endl;
		cout << "GLP_EBOUND\t"<< (err==GLP_EBOUND) << endl;
		cout << "GLP_EFAIL\t"<< (err==GLP_EFAIL) << endl;
		cout << "GLP_EOBJLL\t"<< (err==GLP_EOBJLL) << endl;
		cout << "GLP_EOBJUL\t"<< (err==GLP_EOBJUL) << endl;
		cout << "GLP_EITLIM\t"<< (err==GLP_EITLIM) << endl;
		cout << "GLP_ETMLIM\t"<< (err==GLP_ETMLIM) << endl;
		cout << "GLP_ENOPFS\t"<< (err==GLP_ENOPFS) << endl;
		cout << "GLP_ENODFS\t"<< (err==GLP_ENODFS) << endl;
		cout << "Error: LP simplex method didn't finish successfully.\n";
		ofstream output1("test_lpfailure_exp", ios::out);
		ofstream output2("test_lpfailure_obs", ios::out);
		for (int32_t i = 0; i < exp.rows(); i++){
			for (int32_t j = 0; j < exp.cols()-1; j++)
				output1 << exp(i,j) <<" ";
			output1 << exp(i, exp.cols()-1) << endl;
			output2 << obs(i) << endl;
		}
		output1.close();
		output2.close();
	}
	assert(err == 0);

	// random output for testing
	ofstream output1("test_lppass_exp", ios::out);
	ofstream output2("test_lppass_obs", ios::out);
	for (int32_t i = 0; i < exp.rows(); i++){
		for (int32_t j = 0; j < exp.cols()-1; j++)
			output1 << exp(i,j) <<" ";
		output1 << exp(i, exp.cols()-1) << endl;
		output2 << obs(i) << endl;
	}
	output1.close();
	output2.close();

	vector<double> alpha;
	for (int32_t i = 0; i < exp.cols(); i++)
		alpha.push_back( max(0.0, glp_get_col_prim(lp, i+1)) );

	double obj = glp_get_obj_val(lp);
	return alpha;
};*/


vector<double> LPReassign_t::ReassignReads(vector< vector<double> >& newAssignment, map< string,vector<MovingRead_t> >& MovingMat)
{
	// clear result variable
	newAssignment.clear();
	for (int32_t i = 0; i < AdjustmentList.size(); i++){
		vector<double> tmp;
		newAssignment.push_back(tmp);
	}
	MovingMat.clear();
	vector<double> LPExp;
	for (int32_t i = 0; i < ObservedBin.size(); i++)
		LPExp.push_back(ObservedBin[i].sum());
	// temprary indexes according to AdjustmentList
	map<int32_t,int32_t> TempIndex;
	for(int32_t i = 0; i < AdjustmentList.size(); i++)
		TempIndex[AdjustmentList[i]] = i;
	// involved genes
	vector<string> involved_genes;
	for(int32_t i = 0; i < AdjustmentList.size(); i++)
		involved_genes.push_back(TransGeneMap[AdjustmentList[i]]);
	sort(involved_genes.begin(), involved_genes.end());
	involved_genes.resize( distance(involved_genes.begin(), unique(involved_genes.begin(),involved_genes.end())) );
	// for each involved gene, reassign reads
	for (const string& g : involved_genes){
		// collect intersection between isoforms of the gene and AdjustmentList
		vector<int32_t> tids;
		for(int32_t& t : GeneTransMap[g]){
			map<int32_t,int32_t>::const_iterator ittempidx = TempIndex.find(t);
			if (ittempidx != TempIndex.cend())
				tids.push_back(t);
		}
		assert(tids.size() > 0);

		int32_t npos = ExpectedBinNorm[tids[0]].size();
		Eigen::MatrixXd exp(npos, tids.size());
		for(int32_t i = 0; i < tids.size(); i++)
			exp.col(i) = ExpectedBinNorm[tids[i]];
		Eigen::VectorXd obs = Eigen::VectorXd::Zero(npos);
		for(int32_t i = 0; i < tids.size(); i++)
			obs += ObservedBin[tids[i]];
		// quantify
		vector< Eigen::VectorXd > obs_junction_full = JunctionObserveFull[g];
		vector< Eigen::VectorXd > obs_junction_bin = JunctionObserveBin[g];
		// remove the junction reads from transcripts not involved in reassignment
		for (int32_t& t : GeneTransMap[g]) {
			map<int32_t,int32_t>::const_iterator ittempidx = TempIndex.find(t);
			if (ittempidx == TempIndex.cend()) {
				for (int32_t i = 0; i < obs_junction_full.size(); i++) {
					obs_junction_full[i] -= Original_junction_trans_full[t];
					obs_junction_bin[i] -= Original_junction_trans_bin[t];
					// since Original_junction_trans_full contain the junction reads of all junctions, but obs_junction_full[i] only contains 1 junction, the subtraction may lead to negative values of the other junctions
					// set the negative junction coverage to 0
					for (int32_t pos = 0; pos < obs_junction_full[i].size(); pos++) {
						if (obs_junction_full[i](pos) < 0)
							obs_junction_full[i](pos) = 0;
					}
					for (int32_t pos = 0; pos < obs_junction_bin[i].size(); pos++) {
						if (obs_junction_bin[i](pos) < 0)
							obs_junction_bin[i](pos) = 0;
					}
				}
			}
		}
		vector< Eigen::VectorXi > relevance = JunctionRelevance[g];
		vector< Eigen::VectorXi > existence;
		for (int32_t i = 0; i < tids.size(); i++)
			existence.push_back( JunctionExistence[tids[i]] );
		vector<double> alpha = Quantify_singlecase_junction(exp, obs, relevance, obs_junction_bin, existence);
		// re-assign on base-pair level
		Eigen::MatrixXd expfull(ExpectedFullNorm[tids[0]].size(), tids.size());
		for(int32_t i = 0; i < tids.size(); i++)
			expfull.col(i) = ExpectedFullNorm[tids[i]];
		Eigen::VectorXd obsfull = Eigen::VectorXd::Zero(ExpectedFullNorm[tids[0]].size());
		for(int32_t i = 0; i < tids.size(); i++)
			obsfull += ObservedFull[tids[i]];
		for (int32_t i = 0; i < obs_junction_full.size(); i++)
			obsfull -= obs_junction_full[i];
		// sanity check that obsfull is non-negative
		for (int32_t i = 0; i < obsfull.size(); i++) {
			if (obsfull[i] < -1e-8)
				cout << "watch here\n";
			assert(obsfull(i) > -1e-8);
		}
		// reads that don't span junctions
		Eigen::MatrixXd tmp_assign_full = Eigen::MatrixXd::Zero(expfull.rows(), expfull.cols());
		for (int32_t i = 0; i < obsfull.size(); i++) {
			Eigen::VectorXd tmpalpha = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(alpha.data(), alpha.size());
			Eigen::ArrayXd assign_theo = expfull.row(i).array().transpose() * tmpalpha.array();
			if (assign_theo.sum() != 0)
				tmp_assign_full.row(i) = obsfull(i) * assign_theo / (assign_theo.sum());
		}
		// junction reads
		for (int32_t i = 0; i < obs_junction_full.size(); i++) {
			Eigen::VectorXd tmpalpha = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(alpha.data(), alpha.size());
			Eigen::VectorXd this_exist = Eigen::VectorXd::Zero(expfull.cols());
			for (int32_t k = 0; k < expfull.cols(); k++)
				this_exist(k) = existence[k](i);
			for (int32_t pos = 0; pos < expfull.rows(); pos++) {
				Eigen::ArrayXd assign_theo = expfull.row(pos).array().transpose() * tmpalpha.array() * this_exist.array();
				// assign_theo = assign_theo.transpose() * this_exist.array();
				if (assign_theo.sum() != 0)
					tmp_assign_full.row(pos) += obs_junction_full[i](pos) * assign_theo.matrix() / (assign_theo.sum());
			}
		}
		// sanity check non-negative
		for (int32_t i = 0; i < tmp_assign_full.rows(); i++)
			for (int32_t j = 0; j < tmp_assign_full.cols(); j++) {
				if (tmp_assign_full(i,j) < -1e-8)
					cout << "watch here\n";
				assert(tmp_assign_full(i,j) > -1e-8);
				if (tmp_assign_full(i,j) < 0)
					tmp_assign_full(i,j) = 0;
			}
		// end test

		// Eigen::MatrixXd tmp_assign_full = Eigen::MatrixXd::Zero(expfull.rows(), expfull.cols());
		// for(int32_t i = 0; i < obsfull.size(); i++) {
		// 	Eigen::VectorXd tmpalpha = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(alpha.data(), alpha.size());
		// 	Eigen::ArrayXd assign_theo = expfull.row(i).array().transpose() * tmpalpha.array();
		// 	assert( obsfull[i] >= 0 );
		// 	for (int32_t j = 0; j < assign_theo.size(); j++)
		// 		assert( assign_theo[j] >= 0);
		// 	if (assign_theo.sum() != 0)
		// 		tmp_assign_full.row(i) = obsfull(i) * assign_theo / (assign_theo.sum());
		// }
		// convert to separate transcript index based on PositionExistence matrix
		for(int32_t i = 0; i < tids.size(); i++){
			vector<double> tmp;
			assert(tmp_assign_full.rows() == PositionExistence[tids[i]].size());
			for (int32_t j = 0; j < tmp_assign_full.rows(); j++)
				if (PositionExistence[tids[i]][j])
					tmp.push_back(tmp_assign_full(j,i));
			newAssignment[TempIndex[tids[i]]] = tmp;
			assert(tmp_assign_full.col(i).sum() >= 0);
			LPExp[tids[i]] = tmp_assign_full.col(i).sum();
		}
		// calculate MovingMat
		vector<MovingRead_t> movingreads;
		for (int32_t i = 0; i < tids.size(); i++)
			for (int32_t j = i+1; j < tids.size(); j++) {
				double moved_t1 = 0;
				double moved_t2 = 0;
				for (int32_t k = 0; k < tmp_assign_full.rows(); k++) {
					if (PositionExistence[tids[i]][k] > 0 && PositionExistence[tids[j]][k] > 0) {
						moved_t1 += tmp_assign_full(k, i) - ObservedFull[tids[i]](k);
						moved_t2 += tmp_assign_full(k, j) - ObservedFull[tids[j]](k);
					}
				}
				MovingRead_t tmpmove1(tids[i], tids[j], moved_t1);
				MovingRead_t tmpmove2(tids[j], tids[i], moved_t2);
				movingreads.push_back(tmpmove1);
				movingreads.push_back(tmpmove2);
			}
		sort(movingreads.begin(), movingreads.end());
		MovingMat[g] = movingreads;
	}

	return LPExp;
};


vector<double> LPReassign_t::ReassignReads(vector< vector<double> >& newAssignment, map< string,vector<MovingRead_t> >& MovingMat, const vector<int32_t>& _AdjustmentList)
{
	// clear result variable
	newAssignment.clear();
	for (int32_t i = 0; i < _AdjustmentList.size(); i++){
		vector<double> tmp;
		newAssignment.push_back(tmp);
	}
	MovingMat.clear();
	vector<double> LPExp;
	for (int32_t i = 0; i < ObservedBin.size(); i++)
		LPExp.push_back(ObservedBin[i].sum());
	// temprary indexes according to _AdjustmentList
	map<int32_t,int32_t> TempIndex;
	for(int32_t i = 0; i < _AdjustmentList.size(); i++)
		TempIndex[_AdjustmentList[i]] = i;
	// involved genes
	vector<string> involved_genes;
	for(int32_t i = 0; i < _AdjustmentList.size(); i++)
		involved_genes.push_back(TransGeneMap[_AdjustmentList[i]]);
	sort(involved_genes.begin(), involved_genes.end());
	involved_genes.resize( distance(involved_genes.begin(), unique(involved_genes.begin(),involved_genes.end())) );
	// for each involved gene, reassign reads
	for (const string& g : involved_genes){
		// collect intersection between isoforms of the gene and _AdjustmentList
		vector<int32_t> tids;
		for(int32_t& t : GeneTransMap[g]){
			map<int32_t,int32_t>::const_iterator ittempidx = TempIndex.find(t);
			if (ittempidx != TempIndex.cend())
				tids.push_back(t);
		}
		assert(tids.size() > 0);

		int32_t npos = ExpectedBinNorm[tids[0]].size();
		Eigen::MatrixXd exp(npos, tids.size());
		for(int32_t i = 0; i < tids.size(); i++)
			exp.col(i) = ExpectedBinNorm[tids[i]];
		Eigen::VectorXd obs = Eigen::VectorXd::Zero(npos);
		for(int32_t i = 0; i < tids.size(); i++)
			obs += ObservedBin[tids[i]];
		// quantify
		vector<double> alpha = Quantify_singlecase(exp, obs);
		// re-assign on base-pair level
		Eigen::MatrixXd expfull(ExpectedFullNorm[tids[0]].size(), tids.size());
		for(int32_t i = 0; i < tids.size(); i++)
			expfull.col(i) = ExpectedFullNorm[tids[i]];
		Eigen::VectorXd obsfull = Eigen::VectorXd::Zero(ExpectedFullNorm[tids[0]].size());
		for(int32_t i = 0; i < tids.size(); i++)
			obsfull += ObservedFull[tids[i]];
		Eigen::MatrixXd tmp_assign_full = Eigen::MatrixXd::Zero(expfull.rows(), expfull.cols());
		for(int32_t i = 0; i < obsfull.size(); i++) {
			Eigen::VectorXd tmpalpha = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(alpha.data(), alpha.size());
			Eigen::ArrayXd assign_theo = expfull.row(i).array().transpose() * tmpalpha.array();
			if (assign_theo.sum() != 0)
				tmp_assign_full.row(i) = obsfull(i) * assign_theo / (assign_theo.sum());
		}
		// convert to separate transcript index based on PositionExistence matrix
		for(int32_t i = 0; i < tids.size(); i++){
			vector<double> tmp;
			assert(tmp_assign_full.rows() == PositionExistence[tids[i]].size());
			for (int32_t j = 0; j < tmp_assign_full.rows(); j++)
				if (PositionExistence[tids[i]][j])
					tmp.push_back(tmp_assign_full(j,i));
			newAssignment[TempIndex[tids[i]]] = tmp;
			LPExp[tids[i]] = tmp_assign_full.col(i).sum();
		}
		// calculate MovingMat
		vector<MovingRead_t> movingreads;
		for (int32_t i = 0; i < tids.size(); i++)
			for (int32_t j = i+1; j < tids.size(); j++) {
				double moved_t1 = 0;
				double moved_t2 = 0;
				for (int32_t k = 0; k < tmp_assign_full.rows(); k++) {
					if (PositionExistence[tids[i]][k] > 0 && PositionExistence[tids[j]][k] > 0) {
						moved_t1 += tmp_assign_full(k, i) - ObservedFull[tids[i]](k);
						moved_t2 += tmp_assign_full(k, j) - ObservedFull[tids[j]](k);
					}
				}
				MovingRead_t tmpmove1(tids[i], tids[j], moved_t1);
				MovingRead_t tmpmove2(tids[j], tids[i], moved_t2);
				movingreads.push_back(tmpmove1);
				movingreads.push_back(tmpmove2);
			}
		sort(movingreads.begin(), movingreads.end());
		MovingMat[g] = movingreads;
	}
	
	return LPExp;
};


vector<int32_t> LPReassign_t::RefineLPTrans_singlecase(const vector<int32_t>& shortAdjustmentList, DistTest_t& dt, 
	const vector<PRegion_t>& PValuesPos_salmon, const vector<PRegion_t>& PValuesNeg_salmon, 
	const vector<PRegion_t>& AdjPValuesPos_salmon, const vector<PRegion_t>& AdjPValuesNeg_salmon, 
	const vector<PRegion_t>& PValuesPos_lp, const vector<PRegion_t>& PValuesNeg_lp, 
	const vector<PRegion_t>& AdjPValuesPos_lp, const vector<PRegion_t>& AdjPValuesNeg_lp, double PvalueThresh)
{
	assert(is_sorted(shortAdjustmentList.cbegin(), shortAdjustmentList.cend()));
	// collect significant transcripts, and all Pvalue information
	vector<int32_t> sig_salmon;
	vector<int32_t> sig_lp;
	vector<PRegion_t> shortPValuesPos_salmon;
	vector<PRegion_t> shortPValuesNeg_salmon;
	PRegion_t test;
	for (vector<PRegion_t>::const_iterator it = AdjPValuesPos_salmon.cbegin(); it != AdjPValuesPos_salmon.cend(); it++){
		if (binary_search(shortAdjustmentList.begin(), shortAdjustmentList.end(), it->TID)) {
			shortPValuesPos_salmon.push_back(PValuesPos_salmon[distance(AdjPValuesPos_salmon.cbegin(), it)]);
			if (it->Pvalue < PvalueThresh)
				sig_salmon.push_back(it->TID);
		}
	}
	for (vector<PRegion_t>::const_iterator it = AdjPValuesNeg_salmon.cbegin(); it != AdjPValuesNeg_salmon.cend(); it++){
		if (binary_search(shortAdjustmentList.begin(), shortAdjustmentList.end(), it->TID)) {
			shortPValuesNeg_salmon.push_back(PValuesNeg_salmon[distance(AdjPValuesNeg_salmon.cbegin(), it)]);
			if (it->Pvalue < PvalueThresh)
				sig_salmon.push_back(it->TID);
		}
	}
	sort(sig_salmon.begin(), sig_salmon.end());
	sig_salmon.resize( distance(sig_salmon.begin(), unique(sig_salmon.begin(),sig_salmon.end())) );
	for (vector<PRegion_t>::const_iterator it = AdjPValuesPos_lp.cbegin(); it != AdjPValuesPos_lp.cend(); it++){
		if (binary_search(shortAdjustmentList.begin(), shortAdjustmentList.end(), it->TID)) {
			if (it->Pvalue < PvalueThresh)
				sig_lp.push_back(it->TID);
		}
	}
	for (vector<PRegion_t>::const_iterator it = AdjPValuesNeg_lp.cbegin(); it != AdjPValuesNeg_lp.cend(); it++){
		if (binary_search(shortAdjustmentList.begin(), shortAdjustmentList.end(), it->TID)) {
			if (it->Pvalue < PvalueThresh)
				sig_lp.push_back(it->TID);
		}
	}
	sort(sig_lp.begin(), sig_lp.end());
	sig_lp.resize( distance(sig_lp.begin(), unique(sig_lp.begin(),sig_lp.end())) );
	int32_t num_improved = sig_salmon.size() - sig_lp.size();
	// remove transcripts from shortAdjustmentList one by one
	// create a copy of shortAdjustmentList
	vector<int32_t> tmpAdjList(shortAdjustmentList.cbegin(), shortAdjustmentList.cend());
	while (tmpAdjList.size() > 2) {
		bool whetherchanged = false;
		for (int32_t i = 0; i < tmpAdjList.size(); i++) {
			if (binary_search(sig_salmon.begin(), sig_salmon.end(), tmpAdjList[i]))
				continue;
			vector<int32_t> removedAdjList;
			for (int32_t j = 0; j < tmpAdjList.size(); j++){
				if (j != i)
					removedAdjList.push_back(tmpAdjList[j]);
			}
			assert(removedAdjList.size() + 1 == tmpAdjList.size());
			assert(is_sorted(removedAdjList.begin(), removedAdjList.end()));
			// LP re-assign
			vector< vector<double> > removedAssignment;
			map< string,vector<MovingRead_t> > removedMovingMat;
			vector<double> removedLPExp = ReassignReads(removedAssignment, removedMovingMat, removedAdjList);
			// calculate P value
			vector<PRegion_t> removedPValuePos;
			vector<PRegion_t> removedPValueNeg;
			dt.UpdateObserved(removedAdjList, removedAssignment);
			dt.CalDeletionScore(removedAdjList);
			dt.PValue_regional(removedAdjList, removedPValuePos, removedPValueNeg, false);
			// adding the rest LP pvalues
			for (int32_t j = 0; j < shortPValuesPos_salmon.size(); j++) {
				if (!binary_search(removedAdjList.begin(), removedAdjList.end(), shortPValuesPos_salmon[j].TID))
					removedPValuePos.push_back(shortPValuesPos_salmon[j]);
			}
			for (int32_t j = 0; j < shortPValuesNeg_salmon.size(); j++) {
				if (!binary_search(removedAdjList.begin(), removedAdjList.end(), shortPValuesNeg_salmon[j].TID))
					removedPValueNeg.push_back(shortPValuesNeg_salmon[j]);
			}
			// use the lower bound of adjusted p value in AdjPValuesPos_lp and AdjPValuesNeg_lp to estimate the adjusted p value of removed list
			// collect significant transcripts with this removal
			vector<int32_t> removed_sig_lp;
			for (vector<PRegion_t>::iterator it = removedPValuePos.begin(); it != removedPValuePos.end(); it++){
				vector<PRegion_t>::const_iterator lb = lower_bound(PValuesPos_lp.cbegin(), PValuesPos_lp.cend(), *it, PRegion_t::CompPvalue);
				double adjpvalue = 1;
				if (lb == PValuesPos_lp.cend())
					assert(PValuesPos_lp.back().Pvalue <= it->Pvalue);
				else{
					assert(lb == PValuesPos_lp.cbegin() || (lb-1)->Pvalue <= it->Pvalue);
					assert((lb+1) == PValuesPos_lp.cend() || (lb+1)->Pvalue >= it->Pvalue);
					assert(distance(PValuesPos_lp.cbegin(), lb) <= AdjPValuesPos_lp.size());
					adjpvalue = AdjPValuesPos_lp[distance(PValuesPos_lp.cbegin(), lb)].Pvalue;
				}
				if (adjpvalue < PvalueThresh)
					removed_sig_lp.push_back(it->TID);
			}
			for (vector<PRegion_t>::iterator it = removedPValueNeg.begin(); it != removedPValueNeg.end(); it++){
				vector<PRegion_t>::const_iterator lb = lower_bound(PValuesNeg_lp.cbegin(), PValuesNeg_lp.cend(), *it, PRegion_t::CompPvalue);
				double adjpvalue = 1;
				if (lb == PValuesNeg_lp.cend())
					assert(PValuesNeg_lp.back().Pvalue <= it->Pvalue);
				else{
					assert(lb == PValuesNeg_lp.cbegin() || (lb-1)->Pvalue <= it->Pvalue);
					assert((lb+1) == PValuesNeg_lp.cend() || (lb+1)->Pvalue >= it->Pvalue);
					assert(distance(PValuesNeg_lp.cbegin(), lb) <= AdjPValuesNeg_lp.size());
					adjpvalue = AdjPValuesNeg_lp[distance(PValuesNeg_lp.cbegin(), lb)].Pvalue;
				}
				if (adjpvalue < PvalueThresh)
					removed_sig_lp.push_back(it->TID);
			}
			sort(removed_sig_lp.begin(), removed_sig_lp.end());
			removed_sig_lp.resize( distance(removed_sig_lp.begin(), unique(removed_sig_lp.begin(), removed_sig_lp.end())) );
			// compare the number of improvement
			if (removed_sig_lp.size() <= sig_lp.size()) {
				tmpAdjList = removedAdjList;
				whetherchanged = true;
				break;
			}
		}
		// if removing any of the transcripts will reduce num_improved, then don't remove anything and exit
		if (!whetherchanged)
			break;
	}

	return tmpAdjList;
};


vector<int32_t> LPReassign_t::RefineLPTrans(DistTest_t& dt, 
	const vector<PRegion_t>& PValuesPos_salmon, const vector<PRegion_t>& PValuesNeg_salmon, 
	const vector<PRegion_t>& AdjPValuesPos_salmon, const vector<PRegion_t>& AdjPValuesNeg_salmon, 
	const vector<PRegion_t>& PValuesPos_lp, const vector<PRegion_t>& PValuesNeg_lp, 
	const vector<PRegion_t>& AdjPValuesPos_lp, const vector<PRegion_t>& AdjPValuesNeg_lp, double PvalueThresh)
{
	// copy the original AdjustmentList
	vector<int32_t> AdjustmentList_copy = AdjustmentList;
	assert(is_sorted(AdjustmentList_copy.begin(), AdjustmentList_copy.end()));
	// create result vector
	vector<int32_t> RefinedAdjustmentList;
	// involved genes
	vector<string> involved_genes;
	for(int32_t i = 0; i < AdjustmentList.size(); i++)
		involved_genes.push_back(TransGeneMap[AdjustmentList[i]]);
	sort(involved_genes.begin(), involved_genes.end());
	involved_genes.resize( distance(involved_genes.begin(), unique(involved_genes.begin(),involved_genes.end())) );
	
	// create copies of P values and sort
	// PValuesPos_salmon
	vector<int32_t> indexes(PValuesPos_salmon.size(), 0);
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesPos_salmon(PValuesPos_salmon.cbegin(), PValuesPos_salmon.cend());
	vector<PRegion_t> copyAdjPValuesPos_salmon(AdjPValuesPos_salmon.cbegin(), AdjPValuesPos_salmon.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesPos_salmon](int32_t a, int32_t b){return PValuesPos_salmon[a].Pvalue < PValuesPos_salmon[b].Pvalue;} );
	reorder(copyPValuesPos_salmon, indexes);
	assert( is_sorted(copyPValuesPos_salmon.cbegin(), copyPValuesPos_salmon.cend(), PRegion_t::CompPvalue) );
	reorder(copyAdjPValuesPos_salmon, indexes);
	// PValuesNeg_salmon
	indexes.resize(PValuesNeg_salmon.size());
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesNeg_salmon(PValuesNeg_salmon.cbegin(), PValuesNeg_salmon.cend());
	vector<PRegion_t> copyAdjPValuesNeg_salmon(AdjPValuesNeg_salmon.cbegin(), AdjPValuesNeg_salmon.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesNeg_salmon](int32_t a, int32_t b){return PValuesNeg_salmon[a].Pvalue < PValuesNeg_salmon[b].Pvalue;} );
	reorder(copyPValuesNeg_salmon, indexes);
	assert( is_sorted(copyPValuesNeg_salmon.cbegin(), copyPValuesNeg_salmon.cend(), PRegion_t::CompPvalue) );
	reorder(copyAdjPValuesNeg_salmon, indexes);
	// PValuesPos_lp
	indexes.resize(PValuesPos_lp.size());
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesPos_lp(PValuesPos_lp.cbegin(), PValuesPos_lp.cend());
	vector<PRegion_t> copyAdjPValuesPos_lp(AdjPValuesPos_lp.cbegin(), AdjPValuesPos_lp.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesPos_lp](int32_t a, int32_t b){return PValuesPos_lp[a].Pvalue < PValuesPos_lp[b].Pvalue;} );
	reorder(copyPValuesPos_lp, indexes);
	assert( is_sorted(copyPValuesPos_lp.cbegin(), copyPValuesPos_lp.cend(), PRegion_t::CompPvalue) );
	reorder(copyAdjPValuesPos_lp, indexes);
	// PValuesNeg_lp
	indexes.resize(PValuesNeg_lp.size());
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesNeg_lp(PValuesNeg_lp.cbegin(), PValuesNeg_lp.cend());
	vector<PRegion_t> copyAdjPValuesNeg_lp(AdjPValuesNeg_lp.cbegin(), AdjPValuesNeg_lp.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesNeg_lp](int32_t a, int32_t b){return PValuesNeg_lp[a].Pvalue < PValuesNeg_lp[b].Pvalue;} );
	reorder(copyPValuesNeg_lp, indexes);
	assert( is_sorted(copyPValuesNeg_lp.cbegin(), copyPValuesNeg_lp.cend(), PRegion_t::CompPvalue) );
	reorder(copyAdjPValuesNeg_lp, indexes);

	// refine for each gene
	for (const string& g : involved_genes) {
		// collect intersection between isoforms of the gene and AdjustmentList
		vector<int32_t> tids;
		for(int32_t& t : GeneTransMap[g]){
			if (binary_search(AdjustmentList_copy.begin(), AdjustmentList_copy.end(), t))
				tids.push_back(t);
		}
		assert(tids.size() > 0);
		sort(tids.begin(), tids.end());
		if (tids.size() > 2){
			// calculate refine
			vector<int32_t> tmpAdjList = RefineLPTrans_singlecase(tids, dt, copyPValuesPos_salmon, copyPValuesNeg_salmon, 
				copyAdjPValuesPos_salmon, copyAdjPValuesNeg_salmon, copyPValuesPos_lp, copyPValuesNeg_lp, 
				copyAdjPValuesPos_lp, copyAdjPValuesNeg_lp);
			assert(tmpAdjList.size() > 0);
			assert(tmpAdjList.size() <= tids.size());
			RefinedAdjustmentList.insert(RefinedAdjustmentList.end(), tmpAdjList.begin(), tmpAdjList.end());
		}
		else
			RefinedAdjustmentList.insert(RefinedAdjustmentList.end(), tids.begin(), tids.end());
	}
	sort(RefinedAdjustmentList.begin(), RefinedAdjustmentList.end());
	// sanity check: AdjustmentList not changed
	assert(AdjustmentList_copy.size() == AdjustmentList.size());
	return RefinedAdjustmentList;
};

