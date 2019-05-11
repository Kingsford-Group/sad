/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "LPReassign.hpp"
#include "Transcript.hpp"
#include "DistTest.hpp"

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


LPReassign_t::LPReassign_t(const map<string,int32_t>& _TransIndex, const map<string,string>& _TransGeneMap, 
	const map<string, vector<string> >& _GeneTransMap, const vector<Transcript_t>& Transcripts, 
	const vector< vector<double> >& Expected, const vector< vector<double> >& Observed, int32_t _Bin_Size, int32_t _Num_Threads)
{
	Bin_Size = _Bin_Size;
	Num_Threads = _Num_Threads;
	LPReassign_t(_TransIndex, _TransGeneMap, _GeneTransMap, Transcripts, Expected, Observed);
};

//****************** test *********************
vector<double> LPReassign_t::Quantify_singlecase(Eigen::MatrixXd& exp, Eigen::VectorXd& obs)
{
	assert(exp.rows() == obs.size());
	// create GUROBI environment
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
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
			assert( obsfull[i] >= 0 );
			for (int32_t j = 0; j < assign_theo.size(); j++)
				assert( assign_theo[j] >= 0);
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

	cout << (MovingMat.size()) <<"\t"<< (involved_genes.size()) << endl;
	
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
			dt.PValue_regional(removedAdjList, removedPValuePos, removedPValueNeg);
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
	cout << is_sorted(copyPValuesPos_salmon.cbegin(), copyPValuesPos_salmon.cend(), PRegion_t::CompPvalue) << endl;
	reorder(copyAdjPValuesPos_salmon, indexes);
	// PValuesNeg_salmon
	indexes.resize(PValuesNeg_salmon.size());
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesNeg_salmon(PValuesNeg_salmon.cbegin(), PValuesNeg_salmon.cend());
	vector<PRegion_t> copyAdjPValuesNeg_salmon(AdjPValuesNeg_salmon.cbegin(), AdjPValuesNeg_salmon.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesNeg_salmon](int32_t a, int32_t b){return PValuesNeg_salmon[a].Pvalue < PValuesNeg_salmon[b].Pvalue;} );
	reorder(copyPValuesNeg_salmon, indexes);
	cout << is_sorted(copyPValuesNeg_salmon.cbegin(), copyPValuesNeg_salmon.cend(), PRegion_t::CompPvalue) << endl;
	reorder(copyAdjPValuesNeg_salmon, indexes);
	// PValuesPos_lp
	indexes.resize(PValuesPos_lp.size());
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesPos_lp(PValuesPos_lp.cbegin(), PValuesPos_lp.cend());
	vector<PRegion_t> copyAdjPValuesPos_lp(AdjPValuesPos_lp.cbegin(), AdjPValuesPos_lp.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesPos_lp](int32_t a, int32_t b){return PValuesPos_lp[a].Pvalue < PValuesPos_lp[b].Pvalue;} );
	reorder(copyPValuesPos_lp, indexes);
	cout << is_sorted(copyPValuesPos_lp.cbegin(), copyPValuesPos_lp.cend(), PRegion_t::CompPvalue) << endl;
	reorder(copyAdjPValuesPos_lp, indexes);
	// PValuesNeg_lp
	indexes.resize(PValuesNeg_lp.size());
	iota(indexes.begin(), indexes.end(), 0);
	vector<PRegion_t> copyPValuesNeg_lp(PValuesNeg_lp.cbegin(), PValuesNeg_lp.cend());
	vector<PRegion_t> copyAdjPValuesNeg_lp(AdjPValuesNeg_lp.cbegin(), AdjPValuesNeg_lp.cend());
	sort(indexes.begin(), indexes.end(), [&PValuesNeg_lp](int32_t a, int32_t b){return PValuesNeg_lp[a].Pvalue < PValuesNeg_lp[b].Pvalue;} );
	reorder(copyPValuesNeg_lp, indexes);
	cout << is_sorted(copyPValuesNeg_lp.cbegin(), copyPValuesNeg_lp.cend(), PRegion_t::CompPvalue) << endl;
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

