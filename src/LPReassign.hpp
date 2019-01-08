#ifndef __LPREASSIGN__
#define __LPREASSIGN__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
#include <mutex>
#include "Eigen/Dense"
#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions/normal.hpp"
#include "gurobi_c++.h"
// #include "glpk.h"

using namespace std;

class Transcript_t;

class DistTest_t;

class PRegion_t;

class MovingRead_t
{
public:
	int32_t TID1;
	int32_t TID2;
	double NReads;
public:
	MovingRead_t(){};
	MovingRead_t(int32_t TID1, int32_t TID2, double NReads): TID1(TID1), TID2(TID2), NReads(NReads) {};

	bool operator < (const MovingRead_t& rhs) const {
		if (TID1 != rhs.TID1)
			return TID1 < rhs.TID1;
		else
			return TID2 < rhs.TID2;
	};

	bool operator == (const MovingRead_t& rhs) const {
		return TID1 == rhs.TID1 && TID2 == rhs.TID2;
	};
};


class LPReassign_t
{
public:
	int32_t Bin_Size = 50;
	int32_t Num_Threads = 4;

	vector<int32_t> AdjustmentList;
	map<string,int32_t> TransIndex;
	vector<string> TransGeneMap;
	map<string, vector<int32_t> > GeneTransMap;
private:
	vector< Eigen::VectorXd > ExpectedFullNorm; // full length expected distribution in genomic coordinate
	vector< Eigen::VectorXd > ObservedFull; // full length observed distribution in genomic coordinate, not normalized
	vector< vector<bool> > PositionExistence; // indicate for each transcript, whether it contains the position in full length vector
	vector< Eigen::VectorXd > ExpectedBinNorm;
	vector< Eigen::VectorXd > ObservedBin;

public:
	LPReassign_t(){};
	LPReassign_t(const map<string,int32_t>& _TransIndex, const map<string,string>& _TransGeneMap, 
		const map<string, vector<string> >& _GeneTransMap, const vector<Transcript_t>& Transcripts, 
		const vector< vector<double> >& Expected, const vector< vector<double> >& Observed);
	LPReassign_t(const map<string,int32_t>& _TransIndex, const map<string,string>& _TransGeneMap, 
		const map<string, vector<string> >& _GeneTransMap, const vector<Transcript_t>& Transcripts, 
		const vector< vector<double> >& Expected, const vector< vector<double> >& Observed, int32_t _Bin_Size, int32_t _Num_Threads);

public:
	void SetAdjustmentList(const vector<int32_t>& _AdjustmentList){
		AdjustmentList.clear();
		AdjustmentList.assign(_AdjustmentList.cbegin(), _AdjustmentList.cend());
	};

	// XXX: remove AdjustmentList from member attribute, but add as an argument
	vector<double> Quantify_singlecase(Eigen::MatrixXd& exp, Eigen::VectorXd& obs);
	vector<double> ReassignReads(vector< vector<double> >& newAssignment, map< string,vector<MovingRead_t> >& MovingMat);
	vector<double> ReassignReads(vector< vector<double> >& newAssignment, map< string,vector<MovingRead_t> >& MovingMat, const vector<int32_t>& _AdjustmentList);

	vector<int32_t> RefineLPTrans_singlecase(const vector<int32_t>& shortAdjustmentList, DistTest_t& dt, 
		const vector<PRegion_t>& PValuesPos_salmon, const vector<PRegion_t>& PValuesNeg_salmon, 
		const vector<PRegion_t>& AdjPValuesPos_salmon, const vector<PRegion_t>& AdjPValuesNeg_salmon, 
		const vector<PRegion_t>& PValuesPos_lp, const vector<PRegion_t>& PValuesNeg_lp, 
		const vector<PRegion_t>& AdjPValuesPos_lp, const vector<PRegion_t>& AdjPValuesNeg_lp, double PvalueThresh=0.1);
	vector<int32_t> RefineLPTrans(DistTest_t& dt, 
		const vector<PRegion_t>& PValuesPos_salmon, const vector<PRegion_t>& PValuesNeg_salmon, 
		const vector<PRegion_t>& AdjPValuesPos_salmon, const vector<PRegion_t>& AdjPValuesNeg_salmon, 
		const vector<PRegion_t>& PValuesPos_lp, const vector<PRegion_t>& PValuesNeg_lp, 
		const vector<PRegion_t>& AdjPValuesPos_lp, const vector<PRegion_t>& AdjPValuesNeg_lp, double PvalueThresh=0.1);
};

void reorder(vector<int32_t>& arr, vector<int32_t>& index);

#endif