#ifndef __DISTTEST__
#define __DISTTEST__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
#include <mutex>
#include <random>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions/normal.hpp"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;

class MovingRead_t;


class PRegion_t
{
public:
	int32_t TID;
	int32_t BinStart;
	int32_t BinEnd;
	double AnomalyScore;
	double Pvalue;
public:
	PRegion_t(){};
	PRegion_t(int32_t TID, int32_t BinStart, int32_t BinEnd, double AnomalyScore, double Pvalue): 
		TID(TID), BinStart(BinStart), BinEnd(BinEnd), AnomalyScore(AnomalyScore), Pvalue(Pvalue) {};

	static bool CompPvalue(const PRegion_t& a, const PRegion_t& b) {
		return a.Pvalue < b.Pvalue;
	};
};


class DistTest_t
{
public:
	int32_t Num_Threads;

	vector<int32_t> lenBounds;
	vector<int32_t> nBins;
	map<string,int32_t> TransIndex;
	vector<int32_t> TransLength;
	vector<string> TransNames;
	vector<double> TransCov;
	vector<double> percentile;
	vector<double> Sampling_NSTDS;
	vector<double> Sampling_WEIGHTS;
	
	vector<double> DeletionScore_pos; // exp - obs
	vector<double> DeletionScore_neg; // obs - exp
	vector< pair<int32_t,int32_t> > DeletionRegion_pos;
	vector< pair<int32_t,int32_t> > DeletionRegion_neg;

private:
	vector< Eigen::VectorXd > ExpectedBinNorm;
	vector< Eigen::VectorXd > ObservedBinNorm;
	vector<int32_t> LenClass;
	vector<Eigen::VectorXd> Mean;
	vector<Eigen::MatrixXd> Covariance;
	vector<Eigen::MatrixXd> CholeskyDecom;
	
public:
	DistTest_t();
	DistTest_t(const map<string,int32_t>& _TransIndex, const map<string,int32_t>& _TransLength, const vector< vector<double> >& Expected, const vector< vector<double> >& Observed);
	DistTest_t(const map<string,int32_t>& _TransIndex, const map<string,int32_t>& _TransLength, const vector< vector<double> >& Expected, const vector< vector<double> >& Observed, int32_t _Num_Threads);

	vector< Eigen::VectorXd > GetExpectedBinNorm() {
		return ExpectedBinNorm;
	};
	vector< Eigen::VectorXd > GetObservedBinNorm() {
		return ObservedBinNorm;
	};
	vector<Eigen::MatrixXd> GetCovariance() {
		return Covariance;
	};
	vector<int32_t> GetLenClass() {
		return LenClass;
	};
	
	void AdjustExpected(double nstdthresh = 1);
	void UpdateObserved(const vector< vector<double> >& NewObserved);
	void UpdateObserved(const vector<int32_t>& AdjustmentList, const vector< vector<double> >& NewObserved);

	void CalDeletionScore(double covthresh = 0.01, double numthresh = 20);
	void CalDeletionScore(const vector<int32_t>& AdjustmentList, double covthresh = 0.01, double numthresh = 20);
	static double CalDeletionScore_single(Eigen::VectorXd Dist1, Eigen::VectorXd Dist2, pair<int32_t,int32_t>& region, 
		bool check1 = false, double numthresh = 20);
	void PValue_regional(vector<PRegion_t>& PValuesPos, vector<PRegion_t>& PValuesNeg);
	void PValue_regional(const vector<int32_t>& AdjustmentList, vector<PRegion_t>& PValuesPos, vector<PRegion_t>& PValuesNeg);
	void PValue_overall_empirical(vector<double>& PValuesPos, vector<double>& PValuesNeg, vector<bool>& Choices);
	void PValue_overall_empirical(const vector<int32_t>& AdjustmentList, vector<double>& PValuesPos, vector<double>& PValuesNeg, vector<bool>& Choices);

// private:
	vector<PRegion_t> SinglePvalue_regional_pos(int32_t ind);
	vector<PRegion_t> SinglePvalue_regional_neg(int32_t ind);
	pair<double,double> SinglePvalue_overall(int32_t ind);
	pair<double,double> SinglePvalue_overall_empirical(int32_t ind, int32_t num_sampling, const gsl_rng * r);
};


void BHAdjusting(const vector<PRegion_t>& RawPRegions, vector<PRegion_t>& AdjustedPRegions);

void BHAdjusting(const vector<double>& RawPvalues, vector<double>& AdjustedPvalues);

vector<int32_t> ReduceAssignmentList(const vector<int32_t>& AssignmentList, const map< string,vector<string> >& GeneTransMap, 
	const map<string,int32_t>& TransIndex, const vector<double>& SalmonExp, const vector<double>& LPExp, 
	const map< string,vector<MovingRead_t> >& MovingMat, const vector<PRegion_t>& AdjPValuesPos_salmon, const vector<PRegion_t>& AdjPValuesNeg_salmon, 
	const vector<PRegion_t>& AdjPValuesPos_lp, const vector<PRegion_t>& AdjPValuesNeg_lp, double PvalueThresh=0.1, double MoveReadThresh=0.1);

#endif