#ifndef __IO_H__
#define __IO_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <omp.h>
#include <mutex>
#include "boost/algorithm/string.hpp"
#include "Eigen/Dense"

using namespace std;

class Transcript_t;
class PRegion_t;
class DistTest_t;

void ReadSalmonQuant(string quantfile, map<string,double>& SalmonExp, map<string,double>& TPM, map<string,int32_t>& TransLength);

void ReadCorrection(string correctionfile, const map<string,int32_t>& TransIndex, const map<string,int32_t>& TransLength, vector< vector<double> >& Expected);

void ReadStartpos(string startposfile, const map<string,int32_t>& TransIndex, const map<string,int32_t>& TransLength, vector< vector<double> >& Observed);

void WriteNewAssignment_NumReads(string outputfile, const vector<Transcript_t>& Transcripts, 
	const vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment);

void WriteNewAssignment_Distribution(string outputfile, const vector<Transcript_t>& Transcripts, 
	const vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment);

void WriteAdjustExpectedDistribution(string outputfile, const vector<string>& TransNames, vector< Eigen::VectorXd >& ExpectedBinNorm);

void WriteCovarianceMatrix(string outputfile, const vector<string>& TransNames, vector<int32_t>& LenClass, vector<Eigen::MatrixXd>& Covariance);

void ReadNewAssignment_Distribution(string inputfile, const map<string,int32_t>& TransIndex, 
	vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment);

void ReadNewAssignment_Distribution_old(string numreadsfile, string distfile, const map<string,int32_t>& TransIndex,
	vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment);

void WriteRegionalPvalue(string outputfile, const vector<Transcript_t>& Transcripts,
	const vector<double>& TransCov, const vector<PRegion_t>& RawPRegion, const vector<PRegion_t>& AdjustedPRegion);

void WriteOverallPvalue(string outputfile, const vector<Transcript_t>& Transcripts, const vector<double>& TransCov, 
	const DistTest_t& dt, const vector<double>& PValuesPos, const vector<double>& PValuesNeg, const vector<bool>& Choices);

void WriteOverallPvalue(string outputfile, const vector<int32_t>& AdjustmentList, const vector<Transcript_t>& Transcripts, const vector<double>& TransCov, 
	const DistTest_t& dt, const vector<double>& PValuesPos, const vector<double>& PValuesNeg, const vector<bool>& Choices);

void OldWriteOverallPvalue(string outputfile, const vector<Transcript_t>& Transcripts, const vector<double>& TransCov, 
	const DistTest_t& dt, const vector<double>& PValues, const vector<bool>& Choices);

#endif