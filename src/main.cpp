/*
Part of Salmon Anomaly Detection
(c) 2019 by  Cong Ma, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "Transcript.hpp"
#include "IO.hpp"
#include "DistTest.hpp"
#include "LPReassign.hpp"

using namespace std;

double PvalueThresh=0.1;

int32_t main (int32_t argc, char* argv[])
{
	if (argc < 7 || argc > 8) {
		std::cerr << "SAD <mode: 0 for Salmon, 1 for RSEM> <GTFfile> <QuantFile> <CorrectionFile> <StartposFile> <OutputPrefix> (number_threads)"
			  << std::endl;
		exit(1);
	}

	string mode(argv[1]);
	if (mode != "0" && mode != "1") {
		cerr << "Wrong mode." << endl;
		exit(1);
	}
	string GTFfile(argv[2]);
	string QuantFile(argv[3]);
	string CorrectionFile(argv[4]);
	string StartposFile(argv[5]);
	string OutputPrefix(argv[6]);
	int32_t num_threads = 8;
	if (argc > 7)
		num_threads = atoi(argv[7]);

	map<string,double> ExpMap;
	vector<double> Exp;
	map<string,double> TPM;
	map<string,int32_t> TransLength;
	vector<Transcript_t> Transcripts;
	map<string,int32_t> TransIndex;
	map<string,string> TransGeneMap;
	map< string,vector<string> > GeneTransMap;

	if (mode == "0")
		ReadSalmonQuant(QuantFile, ExpMap, TPM, TransLength);
	else
		ReadRSEMQuant(QuantFile, ExpMap, TPM, TransLength);
	ReadGTF(GTFfile, Transcripts);
	TrimTranscripts(ExpMap, Transcripts);
	Map_Gene_Trans(Transcripts, TransIndex, TransGeneMap, GeneTransMap);
	Exp.assign(TransIndex.size(), 0);
	for (map<string,int32_t>::const_iterator it = TransIndex.cbegin(); it != TransIndex.cend(); it++)
		Exp[it->second] = ExpMap[it->first];

	vector< vector<double> > Expected;
	vector< vector<double> > Observed;
	bool has_junction_file = false;
	vector< vector<double> > junc_obs;
	ReadCorrection(CorrectionFile, TransIndex, TransLength, Expected);
	ReadStartpos(StartposFile, TransIndex, TransLength, Observed);
	{
		// read junction observed read count
		string junction_filename = "";
		size_t s1 = StartposFile.find_last_of("/");
		size_t s2 = StartposFile.find_last_of(".");
		if (s2 != string::npos && s2 > s1)
			junction_filename = StartposFile.substr(0, s2) + "_junction" + StartposFile.substr(s2);
		else
			junction_filename = StartposFile + "_junction";
		if (access( junction_filename.c_str(), F_OK ) != -1)
			has_junction_file = true;
		if (has_junction_file)
			ReadStartpos(junction_filename, TransIndex, TransLength, junc_obs);
	}

	// sanity check
	for (map<string,double>::const_iterator ittpm = TPM.cbegin(); ittpm != TPM.cend(); ittpm++){
		map<string,int32_t>::const_iterator itidx = TransIndex.find(ittpm->first);
		assert(itidx != TransIndex.cend());
		double sumreads = 0;
		for (uint32_t i=0; i<Observed[itidx->second].size(); i++)
			sumreads += Observed[itidx->second][i];
		if (fabs(ittpm->second - 1.0*sumreads/Observed[itidx->second].size()) > 1e-2)
			cout << (ittpm->first) <<"\t"<< (ittpm->second) <<"\t"<< (itidx->second) <<"\t"<< sumreads <<"\t"<< (Observed[itidx->second].size()) <<"\n";
	}

	DistTest_t dt(TransIndex, TransLength, Expected, Observed, num_threads);
	dt.AdjustExpected();
	dt.CalDeletionScore();

	// write adjusted expected distribution and covariance matrix
	vector< Eigen::VectorXd > ExpectedBinNorm = dt.GetExpectedBinNorm();
	WriteAdjustExpectedDistribution(OutputPrefix+"_expectedbinnorm.dat", dt.TransNames, ExpectedBinNorm);
	vector<Eigen::MatrixXd> Covariance = dt.GetCovariance();
	const vector<uint32_t>& LenClass = dt.GetLenClass();
	WriteCovarianceMatrix(OutputPrefix+"_covariance.dat", dt.TransNames, LenClass, Covariance);

	// calculate P value based on salmon assignment
	vector<PRegion_t> PValuesPos_salmon;
	vector<PRegion_t> AdjPValuesPos_salmon;
	vector<PRegion_t> PValuesNeg_salmon;
	vector<PRegion_t> AdjPValuesNeg_salmon;
	vector<PRegion_t> PValuesPos_lp;
	vector<PRegion_t> AdjPValuesPos_lp;
	vector<PRegion_t> PValuesNeg_lp;
	vector<PRegion_t> AdjPValuesNeg_lp;
	dt.PValue_regional(PValuesPos_salmon, PValuesNeg_salmon);
	BHAdjusting(PValuesPos_salmon, AdjPValuesPos_salmon);
	BHAdjusting(PValuesNeg_salmon, AdjPValuesNeg_salmon);

	WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_pos.tsv", Transcripts, dt.TransCov, PValuesPos_salmon, AdjPValuesPos_salmon);
	WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_neg.tsv", Transcripts, dt.TransCov, PValuesNeg_salmon, AdjPValuesNeg_salmon);

	// collect initial adjustment list
	vector<int32_t> AdjustmentList;
	vector<string> RelatedGenes;
	for (uint32_t i = 0; i < AdjPValuesPos_salmon.size(); i++){
		if (AdjPValuesPos_salmon[i].Pvalue < PvalueThresh){
			assert(AdjPValuesPos_salmon[i].TID < (int32_t)Transcripts.size());
			map<string,string>::iterator itmap = TransGeneMap.find(Transcripts[AdjPValuesPos_salmon[i].TID].TransID);
			assert(itmap != TransGeneMap.end());
			RelatedGenes.push_back( itmap->second );
		}
	}
	for (uint32_t i = 0; i < AdjPValuesNeg_salmon.size(); i++){
		if (AdjPValuesNeg_salmon[i].Pvalue < PvalueThresh){
			assert(AdjPValuesNeg_salmon[i].TID < (int32_t)Transcripts.size());
			map<string,string>::iterator itmap = TransGeneMap.find(Transcripts[AdjPValuesNeg_salmon[i].TID].TransID);
			assert(itmap != TransGeneMap.end());
			RelatedGenes.push_back( itmap->second );
		}
	}
	sort(RelatedGenes.begin(), RelatedGenes.end());
	RelatedGenes.resize( distance(RelatedGenes.begin(), unique(RelatedGenes.begin(),RelatedGenes.end())) );
	for (const string& g : RelatedGenes) {
		const vector<string>& trans = GeneTransMap[g];
		vector<string> trans_with_gaussianerror = dt.SelectTransWithGaussianError(trans);
		for (const string& t : trans_with_gaussianerror)
			AdjustmentList.push_back( TransIndex[t] );
	}
	sort(AdjustmentList.begin(), AdjustmentList.end());
	AdjustmentList.resize( distance(AdjustmentList.begin(), unique(AdjustmentList.begin(),AdjustmentList.end())) );

	// {
	// 	// calculate the overall p value of the regional-significant transcripts
	// 	cout <<"Calculating the transcript-level p value corresponding to the regional anomalies\n";
	// 	dt.Num_Threads = 8;
	// 	vector<double> PValuesPos;
	// 	vector<double> PValuesNeg;
	// 	vector<bool> Choices;
	// 	dt.PValue_overall_empirical(AdjustmentList, PValuesPos, PValuesNeg, Choices);
	// 	WriteOverallPvalue(OutputPrefix+"_initial_pvalue_overall.tsv", AdjustmentList, Transcripts, dt.TransCov, dt, PValuesPos, PValuesNeg, Choices);
	// }

	// initial LP
	LPReassign_t LP(TransIndex, TransGeneMap, GeneTransMap, Transcripts, Expected, Observed);
	if (has_junction_file)
		LP.InitializeJunction(Transcripts, junc_obs);
	int32_t round = 1;
	while (true) {
		// LP adjustment
		LP.SetAdjustmentList(AdjustmentList);
		vector< vector<double> > newAssignment;
		map< string,vector<MovingRead_t> > MovingMat;
		vector<double> LPExp = LP.ReassignReads(newAssignment, MovingMat);
		WriteNewAssignment_NumReads(OutputPrefix+"_LP"+to_string(round)+".tsv", Transcripts, AdjustmentList, newAssignment);

		// calculate LP P value
		// reset everything to salmon
		dt.UpdateObserved(Observed);
		// set adjustment list to newAssignment
		dt.UpdateObserved(AdjustmentList, newAssignment);
		dt.CalDeletionScore();
		// sanity check
		assert(dt.TransCov.size() == dt.DeletionScore_pos.size() && dt.TransCov.size() == dt.DeletionScore_neg.size());
		for (uint32_t i = 0; i < dt.TransCov.size(); i++) {
			if (dt.TransCov[i] < 0.01 && (dt.DeletionScore_pos[i] != -1 || dt.DeletionScore_neg[i] != -1))
				cout << "Error: low coverage transcript shouldn't be evaluated.\t"<< i <<"\t"<< dt.TransCov[i] <<"\t"<< dt.DeletionScore_pos[i] <<"\t"<< dt.DeletionScore_neg[i] << endl;
			assert(dt.TransCov[i] >= 0.01 || (dt.DeletionScore_pos[i] == -1 && dt.DeletionScore_neg[i] == -1));
		}
		// calculate p value
		dt.PValue_regional(AdjustmentList, PValuesPos_lp, PValuesNeg_lp);
		for (uint32_t i = 0; i < PValuesPos_salmon.size(); i++) {
			if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesPos_salmon[i].TID))
				PValuesPos_lp.push_back(PValuesPos_salmon[i]);
		}
		for (uint32_t i = 0; i < PValuesNeg_salmon.size(); i++) {
			if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesNeg_salmon[i].TID))
				PValuesNeg_lp.push_back(PValuesNeg_salmon[i]);
		}
		BHAdjusting(PValuesPos_lp, AdjPValuesPos_lp);
		BHAdjusting(PValuesNeg_lp, AdjPValuesNeg_lp);
		WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_pos_LP"+to_string(round)+".tsv", Transcripts, dt.TransCov, PValuesPos_lp, AdjPValuesPos_lp);
		WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_neg_LP"+to_string(round)+".tsv", Transcripts, dt.TransCov, PValuesNeg_lp, AdjPValuesNeg_lp);

		// update adjustment list
		int32_t oldsize = AdjustmentList.size();
		AdjustmentList = ReduceAssignmentList(AdjustmentList, GeneTransMap, TransIndex, Exp, LPExp, MovingMat, 
						      AdjPValuesPos_salmon, AdjPValuesNeg_salmon, AdjPValuesPos_lp, AdjPValuesNeg_lp);
		if (oldsize - AdjustmentList.size() < 100)
			break;
		round++;
	}

	// refine adjustment
	dt.Num_Threads = 1;
	AdjustmentList = LP.RefineLPTrans(dt, PValuesPos_salmon, PValuesNeg_salmon, AdjPValuesPos_salmon, 
					  AdjPValuesNeg_salmon, PValuesPos_lp, PValuesNeg_lp, AdjPValuesPos_lp, AdjPValuesNeg_lp);
	LP.SetAdjustmentList(AdjustmentList);
	vector< vector<double> > newAssignment;
	map< string,vector<MovingRead_t> > MovingMat;
	vector<double> LPExp = LP.ReassignReads(newAssignment, MovingMat);
	WriteNewAssignment_NumReads(OutputPrefix+"_adjusted_quantification.tsv", Transcripts, AdjustmentList, newAssignment);
	WriteNewAssignment_Distribution(OutputPrefix+"_adjusted_observed_distribution.dat", Transcripts, AdjustmentList, newAssignment);

	// update observed distribution and deletion score
	dt.UpdateObserved(Observed);
	dt.UpdateObserved(AdjustmentList, newAssignment);
	dt.CalDeletionScore();
	dt.PValue_regional(AdjustmentList, PValuesPos_lp, PValuesNeg_lp);
	for (uint32_t i = 0; i < PValuesPos_salmon.size(); i++) {
		if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesPos_salmon[i].TID))
			PValuesPos_lp.push_back(PValuesPos_salmon[i]);
	}
	for (uint32_t i = 0; i < PValuesNeg_salmon.size(); i++) {
		if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesNeg_salmon[i].TID))
			PValuesNeg_lp.push_back(PValuesNeg_salmon[i]);
	}
	BHAdjusting(PValuesPos_lp, AdjPValuesPos_lp);
	BHAdjusting(PValuesNeg_lp, AdjPValuesNeg_lp);
	WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_pos_LP_refined.tsv", Transcripts, dt.TransCov, PValuesPos_lp, AdjPValuesPos_lp);
	WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_neg_LP_refined.tsv", Transcripts, dt.TransCov, PValuesNeg_lp, AdjPValuesNeg_lp);

	// sanity check about coverage and deletion score
	assert(dt.TransCov.size() == dt.DeletionScore_pos.size() && dt.TransCov.size() == dt.DeletionScore_neg.size());
	for (uint32_t i = 0; i < dt.TransCov.size(); i++) {
		if (dt.TransCov[i] < 0.01 && (dt.DeletionScore_pos[i] != -1 || dt.DeletionScore_neg[i] != -1))
			cout << "Error: low coverage transcript shouldn't be evaluated.\t"<< i <<"\t"<< dt.TransCov[i] <<"\t"<< dt.DeletionScore_pos[i] <<"\t"<< dt.DeletionScore_neg[i] << endl;
		assert(dt.TransCov[i] >= 0.01 || (dt.DeletionScore_pos[i] == -1 && dt.DeletionScore_neg[i] == -1));
	}

	// calculate overlap p value
	dt.Num_Threads = num_threads;
	vector<double> PValuesPos;
	vector<double> PValuesNeg;
	vector<bool> Choices;
	dt.PValue_overall_empirical(PValuesPos, PValuesNeg, Choices);
	// OldWriteOverallPvalue(OutputPrefix+"_old_pvalue_overall", Transcripts, dt.TransCov, dt, PValues, Choices);
	WriteOverallPvalue(OutputPrefix+"_unadjustable_pvalue.tsv", Transcripts, dt.TransCov, dt, PValuesPos, PValuesNeg, Choices);
}

