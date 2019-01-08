#include "Transcript.hpp"
#include "IO.hpp"
#include "DistTest.hpp"
#include "LPReassign.hpp"

using namespace std;

double PvalueThresh=0.1;

int32_t main (int32_t argc, char* argv[])
{
	if (argc==1){
		printf("SAD <GTFfile> <SalmonQuantFile> <CorrectionFile> <StartposFile> <OutputPrefix>\n");
	}
	else{
		string GTFfile(argv[1]);
		string SalmonQuantFile(argv[2]);
		string CorrectionFile(argv[3]);
		string StartposFile(argv[4]);
		string OutputPrefix(argv[5]);

		map<string,double> SalmonExpMap;
		vector<double> SalmonExp;
		map<string,double> TPM;
		map<string,int32_t> TransLength;
		vector<Transcript_t> Transcripts;
		map<string,int32_t> TransIndex;
		map<string,string> TransGeneMap;
		map< string,vector<string> > GeneTransMap;

		ReadSalmonQuant(SalmonQuantFile, SalmonExpMap, TPM, TransLength);
		ReadGTF(GTFfile, Transcripts);
		TrimTranscripts(SalmonExpMap, Transcripts);
		Map_Gene_Trans(Transcripts, TransIndex, TransGeneMap, GeneTransMap);
		SalmonExp.assign(TransIndex.size(), 0);
		for (map<string,int32_t>::const_iterator it = TransIndex.cbegin(); it != TransIndex.cend(); it++)
			SalmonExp[it->second] = SalmonExpMap[it->first];

		vector< vector<double> > Expected;
		vector< vector<double> > Observed;
		ReadCorrection(CorrectionFile, TransIndex, TransLength, Expected);
		ReadStartpos(StartposFile, TransIndex, TransLength, Observed);

		// sanity check
		for (map<string,double>::const_iterator ittpm = TPM.cbegin(); ittpm != TPM.cend(); ittpm++){
			map<string,int32_t>::const_iterator itidx = TransIndex.find(ittpm->first);
			assert(itidx != TransIndex.cend());
			double sumreads = 0;
			for (int32_t i=0; i<Observed[itidx->second].size(); i++)
				sumreads += Observed[itidx->second][i];
			if (fabs(ittpm->second - 1.0*sumreads/Observed[itidx->second].size()) > 1e-2)
				cout << (ittpm->first) <<"\t"<< (ittpm->second) <<"\t"<< (itidx->second) <<"\t"<< sumreads <<"\t"<< (Observed[itidx->second].size()) <<"\n";
		}

		DistTest_t dt(TransIndex, TransLength, Expected, Observed);
		dt.AdjustExpected();
		dt.CalDeletionScore();

		// write adjusted expected distribution and covariance matrix
		vector< Eigen::VectorXd > ExpectedBinNorm = dt.GetExpectedBinNorm();
		WriteAdjustExpectedDistribution(OutputPrefix+"_expectedbinnorm.dat", dt.TransNames, ExpectedBinNorm);
		vector<Eigen::MatrixXd> Covariance = dt.GetCovariance();
		vector<int32_t> LenClass = dt.GetLenClass();
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

		WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_pos", Transcripts, dt.TransCov, PValuesPos_salmon, AdjPValuesPos_salmon);
		WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_neg", Transcripts, dt.TransCov, PValuesNeg_salmon, AdjPValuesNeg_salmon);

		// collect initial adjustment list
		vector<int32_t> AdjustmentList;
		vector<string> RelatedGenes;
		for (int32_t i = 0; i < AdjPValuesPos_salmon.size(); i++){
			if (AdjPValuesPos_salmon[i].Pvalue < PvalueThresh){
				assert(AdjPValuesPos_salmon[i].TID < Transcripts.size());
				map<string,string>::iterator itmap = TransGeneMap.find(Transcripts[AdjPValuesPos_salmon[i].TID].TransID);
				assert(itmap != TransGeneMap.end());
				RelatedGenes.push_back( itmap->second );
			}
		}
		for (int32_t i = 0; i < AdjPValuesNeg_salmon.size(); i++){
			if (AdjPValuesNeg_salmon[i].Pvalue < PvalueThresh){
				assert(AdjPValuesNeg_salmon[i].TID < Transcripts.size());
				map<string,string>::iterator itmap = TransGeneMap.find(Transcripts[AdjPValuesNeg_salmon[i].TID].TransID);
				assert(itmap != TransGeneMap.end());
				RelatedGenes.push_back( itmap->second );
			}
		}
		sort(RelatedGenes.begin(), RelatedGenes.end());
		RelatedGenes.resize( distance(RelatedGenes.begin(), unique(RelatedGenes.begin(),RelatedGenes.end())) );
		for (const string& g : RelatedGenes) {
			const vector<string>& trans = GeneTransMap[g];
			for (const string& t : trans)
				AdjustmentList.push_back( TransIndex[t] );
		}
		sort(AdjustmentList.begin(), AdjustmentList.end());
		AdjustmentList.resize( distance(AdjustmentList.begin(), unique(AdjustmentList.begin(),AdjustmentList.end())) );

		{
			// calculate the overall p value of the regional-significant transcripts
			cout <<"Calculating the transcript-level p value corresponding to the regional anomalies\n";
			dt.Num_Threads = 8;
			vector<double> PValuesPos;
			vector<double> PValuesNeg;
			vector<bool> Choices;
			dt.PValue_overall_empirical(AdjustmentList, PValuesPos, PValuesNeg, Choices);
			WriteOverallPvalue(OutputPrefix+"_initial_pvalue_overall", AdjustmentList, Transcripts, dt.TransCov, dt, PValuesPos, PValuesNeg, Choices);
		}

		// initial LP
		LPReassign_t LP(TransIndex, TransGeneMap, GeneTransMap, Transcripts, Expected, Observed);
		int32_t round = 1;
		while (true) {
			// LP adjustment
			LP.SetAdjustmentList(AdjustmentList);
			vector< vector<double> > newAssignment;
			map< string,vector<MovingRead_t> > MovingMat;
			vector<double> LPExp = LP.ReassignReads(newAssignment, MovingMat);
			WriteNewAssignment_NumReads(OutputPrefix+"_LP"+to_string(round), Transcripts, AdjustmentList, newAssignment);

			// calculate LP P value
			// reset everything to salmon
			dt.UpdateObserved(Observed);
			// set adjustment list to newAssignment
			dt.UpdateObserved(AdjustmentList, newAssignment);
			dt.CalDeletionScore();
			// sanity check
			assert(dt.TransCov.size() == dt.DeletionScore_pos.size() && dt.TransCov.size() == dt.DeletionScore_neg.size());
			for (int32_t i = 0; i < dt.TransCov.size(); i++) {
				if (dt.TransCov[i] < 0.01 && (dt.DeletionScore_pos[i] != -1 || dt.DeletionScore_neg[i] != -1))
					cout << "Error: low coverage transcript shouldn't be evaluated.\t"<< i <<"\t"<< dt.TransCov[i] <<"\t"<< dt.DeletionScore_pos[i] <<"\t"<< dt.DeletionScore_neg[i] << endl;
				assert(dt.TransCov[i] >= 0.01 || (dt.DeletionScore_pos[i] == -1 && dt.DeletionScore_neg[i] == -1));
			}
			// calculate p value
			dt.PValue_regional(AdjustmentList, PValuesPos_lp, PValuesNeg_lp);
			for (int32_t i = 0; i < PValuesPos_salmon.size(); i++) {
				if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesPos_salmon[i].TID))
					PValuesPos_lp.push_back(PValuesPos_salmon[i]);
			}
			for (int32_t i = 0; i < PValuesNeg_salmon.size(); i++) {
				if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesNeg_salmon[i].TID))
					PValuesNeg_lp.push_back(PValuesNeg_salmon[i]);
			}
			BHAdjusting(PValuesPos_lp, AdjPValuesPos_lp);
			BHAdjusting(PValuesNeg_lp, AdjPValuesNeg_lp);
			WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_pos_LP"+to_string(round), Transcripts, dt.TransCov, PValuesPos_lp, AdjPValuesPos_lp);
			WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_neg_LP"+to_string(round), Transcripts, dt.TransCov, PValuesNeg_lp, AdjPValuesNeg_lp);

			// update adjustment list
			int32_t oldsize = AdjustmentList.size();
			AdjustmentList = ReduceAssignmentList(AdjustmentList, GeneTransMap, TransIndex, SalmonExp, LPExp, MovingMat, 
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
		WriteNewAssignment_NumReads(OutputPrefix+"_LP_refined", Transcripts, AdjustmentList, newAssignment);
		WriteNewAssignment_Distribution(OutputPrefix+"_LP_refined_dist", Transcripts, AdjustmentList, newAssignment);

		// update observed distribution and deletion score
		dt.UpdateObserved(Observed);
		dt.UpdateObserved(AdjustmentList, newAssignment);
		dt.CalDeletionScore();
		dt.PValue_regional(AdjustmentList, PValuesPos_lp, PValuesNeg_lp);
		for (int32_t i = 0; i < PValuesPos_salmon.size(); i++) {
			if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesPos_salmon[i].TID))
				PValuesPos_lp.push_back(PValuesPos_salmon[i]);
		}
		for (int32_t i = 0; i < PValuesNeg_salmon.size(); i++) {
			if (!binary_search(AdjustmentList.begin(), AdjustmentList.end(), PValuesNeg_salmon[i].TID))
				PValuesNeg_lp.push_back(PValuesNeg_salmon[i]);
		}
		BHAdjusting(PValuesPos_lp, AdjPValuesPos_lp);
		BHAdjusting(PValuesNeg_lp, AdjPValuesNeg_lp);
		WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_pos_LP_refined", Transcripts, dt.TransCov, PValuesPos_lp, AdjPValuesPos_lp);
		WriteRegionalPvalue(OutputPrefix+"_pvalue_regional_neg_LP_refined", Transcripts, dt.TransCov, PValuesNeg_lp, AdjPValuesNeg_lp);

		// sanity check about coverage and deletion score
		assert(dt.TransCov.size() == dt.DeletionScore_pos.size() && dt.TransCov.size() == dt.DeletionScore_neg.size());
		for (int32_t i = 0; i < dt.TransCov.size(); i++) {
			if (dt.TransCov[i] < 0.01 && (dt.DeletionScore_pos[i] != -1 || dt.DeletionScore_neg[i] != -1))
				cout << "Error: low coverage transcript shouldn't be evaluated.\t"<< i <<"\t"<< dt.TransCov[i] <<"\t"<< dt.DeletionScore_pos[i] <<"\t"<< dt.DeletionScore_neg[i] << endl;
			assert(dt.TransCov[i] >= 0.01 || (dt.DeletionScore_pos[i] == -1 && dt.DeletionScore_neg[i] == -1));
		}

		// calculate overlap p value
		dt.Num_Threads = 8;
		vector<double> PValuesPos;
		vector<double> PValuesNeg;
		vector<bool> Choices;
		dt.PValue_overall_empirical(PValuesPos, PValuesNeg, Choices);
		// OldWriteOverallPvalue(OutputPrefix+"_old_pvalue_overall", Transcripts, dt.TransCov, dt, PValues, Choices);
		WriteOverallPvalue(OutputPrefix+"_pvalue_overall", Transcripts, dt.TransCov, dt, PValuesPos, PValuesNeg, Choices);
	}
};