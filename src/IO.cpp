#include "IO.hpp"
#include "Transcript.hpp"
#include "DistTest.hpp"

using namespace std;

void ReadSalmonQuant(string quantfile, map<string,double>& SalmonExp, map<string,double>& TPM, 
	map<string,int32_t>& TransLength)
{
	// clear variables
	SalmonExp.clear();
	TPM.clear();
	TransLength.clear();
	// reading file
	ifstream input(quantfile);
	string line;
	int32_t linecount = 0;
	while (getline(input, line)){
		linecount++;
		if (linecount == 1)
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		SalmonExp[strs[0]] = stod(strs[4]);
		TPM[strs[0]] = stod(strs[4]) / stoi(strs[1]);
		TransLength[strs[0]] = stoi(strs[1]);
	}
	input.close();
};

void ReadCorrection(string correctionfile, const map<string,int32_t>& TransIndex, const map<string,int32_t>& TransLength, 
	vector< vector<double> >& Expected)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Start reading "<<correctionfile<<endl;
	// initialize variables
	Expected.clear();
	for(int32_t i = 0; i < TransIndex.size(); i++){
		vector<double> tmp;
		Expected.push_back(tmp);
	}
	for(map<string,int32_t>::const_iterator it = TransIndex.cbegin(); it != TransIndex.cend(); it++){
		map<string,int32_t>::const_iterator itlen = TransLength.find(it->first);
		assert(itlen != TransLength.cend());
		vector<double> tmp(itlen->second, 0);
		Expected[it->second] = tmp;
	}
	// reading binary file
	ifstream input(correctionfile, ios::binary);
	int32_t numtrans;
	input.read((char*)(&numtrans), sizeof(int32_t));
	for (int32_t i = 0; i < numtrans; i++){
		// reading length info
		int32_t namelen;
		int32_t seqlen;
		input.read((char*)&namelen, sizeof(int32_t));
		input.read((char*)&seqlen, sizeof(int32_t));
		// reading actual values
		char* buf=new char[namelen];
		string name;
		vector<double> mydata(seqlen, 0);
		input.read(buf, namelen*sizeof(char));
		input.read((char*)(mydata.data()), seqlen*sizeof(double));
		name.assign(buf, namelen);
		// find corresponding position
		map<string,int32_t>::const_iterator itidx = TransIndex.find(name);
		if (itidx != TransIndex.cend()){
			assert(Expected[itidx->second].size() == seqlen);
			Expected[itidx->second] = mydata;
		}
	}
	input.close();
	// time info
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish reading"<<endl;
};

void ReadStartpos(string startposfile, const map<string,int32_t>& TransIndex, const map<string,int32_t>& TransLength, 
	vector< vector<double> >& Observed)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Start reading "<<startposfile<<endl;
	// initialize variables
	Observed.clear();
	for(int32_t i = 0; i < TransIndex.size(); i++){
		vector<double> tmp;
		Observed.push_back(tmp);
	}
	for(map<string,int32_t>::const_iterator it = TransIndex.cbegin(); it != TransIndex.cend(); it++){
		map<string,int32_t>::const_iterator itlen = TransLength.find(it->first);
		assert(itlen != TransLength.cend());
		vector<double> tmp(itlen->second, 0);
		Observed[it->second] = tmp;
	}
	// reading binary file
	ifstream input(startposfile, ios::binary);
	int32_t numtrans;
	input.read((char*)(&numtrans), sizeof(int32_t));
	for (int32_t i = 0; i < numtrans; i++){
		// reading length info
		int32_t namelen;
		int32_t seqlen;
		input.read((char*)(&namelen), sizeof(int32_t));
		input.read((char*)(&seqlen), sizeof(int32_t));
		// reading actual values
		char* buf = new char[namelen];
		string name;
		vector<int32_t> poses(seqlen, 0);
		vector<double> counts(seqlen, 0);
		input.read(buf, namelen*sizeof(char));
		name.assign(buf, namelen);
		input.read((char*)(poses.data()), seqlen*sizeof(int32_t));
		input.read((char*)(counts.data()), seqlen*sizeof(double));
		// generate distribution from positions and counts
		vector<double> mydata(poses.back() + 1, 0);
		for(int32_t j = 0; j < counts.size(); j++)
			mydata[poses[j]] += counts[j];
		// find corresponding position
		map<string,int32_t>::const_iterator itidx = TransIndex.find(name);
		if (itidx != TransIndex.cend()){
			assert(Observed[itidx->second].size() == mydata.size());
			Observed[itidx->second] = mydata;
		}
	}
	input.close();
	// time info
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish reading"<<endl;
};


void WriteNewAssignment_NumReads(string outputfile, const vector<Transcript_t>& Transcripts, 
	const vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment)
{
	assert(AdjustmentList.size() == newAssignment.size());
	ofstream output(outputfile, ios::out);
	output << "# Name\tLength\tNumReads\n";
	for (int32_t i = 0; i < AdjustmentList.size(); i++) {
		Eigen::VectorXd newobs = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(newAssignment[i].data(), newAssignment[i].size());
		output << (Transcripts[AdjustmentList[i]].TransID) <<"\t"<< (newobs.size()) <<"\t"<< (newobs.sum()) << endl;
	}
	output.close();
};


void WriteNewAssignment_Distribution(string outputfile, const vector<Transcript_t>& Transcripts, 
	const vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment)
{
	assert(AdjustmentList.size() == newAssignment.size());
	ofstream output(outputfile, ios::out | ios::binary);
	int32_t n = newAssignment.size();
	output.write(reinterpret_cast<char*>(&n), sizeof(int32_t));
	for (int32_t i = 0; i < newAssignment.size(); i++) {
		string tname = Transcripts[AdjustmentList[i]].TransID;
		int32_t namelen = tname.size();
		int32_t seqlen = newAssignment[i].size();
		output.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
		output.write(reinterpret_cast<char*>(&seqlen), sizeof(int32_t));
		output.write(tname.c_str(), namelen*sizeof(char));
		output.write(reinterpret_cast<char*>(newAssignment[i].data()), seqlen*sizeof(double));
	}
	output.close();
};


void WriteAdjustExpectedDistribution(string outputfile, const vector<string>& TransNames, vector< Eigen::VectorXd >& ExpectedBinNorm)
{
	assert(TransNames.size() == ExpectedBinNorm.size());
	ofstream output(outputfile, ios::out | ios::binary);
	int32_t n = ExpectedBinNorm.size();
	output.write(reinterpret_cast<char*>(&n), sizeof(int32_t));
	for (int32_t i = 0; i < ExpectedBinNorm.size(); i++) {
		string tname = TransNames[i];
		int32_t namelen = tname.size();
		int32_t seqlen = ExpectedBinNorm[i].size();
		output.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
		output.write(reinterpret_cast<char*>(&seqlen), sizeof(int32_t));
		output.write(tname.c_str(), namelen*sizeof(char));
		for (int32_t j = 0; j < seqlen; j++) {
			double value = ExpectedBinNorm[i](j);
			output.write(reinterpret_cast<char*>(&value), sizeof(double));
		}
	}
	output.close();
};


void WriteCovarianceMatrix(string outputfile, const vector<string>& TransNames, vector<int32_t>& LenClass, vector<Eigen::MatrixXd>& Covariance)
{
	assert(TransNames.size() == LenClass.size());
	ofstream output(outputfile, ios::out | ios::binary);
	// write the matrix of each class
	int32_t nclass = Covariance.size();
	output.write(reinterpret_cast<char*>(&nclass), sizeof(int32_t));
	for (int32_t i = 0; i < Covariance.size(); i++) {
		int32_t matrixsize = Covariance[i].cols();
		output.write(reinterpret_cast<char*>(&matrixsize), sizeof(int32_t));
		for (int32_t j = 0; j < matrixsize; j++)
			for (int32_t k = 0; k < matrixsize; k++) {
				double value = Covariance[i](j,k);
				output.write(reinterpret_cast<char*>(&value), sizeof(double));
			}
	}
	// write the class label of each transcript
	int32_t ntrans = TransNames.size();
	output.write(reinterpret_cast<char*>(&ntrans), sizeof(int32_t));
	for (int32_t i = 0; i < ntrans; i++) {
		string tname = TransNames[i];
		int32_t namelen = tname.size();
		output.write(reinterpret_cast<char*>(&namelen), sizeof(int32_t));
		output.write(tname.c_str(), namelen*sizeof(char));
		output.write(reinterpret_cast<char*>(&(LenClass[i])), sizeof(int32_t));
	}
	output.close();
};


void ReadNewAssignment_Distribution(string inputfile, const map<string,int32_t>& TransIndex, 
	vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment)
{
	// clear variable
	AdjustmentList.clear();
	newAssignment.clear();
	// read binary file
	ifstream input(inputfile, ios::binary);
	int32_t n;
	input.read((char*)&n, sizeof(int32_t));
	for (int32_t i = 0; i < n; i++) {
		int32_t namelen;
		int32_t seqlen;
		input.read((char*)&namelen, sizeof(int32_t));
		input.read((char*)&seqlen, sizeof(int32_t));

		char* buf=new char[namelen];
		string name;
		vector<double> mydata(seqlen, 0);
		input.read(buf, namelen*sizeof(char));
		input.read((char*)(mydata.data()), seqlen*sizeof(double));
		for (int32_t j = 0; j < namelen; j++)
			name += buf[j];

		// add to variable
		map<string,int32_t>::const_iterator itmap = TransIndex.find(name);
		if (itmap == TransIndex.cend())
			cout << name << endl;
		assert(itmap != TransIndex.cend());
		AdjustmentList.push_back( itmap->second );
		newAssignment.push_back(mydata);
	}
	input.close();
};


void ReadNewAssignment_Distribution_old(string numreadsfile, string distfile, const map<string,int32_t>& TransIndex,
	vector<int32_t>& AdjustmentList, vector< vector<double> >& newAssignment)
{
	// clear variables
	AdjustmentList.clear();
	newAssignment.clear();
	vector<double> NewQuant;
	// read new quant file
	ifstream input1(numreadsfile);
	string line;
	while (getline(input1, line)) {
		if (line[0] == '#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		map<string,int32_t>::const_iterator itmap = TransIndex.find(strs[0]);
		assert(itmap != TransIndex.cend());
		AdjustmentList.push_back(itmap->second);
		NewQuant.push_back(stod(strs[2]));
	}
	input1.close();
	// read binary distribution file
	ifstream input2(distfile, ios::binary);
	int32_t numtrans;
	input2.read((char*)&numtrans, sizeof(int32_t));
	assert(numtrans == AdjustmentList.size());
	for (int32_t i = 0; i < numtrans; i++) {
		int32_t seqlen;
		input2.read((char*)&seqlen, sizeof(int32_t));
		vector<double> mydata(seqlen, 0);
		input2.read((char*)(mydata.data()), seqlen*sizeof(double));
		newAssignment.push_back(mydata);
		// sanity check: number of reads should be equal to quant
		double sum = 0;
		for (int32_t j = 0; j < mydata.size(); j++)
			sum += mydata[j];
		if (fabs(sum - NewQuant[i]) > 0.1)
			cout << "watch here\t"<< sum <<"\t"<< NewQuant[i] << endl;
		assert(fabs(sum - NewQuant[i]) < 0.1 || fabs(sum - NewQuant[i])/sum < 1e-3);
	}
	input2.close();
};


void WriteRegionalPvalue(string outputfile, const vector<Transcript_t>& Transcripts,
	const vector<double>& TransCov, const vector<PRegion_t>& RawPRegion, const vector<PRegion_t>& AdjustedPRegion)
{
	assert(RawPRegion.size() == AdjustedPRegion.size());
	ofstream output(outputfile, ios::out);
	output << "# Name\tBinStart\tBinEnd\tCoverage\tAnomalyScore\tRawPvalue\tAdjustedPvalue\n";
	for (int32_t i = 0; i < RawPRegion.size(); i++){
		const PRegion_t& rawp = RawPRegion[i];
		const PRegion_t& adjp = AdjustedPRegion[i];
		if (rawp.TID != adjp.TID)
			cout <<"watch here\n";
		assert(rawp.TID == adjp.TID);
		assert(rawp.BinStart == adjp.BinStart && rawp.BinEnd == adjp.BinEnd);
		output << (Transcripts[rawp.TID].TransID) <<"\t"<< (rawp.BinStart) <<"\t"<< (rawp.BinEnd) <<"\t";
		output << (TransCov[rawp.TID]) <<"\t"<< (rawp.AnomalyScore) <<"\t"<< (rawp.Pvalue) <<"\t"<< (adjp.Pvalue) <<"\n";
	}
};


void WriteOverallPvalue(string outputfile, const vector<Transcript_t>& Transcripts, const vector<double>& TransCov, 
	const DistTest_t& dt, const vector<double>& PValuesPos, const vector<double>& PValuesNeg, const vector<bool>& Choices)
{
	// sort p values and Transcript name
	assert(PValuesPos.size() == Choices.size());
	assert(PValuesNeg.size() == Choices.size());

	// adjust p value only for non-negative ones
	vector<double> NNPvaluesPos;
	vector<double> NNAdjPvaluePos;
	vector<double> NNPvaluesNeg;
	vector<double> NNAdjPvalueNeg;

	for (int32_t i = 0; i < PValuesPos.size(); i++){
		if (PValuesPos[i] > -0.5 && PValuesNeg[i] > -0.5) {
			NNPvaluesPos.push_back(PValuesPos[i]);
			NNPvaluesNeg.push_back(PValuesNeg[i]);
		}
	}
	BHAdjusting(NNPvaluesPos, NNAdjPvaluePos);
	BHAdjusting(NNPvaluesNeg, NNAdjPvalueNeg);
	// write non-negative p value to file
	ofstream output(outputfile, ios::out);
	int32_t count = 0;
	output << "#Name\tCoverage\tAnomalyScorePos\tRegionStartPos\tRegionEndPos\tPValue_Pos\tAdjPValue_Pos\tAnomalyScoreNeg\tRegionStartNeg\tRegionEndNeg\tPValue_Neg\tAdjPValue_Neg\tMinAdjPValue\tChoice\n";
	for (int32_t i = 0; i < PValuesPos.size(); i++) {
		if (PValuesPos[i] < -0.5 || PValuesNeg[i] < -0.5)
			continue;
		// calculate anomaly region in original coordinate
		int32_t length = dt.TransLength[i];
		vector<int32_t>::const_iterator ub = upper_bound(dt.lenBounds.cbegin(), dt.lenBounds.cend(), length, [](int a, int b){return a<=b;} ); // XXX whether upper_bound is ccorrect?
		int32_t len_class = distance(dt.lenBounds.cbegin(), ub);
		pair<int32_t, int32_t> region_pos = dt.DeletionRegion_pos[i];
		pair<int32_t, int32_t> region_neg = dt.DeletionRegion_neg[i];
		// sanity check about the binned region
		if (region_pos.first < 0 || region_pos.second < 0)
			cout <<"invalid binned region: "<< region_pos.first <<" "<< region_pos.second;
		if (region_neg.first < 0 || region_neg.second < 0)
			cout <<"invalid binned region: "<< region_neg.first <<" "<< region_neg.second;
		assert(region_pos.first >= 0 && region_pos.second >= 0);
		assert(region_neg.first >= 0 && region_neg.second >= 0);
		// convert to transcript coordinate
		pair<int32_t, int32_t> oriregion_pos = make_pair(int32_t(1.0*region_pos.first/dt.nBins[len_class]*length), int32_t(1.0*region_pos.second/dt.nBins[len_class]*length));
		pair<int32_t, int32_t> oriregion_neg = make_pair(int32_t(1.0*region_neg.first/dt.nBins[len_class]*length), int32_t(1.0*region_neg.second/dt.nBins[len_class]*length));
		// sanity check about the transcript-coordinate region
		if (oriregion_pos.first < 0 || oriregion_pos.second < 0)
			cout <<"invalid binned region: "<< oriregion_pos.first <<" "<< oriregion_pos.second;
		if (oriregion_neg.first < 0 || oriregion_neg.second < 0)
			cout <<"invalid binneds region: "<< oriregion_neg.first <<" "<< oriregion_neg.second;
		assert(oriregion_pos.first >= 0 && oriregion_pos.second >= 0);
		assert(oriregion_neg.first >= 0 && oriregion_neg.second >= 0);
		// write to file
		assert(count < NNPvaluesPos.size());
		output << Transcripts[i].TransID <<"\t"<< TransCov[i] <<"\t"<< (dt.DeletionScore_pos[i]) <<"\t"<< oriregion_pos.first <<"\t"<< oriregion_pos.second <<"\t";
		output << NNPvaluesPos[count] <<"\t"<< NNAdjPvaluePos[count] <<"\t";
		output << (dt.DeletionScore_neg[i]) <<"\t"<< oriregion_neg.first <<"\t"<< oriregion_neg.second <<"\t";
		output << NNPvaluesNeg[count] <<"\t"<< NNAdjPvalueNeg[count] <<"\t"<< min(NNAdjPvaluePos[count], NNAdjPvalueNeg[count]) <<"\t"<< Choices[i] << endl;
		// update counter
		count ++;
	}
	assert(count == NNPvaluesPos.size());
	output.close();
};


void WriteOverallPvalue(string outputfile, const vector<int32_t>& AdjustmentList, const vector<Transcript_t>& Transcripts, const vector<double>& TransCov, 
	const DistTest_t& dt, const vector<double>& PValuesPos, const vector<double>& PValuesNeg, const vector<bool>& Choices)
{
	// sort p values and Transcript name
	assert(AdjustmentList.size() == PValuesPos.size());
	assert(AdjustmentList.size() == PValuesNeg.size());
	assert(PValuesPos.size() == Choices.size());
	assert(PValuesNeg.size() == Choices.size());

	// adjust p value only for non-negative ones
	vector<double> NNPvaluesPos;
	vector<double> NNAdjPvaluePos;
	vector<double> NNPvaluesNeg;
	vector<double> NNAdjPvalueNeg;

	for (int32_t i = 0; i < PValuesPos.size(); i++){
		if (PValuesPos[i] > -0.5 && PValuesNeg[i] > -0.5) {
			NNPvaluesPos.push_back(PValuesPos[i]);
			NNPvaluesNeg.push_back(PValuesNeg[i]);
		}
	}
	BHAdjusting(NNPvaluesPos, NNAdjPvaluePos);
	BHAdjusting(NNPvaluesNeg, NNAdjPvalueNeg);
	// write non-negative p value to file
	ofstream output(outputfile, ios::out);
	int32_t count = 0;
	output << "#Name\tCoverage\tAnomalyScorePos\tRegionStartPos\tRegionEndPos\tPValue_Pos\tAdjPValue_Pos\tAnomalyScoreNeg\tRegionStartNeg\tRegionEndNeg\tPValue_Neg\tAdjPValue_Neg\tMinAdjPValue\tChoice\n";
	for (int32_t i = 0; i < PValuesPos.size(); i++) {
		if (PValuesPos[i] < -0.5 || PValuesNeg[i] < -0.5)
			continue;
		// calculate anomaly region in original coordinate
		int32_t ind = AdjustmentList[i];
		int32_t length = dt.TransLength[ind];
		vector<int32_t>::const_iterator ub = upper_bound(dt.lenBounds.cbegin(), dt.lenBounds.cend(), length, [](int a, int b){return a<=b;} ); // XXX whether upper_bound is ccorrect?
		int32_t len_class = distance(dt.lenBounds.cbegin(), ub);
		pair<int32_t, int32_t> region_pos = dt.DeletionRegion_pos[ind];
		pair<int32_t, int32_t> region_neg = dt.DeletionRegion_neg[ind];
		// sanity check about the binned region
		if (region_pos.first < 0 || region_pos.second < 0)
			cout <<"invalid binned region: "<< region_pos.first <<" "<< region_pos.second;
		if (region_neg.first < 0 || region_neg.second < 0)
			cout <<"invalid binned region: "<< region_neg.first <<" "<< region_neg.second;
		assert(region_pos.first >= 0 && region_pos.second >= 0);
		assert(region_neg.first >= 0 && region_neg.second >= 0);
		// convert to transcript coordinate
		pair<int32_t, int32_t> oriregion_pos = make_pair(int32_t(1.0*region_pos.first/dt.nBins[len_class]*length), int32_t(1.0*region_pos.second/dt.nBins[len_class]*length));
		pair<int32_t, int32_t> oriregion_neg = make_pair(int32_t(1.0*region_neg.first/dt.nBins[len_class]*length), int32_t(1.0*region_neg.second/dt.nBins[len_class]*length));
		// sanity check about the transcript-coordinate region
		if (oriregion_pos.first < 0 || oriregion_pos.second < 0)
			cout <<"invalid binned region: "<< oriregion_pos.first <<" "<< oriregion_pos.second;
		if (oriregion_neg.first < 0 || oriregion_neg.second < 0)
			cout <<"invalid binneds region: "<< oriregion_neg.first <<" "<< oriregion_neg.second;
		assert(oriregion_pos.first >= 0 && oriregion_pos.second >= 0);
		assert(oriregion_neg.first >= 0 && oriregion_neg.second >= 0);
		// write to file
		assert(count < NNPvaluesPos.size());
		output << Transcripts[ind].TransID <<"\t"<< TransCov[ind] <<"\t"<< (dt.DeletionScore_pos[ind]) <<"\t"<< oriregion_pos.first <<"\t"<< oriregion_pos.second <<"\t";
		output << NNPvaluesPos[count] <<"\t"<< NNAdjPvaluePos[count] <<"\t";
		output << (dt.DeletionScore_neg[ind]) <<"\t"<< oriregion_neg.first <<"\t"<< oriregion_neg.second <<"\t";
		output << NNPvaluesNeg[count] <<"\t"<< NNAdjPvalueNeg[count] <<"\t"<< min(NNAdjPvaluePos[count], NNAdjPvalueNeg[count]) <<"\t"<< Choices[i] << endl;
		// update counter
		count ++;
	}
	assert(count == NNPvaluesPos.size());
	output.close();
};


void OldWriteOverallPvalue(string outputfile, const vector<Transcript_t>& Transcripts, const vector<double>& TransCov, 
	const DistTest_t& dt, const vector<double>& PValues, const vector<bool>& Choices)
{
	assert(PValues.size() == Choices.size());
	// write non-negative p value to file
	ofstream output(outputfile, ios::out);
	output << "#Name\tCoverage\tDeletionScorePos\tDeletioScoreNeg\tRawPvalue\tAdjustedPvalue\tChoice\n";
	for (int32_t i = 0; i < PValues.size(); i++) {
		output << Transcripts[i].TransID <<"\t"<< TransCov[i] <<"\t"<< (dt.DeletionScore_pos[i]) <<"\t";
		output << (dt.DeletionScore_neg[i]) <<"\t"<< PValues[i] <<"\t.\t"<< Choices[i] << endl;
	}
	output.close();
};
