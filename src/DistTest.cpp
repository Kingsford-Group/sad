#include "DistTest.hpp"
#include "LPReassign.hpp"

using namespace std;

DistTest_t::DistTest_t(const map<string,int32_t>& _TransIndex, const map<string,int32_t>& _TransLength, 
	const vector< vector<double> >& Expected, const vector< vector<double> >& Observed)
{
	TransIndex = _TransIndex;
	percentile = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1};
	TransLength.resize(TransIndex.size(), 0);
	for (map<string,int32_t>::const_iterator it = TransIndex.cbegin(); it != TransIndex.cend(); it++){
		map<string,int32_t>::const_iterator itlen = _TransLength.find(it->first);
		assert(itlen != _TransLength.cend());
		TransLength[it->second] = itlen->second;
	}
	TransNames.resize(TransIndex.size(), "");
	for (map<string,int32_t>::const_iterator it = TransIndex.cbegin(); it != TransIndex.cend(); it++)
		TransNames[ it->second ] = it->first;
	TransCov.clear();
	// sanity check
	for (int32_t i = 0; i < TransIndex.size(); i++){
		assert(TransNames[i] != "");
		assert(TransLength[i] > 0);
	}
	// calculate lenBounds and nBins
	lenBounds.clear();
	nBins.clear();
	vector<int32_t> sortedTransLength(TransLength.begin(), TransLength.end());
	sort(sortedTransLength.begin(), sortedTransLength.end());
	for (int32_t i = 0; i < percentile.size(); i++) {
		// lenBounds is the largest length of each group
		int32_t ind = sortedTransLength.size()*percentile[i] - 1;
		lenBounds.push_back(sortedTransLength[ind]);
		// binning size is calculated based on middle length of each group
		double midper = percentile[i]/2;
		if (i != 0)
			midper = (percentile[i-1]+percentile[i])/2;
		int32_t ind_mid = sortedTransLength.size()*midper;
		nBins.push_back(int32_t(round(1.0*sortedTransLength[ind_mid]/25)));
		if (nBins.back() == 0)
			nBins.back() = 1;
	}
	// binning and normalizeing Expected and Observed
	ExpectedBinNorm.clear();
	ObservedBinNorm.clear();
	LenClass.clear();
	assert(Expected.size() == Observed.size());
	for (int32_t i = 0; i < Expected.size(); i++) {
		// identify len_class
		assert(Expected[i].size() == Observed[i].size());
		int32_t length = Expected[i].size();
		vector<int32_t>::iterator ub = upper_bound(lenBounds.begin(), lenBounds.end(), length, [](int a, int b){return a<=b;} ); // XXX whether upper_bound is ccorrect?
		int32_t len_class = distance(lenBounds.begin(), ub);
		// binsize based on corresponding len_class and nBins
		vector<double> tmpexp(nBins[len_class], 0);
		vector<double> tmpobs(nBins[len_class], 0);
		int32_t current_bin = 0;
		double current_bin_end = 1.0*(current_bin+1)*length/nBins[len_class];
		for (int32_t j = 0; j < Expected[i].size(); j++){
			if (j > current_bin_end){
				current_bin++;
				current_bin_end = 1.0*(current_bin+1)*length/nBins[len_class];
				assert(current_bin < nBins[len_class]);
			}
			tmpexp[current_bin] += Expected[i][j];
			tmpobs[current_bin] += Observed[i][j];
		}
		// normalize
		double sumexp = 0;
		double sumobs = 0;
		for (int32_t j = 0; j < tmpexp.size(); j++){
			sumexp += tmpexp[j];
			sumobs += tmpobs[j];
		}
		if (sumexp > 0){
			for (int32_t j = 0; j < tmpexp.size(); j++)
				tmpexp[j] /= sumexp;
		}
		if (sumobs > 0){
			for (int32_t j = 0; j < tmpobs.size(); j++)
				tmpobs[j] /= sumobs;
		}
		// push back
		Eigen::VectorXd tmpexp_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(tmpexp.data(), tmpexp.size());
		Eigen::VectorXd tmpobs_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(tmpobs.data(), tmpobs.size());
		ExpectedBinNorm.push_back(tmpexp_eigen);
		ObservedBinNorm.push_back(tmpobs_eigen);
		LenClass.push_back(len_class);
		TransCov.push_back(sumobs / TransLength[i]);
		// if (sumexp == 0)
		// 	TransCov.push_back(0);
		// else
		// 	TransCov.push_back(sumobs / sumexp);
	}
	assert(ExpectedBinNorm.size() == LenClass.size() && ObservedBinNorm.size() == LenClass.size());

	// initialize Sampling_NSTDS
	Sampling_NSTDS.clear();
	for (int32_t i = -25; i < -15; i+=2)
		Sampling_NSTDS.push_back(1.0*i/10);
	for (int32_t i = -15; i < 15; i++)
		Sampling_NSTDS.push_back(1.0*i/10);
	for (int32_t i = 15; i < 26; i+=2)
		Sampling_NSTDS.push_back(1.0*i/10);
	// initialize Sampling_WEIGHTS
	Sampling_WEIGHTS.clear();
	boost::math::normal s;
	double sumweights = 0;
	for(int32_t i = 0; i < Sampling_NSTDS.size(); i++){
		Sampling_WEIGHTS.push_back( boost::math::pdf(s, Sampling_NSTDS[i]) );
		sumweights += Sampling_WEIGHTS.back();
	}
	for(int32_t i = 0; i < Sampling_WEIGHTS.size(); i++)
		Sampling_WEIGHTS[i] /= sumweights;

	// set default number of threads
	Num_Threads = 4;
};


DistTest_t::DistTest_t(const map<string,int32_t>& _TransIndex, const map<string,int32_t>& _TransLength, 
	const vector< vector<double> >& Expected, const vector< vector<double> >& Observed, int32_t _Num_Threads)
{
	DistTest_t(_TransIndex, _TransLength, Expected, Observed);
	Num_Threads = _Num_Threads;
};


void DistTest_t::AdjustExpected(double nstdthresh)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Adjust expected distribution."<<endl;

	// group expressed transcripts based on LenClass
	vector< vector<int32_t> > lenTransIndexes;
	for (int32_t i = 0; i < lenBounds.size(); i++){
		vector<int32_t> tmp;
		tmp.reserve(int32_t(TransIndex.size()*0.1));
		lenTransIndexes.push_back(tmp);
	}
	for (int32_t i = 0; i < ExpectedBinNorm.size(); i++){
		if (TransCov[i] > 0)
			lenTransIndexes[ LenClass[i] ].push_back(i);
	}
	// clear Mean and covariance variable
	// for each len_class, len(Mean) = nBins[len_class]; Covariance.shape = nBins[len_class]; CholeskyDecom.shape = nBins[len_class]-1
	Mean.clear();
	Covariance.clear();
	CholeskyDecom.clear();
	// for each LenClass, calculate the mean / std of expression, select the highly expressed to calculate shift
	for (int32_t i = 0; i < lenBounds.size(); i++){
		// calculate logTPM threshold for high expression
		Eigen::VectorXd logTPM(lenTransIndexes[i].size());
		for (int32_t j = 0; j < lenTransIndexes[i].size(); j++){
			int32_t ind = lenTransIndexes[i][j];
			double tpm = TransCov[ind];
			logTPM(j) = std::log(tpm);
		}
		double logTPM_mean = logTPM.mean();
		double logTPM_std = std::sqrt( (logTPM.array()-logTPM.mean()).square().sum() / (logTPM.size()-1) );
		double logTPM_thresh = logTPM_mean + nstdthresh * logTPM_std;
		// make a matrix of shifts for high expression transcripts
		vector<int32_t> passedIndex;
		for (int32_t j = 0; j < lenTransIndexes[i].size(); j++){
			if (logTPM(j) > logTPM_thresh)
				passedIndex.push_back(lenTransIndexes[i][j]);
		}
		assert(passedIndex.size() > 10);
		assert(ObservedBinNorm[passedIndex[0]].size() == nBins[i]);
		Eigen::MatrixXd Shifts(passedIndex.size(), nBins[i] - 1);
		for (int32_t j = 0; j < passedIndex.size(); j++){
			assert(fabs(ObservedBinNorm[passedIndex[j]].sum() - 1) < 1e-8);
			assert(fabs(ExpectedBinNorm[passedIndex[j]].sum() - 1) < 1e-8);
			Shifts.row(j) = (ObservedBinNorm[passedIndex[j]] - ExpectedBinNorm[passedIndex[j]]).head(nBins[i] - 1);
			assert(fabs(Shifts.row(j).sum()) < 2);
		}
		// calculate mean shifts of this LenClass
		Eigen::VectorXd mean_shift = Shifts.colwise().mean();
		mean_shift.conservativeResize(nBins[i]);
		mean_shift(Shifts.cols()) = 0;
		mean_shift(Shifts.cols()) = -mean_shift.sum();
		assert(fabs(mean_shift.sum()) <= 1e-8);
		Mean.push_back(mean_shift);
		assert(Mean.back().size() == nBins[i]);
		// calculate covariance shifts
		Eigen::MatrixXd Shifts_meancenter = Shifts.rowwise() - Shifts.colwise().mean();
		Eigen::MatrixXd cov_shift = Shifts_meancenter.transpose() * Shifts_meancenter / (passedIndex.size() - 1);
		// add a small diagonal matrix to increase matrix stability
		cov_shift += Eigen::MatrixXd::Identity(cov_shift.rows(), cov_shift.cols()) * 1e-18;
		// calculate cholesky decomposition of covariance shifts
		Eigen::LLT<Eigen::MatrixXd> lltOfA(cov_shift);
		Eigen::MatrixXd chol_L = lltOfA.matrixL();
		CholeskyDecom.push_back(chol_L);
		// resize covariance shifts and add to member attributes
		cov_shift.conservativeResize(nBins[i], nBins[i]);
		for(int32_t j = 0; j < nBins[i]; j++){
			cov_shift(j, nBins[i]-1) = 0;
			cov_shift(nBins[i]-1, j) = 0;
		}
		Covariance.push_back(cov_shift);
		// adjust all ExpectedBinNorm in this LenClass
		for (int32_t j = 0; j < ExpectedBinNorm.size(); j++){
			if (LenClass[j] == i){
				ExpectedBinNorm[j] += Mean.back();
				bool needrenormal = false;
				for (int32_t k = 0; k < nBins[i]; k++)
					if (ExpectedBinNorm[j](k) < 0){
						ExpectedBinNorm[j](k) = 1e-8;
						needrenormal = true;
					}
				if (needrenormal)
					ExpectedBinNorm[j] /= (ExpectedBinNorm[j].sum());
				assert(fabs(ExpectedBinNorm[j].sum() - 1) < 1e-8);
			}
		}
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish Adjustment.\n";
};


void DistTest_t::UpdateObserved(const vector< vector<double> >& NewObserved)
{
	assert(NewObserved.size() == ObservedBinNorm.size());
	for (int32_t i = 0; i < NewObserved.size(); i++){
		assert(NewObserved[i].size() == TransLength[i]);
		int32_t len_class = LenClass[i];
		// binning and calculate sum
		int32_t length = NewObserved[i].size();
		double sumobs = 0;
		vector<double> tmpobs(nBins[len_class], 0);
		int32_t current_bin = 0;
		double current_bin_end = 1.0*(current_bin+1)*length/nBins[len_class];
		for (int32_t j = 0; j < NewObserved[i].size(); j++){
			if (j > current_bin_end){
				current_bin++;
				current_bin_end = 1.0*(current_bin+1)*length/nBins[len_class];
				assert(current_bin < nBins[len_class]);
			}
			tmpobs[current_bin] += NewObserved[i][j];
			sumobs += NewObserved[i][j];
		}
		// normalize
		if (sumobs > 0){
			for (int32_t j = 0; j < tmpobs.size(); j++)
				tmpobs[j] /= sumobs;
		}
		// update
		Eigen::VectorXd tmpobs_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(tmpobs.data(), tmpobs.size());
		ObservedBinNorm[i] = tmpobs_eigen;
		TransCov[i] = sumobs / TransLength[i];
	}
};


void DistTest_t::UpdateObserved(const vector<int32_t>& AdjustmentList, const vector< vector<double> >& NewObserved)
{
	assert(NewObserved.size() == AdjustmentList.size());
	for (int32_t i = 0; i < AdjustmentList.size(); i++){
		assert(NewObserved[i].size() == TransLength[AdjustmentList[i]]);
		int32_t len_class = LenClass[AdjustmentList[i]];
		// binning and calculate sum
		int32_t length = NewObserved[i].size();
		double sumobs = 0;
		vector<double> tmpobs(nBins[len_class], 0);
		int32_t current_bin = 0;
		double current_bin_end = 1.0*(current_bin+1)*length/nBins[len_class];
		for (int32_t j = 0; j < NewObserved[i].size(); j++){
			if (j > current_bin_end){
				current_bin++;
				current_bin_end = 1.0*(current_bin+1)*length/nBins[len_class];
				assert(current_bin < nBins[len_class]);
			}
			tmpobs[current_bin] += NewObserved[i][j];
			sumobs += NewObserved[i][j];
		}
		// normalize
		if (sumobs > 0){
			for (int32_t j = 0; j < tmpobs.size(); j++)
				tmpobs[j] /= sumobs;
		}
		// update
		Eigen::VectorXd tmpobs_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(tmpobs.data(), tmpobs.size());
		ObservedBinNorm[AdjustmentList[i]] = tmpobs_eigen;
		TransCov[AdjustmentList[i]] = sumobs / TransLength[AdjustmentList[i]];
	}
};


void DistTest_t::CalDeletionScore(double covthresh, double numthresh)
{
	// clear variables: pos = exp - obs; neg = obs - exp
	DeletionScore_pos.assign(ExpectedBinNorm.size(), -1);
	DeletionRegion_pos.assign(ExpectedBinNorm.size(), make_pair(-1,-1));
	DeletionScore_neg.assign(ExpectedBinNorm.size(), -1);
	DeletionRegion_neg.assign(ExpectedBinNorm.size(), make_pair(-1,-1));
	// calculate exp - obs
	for (int32_t i = 0; i < ExpectedBinNorm.size(); i++){
		if (TransCov[i] < covthresh || (TransCov[i]*TransLength[i]) < numthresh)
			continue;
		Eigen::VectorXd Diff = ExpectedBinNorm[i] - ObservedBinNorm[i];
		double maxdiff = 0;
		double currentdiff = 0;
		int32_t maxregionstart = -1, maxregionend = -1;
		int32_t currentregionstart = 0, currentregionend = 0;
		for (int32_t cursor = 0; cursor < Diff.size(); cursor++){
			currentdiff += Diff(cursor);
			if (currentdiff < 0){
				currentregionstart = cursor + 1;
				currentdiff = 0;
			}
			else{
				currentregionend = cursor + 1;
				if (currentdiff > maxdiff){
					maxdiff = currentdiff;
					maxregionstart = currentregionstart;
					maxregionend = currentregionend;
				}
			}
		}
		if (maxdiff > 0){
			assert(maxregionstart >= 0 && maxregionend <= Diff.size() && maxregionstart < maxregionend);
			assert(maxdiff < 1);
			DeletionScore_pos[i] = maxdiff;
			DeletionRegion_pos[i] = make_pair(maxregionstart, maxregionend);
		}
	}
	// calculate obs - exp
	for (int32_t i = 0; i < ObservedBinNorm.size(); i++){
		if (TransCov[i] < covthresh || (TransCov[i]*TransLength[i]) < numthresh)
			continue;
		Eigen::VectorXd Diff = ObservedBinNorm[i] - ExpectedBinNorm[i];
		double maxdiff = 0;
		double currentdiff = 0;
		int32_t maxregionstart = -1, maxregionend = -1;
		int32_t currentregionstart = 0, currentregionend = 0;
		for (int32_t cursor = 0; cursor < Diff.size(); cursor++){
			currentdiff += Diff(cursor);
			if (currentdiff < 0){
				currentregionstart = cursor + 1;
				currentdiff = 0;
			}
			else{
				currentregionend = cursor + 1;
				if (currentdiff > maxdiff){
					maxdiff = currentdiff;
					maxregionstart = currentregionstart;
					maxregionend = currentregionend;
				}
			}
		}
		if (maxdiff > 0){
			assert(maxregionstart >= 0 && maxregionend <= Diff.size() && maxregionstart < maxregionend);
			assert(maxdiff < 1);
			DeletionScore_neg[i] = maxdiff;
			DeletionRegion_neg[i] = make_pair(maxregionstart, maxregionend);
		}
	}
};


void DistTest_t::CalDeletionScore(const vector<int32_t>& AdjustmentList, double covthresh, double numthresh)
{
	// sanity check: DeletionScore of salmon is calculated
	assert(DeletionScore_pos.size() == ExpectedBinNorm.size());
	assert(DeletionRegion_pos.size() == ExpectedBinNorm.size());
	assert(DeletionScore_neg.size() == ExpectedBinNorm.size());
	assert(DeletionRegion_neg.size() == ExpectedBinNorm.size());
	// calculate exp - obs
	for (int32_t tmpi = 0; tmpi < AdjustmentList.size(); tmpi++){
		int32_t i = AdjustmentList[tmpi];
		if (TransCov[i] < covthresh || (TransCov[i]*TransLength[i]) < numthresh) {
			DeletionScore_pos[i] = -1;
			DeletionRegion_pos[i] = make_pair(-1, -1);
			DeletionScore_neg[i] = -1;
			DeletionRegion_neg[i] = make_pair(-1, -1);
			continue;
		}
		Eigen::VectorXd Diff = ExpectedBinNorm[i] - ObservedBinNorm[i];
		double maxdiff = 0;
		double currentdiff = 0;
		int32_t maxregionstart = -1, maxregionend = -1;
		int32_t currentregionstart = 0, currentregionend = 0;
		for (int32_t cursor = 0; cursor < Diff.size(); cursor++){
			currentdiff += Diff(cursor);
			if (currentdiff < 0){
				currentregionstart = cursor + 1;
				currentdiff = 0;
			}
			else{
				currentregionend = cursor + 1;
				if (currentdiff > maxdiff){
					maxdiff = currentdiff;
					maxregionstart = currentregionstart;
					maxregionend = currentregionend;
				}
			}
		}
		if (maxdiff > 0){
			assert(maxregionstart >= 0 && maxregionend <= Diff.size() && maxregionstart < maxregionend);
			assert(maxdiff < 1);
			DeletionScore_pos[i] = maxdiff;
			DeletionRegion_pos[i] = make_pair(maxregionstart, maxregionend);
		}
	}
	// calculate obs - exp
	for (int32_t tmpi = 0; tmpi < AdjustmentList.size(); tmpi++){
		int32_t i = AdjustmentList[tmpi];
		if (TransCov[i] < covthresh || (TransCov[i]*TransLength[i]) < numthresh)
			continue;
		Eigen::VectorXd Diff = ObservedBinNorm[i] - ExpectedBinNorm[i];
		double maxdiff = 0;
		double currentdiff = 0;
		int32_t maxregionstart = -1, maxregionend = -1;
		int32_t currentregionstart = 0, currentregionend = 0;
		for (int32_t cursor = 0; cursor < Diff.size(); cursor++){
			currentdiff += Diff(cursor);
			if (currentdiff < 0){
				currentregionstart = cursor + 1;
				currentdiff = 0;
			}
			else{
				currentregionend = cursor + 1;
				if (currentdiff > maxdiff){
					maxdiff = currentdiff;
					maxregionstart = currentregionstart;
					maxregionend = currentregionend;
				}
			}
		}
		if (maxdiff > 0){
			assert(maxregionstart >= 0 && maxregionend <= Diff.size() && maxregionstart < maxregionend);
			assert(maxdiff < 1);
			DeletionScore_neg[i] = maxdiff;
			DeletionRegion_neg[i] = make_pair(maxregionstart, maxregionend);
		}
	}
};

double DistTest_t::CalDeletionScore_single(Eigen::VectorXd Dist1, Eigen::VectorXd Dist2, pair<int32_t,int32_t>& region, 
	bool check1, double numthresh)
{
	// check numreads threshold
	double DelScore = -1;
	region = make_pair(-1, -1);
	if (check1 && Dist1.sum() < numthresh)
		return DelScore;
	else if (!check1 && Dist2.sum() < numthresh)
		return DelScore;

	Dist1 /= Dist1.sum();
	Dist2 /= Dist2.sum();
	Eigen::VectorXd Diff = Dist1 - Dist2;
	// calculate region with the largest difference
	double maxdiff = 0;
	double currentdiff = 0;
	int32_t maxregionstart = -1, maxregionend = -1;
	int32_t currentregionstart = 0, currentregionend = 0;
	for (int32_t cursor = 0; cursor < Diff.size(); cursor++){
		currentdiff += Diff(cursor);
		if (currentdiff < 0){
			currentregionstart = cursor + 1;
			currentdiff = 0;
		}
		else{
			currentregionend = cursor + 1;
			if (currentdiff > maxdiff){
				maxdiff = currentdiff;
				maxregionstart = currentregionstart;
				maxregionend = currentregionend;
			}
		}
	}
	if (maxdiff > 0){
		assert(maxregionstart >= 0 && maxregionend <= Diff.size() && maxregionstart < maxregionend);
		DelScore = maxdiff;
		region = make_pair(maxregionstart, maxregionend);
	}
	return DelScore;
};


vector<PRegion_t> DistTest_t::SinglePvalue_regional_pos(int32_t ind)
{
	vector<PRegion_t> PVs;
	// for pos = exp - obs
	if (DeletionScore_pos[ind] > 0){
		int32_t Spacing = DeletionRegion_pos[ind].second - DeletionRegion_pos[ind].first;
		int32_t len_class = LenClass[ind];
		for (int32_t binstart = 0; binstart < ExpectedBinNorm[ind].size(); binstart++) {
			int32_t binend = binstart + Spacing;
			if (binend >= ExpectedBinNorm[ind].size())
				break;
			// theoretical probability
			double prob_theo = ExpectedBinNorm[ind].segment(binstart, binend-binstart).sum();
			double rawdiff = prob_theo - ObservedBinNorm[ind].segment(binstart, binend-binstart).sum();
			double diff = fabs(rawdiff);
			double variance = Covariance[len_class].block(binstart, binstart, binend-binstart, binend-binstart).sum();
			assert(variance > -1e-4);
			// total number of reads and extreme number rod reads
			double numreads_total = TransCov[ind] * TransLength[ind];
			int32_t numreads_expected = numreads_total * prob_theo;
			int32_t numreads_extreme_lb = numreads_total * (prob_theo - diff);
			int32_t numreads_extreme_ub = ceil(numreads_total * (prob_theo + diff));
			assert(numreads_total >= 0);
			// calculate 2-sided probability
			double tmppvalue_lb = 0;
			double tmppvalue_ub = 0;
			if (variance > 0){
				double sumeffect_weight_lb = 0;
				double sumeffect_weight_ub = 0;
				for(int32_t i = 0; i < Sampling_NSTDS.size(); i++) {
					double prob_sample = prob_theo + Sampling_NSTDS[i] * std::sqrt(variance);
					if (prob_sample <= 0 || prob_sample >= 1)
						continue;
					assert(int32_t(round(numreads_total)) >= 0);
					boost::math::binomial b(int32_t(round(numreads_total)), prob_sample);
					// lower prob: P(reads <= numreads_extreme_lb)
					if (numreads_extreme_lb >= 0){
						tmppvalue_lb += boost::math::cdf(b, numreads_extreme_lb) * Sampling_WEIGHTS[i];
						// tmppvalue_lb += boost::math::cdf(b, numreads_extreme_lb)/boost::math::cdf(b, numreads_expected) * Sampling_WEIGHTS[i];
						sumeffect_weight_lb += Sampling_WEIGHTS[i];
					}
					// upper prob: P(reads >= numreads_extreme_ub)
					// complement calculated P(reads > numreads_extreme_ub), so we need to reduce numreads_extreme_ub by 1 to include equality
					if (numreads_extreme_ub <= numreads_total){
						tmppvalue_ub += boost::math::cdf(complement(b, max(0,numreads_extreme_ub-1))) * Sampling_WEIGHTS[i];
						sumeffect_weight_ub += Sampling_WEIGHTS[i];
					}
				}
				if (sumeffect_weight_lb > 0)
					tmppvalue_lb /= sumeffect_weight_lb;
				if (sumeffect_weight_ub > 0)
					tmppvalue_ub /= sumeffect_weight_ub;
			}
			PRegion_t tmp(ind, binstart, binend, rawdiff, tmppvalue_lb+tmppvalue_ub);
			PVs.push_back(tmp);
		}
	}
	return PVs;
};


vector<PRegion_t> DistTest_t::SinglePvalue_regional_neg(int32_t ind)
{
	vector<PRegion_t> PVs;
	// for neg = obs - exp
	if (DeletionScore_neg[ind] > 0){
		int32_t Spacing = DeletionRegion_neg[ind].second - DeletionRegion_neg[ind].first;
		int32_t len_class = LenClass[ind];
		for (int32_t binstart = 0; binstart < ExpectedBinNorm[ind].size(); binstart++) {
			int32_t binend = binstart + Spacing;
			if (binend >= ExpectedBinNorm[ind].size())
				break;
			// theoretical probability
			double prob_theo = ExpectedBinNorm[ind].segment(binstart, binend-binstart).sum();
			double rawdiff = ObservedBinNorm[ind].segment(binstart, binend-binstart).sum() - prob_theo;
			double diff = fabs(rawdiff);
			double variance = Covariance[len_class].block(binstart, binstart, binend-binstart, binend-binstart).sum();
			assert(variance > -1e-4);
			// total number of reads and extreme number rod reads
			double numreads_total = TransCov[ind] * TransLength[ind];
			int32_t numreads_extreme_lb = numreads_total * (prob_theo - diff);
			int32_t numreads_extreme_ub = ceil(numreads_total * (prob_theo + diff));
			assert(numreads_total >= 0);
			// calculate 2-sided probability
			double tmppvalue_lb = 0;
			double tmppvalue_ub = 0;
			if (variance > 0){
				double sumeffect_weight_lb = 0;
				double sumeffect_weight_ub = 0;
				for(int32_t i = 0; i < Sampling_NSTDS.size(); i++) {
					double prob_sample = prob_theo + Sampling_NSTDS[i] * std::sqrt(variance);
					if (prob_sample <= 0 || prob_sample >= 1)
						continue;
					assert(int32_t(round(numreads_total)) >= 0);
					boost::math::binomial b(int32_t(round(numreads_total)), prob_sample);
					// lower prob: P(reads <= numreads_extreme_lb)
					if (numreads_extreme_lb >= 0){
						tmppvalue_lb += boost::math::cdf(b, numreads_extreme_lb) * Sampling_WEIGHTS[i];
						sumeffect_weight_lb += Sampling_WEIGHTS[i];
					}
					// upper prob: P(reads >= numreads_extreme_ub)
					// complement calculated P(reads > numreads_extreme_ub), so we need to reduce numreads_extreme_ub by 1 to include equality
					if (numreads_extreme_ub <= numreads_total){
						tmppvalue_ub += boost::math::cdf(complement(b, max(0,numreads_extreme_ub-1))) * Sampling_WEIGHTS[i];
						sumeffect_weight_ub += Sampling_WEIGHTS[i];
					}
				}
				if (sumeffect_weight_lb > 0)
					tmppvalue_lb /= sumeffect_weight_lb;
				if (sumeffect_weight_ub > 0)
					tmppvalue_ub /= sumeffect_weight_ub;
			}
			PRegion_t tmp(ind, binstart, binend, rawdiff, tmppvalue_lb+tmppvalue_ub);
			PVs.push_back(tmp);
		}
	}
	return PVs;
};


pair<double,double> DistTest_t::SinglePvalue_overall(int32_t ind)
{
	double pvalue_pos = -1;
	double pvalue_neg = -1;
	bool choice = 0;
	// for pos = exp - obs
	if (DeletionScore_pos[ind] > 0){
		pvalue_pos = 0;
		double diff = DeletionScore_pos[ind];
		pair<int32_t,int32_t> region = DeletionRegion_pos[ind];
		// theoretical probability
		int32_t len_class = LenClass[ind];
		double numreads_total = TransCov[ind] * TransLength[ind];
		bool overbound = false;
		for (int32_t i = 0; i < nBins[len_class]; i++){
			for (int32_t j = i+max(1, int32_t((region.second-region.first)/2)); j < nBins[len_class]; j++){
				double prob_theo = ExpectedBinNorm[ind].segment(i, j-i).sum();
				double variance = Covariance[len_class].block(i, i, j-i, j-i).sum();
				assert(variance > -1e-4);
				// extreme number of reads
				int32_t numreads_extreme = numreads_total * (prob_theo - diff);
				if (variance > 0 && numreads_extreme >= 0){
					for (int32_t k = 0; k < Sampling_NSTDS.size(); k++){
						double prob_sample = prob_theo + Sampling_NSTDS[k] * std::sqrt(variance);
						if (prob_sample <= 0 || prob_sample >= 1)
							continue;
						boost::math::binomial b(int32_t(round(numreads_total)), prob_sample);
						double tmppvalue_pos = boost::math::cdf(b, numreads_extreme) * Sampling_WEIGHTS[k];
						pvalue_pos = (tmppvalue_pos > pvalue_pos) ? tmppvalue_pos : pvalue_pos;
						if (pvalue_pos >= 0.999){
							overbound = true;
							break;
						}
					}
				}
				if (overbound)
					break;
			}
			if (overbound)
				break;
		}
	}
	// for neg = obs - exp
	if (DeletionScore_neg[ind] > 0){
		pvalue_neg = 0;
		double diff = DeletionScore_neg[ind];
		pair<int32_t,int32_t> region = DeletionRegion_neg[ind];
		// theoretical probability
		int32_t len_class = LenClass[ind];
		double numreads_total = TransCov[ind] * TransLength[ind];
		bool overbound = false;
		for (int32_t i = 0; i < nBins[len_class]; i++){
			for (int32_t j = i+max(1, int32_t(region.second-region.first)/2); j < nBins[len_class]; j++){
				double prob_theo = ExpectedBinNorm[ind].segment(i, j-i).sum();
				double variance = Covariance[len_class].block(i, i, j-i, j-i).sum();
				assert(variance > -1e-4);
				// extreme number of reads
				int32_t numreads_extreme = std::min(numreads_total, ceil(numreads_total * (prob_theo + diff)));
				if (variance > 0){
					for (int32_t k = 0; k < Sampling_NSTDS.size(); k++){
						double prob_sample = prob_theo + Sampling_NSTDS[k] * std::sqrt(variance);
						if (prob_sample <= 0 || prob_sample >= 1)
							continue;
						boost::math::binomial b(int32_t(round(numreads_total)), prob_sample);
						double tmppvalue_neg = boost::math::cdf(complement(b, numreads_extreme)) * Sampling_WEIGHTS[k];
						pvalue_neg = (tmppvalue_neg > pvalue_neg) ? tmppvalue_neg : pvalue_neg;
						if (pvalue_neg >= 0.999){
							overbound = true;
							break;
						}
					}
				}
				if (overbound)
					break;
			}
			if (overbound)
				break;
		}
	}

	return make_pair(pvalue_pos, pvalue_neg);
};


pair<double, double> DistTest_t::SinglePvalue_overall_empirical(int32_t ind, int32_t num_sampling, const gsl_rng * r)
{
	double pvalue_pos = 1;
	double pvalue_neg = 1;
	int32_t len_class = LenClass[ind];
	double diff_pos = DeletionScore_pos[ind];
	double diff_neg = DeletionScore_neg[ind];
	int32_t numreads = int32_t(round(TransCov[ind] * TransLength[ind]));
	// check threshold
	if (DeletionScore_pos[ind] < 0 || DeletionScore_neg[ind] < 0)
		return make_pair(-1, -1);
	// // tmp: write simulated files
	// ofstream output1("tmp_choL", ios::out);
	// ofstream output2("tmp_expected", ios::out);
	// ofstream output3("tmp_norm", ios::out);
	// ofstream output4("tmp_sampledobs", ios::out);
	// for (int32_t i = 0; i < CholeskyDecom[len_class].rows(); i++){
	// 	for (int32_t j = 0; j < CholeskyDecom[len_class].cols()-1; j++)
	// 		output1 << CholeskyDecom[len_class](i,j) <<" ";
	// 	output1 << CholeskyDecom[len_class](i, CholeskyDecom[len_class].cols()-1) << endl;
	// }
	// for (int32_t i = 0; i < ExpectedBinNorm[ind].size()-1; i++)
	// 	output2 << ExpectedBinNorm[ind](i) <<" ";
	// output2 << ExpectedBinNorm[ind](ExpectedBinNorm[ind].size()-1) << endl;
	// sample from standard gaussian distribution
	Eigen::MatrixXd norm = Eigen::MatrixXd::Zero(nBins[len_class]-1, num_sampling); 
	for (int32_t i = 0; i < nBins[len_class]-1; i++) {
		for (int32_t j = 0; j < num_sampling; j++) {
			norm(i,j) = gsl_ran_gaussian(r, 1.0);
		}
	}
	Eigen::MatrixXd rand = (CholeskyDecom[len_class]*norm).colwise() + ExpectedBinNorm[ind].segment(0, nBins[len_class]-1);
	rand.transposeInPlace();
	// // tmp: write simulated files
	// for (int32_t i = 0; i < norm.rows(); i++) {
	// 	for (int32_t j = 0; j < norm.cols()-1; j++)
	// 		output3 << norm(i,j) << " ";
	// 	output3 << norm(i, norm.cols()-1) << endl;
	// }
	// sample multinomial distribution
	vector<double> DelScore_pos(5*num_sampling, 0);
	vector<double> DelScore_neg(5*num_sampling, 0);
	pair<int32_t,int32_t> region;
	Eigen::MatrixXd sampled_obs = Eigen::MatrixXd::Zero(5*num_sampling, nBins[len_class]);
	for (int32_t i = 0; i < num_sampling; i++){
		Eigen::VectorXd tmp_theo(nBins[len_class]);
		tmp_theo.segment(0, nBins[len_class]-1) = rand.row(i);
		tmp_theo(nBins[len_class]-1) = 1 - rand.row(i).sum();
		for (int32_t j = 0; j < nBins[len_class]; j++)
			tmp_theo(j) = max(0.0, tmp_theo(j));
		tmp_theo /= tmp_theo.sum();
		for (int32_t j = 0; j < 5; j++) {
			uint32_t* tmp_obs = new uint32_t[nBins[len_class]];
			gsl_ran_multinomial(r, nBins[len_class], numreads, tmp_theo.data(), tmp_obs);
			for (int32_t k = 0; k < nBins[len_class]; k++)
				sampled_obs(i*5+j, k) = tmp_obs[k];
			assert(sampled_obs.row(i*5+j).sum() == numreads);
			// // tmp: write simulated files
			// for (int32_t k = 0; k < sampled_obs.cols()-1; k++)
			// 	output4 << sampled_obs(i*5+j, k) << " ";
			// output4 << sampled_obs(i*5+j, sampled_obs.cols()-1) << endl;
			DelScore_pos[i*5+j] = CalDeletionScore_single(ExpectedBinNorm[ind], sampled_obs.row(i*5+j), region, false);
			DelScore_neg[i*5+j] = CalDeletionScore_single(sampled_obs.row(i*5+j), ExpectedBinNorm[ind], region, true);
			delete [] tmp_obs;
		}
	}
	// // tmp: write simulated files
	// output1.close();
	// output2.close();
	// output3.close();
	// output4.close();
	// count number of samples that have a more extreme deletion score
	int32_t count_pos = 0;
	int32_t count_neg = 0;
	for (int32_t i = 0; i < DelScore_pos.size(); i++) {
		if (DelScore_pos[i] >= diff_pos)
			count_pos ++;
	}
	for (int32_t i = 0; i < DelScore_neg.size(); i++) {
		if (DelScore_neg[i] >= diff_neg)
			count_neg ++;
	}
	pvalue_pos = 1.0*count_pos / DelScore_pos.size();
	pvalue_neg = 1.0*count_neg / DelScore_neg.size();
	// return a pair of pvalues for positive deletion score and negative deletion score
	return make_pair(pvalue_pos, pvalue_neg);
};


void DistTest_t::PValue_regional(vector<PRegion_t>& PValuesPos, vector<PRegion_t>& PValuesNeg)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating regional P-value."<<endl;

	// clear variables
	PValuesPos.clear();
	PValuesNeg.clear();
	// calculate p value
	mutex PValues_mutex;
	omp_set_num_threads(Num_Threads);
	#pragma omp parallel for
	for (int32_t i = 0; i < TransNames.size(); i++){
		vector<PRegion_t> tmp_pos = SinglePvalue_regional_pos(i);
		{
			lock_guard<std::mutex> guard(PValues_mutex);
			PValuesPos.insert(PValuesPos.end(), tmp_pos.begin(), tmp_pos.end());
		}
		vector<PRegion_t> tmp_neg = SinglePvalue_regional_neg(i);
		{
			lock_guard<std::mutex> guard(PValues_mutex);
			PValuesNeg.insert(PValuesNeg.end(), tmp_neg.begin(), tmp_neg.end());
		}
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish regional P-value calculation. There are ";
	cout << (PValuesPos.size()) << " positive bins and " << (PValuesNeg.size()) << " negative bins." << "\n";
};


void DistTest_t::PValue_regional(const vector<int32_t>& AdjustmentList, vector<PRegion_t>& PValuesPos, vector<PRegion_t>& PValuesNeg)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating regional P-value."<<endl;

	// clear variables
	PValuesPos.clear();
	PValuesNeg.clear();
	// calculate p value
	mutex PValues_mutex;
	omp_set_num_threads(Num_Threads);
	#pragma omp parallel for
	for (int32_t tmpi = 0; tmpi < AdjustmentList.size(); tmpi++){
		int32_t i = AdjustmentList[tmpi];
		vector<PRegion_t> tmp_pos = SinglePvalue_regional_pos(i);
		{
			lock_guard<std::mutex> guard(PValues_mutex);
			PValuesPos.insert(PValuesPos.end(), tmp_pos.begin(), tmp_pos.end());
		}
		vector<PRegion_t> tmp_neg = SinglePvalue_regional_neg(i);
		{
			lock_guard<std::mutex> guard(PValues_mutex);
			PValuesNeg.insert(PValuesNeg.end(), tmp_neg.begin(), tmp_neg.end());
		}
	}

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish regional P-value calculation. There are ";
	cout << (PValuesPos.size()) << " positive bins and " << (PValuesNeg.size()) << " negative bins." << "\n";
};


void DistTest_t::PValue_overall_empirical(vector<double>& PValuesPos, vector<double>& PValuesNeg, vector<bool>& Choices)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating empirical P-value."<<endl;

	// clear variables
	PValuesPos.resize(TransNames.size(), 1);
	PValuesNeg.resize(TransNames.size(), 1);
	Choices.resize(TransNames.size(), 0);
	// sampling from empirical p value
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	omp_set_num_threads(Num_Threads);
	#pragma omp parallel for
	for (int32_t i = 0; i < TransNames.size(); i++){
		pair<double,double> pvs = SinglePvalue_overall_empirical(i, 100, r);
		PValuesPos[i] = pvs.first;
		PValuesNeg[i] = pvs.second;
		Choices[i] = (pvs.first > pvs.second);
		if (i%5000 == 0){
			time(&CurrentTime);
			CurrentTimeStr=ctime(&CurrentTime);
			cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << i <<endl;
		}
	}
	// count small p values due to lack of sampling
	int32_t count = 0;
	for (int32_t i = 0; i < TransNames.size(); i++) {
		if (min(PValuesPos[i], PValuesNeg[i]) < 0.01 && min(PValuesPos[i], PValuesNeg[i]) >= 0)
			count++;
	}
	cout << "With sampling size = 100, there are " << count << " transcripts with <0.01 pvalue.\n";
	// second round sampling with enlarged sampling size
	omp_set_num_threads(Num_Threads);
	#pragma omp parallel for
	for (int32_t i = 0; i < TransNames.size(); i++){
		if (min(PValuesPos[i], PValuesNeg[i]) > 0.01 || min(PValuesPos[i], PValuesNeg[i]) < -0.5)
			continue;
		// calculate approximate p value
		pair<double,double> pvs1 = SinglePvalue_overall(i);
		pair<double,double> pvs2 = make_pair(pvs1.first, pvs1.second);
		if (min(pvs1.first, pvs2.second) < 1e-4)
			pvs2 = SinglePvalue_overall_empirical(i, 5000, r);
		// comparing approximation and larger sampling, select the one with a larger p value
		// since approximation is a lower bound, and not enough sampling gives a 0 p value. Selecting the larger one makes sense.
		double pvmin1 = min(pvs1.first, pvs1.second);
		double pvmin2 = min(pvs2.first, pvs2.second);
		bool choice1 = (pvs1.first > pvs1.second);
		bool choice2 = (pvs2.first > pvs2.second);
		if (max(pvmin1, pvmin2) > min(PValuesPos[i], PValuesNeg[i])) {
			if (pvmin1 < pvmin2) {
				PValuesPos[i] = pvs2.first;
				PValuesNeg[i] = pvs2.second;
				Choices[i] = choice2;
			}
			else {
				PValuesPos[i] = pvs1.first;
				PValuesNeg[i] = pvs1.second;
				Choices[i] = choice1;
			}
		}
	}
	// count small p values due to lack of sampling
	count = 0;
	for (int32_t i = 0; i < TransNames.size(); i++) {
		if (min(PValuesPos[i], PValuesNeg[i]) < 0.01 && min(PValuesPos[i], PValuesNeg[i]) >= 0)
			count++;
	}
	cout << "With sampling size = 5000, there are " << count << " transcripts with <0.01 pvalue.\n";


	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish empirical P-value calculation.\n";
};


void DistTest_t::PValue_overall_empirical(const vector<int32_t>& AdjustmentList, vector<double>& PValuesPos, vector<double>& PValuesNeg, vector<bool>& Choices)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Calculating empirical P-value."<<endl;

	// clear variables
	PValuesPos.resize(AdjustmentList.size(), 1);
	PValuesNeg.resize(AdjustmentList.size(), 1);
	Choices.resize(AdjustmentList.size(), 0);
	// sampling from empirical p value
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	omp_set_num_threads(Num_Threads);
	#pragma omp parallel for
	for (int32_t i = 0; i < AdjustmentList.size(); i++){
		pair<double,double> pvs = SinglePvalue_overall_empirical(AdjustmentList[i], 100, r);
		PValuesPos[i] = pvs.first;
		PValuesNeg[i] = pvs.second;
		Choices[i] = (pvs.first > pvs.second);
		if (i%5000 == 0){
			time(&CurrentTime);
			CurrentTimeStr=ctime(&CurrentTime);
			cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << i <<endl;
		}
	}
	// count small p values due to lack of sampling
	int32_t count = 0;
	for (const int32_t& i : AdjustmentList) {
		if (min(PValuesPos[i], PValuesNeg[i]) < 0.01 && min(PValuesPos[i], PValuesNeg[i]) >= 0)
			count++;
	}
	cout << "With sampling size = 100, there are " << count << " transcripts with <0.01 pvalue.\n";
	// second round sampling with enlarged sampling size
	omp_set_num_threads(Num_Threads);
	#pragma omp parallel for
	for (int32_t i = 0; i < AdjustmentList.size(); i++){
		if (min(PValuesPos[i], PValuesNeg[i]) > 0.01 || min(PValuesPos[i], PValuesNeg[i]) < -0.5)
			continue;
		// calculate approximate p value
		pair<double,double> pvs1 = SinglePvalue_overall(AdjustmentList[i]);
		pair<double,double> pvs2 = make_pair(pvs1.first, pvs1.second);
		if (min(pvs1.first, pvs2.second) < 1e-4)
			pvs2 = SinglePvalue_overall_empirical(AdjustmentList[i], 5000, r);
		// comparing approximation and larger sampling, select the one with a larger p value
		// since approximation is a lower bound, and not enough sampling gives a 0 p value. Selecting the larger one makes sense.
		double pvmin1 = min(pvs1.first, pvs1.second);
		double pvmin2 = min(pvs2.first, pvs2.second);
		bool choice1 = (pvs1.first > pvs1.second);
		bool choice2 = (pvs2.first > pvs2.second);
		if (max(pvmin1, pvmin2) > min(PValuesPos[i], PValuesNeg[i])) {
			if (pvmin1 < pvmin2) {
				PValuesPos[i] = pvs2.first;
				PValuesNeg[i] = pvs2.second;
				Choices[i] = choice2;
			}
			else {
				PValuesPos[i] = pvs1.first;
				PValuesNeg[i] = pvs1.second;
				Choices[i] = choice1;
			}
		}
	}
	// count small p values due to lack of sampling
	count = 0;
	for (int32_t i = 0; i < AdjustmentList.size(); i++) {
		if (min(PValuesPos[i], PValuesNeg[i]) > 0.01 || min(PValuesPos[i], PValuesNeg[i]) < -0.5)
			count++;
	}
	cout << "With sampling size = 5000, there are " << count << " transcripts with <0.01 pvalue.\n";


	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish empirical P-value calculation.\n";
};


void BHAdjusting(const vector<PRegion_t>& RawPRegions, vector<PRegion_t>& AdjustedPRegions)
{
	// clear variable
	AdjustedPRegions.clear();
	for (int32_t i = 0; i < RawPRegions.size(); i++) {
		PRegion_t tmp;
		AdjustedPRegions.push_back(tmp);
	}
	// sort raw pvalues from the largest to the smallest
	vector<int32_t> Indexes(RawPRegions.size(), 0);
	std::iota(Indexes.begin(), Indexes.end(), 0);
	vector<PRegion_t> RawPRegions_sort(RawPRegions.cbegin(), RawPRegions.cend());
	sort(Indexes.begin(), Indexes.end(), [&RawPRegions](int32_t a, int32_t b){return RawPRegions[a].Pvalue > RawPRegions[b].Pvalue;} );
	for (int32_t i = 0; i < RawPRegions.size(); i++)
		RawPRegions_sort[i] = RawPRegions[Indexes[i]];
	// BH adjustment
	double cumulative_min = 1;
	int32_t n = RawPRegions_sort.size();
	for(int32_t i = 0; i < RawPRegions_sort.size(); i++){
		double adjp = min(cumulative_min, n / (n-i) * RawPRegions_sort[i].Pvalue);
		RawPRegions_sort[i].Pvalue = adjp;
		cumulative_min = min(cumulative_min, adjp);
	}
	// add to result
	for(int32_t i = 0; i < RawPRegions_sort.size(); i++){
		AdjustedPRegions[Indexes[i]] = RawPRegions_sort[i];
	}
};


void BHAdjusting(const vector<double>& RawPvalues, vector<double>& AdjustedPvalues)
{
	// clear variable
	AdjustedPvalues.resize(RawPvalues.size(), 0);
	// sort raw pvalues in decreasing order
	vector<int32_t> Indexes(RawPvalues.size(), 0);
	std::iota(Indexes.begin(), Indexes.end(), 0);
	vector<double> RawPvalues_sort(RawPvalues.cbegin(), RawPvalues.cend());
	sort(Indexes.begin(), Indexes.end(), [&RawPvalues_sort](int32_t a, int32_t b){return RawPvalues_sort[a] > RawPvalues_sort[b];} );
	for (int32_t i = 0; i < RawPvalues.size(); i++)
		RawPvalues_sort[i] = RawPvalues[Indexes[i]];
	// BH adjustment
	for(int32_t i = 0; i < RawPvalues_sort.size(); i++)
		RawPvalues_sort[i] = RawPvalues_sort.size() / ((int32_t)RawPvalues_sort.size() - i) * RawPvalues_sort[i];
	// cumulative min
	for (int32_t i = 1; i < RawPvalues_sort.size(); i++)
		RawPvalues_sort[i] = min(RawPvalues_sort[i-1], RawPvalues_sort[i]);
	// add to result
	for(int32_t i = 0; i < RawPvalues_sort.size(); i++)
		AdjustedPvalues[Indexes[i]] = RawPvalues_sort[i];
};


// XXX 
vector<int32_t> ReduceAssignmentList(const vector<int32_t>& AssignmentList, const map< string,vector<string> >& GeneTransMap, 
	const map<string,int32_t>& TransIndex, const vector<double>& SalmonExp, const vector<double>& LPExp, 
	const map< string,vector<MovingRead_t> >& MovingMat, const vector<PRegion_t>& AdjPValuesPos_salmon, const vector<PRegion_t>& AdjPValuesNeg_salmon, 
	const vector<PRegion_t>& AdjPValuesPos_lp, const vector<PRegion_t>& AdjPValuesNeg_lp, double PvalueThresh, double MoveReadThresh)
{
	// require old AssignmentList is sorted
	assert(is_sorted(AssignmentList.begin(), AssignmentList.end()));
	vector<int32_t> newAssignmentList;
	// collect significant salmon and LP
	vector<int32_t> Significant_salmon;
	vector<int32_t> Significant_lp;
	for (int32_t i = 0; i < AdjPValuesPos_salmon.size(); i++)
		if (AdjPValuesPos_salmon[i].Pvalue < PvalueThresh)
			Significant_salmon.push_back(AdjPValuesPos_salmon[i].TID);
	for (int32_t i = 0; i < AdjPValuesNeg_salmon.size(); i++)
		if (AdjPValuesNeg_salmon[i].Pvalue < PvalueThresh)
			Significant_salmon.push_back(AdjPValuesNeg_salmon[i].TID);
	for (int32_t i = 0; i < AdjPValuesPos_lp.size(); i++)
		if (AdjPValuesPos_lp[i].Pvalue < PvalueThresh)
			Significant_lp.push_back(AdjPValuesPos_lp[i].TID);
	for (int32_t i = 0; i < AdjPValuesNeg_lp.size(); i++)
		if (AdjPValuesNeg_lp[i].Pvalue < PvalueThresh)
			Significant_lp.push_back(AdjPValuesNeg_lp[i].TID);
	sort(Significant_salmon.begin(), Significant_salmon.end());
	sort(Significant_lp.begin(), Significant_lp.end());
	Significant_salmon.resize( distance(Significant_salmon.begin(), unique(Significant_salmon.begin(),Significant_salmon.end())) );
	Significant_lp.resize( distance(Significant_lp.begin(), unique(Significant_lp.begin(),Significant_lp.end())) );
	// filter transcript for each gene
	for (map< string,vector<string> >::const_iterator it = GeneTransMap.cbegin(); it != GeneTransMap.cend(); it++) {
		vector<int32_t> tids;
		for (const string& tname : it->second){
			map<string,int32_t>::const_iterator itidx = TransIndex.find(tname);
			assert(itidx != TransIndex.cend());
			if (binary_search(AssignmentList.begin(), AssignmentList.end(), itidx->second))
				tids.push_back(itidx->second);
		}
		if (tids.size() == 0)
			continue;
		map< string,vector<MovingRead_t> >::const_iterator itmoving = MovingMat.find(it->first);
		assert(itmoving != MovingMat.cend());
		const vector<MovingRead_t>& g_moving = itmoving->second;
		vector<int32_t> trans_better;
		vector<int32_t> trans_worse;
		vector<int32_t> trans_helper;
		for (const int32_t& t : tids){
			bool is_old_sig = binary_search(Significant_salmon.begin(), Significant_salmon.end(), t);
			bool is_new_sig = binary_search(Significant_lp.begin(), Significant_lp.end(), t);
			if (is_old_sig && ! is_new_sig)
				trans_better.push_back(t);
			else if(!is_old_sig && is_new_sig)
				trans_worse.push_back(t);
		}
		sort(trans_better.begin(), trans_better.end());
		trans_better.resize( distance(trans_better.begin(), unique(trans_better.begin(),trans_better.end())) );
		sort(trans_worse.begin(), trans_worse.end());
		trans_worse.resize( distance(trans_worse.begin(), unique(trans_worse.begin(),trans_worse.end())) );
		for (const int32_t& t1 : trans_better) {
			for (const int32_t& t2 : tids) {
				if (t2 != t1 && !binary_search(trans_worse.begin(), trans_worse.end(), t2)) {
					MovingRead_t tmp1 (t1, t2, 0);
					MovingRead_t tmp2 (t2, t1, 0);
					vector<MovingRead_t>::const_iterator lb1 = lower_bound(g_moving.cbegin(), g_moving.cend(), tmp1);
					vector<MovingRead_t>::const_iterator lb2 = lower_bound(g_moving.cbegin(), g_moving.cend(), tmp2);
					assert(lb1 != g_moving.cend() && lb1->TID1 == t1 && lb1->TID2 == t2);
					assert(lb2 != g_moving.cend() && lb2->TID1 == t2 && lb2->TID2 == t1);
					if ((lb1->NReads)*(lb2->NReads) < 0 && fabs(lb2->NReads) >= MoveReadThresh*fabs(lb1->NReads))
						trans_helper.push_back(t2);
				}
			}
		}
		sort(trans_helper.begin(), trans_helper.end());
		trans_helper.resize( distance(trans_helper.begin(), unique(trans_helper.begin(), trans_helper.end())) );
		// new Adjustment list defined as trans_better | trans_helper - trans_worse
		vector<int32_t> tmpAdjustmentList;
		set_difference(trans_better.begin(), trans_better.end(), trans_worse.begin(), trans_worse.end(), 
			inserter(tmpAdjustmentList, tmpAdjustmentList.begin()));
		set_difference(trans_helper.begin(), trans_helper.end(), trans_worse.begin(), trans_worse.end(), 
			inserter(tmpAdjustmentList, tmpAdjustmentList.end()));
		sort(tmpAdjustmentList.begin(), tmpAdjustmentList.end());
		tmpAdjustmentList.resize( distance(tmpAdjustmentList.begin(), unique(tmpAdjustmentList.begin(),tmpAdjustmentList.end())) );
		newAssignmentList.insert(newAssignmentList.end(), tmpAdjustmentList.begin(), tmpAdjustmentList.end());
	}
	// sort new AssignmentList
	sort(newAssignmentList.begin(), newAssignmentList.end());
	return newAssignmentList;
};


