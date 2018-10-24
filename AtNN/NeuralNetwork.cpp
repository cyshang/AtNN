#include "NeuralNetwork.h"
#include "SymFunction.h"
#include "NetworkInfo.h"
#include "Group.h"
#include <algorithm>
#include <iomanip>

using std::left;
using std::setw;
using std::setprecision;
using namespace Eigen;

NeuralNetwork::NeuralNetwork(SymFunction *_pSymfunc) :
	networkinfo(), pSymfunc(NULL), rawX(NULL, 0, 0), rawEnergy(NULL, 0)
{
	pSymfunc = _pSymfunc;
	Group::pNetwork = this;
	Group::pInfo = &networkinfo;

	networkinfo.GetInfo();
	for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		netGroup.push_back(iGroup);
	}
}

NeuralNetwork::~NeuralNetwork() {}

void NeuralNetwork::Construct()
{
#ifdef DEBUG_NEURALNETWORK
	debug << LMARK << "NeuralNetwork::ConstructNetwork" << RMARK << endl;	
#endif

	networkinfo.Construct(pSymfunc);

	//---------- Calc dimX ----------	
	dimX = pSymfunc->dimX;

	//---------- Init rawX & inputX ----------
	new (&rawX) Map<MatrixXd>(pSymfunc->X, dimX, parameter.nSample);
	inputX.resize(dimX, parameter.nSample);

	maxX.resize(dimX);
	minX.resize(dimX);
	avgX.resize(dimX);

	//---------- Init rawEnergy & targetEnergy ----------
	new (&rawEnergy) Map<RowVectorXd>(pSymfunc->Energy, parameter.nSample);	
	targetEnergy.resize(parameter.nSample);

	tEnergy.resize(networkinfo.tSample);
	vEnergy.resize(networkinfo.vSample);
	tErr.resize(networkinfo.tSample);
	vErr.resize(networkinfo.vSample);

	// ---------- Construct Groups ----------
	nWeight = 0;
	for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		netGroup[iGroup].Construct();
		nWeight += netGroup[iGroup].nWeight;
	}

	Jac.resize(networkinfo.tSample, nWeight);
	JtJ.resize(nWeight, nWeight);
	JtErr.resize(nWeight);
	dWeight.resize(nWeight);

	//---------- Init RandomList ----------
	random_sequence.resize(parameter.nSample);
	for (long i = 0; i < parameter.nSample; ++i)
		random_sequence[i] = i;
		
#ifdef DEBUG_NEURALNETWORK

#endif
}

void NeuralNetwork::PreprocessData()
{

	maxX = rawX.rowwise().maxCoeff();
	minX = rawX.rowwise().minCoeff();
	avgX = rawX.rowwise().mean();

	maxEnergy = rawEnergy.maxCoeff();
	minEnergy = rawEnergy.minCoeff();
	avgEnergy = rawEnergy.mean();

	rawX = (rawX.colwise() - avgX).array().colwise() / (maxX - minX).array();
	rawEnergy = (rawEnergy.array() - avgEnergy) / (maxEnergy - minEnergy);
}

void NeuralNetwork::ShuffleData()
{
	std::random_shuffle(random_sequence.begin(), random_sequence.end());

	for (long i = 0; i < parameter.nSample; ++i) {
		inputX.col(i) = rawX.col(random_sequence[i]);
		targetEnergy.col(i) = rawEnergy.col(random_sequence[i]);
	}
}

void NeuralNetwork::TrainNetwork()
{
	string FileName;	

	FileName = parameter.output_folder +
		"report_" + parameter.projectName + "_F" + parameter.funcinfo_id + ".dat";
	ofstream fout(FileName.c_str(), ofstream::out);

	PreprocessData();
	for (int iFit = 1; iFit <= networkinfo.nFitting; ++iFit) {

		tRMSE_last = 9e9;
		vRMSE_last = 9e9;

		for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
			netGroup[iGroup].RandWeight();
		}
		ShuffleData();

		FittingControl(iFit);
		fout << left << setw(12) << tRMSE;
		fout << left << setw(12) << vRMSE << endl;
		SaveNetwork(iFit);
	}

	fout.close();	
	return;
}

void NeuralNetwork::FittingControl(const int & iFit)
{
	ofstream fout;

	char now_time[20];
	double runtime;
	time_t rawtime, start_time, end_time;

	time(&rawtime);
	strftime(now_time, 20, "%Y.%m.%d %X", localtime(&rawtime));

	string LogName;
	ostringstream toStr;
	
	if (networkinfo.nFitting < 10)
		toStr << iFit;
	else if (networkinfo.nFitting < 100)
		toStr << std::setw(2) << std::setfill('0') << std::right << iFit;
	else
		toStr << std::setw(3) << std::setfill('0') << std::right << iFit;

	LogName = parameter.output_folder +
		"log_" + parameter.projectName + "_F" + parameter.funcinfo_id + "." + toStr.str();
	fout.open(LogName.c_str(), ofstream::out);

	fout << now_time << endl << endl;
	fout << "Fit" << iFit << endl;
	fout << left << setw(8) << "Epoch";
	fout << left << setw(12) << "tRMSE(meV)";
	fout << left << setw(12) << "vRMSE(meV)";
	fout << left << setw(8) << "MU" << endl;

	bool IfBreak = false;
	double now_mu = networkinfo.mu;

	time(&start_time);

	for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		netGroup[iGroup].DataInput();
	}

	ForwardProp();
	TrainPerf();

	for (int iEpoch = 1; iEpoch <= networkinfo.maxEpoch; ++iEpoch) {

		for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
			netGroup[iGroup].BackupWeight();
		}
		
		Jac.setZero();
		for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
			netGroup[iGroup].BackProp();
		}
		JtJ = Jac.transpose().eval() * Jac;
		JtErr = Jac.transpose() * tErr.transpose();
		dWeight = (JtJ + now_mu * MatrixXd::Identity(nWeight, nWeight)).llt().solve(JtErr);
		for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
			netGroup[iGroup].UpdateWeight();
		}

		ForwardProp();
		while (TrainPerf()) {
			now_mu *= 2;
			if (now_mu > 1e10) {
				IfBreak = true;
				break;
			}
			for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
				netGroup[iGroup].RestoreWeight();
			}
			dWeight = (JtJ + now_mu * MatrixXd::Identity(nWeight, nWeight)).llt().solve(JtErr);
			for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
				netGroup[iGroup].UpdateWeight();
			}
			ForwardProp();
		}
		if (IfBreak)
			break;

		now_mu /= 2;

		CalDevSet();
		DevPerf();

		fout << left << setw(8) << iEpoch;
		fout << left << setw(12) << tRMSE;
		fout << left << setw(12) << vRMSE;
		fout << left << setw(8) << now_mu << endl;

#ifdef OUTPUT_TO_SCREEN
		cout << left << setw(8) << iEpoch;
		cout << left << setw(12) << tRMSE;
		cout << left << setw(12) << vRMSE;
		cout << left << setw(8) << now_mu << endl;
#endif

		if (networkinfo.IfEarly) {
			if (inc_step > networkinfo.EarlySteps)
				break;
		}
	}

	time(&end_time);
	runtime = difftime(end_time, start_time);
	time(&rawtime);
	strftime(now_time, 20, "%Y.%m.%d %X", localtime(&rawtime));
	fout << endl << now_time << endl;
	if (runtime < 60)
		fout << "Runtime: " << runtime << "s" << endl;
	else if (runtime >= 60 && runtime < 3600)
		fout << "Runtime: " << static_cast<int>(runtime / 60) << "min" << runtime - 60 * static_cast<int>(runtime / 60) << "s" << endl;
	else {
		fout << "Runtime: " << static_cast<int>(runtime / 3600) << "hour, ";
		fout << static_cast<int>((runtime - static_cast<int>(runtime / 3600) * 3600) / 60) << "min, ";
		fout << runtime - 60 * static_cast<int>(runtime / 60) << "s" << endl;
	}

	fout.close();
}

void NeuralNetwork::ForwardProp()
{
	tEnergy.setZero();
	for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		netGroup[iGroup].ForwardProp();
	}
#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::ForwardProp" << RMARK << endl;
	debug << "trainEnergy" << endl << trainEnergy << endl;
#endif // DEBUG_MODE
}



bool NeuralNetwork::TrainPerf()
{
	bool IfBad;

	tErr = targetEnergy.segment(0, networkinfo.tSample) - tEnergy;
	tRMSE = sqrt(tErr.squaredNorm() / networkinfo.tSample) * (maxEnergy - minEnergy) * parameter.energy_correction;

	if(tRMSE > tRMSE_last)
		IfBad = true;
	else {
		IfBad = false;
		tRMSE_last = tRMSE;
	}

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::TrainPerf" << RMARK << endl;
	debug << "tErr" << endl << tErr << endl;
#endif // DEBUG_MODE

	return IfBad;
}

void NeuralNetwork::CalDevSet()
{
	vEnergy.setZero();
	for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		netGroup[iGroup].CalDevEnergy();
	}
}

bool NeuralNetwork::DevPerf()
{
	bool IfBad = false;

	vErr = targetEnergy.segment(networkinfo.tSample, networkinfo.vSample) - vEnergy;	
	vRMSE = sqrt(vErr.squaredNorm() / networkinfo.vSample) * (maxEnergy - minEnergy) * parameter.energy_correction;

	if (vRMSE > vRMSE_last) {
		IfBad = true;
		inc_step++;
	}
	else {
		inc_step = 0;
	}
    vRMSE_last = vRMSE;

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::DevPerf" << RMARK << endl;
	debug << "vErr" << endl << vErr << endl;
#endif // DEBUG_MODE

	return IfBad;
}

void NeuralNetwork::SaveNetwork(const int & iFit)
{
	string OutputName;
	ofstream fout;
	ostringstream toStr;

	if (networkinfo.nFitting < 10)
		toStr << iFit;
	else if (networkinfo.nFitting < 100)
		toStr << std::setw(2) << std::setfill('0') << std::right << iFit;
	else
		toStr << std::setw(3) << std::setfill('0') << std::right << iFit;

	OutputName = parameter.output_folder +
		"net_" + parameter.projectName + "_F" + parameter.funcinfo_id + "." + toStr.str();
	fout.open(OutputName.c_str(), ofstream::out);

	// ========== funcinfo ==========
	fout << "<funcinfo>" << endl << endl;
	pSymfunc->OutputFuncInfo(fout);
	fout << "<end>" << endl << endl;

	// ========== networkinfo =========
	fout << "<networkinfo>" << endl << endl;
	
	networkinfo.SaveInfo(fout);

	fout << "<end>" << endl << endl;

	//--------------- Weight ---------------
	fout << "<network>" << endl << endl;

	fout << "weight" << endl;
	for (int iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		netGroup[iGroup].SaveWeight(fout);
	}
	fout << endl;

	//--------------- minX ---------------
	fout << "minX" << endl;
	for (int i = 0; i < dimX; ++i) {
		fout << setprecision(16) << minX[i] << " ";
	}
	fout << endl;

	//--------------- avgX ---------------
	fout << "avgX" << endl;
	for (int i = 0; i < dimX; ++i) {
		fout << setprecision(16) << avgX[i] << " ";
	}
	fout << endl;

	//--------------- maxX ---------------
	fout << "maxX" << endl;
	for (int i = 0; i < dimX; ++i) {
		fout << setprecision(16) << maxX[i] << " ";
	}
	fout << endl;

	//--------------- minEnergy ---------------
	fout << "minEnergy" << endl;
	fout << minEnergy << endl;

	//--------------- avgEnergy ---------------
	fout << "avgEnergy" << endl;
	fout << avgEnergy << endl;

	//--------------- maxEnergy ---------------
	fout << "maxEnergy" << endl;
	fout << maxEnergy << endl;

	fout << "<end>" << endl;

	fout.close();
}

void NeuralNetwork::OutputDebug(std::ostream & out)
{
	out << LMARK << "NeuralNetwork" << RMARK << endl;
	out << "<Matrix Size>" << endl;
	out << "rawX: " << rawX.rows() << "x" << rawX.cols() << endl;
	out << "inputX: " << inputX.rows() << "x" << inputX.cols() << endl;
	out << "maxX: " << maxX.size() << endl;
	out << "minX: " << minX.size() << endl;
	out << "rawEnergy: " << rawEnergy.size() << endl;
	out << "targetEnergy: " << targetEnergy.size() << endl;
	out << endl;
	out << "<Matrix Value>" << endl;
	out << LMARK << "rawX" << RMARK << endl << rawX << endl << endl;
	out << LMARK << "inputX" << RMARK << endl << inputX << endl << endl;
	out << LMARK << "maxX" << RMARK << endl << maxX << endl << endl;
	out << LMARK << "minX" << RMARK << endl << minX << endl << endl;
	out << LMARK << "rawEnergy" << RMARK << endl << rawEnergy << endl << endl;
	out << LMARK << "targetEnergy" << RMARK << endl << targetEnergy << endl << endl;
	out << "maxEnergy: " << maxEnergy << endl;
	out << "minEnergy: " << minEnergy << endl;
	out << endl;
}
