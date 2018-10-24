#ifndef NEURALNETWORK_H_
#define NEURALNETWORK_H_

#include "global.h"
#include <Eigen/Dense>
#include "NetworkInfo.h"
#include "Group.h"

#ifdef DEBUG_MODE
//#define DEBUG_NEURALNETWORK
#endif

struct SymFunction;

struct NeuralNetwork
{
	SymFunction *pSymfunc;
	NetworkInfo networkinfo;

	int dimX;
	Eigen::Map<Eigen::MatrixXd> rawX;
	Eigen::MatrixXd inputX;
	Eigen::VectorXd maxX;
	Eigen::VectorXd minX;
	Eigen::VectorXd avgX;

	Eigen::Map<Eigen::RowVectorXd> rawEnergy;
	Eigen::RowVectorXd targetEnergy;
	double maxEnergy;
	double minEnergy;
	double avgEnergy;

	std::vector<int> random_sequence;

	std::vector<Group> netGroup;

	Eigen::RowVectorXd tEnergy;
	Eigen::RowVectorXd vEnergy;
	Eigen::RowVectorXd tErr;
	Eigen::RowVectorXd vErr;
	int inc_step;

	int nWeight;
	Eigen::MatrixXd Jac;
	Eigen::MatrixXd JtJ;
	Eigen::VectorXd JtErr;
	Eigen::VectorXd dWeight;

	double tRMSE;	
	double vRMSE;
	double tRMSE_last;
	double vRMSE_last;

	NeuralNetwork(SymFunction *_pSymfunc);
	~NeuralNetwork();

//	void Init(SymFunction *_pSymfunc) {}
	void Construct();
	void TrainNetwork();

	void PreprocessData();
	void ShuffleData(); //	rearrange X and Energy according to random_list
	void FittingControl(const int & iFit);

	void ForwardProp();
	bool TrainPerf();
	void CalDevSet();
	bool DevPerf();

	void OutputDebug(std::ostream & out = std::cout);
	void SaveNetwork(const int & iFit);

	//========Intrinsic Function========
	
};
#endif