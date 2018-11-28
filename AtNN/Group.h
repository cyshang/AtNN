#ifndef GROUP_H_
#define GROUP_H_

#include "global.h"
#include <Eigen/Core>

#ifdef DEBUG_MODE
//#define DEBUG_GROUP
#endif

struct NeuralNetwork;
struct NetworkInfo;

struct Group
{
	static const NetworkInfo *pInfo;
	static NeuralNetwork *pNetwork;
	int iGroup; //id of this group
	int nNet;	//number of sub-nets in this group
	int nLayer;	//number of layers in this group
	std::vector<int> nNeuron;					//An array contrain dims information for each sub-net

	int nWeight;
	int WeightStart;
	int inputStart;

	std::vector<std::vector<Eigen::MatrixXd> > A;
	std::vector<std::vector<Eigen::MatrixXd> > devA;
	std::vector<std::vector<Eigen::MatrixXd> > dFdZ;
	std::vector<Eigen::MatrixXd> W;
	std::vector<Eigen::VectorXd> b;
	std::vector<Eigen::MatrixXd> W_copy;
	std::vector<Eigen::VectorXd> b_copy;

	Group(const int & _Group);
	~Group();
	void Construct();
	void DataInput();
	void ForwardProp();
	void BackProp();
	void UpdateWeight();
	void BackupWeight();
	void RestoreWeight();
	void CalDevEnergy();
	void RandWeight();
	void OutputWeight(std::ostream & outW);
	void SaveWeight(std::ostream & outW);
};

#endif // !GROUP_H

