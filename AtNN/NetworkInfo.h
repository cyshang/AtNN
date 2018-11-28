#ifndef NETWORKINFO_H_
#define NETWORKINFO_H_

#include "global.h"

#ifdef DEBUG_MODE
//#define DEBUG_NETWORKINFO
#endif

struct SymFunction;

struct NetworkInfo {

	std::vector<int> nHidden;
	std::vector<std::vector<int> > HiddenLayer;

	int nGroup;		//ith group in NN
	std::vector<int> nNet;	//number of subnets in this group
	std::vector<int> nLayer;	//number of layers of every subnet in this group
	std::vector<std::vector<int> > nNeuron;	//number of neurons in each layer

	double train_ratio;
	long tSample;
	long vSample;

	int nFitting;
	int maxEpoch;
	bool IfEarly;
	int EarlySteps;

	double mu;

	NetworkInfo();
	~NetworkInfo();

	void GetInfo();
	void Construct(const SymFunction *pSymfunc);
	void SaveInfo(std::ostream & fout);
};

#endif // !NETWORKINFO_H

