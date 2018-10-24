#ifndef FUNCTYPE_H_
#define FUNCTYPE_H_

#include "global.h"
#ifdef DEBUG_MODE
//#define DEBUG_FUNCINFO
#endif // DEBUG_MODE



#define FUNC_PARAMETER 4

struct FuncInfo {

	double funcParam[FUNC_PARAMETER];
	double copyParam[FUNC_PARAMETER];
	//	G1 : {Rc}
	//	G2 : {Rc, Rs, eta}
	//	G3 : {Rc, lambda, eta, xi}
	//	G4 : {Rc, lambda, eta, xi}
	int cutoff;
	int symfunc;
	int elements[2]; 

	FuncInfo(const std::string & str);
	void GetFuncInfo(const std::string & str);
	void Backup();
	void Restore();
	std::string PrintFuncInfo() const;
	std::string OutputPES();
};

#endif // !FUNCTYPE_H

