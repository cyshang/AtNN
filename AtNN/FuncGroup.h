#ifndef FUNCGROUP_H_
#define FUNCGROUP_H_

#include "global.h"

#ifdef DEBUG_MODE
//#define DEBUG_FUNCGROUP
#endif // DEBUG_MODE

struct FuncInfo;

struct FuncGroup {
	int nFunc;
	int element;
	vector<FuncInfo> funcinfo;

	FuncGroup(const int & _element);
	~FuncGroup();
	void Construct(std::istream & fin, const int & _nFunc);	
	void OutputInfo(std::ostream & fout);
};

#endif // !FUNCTIONINFO_H