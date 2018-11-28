#ifndef SYMFUNCTION_H_
#define SYMFUNCTION_H_

#include "global.h"
#include "Molecule.h"

#ifdef DEBUG_MODE
//#define DEBUG_SYMFUNCTION
#endif // DEBUG_MODE

struct FuncGroup;

struct SymFunction {

	struct Array2;

	vector<FuncGroup> funcgroup;

	int dimX;
	double *X;
	double *Energy;

	vector<vector<Array2> > tlist;

	SymFunction();
	~SymFunction();

	void Construct(istream & fin);
	void Clear();	
	void CalX(const long & iSample, const Molecule & molecule);
	void PES_Funcinfo();
	void OutputFuncInfo(ostream & fout);
	
};

struct SymFunction::Array2
{
	int data[2];

	Array2(const int & a, const int & b)
	{
		data[0] = a;
		data[1] = b;
	}

	Array2(const Array2 & arr)
	{
		data[0] = arr.data[0];
		data[1] = arr.data[1];
	}

	Array2 & operator = (const Array2 & arr)
	{
		data[0] = arr.data[0];
		data[1] = arr.data[1];
		
		return *this;
	}

	int & operator [] (const int & n)
	{
		return data[n];
	}
};

#endif // !SYMFUNCTION_H_
