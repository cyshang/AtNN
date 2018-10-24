#include "FuncGroup.h"
#include "FuncInfo.h"

FuncGroup::FuncGroup(const int & _element) :element(_element), nFunc(0) {}

FuncGroup::~FuncGroup() {}

void FuncGroup::Construct(istream & fin, const int & _nFunc)
{
	nFunc = _nFunc;
	
	string line;
	funcinfo.clear();
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		std::getline(fin, line);
		funcinfo.push_back(line);
	}

#ifdef DEBUG_FUNCGROUP

	debug << parameter.num_to_element[element] << " " << nFunc << endl;
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		debug << funcinfo[iFunc].PrintFuncInfo() << endl;
	}
	debug << endl;

#endif // DEBUG_FUNCGROUP

}

void FuncGroup::OutputInfo(ostream & fout)
{
	fout << parameter.num_to_element[element] << " " << nFunc << endl;
	
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		fout << funcinfo[iFunc].PrintFuncInfo() << endl;
	}
}
