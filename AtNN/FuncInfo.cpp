#include "FuncInfo.h"

using std::left;
using std::setw;
using std::scientific;
using std::fixed;

FuncInfo::FuncInfo(const std::string & str)
{
	GetFuncInfo(str);
}

void FuncInfo::GetFuncInfo(const std::string & str)
{
	istringstream sin(str);
	string elementStr, cutoffStr, symStr;

	sin >> cutoffStr >> symStr;
	// Init cutoff_func
	if (cutoffStr == "C1") {
		cutoff = 0;
	}
	else if (cutoffStr == "C2") {
		cutoff = 1;
	}
	//--------------------------------

	// Init sym_func
	if (symStr == "G1") {
		symfunc = 0;
	}
	else if (symStr == "G2") {
		symfunc = 1;
	}
	else if (symStr == "G3") {
		symfunc = 2;
	}
	else if (symStr == "G4") {
		symfunc = 3;
	}
	//---------------------------------

	switch (symfunc) {
	case 0: {

		sin >> elementStr;
		elements[0] = parameter.element_to_num[elementStr];		

		sin >> funcParam[0];
		break;
	}
	case 1: {

		sin >> elementStr;
		elements[0] = parameter.element_to_num[elementStr];

		sin >> funcParam[0] >> funcParam[1] >> funcParam[2];
		break;
	}
	case 2: {

		sin >> elementStr;
		elements[0] = parameter.element_to_num[elementStr];
		sin >> elementStr;
		elements[1] = parameter.element_to_num[elementStr];
		
		sin >> funcParam[0] >> funcParam[1] >> funcParam[2] >> funcParam[3];
		break;
	}
	case 3: {

		sin >> elementStr;
		elements[0] = parameter.element_to_num[elementStr];
		sin >> elementStr;
		elements[1] = parameter.element_to_num[elementStr];

		sin >> funcParam[0] >> funcParam[1] >> funcParam[2] >> funcParam[3];
		break;
	}
	}
}

void FuncInfo::Backup()
{
	for (int i = 0; i < FUNC_PARAMETER; ++i) {
		copyParam[i] = funcParam[i];
	}
}

void FuncInfo::Restore()
{
	for (int i = 0; i < FUNC_PARAMETER; ++i) {
		funcParam[i] = copyParam[i];
	}
}

string FuncInfo::PrintFuncInfo() const
{	
	ostringstream sout;
	int element_size, parameter_size;

	switch (cutoff)
	{
	case 0: sout << setw(4) << left << "C1"; break;
	case 1: sout << setw(4) << left << "C2"; break;
	}

	switch (symfunc)
	{
	case 0:	sout << setw(4) << "G1"; element_size = 1; parameter_size = 1; break;
	case 1:	sout << setw(4) << "G2"; element_size = 1; parameter_size = 3; break;
	case 2:	sout << setw(4) << "G3"; element_size = 2; parameter_size = 4; break;
	case 3:	sout << setw(4) << "G4"; element_size = 2; parameter_size = 4; break;
	}

	for (int i = 0; i < 2; ++i) {
		if (i < element_size)
			sout << setw(4) << left << parameter.num_to_element[elements[i]];
		else
			sout << setw(4) << left << ' ';
	}

	sout << setw(6) << left << funcParam[0];

	if (parameter_size == 4) {
		sout << std::showpos;
		sout << setw(25) << left << funcParam[1];
		sout << std::noshowpos;
		sout << std::setprecision(15) << scientific;
		sout << setw(25) << left << funcParam[2];
		sout << setw(25) << left << funcParam[3];
	}
	else {
		for (int i = 1; i < parameter_size; ++i) {
			sout << std::setprecision(15) << scientific;
			sout << setw(25) << left << funcParam[i];
		}
	}

	return sout.str();
}


std::string FuncInfo::OutputPES()
{
	ostringstream sout;

	sout << setw(4) << left << symfunc;
	
	if (symfunc <= 1) {
		sout << setw(4) << left << elements[0];
		sout << setw(4) << left << 0;
	}
	else {
		sout << setw(4) << left << elements[0];
		sout << setw(4) << left << elements[1];
	}

	sout << std::setprecision(16) << std::scientific;
	if (symfunc == 1) {
		for (int i = 1; i < 3; ++i) {
			sout << setw(25) << left << funcParam[i];
		}
		sout << setw(25) << left << 0;
	}
	else {
		for (int i = 1; i < 4; ++i) {
			sout << setw(25) << left << funcParam[i];
		}
	}

	return sout.str();
}