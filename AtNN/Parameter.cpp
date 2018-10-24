#include "Parameter.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using std::string;
using std::istringstream;
using std::getline;
using std::cout;
using std::endl;


Parameter::Parameter() :
	random_seed(-1), 
	nSample(-1), energy_correction(-1),
	nElement(-1), nAtom(-1),
	nFuncInfo(-1)
{
	function_table.insert(Str_to_Func::value_type("random_seed", &Parameter::get_random_seed));
	//
	function_table.insert(Str_to_Func::value_type("nSample", &Parameter::get_nSample));
	function_table.insert(Str_to_Func::value_type("energy_correction", &Parameter::get_energy_correction));
	function_table.insert(Str_to_Func::value_type("nElement", &Parameter::get_nElement));
	function_table.insert(Str_to_Func::value_type("nAtom", &Parameter::get_nAtom));
	function_table.insert(Str_to_Func::value_type("element_list", &Parameter::get_element_list));
	function_table.insert(Str_to_Func::value_type("atom_list", &Parameter::get_atom_list));
	function_table.insert(Str_to_Func::value_type("run_mode", &Parameter::get_run_mode));
	//
	function_table.insert(Str_to_Func::value_type("nFuncInfo", &Parameter::get_nFuncInfo));
	function_table.insert(Str_to_Func::value_type("projectName", &Parameter::getProjectName));
	function_table.insert(Str_to_Func::value_type("fNetworkInfo", &Parameter::getNetworkInfo));
	function_table.insert(Str_to_Func::value_type("DataType", &Parameter::getDataType));
	function_table.insert(Str_to_Func::value_type("fCartesianData", &Parameter::getCartesianData));
	function_table.insert(Str_to_Func::value_type("fDistanceData", &Parameter::getDistanceData));
	//
	function_table.insert(Str_to_Func::value_type("/symfunc_folder/", &Parameter::get_symfunc_folder));
	function_table.insert(Str_to_Func::value_type("/input_folder/", &Parameter::get_input_folder));
	function_table.insert(Str_to_Func::value_type("/output_folder/", &Parameter::get_output_folder));
}


void Parameter::GetParameter(std::istream & in)
{	
	int pk;
	string var;

	while ((pk = in.peek()) != EOF) {

		if (pk == ' ' || pk == '\n' || pk == '#' || pk == '\r') {
			in.ignore(1024, '\n');
			continue;
		}

		in >> var;
		if (function_table.find(var) == function_table.end()) {
			cout << "Counld not recognize key word: " << var << endl;
		}
		else {
			(this->*function_table[var])(in);
		}
		in.ignore(1024, '\n');
	}
}

bool Parameter::InitAllParameter()
{
	bool IfBad = false;

	if (random_seed == -1) {
		cout << "random_seed uninitialized!" << endl;
		IfBad = true;
	}
	if (nSample == -1) {
		cout << "nSample uninitialized!" << endl;
		IfBad = true;
	}
	if (energy_correction < 0) {
		cout << "energy_correction uninitialized!" << endl;
		IfBad = true;
	}
	if (nElement == -1) {
		cout << "nElement uninitialized!" << endl;
		IfBad = true;
	}
	if (nAtom == -1) {
		cout << "nAtom uninitialized!" << endl;
		IfBad = true;
	}
	if (nFuncInfo == -1) {
		cout << "nFuncInfo uninitialized!" << endl;
		IfBad = true;
	}

	return IfBad;
}

void Parameter::Debug(std::ostream & out)
{
	out << LMARK << "struct Parameter" << RMARK << endl;
	out << "random_seed: " << random_seed << endl;
	out << "nSample: " << nSample << endl;
	out << "energy_correction: " << energy_correction << endl;
	out << "nElement: " << nElement << endl;
	out << "nAtom: " << nAtom << endl;
	out << "element list:" << endl;
	for (int i = 0; i < nElement; ++i) {
		out << num_to_element[i] << "<->" << element_to_num[num_to_element[i]] << endl;
	}
	out << "atom list:" << endl;
	for (int i = 0; i < nAtom; ++i) {
		out << atom_list[i] << " ";
	}
	out << endl;
	out << "nFuncInfo: " << nFuncInfo << endl;
	out << "projectName: " << projectName << endl;
	out << "fNetworkInfo: " << fNetworkInfo << endl;
	out << "DataType: " << DataType << endl;
	out << "fCartesianData: " << fCartesianData << endl;
	out << "fDistanceData: " << fDistanceData << endl;
	out << "symfunc_folder: " << symfunc_folder << endl;
	out << "input_folder: " << input_folder << endl;
	out << "output_folder: " << output_folder << endl;
	out << endl;
}
