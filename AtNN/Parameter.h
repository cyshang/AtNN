#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <iostream>
#include <vector>
#include <string>
#include <map>

#ifndef LMARK
#define LMARK "<=========="
#endif

#ifndef RMARK
#define RMARK "==========>"
#endif 

//--------Typedef--------
typedef std::map<std::string, int> Str_to_Int;
typedef std::map<int, std::string> Int_to_Str;

struct Parameter {

	typedef void(Parameter::*pFunc)(std::istream & in);
	typedef std::map<std::string, pFunc> Str_to_Func;

	Str_to_Func function_table;

	int random_seed;
	void get_random_seed(std::istream & in) { in >> random_seed; }

	long nSample; //Number of Samples in Data
	void get_nSample(std::istream & in) { in >> nSample;}

	double energy_correction;	//	Conversion coefficient between different units of Energy
	void get_energy_correction(std::istream & in) { in >> energy_correction; }

	int nElement;	//	The number of elements in this system
	void get_nElement(std::istream & in) { in >> nElement; }

	int nAtom;		//	The total number of atoms in this system
	void get_nAtom(std::istream & in) { in >> nAtom; }

	Str_to_Int element_to_num;	//	A map convert an element's name into it's number
	Int_to_Str num_to_element;	//	A map convert an element's number into it's name
	void get_element_list(std::istream & in) {
		std::string elementStr;
		for (int i = 0; i < nElement; ++i) {
			in >> elementStr;
			element_to_num.insert(Str_to_Int::value_type(elementStr, i));
			num_to_element.insert(Int_to_Str::value_type(i, elementStr));
		}
	}

	std::vector<int> atom_list;	// A list recording each atom's name in this system
	void get_atom_list(std::istream & in) {
		std::string element;
		atom_list.resize(nAtom);
		for (int i = 0; i < nAtom; ++i) {
			in >> element;
			atom_list[i] = element_to_num[element];
		}
	}

	std::string run_mode;
	void get_run_mode(std::istream & in) { in >> run_mode; }

	int nFuncInfo;
	std::string fFunctionInfo;
	void get_nFuncInfo(std::istream & in) { in >> nFuncInfo; }	

	std::string funcinfo_id;
	void get_funcinfo_id()
	{
		std::string::size_type i = fFunctionInfo.rfind('.');
		if (i + 1 < fFunctionInfo.size())
			funcinfo_id = fFunctionInfo.substr(i + 1);
		else
			funcinfo_id = "";
	}

	std::string projectName;
	void getProjectName(std::istream & in) { in >> projectName; }

	std::string fNetworkInfo;
	void getNetworkInfo(std::istream & in) { in >> fNetworkInfo; }

	std::string DataType;
	void getDataType(std::istream & in) { in >> DataType; }

	std::string fCartesianData;
	void getCartesianData(std::istream & in) {
		in >> fCartesianData;
		if (fCartesianData == "NULL")
			fCartesianData = "";
	}
	
	std::string fDistanceData;
	void getDistanceData(std::istream & in) {
		in >> fDistanceData;
		if (fDistanceData == "NULL")
			fDistanceData = "";
	}

	//--------Folder Address--------
	std::string symfunc_folder;
	void get_symfunc_folder(std::istream & in) {
		in >> symfunc_folder;
		if (symfunc_folder == "NULL")
			symfunc_folder = "";
	}

	std::string input_folder;
	void get_input_folder(std::istream & in) {
		in >> input_folder;
		if (input_folder == "NULL")
			input_folder = "";
	}

	std::string output_folder;
	void get_output_folder(std::istream & in) {
		in >> output_folder;
		if (output_folder == "NULL")
			output_folder = "";
	}

	Parameter();
	void GetParameter(std::istream & in);
	bool InitAllParameter();
	void Debug(std::ostream & out = std::cout);
};

#endif // !PARAMETER_H
