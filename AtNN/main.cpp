#include "global.h"
#include "NeuralNetwork.h"
#include "SymFunction.h"
#include "Molecule.h"

Parameter parameter;
ofstream debug;

int main(void)
{

#ifdef DEBUG_MODE
	debug.open("debug.txt", ofstream::out);
#endif // DEBUG_MODE

	ifstream setup;
	NeuralNetwork *pNetwork = NULL;
	SymFunction *pSymfunc = NULL;
	vector<Molecule> molecule;

	setup.open("../input/setup.cfg", ifstream::in);
	if (setup) {
		parameter.GetParameter(setup);
		setup.close();
	}
	else {
		std::cout << "Failed to open setup.cfg!" << endl;
		return 1;
	}
#ifdef DEBUG_MODE
	parameter.Debug(debug);
#endif

	if (parameter.random_seed == 0) {
		srand(unsigned(time(NULL)));
	}
	else {
		srand(unsigned(parameter.random_seed));
	}
	try {
		string File_FuncInfo;
		string File_Molecule;
		ifstream fin;

		if (parameter.InitAllParameter()) throw "struct parameter::InitAllParameter";

		if (parameter.DataType == "XYZ")
			File_Molecule = parameter.input_folder + parameter.fCartesianData;
		else
			File_Molecule = parameter.input_folder + parameter.fDistanceData;

		fin.open(File_Molecule.c_str(), ifstream::in);
		if (!fin) {
			cout << "open molecule file failed!" << endl;
			throw "File_Molecule";
		}
		molecule.resize(parameter.nSample, Molecule());
		for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
			if (parameter.DataType == "XYZ")
				molecule[iSample].GetCartesian(fin);
			else if (parameter.DataType == "R")
				molecule[iSample].GetDistance(fin);

			molecule[iSample].CalMatrix();
		}
		fin.close();

		pSymfunc = new SymFunction;
		pNetwork = new NeuralNetwork(pSymfunc);

        string ls_order;
        string FuncInfoFiles = "tmp";
        ifstream finfo;

        ls_order = "ls " + parameter.symfunc_folder + " > " +  FuncInfoFiles;
        system(ls_order.c_str());
        finfo.open(FuncInfoFiles.c_str(), ifstream::in);

		for (int iFuncInfo = 1; iFuncInfo <= parameter.nFuncInfo; iFuncInfo++) {
            finfo >> parameter.fFunctionInfo;
			parameter.get_funcinfo_id();
			File_FuncInfo = parameter.symfunc_folder + parameter.fFunctionInfo;
			fin.open(File_FuncInfo.c_str(), ifstream::in);
            if(!fin) {
                cout << "open FuncInfo file: " << parameter.fFunctionInfo << " failed!" << endl;
                throw "open FuncInfo file";
            }

			pSymfunc->Construct(fin);
			fin.close();

			for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
				pSymfunc->CalX(iSample, molecule[iSample]);
			}
#ifdef OUTPUT_TO_SCREEN
			cout << "SymFunction::CalX()" << endl;
#endif // 

			pNetwork->Construct();
			pNetwork->TrainNetwork();
			pSymfunc->Clear();
		}

		system(("rm " + FuncInfoFiles).c_str());

        finfo.close();
        
	}
	catch (const char * error_pos) {
		std::cout << "Error in " << error_pos << std::endl;
	}

#ifdef DEBUG_MODE
	debug.close();
#endif

	return 0;
}
