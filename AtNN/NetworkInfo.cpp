#include "NetworkInfo.h"
#include "SymFunction.h"
#include "FuncGroup.h"

using std::string;
using std::ifstream;
using std::istringstream;
using std::cout;
using std::setw;
using std::left;

using std::getline;

NetworkInfo::NetworkInfo() :
	train_ratio(0), tSample(0), vSample(0),
	nFitting(0), maxEpoch(0), IfEarly(false), EarlySteps(0)
{
	nGroup = parameter.nElement;
	nNet.resize(nGroup, 0);
    for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
        nNet[parameter.atom_list[iAtom]]++;
    }
	nHidden.resize(nGroup);
	HiddenLayer.resize(nGroup);
	nLayer.resize(nGroup);
	nNeuron.resize(nGroup);
}

NetworkInfo::~NetworkInfo() {}

void NetworkInfo::GetInfo()
{
	string FileName;
	ifstream fin;

	FileName = parameter.input_folder + parameter.fNetworkInfo;	
	fin.open(FileName.c_str(), ifstream::in);
	
	int pk;
	string var;
	vector<int>::iterator ii;

	while ((pk = fin.peek()) != EOF) {

		if (pk == ' ' || pk == '\n' || pk == '#' || pk == '\r') {
			fin.ignore(1024, '\n');
			continue;
		}
		
		fin >> var;
		if (var == "nHidden") {
			for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
				fin >> nHidden[iGroup];
			}
		}
		else if (var == "HiddenLayer") {
			for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
				HiddenLayer[iGroup].resize(nHidden[iGroup]);
				for (int iHidden = 0; iHidden < nHidden[iGroup]; ++iHidden) {
					fin >> HiddenLayer[iGroup][iHidden];
				}
			}
		}
		else if (var == "train_ratio") {
			fin >> train_ratio;
			tSample = static_cast<long>(train_ratio * parameter.nSample);
			vSample = parameter.nSample - tSample;
		}
		else if (var == "nFitting") {
			fin >> nFitting;
		}
		else if (var == "maxEpoch") {
			fin >> maxEpoch;
		}
		else if (var == "mu") {
			fin >> mu;
		}
		else if (var == "IfEarly") {
			string tmp;

			fin >> tmp;
			if (tmp == "true")
				IfEarly = true;
			else if (tmp == "false")
				IfEarly = false;			
		}
		else if (var == "EarlySteps") {
			fin >> EarlySteps;
		}

		fin.ignore(1024, '\n');
	}
	fin.close();

#ifdef DEBUG_NETWORKINFO

	debug << "nHidden:" << endl;
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		debug << nHidden[iGroup] << " ";
	}
	debug << endl << endl;

	debug << "HiddenLayer: " << endl;
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		for (int iHidden = 0; iHidden < nHidden[iGroup]; ++iHidden) {
			debug << HiddenLayer[iGroup][iHidden] << " ";
		}
		debug << endl;
	}
	debug << endl;

	debug << "train_ratio: " << train_ratio << endl << endl;
	debug << "nFitting: " << nFitting << endl << endl;
	debug << "maxEpoch: " << maxEpoch << endl << endl;
	debug << "mu: " << mu << endl << endl;
	debug << "IfEarly: " << IfEarly << endl << endl;
	debug << "EarlySteps: " << EarlySteps << endl << endl;

#endif // DEBUG_NETWORKINFO

}

void NetworkInfo::Construct(const SymFunction *pSymfunc)
{
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		nLayer[iGroup] = nHidden[iGroup] + 2;
		nNeuron[iGroup].resize(nLayer[iGroup]);
		nNeuron[iGroup][0] = pSymfunc->funcgroup[iGroup].nFunc;
		
		int iHidden;
		for (iHidden = 0; iHidden < nHidden[iGroup]; ++iHidden) {
			nNeuron[iGroup][iHidden + 1] = HiddenLayer[iGroup][iHidden];
		}
		nNeuron[iGroup][iHidden + 1] = 1;
	}

#ifdef DEBUG_NETWORKINFO

	debug << "nGroup: " << nGroup << endl << endl;
	debug << "nNet: ";
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		debug << nNet[iGroup] << " ";
	}
	debug << endl << endl;
	debug << "nLayer: ";
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		debug << nLayer[iGroup] << " ";
	}
	debug << endl << endl;
	debug << "nNeuron: " << endl;
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		for (int iLayer = 0; iLayer < nLayer[iGroup]; ++iLayer) {
			debug << nNeuron[iGroup][iLayer] << " ";
		}
		debug << endl;
	}
	debug << endl;

#endif

#ifdef OUTPUT_TO_SCREEN
	cout << "NetworkInfo::Construct(const SymFunction *pSymfunc)" << endl;
#endif // OUTPUT_TO_SCREEN
}

void NetworkInfo::SaveInfo(std::ostream & fout)
{
	fout << setw(10) << left << "nGroup" << nGroup << endl;
	fout << setw(10) << left << "nNet";
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		fout << setw(5) << left << nNet[iGroup];
	}
	fout << endl;
	fout << setw(10) << left << "nLayer";
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		fout << setw(5) << left << nLayer[iGroup];
	}
	fout << endl;
	fout << setw(10) << left << "nNeuron" << endl;
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		for (int iLayer = 0; iLayer < nLayer[iGroup]; ++iLayer) {
			fout << setw(5) << left << nNeuron[iGroup][iLayer];
		}
		fout << endl;
	}
}
