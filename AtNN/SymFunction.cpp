#include "SymFunction.h"
#include "FuncGroup.h"
#include "FuncInfo.h"
#include "Molecule.h"
#include "NeuralNetwork.h"
#include <Eigen/Core>

#define DECAY 0.93548387

using std::getline;
using std::setw;
using std::setprecision;
using std::left;

using namespace Eigen;

SymFunction::SymFunction() :dimX(0), X(NULL), Energy(NULL)
{
	for (int i = 0; i < parameter.nElement; ++i) {
		funcgroup.push_back(i);
	}
}

SymFunction::~SymFunction() {}

void SymFunction::Construct(istream & fin)
{
	istringstream getVar;
	int pk;
	string element;
	int nFunc;

	while ((pk = fin.peek()) != EOF) {

		if (pk == ' ' || pk == '\n' || pk == '#' || pk == '\r') {
			fin.ignore(1024, '\n');
			continue;
		}

		fin >> element >> nFunc;
		fin.ignore(1024, '\n');
		funcgroup[parameter.element_to_num[element]].Construct(fin, nFunc);
	}

	dimX = 0;
    for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
        dimX += funcgroup[parameter.atom_list[iAtom]].nFunc;
    }


#ifdef DEBUG_SYMFUNCTION

	//debug << "dimX: " << dimX << endl << endl;

#endif

	//---------- Construct X ----------
	X = new double[parameter.nSample * dimX];
	Energy = new double[parameter.nSample];

	tlist.clear();
	tlist.resize(dimX);
	int iX = 0;
	for (int i = 0; i < parameter.nAtom; ++i) {
		const FuncGroup & func = funcgroup[parameter.atom_list[i]];
		 
		for (int iFunc = 0; iFunc < func.nFunc; ++iFunc) {
			if (func.funcinfo[iFunc].symfunc <= 1) {
				const int element = func.funcinfo[iFunc].elements[0];
				for (int j = 0; j < parameter.nAtom; ++j) {
					if (i != j && element == parameter.atom_list[j]) {
						tlist[iX].push_back(Array2(j, 0));
					}
				}
			}
			else {
				const int element1 = func.funcinfo[iFunc].elements[0];
				const int element2 = func.funcinfo[iFunc].elements[1];
				for (int j = 0; j < parameter.nAtom - 1; ++j) {
					for (int k = j + 1; k < parameter.nAtom; ++k) {
						if (i != j && i != k && (element1 == parameter.atom_list[j] && element2 == parameter.atom_list[k] || element1 == parameter.atom_list[k] && element2 == parameter.atom_list[j])) {
							tlist[iX].push_back(Array2(j, k));
						}
					}
				}
			}			
			++iX;
		}
	}

#ifdef DEBUG_SYMFUNCTION

	//debug << "tlist:" << endl;
	//for (iX = 0; iX < dimX; ++iX) {
	//	vector<std::array<int, 2> >::iterator ii;
	//	
	//	for (ii = tlist[iX].begin(); ii != tlist[iX].end(); ++ii) {
	//		debug << "(" << (*ii)[0] << "," << (*ii)[1] << ") ";
	//	}
	//	debug << endl;
	//}
	//debug << endl;

#endif

}

void SymFunction::Clear()
{
	if (X)
		delete[] X;
	X = NULL;

	if (Energy)
		delete[] Energy;
	Energy = NULL;
}

void SymFunction::CalX(const long & iSample, const Molecule & molecule)
{
	double * const data = &X[iSample * dimX];
	
	Energy[iSample] = molecule.energy;

	int i, j, k;
	MatrixXd cutoff(parameter.nAtom, parameter.nAtom);
	cutoff.setZero();
	double R_ratio, tmp;

	int iX = 0;
	for (i = 0; i < parameter.nAtom; ++i) {
		const FuncGroup & func = funcgroup[parameter.atom_list[i]];

		for (int iFunc = 0; iFunc < func.nFunc; ++iFunc) {
			const FuncInfo & info = func.funcinfo[iFunc];
			const double *param = info.funcParam;
			//	G1 : {Rc}
			//	G2 : {Rc, Rs, eta}
			//	G3 : {Rc, lambda, eta, xi}
			//	G4 : {Rc, lambda, eta, xi}
			//---------- Calculate cutoff matrix ----------
			const double & Rc = param[0];
			if (info.cutoff == 0) {
				for (j = 0; j < parameter.nAtom - 1; ++j) {
					for (k = j + 1; k < parameter.nAtom; ++k) {
						R_ratio = molecule.R(j, k) / Rc;

						if (R_ratio > 1)
							cutoff(j, k) = 0;
						else {
							tmp = 0.5 * (std::cos(PI * R_ratio) + 1);
							cutoff(j, k) = (tmp < 0) ? 0 : ((tmp > 1) ? 1 : tmp);
						}
						cutoff(k, j) = cutoff(j, k);
					}
				}
			}
			else {
				for (j = 0; j < parameter.nAtom - 1; ++j) {
					for (k = j + 1; k < parameter.nAtom; ++k) {
						R_ratio = molecule.R(j, k) / Rc;

						if (R_ratio > 1)
							cutoff(j, k) = 0;
						else {
							tmp = (std::tanh(1 - R_ratio) * std::tanh(1 - R_ratio) * std::tanh(1 - R_ratio));
							cutoff(j, k) = (tmp < 0) ? 0 : ((tmp > 1) ? 1 : tmp);
						}
						cutoff(k, j) = cutoff(j, k);
					}
				}
			}
			// ---------- Calculate data ----------
			data[iX] = 0;
			vector<Array2>::iterator ii;

			switch (info.symfunc) {
			case 0: {
				for (ii = tlist[iX].begin(); ii != tlist[iX].end(); ++ii) {
					j = (*ii)[0];
					data[iX] += cutoff(i, j);
				}
				break;
			}
			case 1: {
				for (ii = tlist[iX].begin(); ii != tlist[iX].end(); ++ii) {
					j = (*ii)[0];
					data[iX] += std::exp(-1 * param[2] * std::pow(molecule.R(i, j) - param[1], 2)) * cutoff(i, j);
				}
				break;
			}
			case 2: {
				for (ii = tlist[iX].begin(); ii != tlist[iX].end(); ++ii) {
					j = (*ii)[0];
					k = (*ii)[1];
					
					data[iX] += 2 * std::pow((1 + param[1] * molecule.cos0[i](j, k)) / 2, param[3])
						* std::exp(-1 * param[2] * (molecule.R2(i, j) + molecule.R2(i, k) + molecule.R2(j, k)))
						* cutoff(i, j) * cutoff(i, k) * cutoff(j, k);
				}
				break;
			}
			case 3: {
				for (ii = tlist[iX].begin(); ii != tlist[iX].end(); ++ii) {
					j = (*ii)[0];
					k = (*ii)[1];

					data[iX] += 2 * std::pow((1 + param[1] * molecule.cos0[i](j, k)) / 2, param[3])
						* std::exp(-1 * param[2] * (molecule.R2(i, j) + molecule.R2(i, k)))
						* cutoff(i, j) * cutoff(i, k);
				}
				break;
			}
			}
			++iX;
		}
	}

#ifdef DEBUG_SYMFUNCTION
	//debug << "data:" << endl;
	//for (iX = 0; iX < dimX; ++iX) {
	//	debug << data[iX] << " ";
	//}
	//debug << Energy[iSample] << endl;
#endif

}

void SymFunction::PES_Funcinfo()
{
	string FileFuncInfo;
	ofstream fout;

	FileFuncInfo = parameter.output_folder + "PES_" + parameter.fFunctionInfo;
	fout.open(FileFuncInfo.c_str(), ofstream::out);

	for (int i = 0; i < parameter.nElement; ++i) {
		fout << funcgroup[i].nFunc << " ";
	}
	cout << endl;

	for (int i = 0; i < parameter.nElement; ++i) {
		for (int iFunc = 0; iFunc < funcgroup[i].nFunc; ++iFunc) {
			fout << funcgroup[i].funcinfo[iFunc].OutputPES() << endl;
		}
	}
	
	fout.close();
}

void SymFunction::OutputFuncInfo(ostream & fout)
{
	for (int iE = 0; iE < parameter.nElement; ++iE) {
		fout << parameter.num_to_element[iE] << " " << funcgroup[iE].nFunc << endl;

		for (int iFunc = 0; iFunc < funcgroup[iE].nFunc; ++iFunc) {
			fout << funcgroup[iE].funcinfo[iFunc].PrintFuncInfo() << endl;
		}
	}
}