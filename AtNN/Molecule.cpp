#include "Molecule.h"
#include "SymFunction.h"
#include "FuncGroup.h"
#include "FuncInfo.h"
#include <cmath>
#include <Eigen/Core>

using namespace Eigen;


using std::setw;
using std::setprecision;
using std::left;

Molecule::Molecule() :
	energy(0), 
	Cartesian(3, parameter.nAtom),
	R(parameter.nAtom, parameter.nAtom),
	R2(parameter.nAtom, parameter.nAtom)
{
	//int length_radial = parameter.nAtom * parameter.nAtom;
	//int length_angular = parameter.nAtom * ((parameter.nAtom * (parameter.nAtom - 1)) / 2);
	R.setZero();
	R2.setZero();
    cos0.resize(parameter.nAtom);
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {	
		cos0[iAtom].resize(parameter.nAtom, parameter.nAtom);
		cos0[iAtom].setZero();		
	}
}

Molecule::~Molecule() {}

void Molecule::GetCartesian(istream & fin) {
	//Ignore first line, which contains the number of atoms
	fin.ignore(1024, '\n');

	fin >> energy;	
	fin.ignore(1024, '\n');

	string element;	
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		fin >> element;
		for (int x = 0; x < 3; ++x) {
			fin >> Cartesian(x, iAtom);
		}
		fin.ignore(1024, '\n');
	}
    

#ifdef DEBUG_MOLECULE
	debug << "Cartesian:" << endl;
	debug << Cartesian << endl << endl;
#endif
}

void Molecule::GetDistance(std::istream & fin)
{
	for (int i = 0; i < parameter.nAtom - 1; ++i) {
		for (int j = i + 1; j < parameter.nAtom; ++j) {
			fin >> R(i, j);
			R(j, i) = R(i, j);
		}
	}
	fin >> energy;
	fin.ignore(1024, '\n');
}

void Molecule::CalMatrix()
{
	if (parameter.DataType == "XYZ") {
		for (int i = 0; i < parameter.nAtom - 1; ++i) {
			for (int j = i + 1; j < parameter.nAtom; ++j) {

				R2(i, j) = (Cartesian.col(i) - Cartesian.col(j)).squaredNorm();
				R2(j, i) = R2(i, j);

				R(i, j) = std::sqrt(R2(i, j));
				R(j, i) = R(i, j);
			}
		}
	}
	else if (parameter.DataType == "R") {
		R2 = R.array().square();
	}

	for (int i = 0; i < parameter.nAtom - 2; ++i) {
		for (int j = i + 1; j < parameter.nAtom - 1; ++j) {
			for (int k = j + 1; k < parameter.nAtom; ++k) {
				cos0[i](j, k) = (R2(i, j) + R2(i, k) - R2(j, k)) / (2 * R(i, k) * R(i, j));
				cos0[j](i, k) = (R2(i, j) + R2(j, k) - R2(i, k)) / (2 * R(i, j) * R(j, k));
				cos0[k](i, j) = (R2(i, k) + R2(j, k) - R2(i, j)) / (2 * R(i, k) * R(j, k));
			}
		}
	}


#ifdef DEBUG_MOLECULE

	debug << "R2:" << endl;
	debug << R2 << endl << endl;

	debug << "R:" << endl;
	debug << R << endl << endl;

	for (int i = 0; i < parameter.nAtom; ++i) {
		debug << "cos0[" << i << "]:" << endl;
		debug << cos0[i] << endl << endl;
	}

#endif
}

void Molecule::OutputDistance(std::ostream * fout) {}
void Molecule::OutputAngle(std::ostream * fout) {}
