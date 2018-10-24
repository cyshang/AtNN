#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "global.h"
#include <Eigen/Core>

#ifdef DEBUG_MODE
//#define DEBUG_MOLECULE
#endif

struct SymFunction;

struct Molecule {

	double energy;
	Eigen::MatrixXd Cartesian;
	Eigen::MatrixXd R;
	Eigen::MatrixXd R2;
    std::vector<Eigen::MatrixXd> cos0;

	Molecule();
	~Molecule();
	void GetCartesian(std::istream & fin);
	void GetDistance(std::istream & fin);
	void CalMatrix();
	void OutputDistance(std::ostream * fout);
	void OutputAngle(std::ostream * fout);
};

#endif // !MOLECULE_H_
