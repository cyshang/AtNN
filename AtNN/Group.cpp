#include "Group.h"
#include "NeuralNetwork.h"
#include "NetworkInfo.h"

using namespace Eigen;
using std::vector;
using std::cout;
using std::setw;
using std::left;
using std::setprecision;

const NetworkInfo * Group::pInfo = NULL;
NeuralNetwork * Group::pNetwork = NULL;

Group::Group(const int & _group) :
	iGroup(_group), nNet(0), nLayer(0), inputStart(0), nWeight(0), WeightStart(0) {}

Group::~Group() {}

void Group::Construct() {

	nNet = pInfo->nNet[iGroup];
	nLayer = pInfo->nLayer[iGroup];
	nNeuron = pInfo->nNeuron[iGroup];

	A.clear();
	devA.clear();
	dFdZ.clear();
	W.clear();
	b.clear();
	W_copy.clear();
	b_copy.clear();

	A.resize(nLayer);
	devA.resize(nLayer);
	dFdZ.resize(nLayer);
	W.resize(nLayer);
	b.resize(nLayer);
	W_copy.resize(nLayer);
	b_copy.resize(nLayer);
	

	for (int iLayer = 0; iLayer < nLayer; ++iLayer) {
		A[iLayer].resize(nNet, MatrixXd(nNeuron[iLayer], pInfo->tSample));
		dFdZ[iLayer].resize(nNet, MatrixXd(nNeuron[iLayer], pInfo->tSample));
		devA[iLayer].resize(nNet, MatrixXd(nNeuron[iLayer], pInfo->vSample));
	}

	for (int iNet = 0; iNet < nNet; ++iNet) {
		dFdZ[nLayer - 1][iNet].setOnes();
	}

	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		W[iLayer].resize(nNeuron[iLayer], nNeuron[iLayer - 1]);
		b[iLayer].resize(nNeuron[iLayer]);
		W_copy[iLayer].resize(nNeuron[iLayer], nNeuron[iLayer - 1]);
		b_copy[iLayer].resize(nNeuron[iLayer]);
	}

	nWeight = 0;
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		nWeight += nNeuron[iLayer] * (nNeuron[iLayer - 1] + 1);
	}

	WeightStart = 0;
	for (int group = 0; group < iGroup; ++group) {
		WeightStart += pNetwork->netGroup[group].nWeight;
	}

	//------------Calculate the start point to get InputData
	inputStart = 0;
	for (int group = 0; group < iGroup; ++group) {
		inputStart += pInfo->nNet[group] * pInfo->nNeuron[group][0];
	}

#ifdef DEBUG_GROUP
	debug << "Group" << iGroup << "::Construct()" << endl;
	debug << "nWeight: " << nWeight << endl;
	debug << "WeightStart: " << WeightStart << endl;
	debug << "inputStart: " << inputStart << endl;
	debug << endl;
#endif	
}

void Group::DataInput()
{
	int input = inputStart;

	for (int iNet = 0; iNet < nNet; ++iNet) {
		A[0][iNet] = pNetwork->inputX.block(input, 0, nNeuron[0], pInfo->tSample);
		devA[0][iNet] = pNetwork->inputX.block(input, pInfo->tSample, nNeuron[0], pInfo->vSample);
		input += nNeuron[0];
	}

#ifdef DEBUG_GROUP
	debug << "Group" << iGroup << "::DataInput()" << endl;
	for (int iNet = 0; iNet < nNet; ++iNet) {
		debug << "A[0][" << iNet << "]:" << endl;
		debug << A[0][iNet] << endl << endl;
	}

	for (int iNet = 0; iNet < nNet; ++iNet) {
		debug << "devA[0][" << iNet << "]:" << endl;
		debug << devA[0][iNet] << endl << endl;
	}
	debug << endl;
#endif	
}

void Group::ForwardProp()
{
	for (int i = 1; i < nLayer - 1; ++i) {
		for (int j = 0; j < nNet; ++j) {
			A[i][j] = ((W[i] * A[i - 1][j]).colwise() + b[i]).array() / 
				sqrt(((W[i] * A[i - 1][j]).colwise() + b[i]).array().square() + 1);
		}
	}
	
	for (int iNet = 0; iNet < nNet; ++iNet) {
		A[nLayer - 1][iNet] = (W[nLayer - 1] * A[nLayer - 2][iNet]).colwise() + b[nLayer - 1];
		pNetwork->tEnergy += A[nLayer - 1][iNet];
	}
}

void Group::BackProp()
{
	for (int iLayer = nLayer - 2; iLayer > 0; --iLayer) {
		for (int iNet = 0; iNet < nNet; ++iNet) {
			dFdZ[iLayer][iNet] = (W[iLayer + 1].transpose() * dFdZ[iLayer + 1][iNet]).array() *
				(1 - A[iLayer][iNet].array().square()).array().pow(1.5);
		}
	}

	int iCol = WeightStart;
	int pre_dim, dim;

	for (int iLayer = 1; iLayer < nLayer - 1; ++iLayer) {
		pre_dim = nNeuron[iLayer - 1];
		dim = nNeuron[iLayer];
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < pre_dim; ++j) {
				for (int iNet = 0; iNet < nNet; ++iNet) {
					pNetwork->Jac.col(iCol + i * pre_dim + j).array() +=
						(dFdZ[iLayer][iNet].array().row(i) * A[iLayer - 1][iNet].array().row(j)).transpose();
				}
			}
		}
		iCol += pre_dim * dim;
		for (int iNet = 0; iNet < nNet; ++iNet) {
			pNetwork->Jac.block(0, iCol, pInfo->tSample, dim) += dFdZ[iLayer][iNet].transpose();
		}
		iCol += dim;
	}
	dim = nNeuron[nLayer - 2];
	for (int iNet = 0; iNet < nNet; ++iNet) {
		pNetwork->Jac.block(0, iCol, pInfo->tSample, dim) += A[nLayer - 2][iNet].transpose();
	}
	iCol += dim;

	pNetwork->Jac.col(iCol) = VectorXd::Ones(pInfo->tSample) * nNet;

#ifdef DEBUG_GROUP
	debug << "Group" << iGroup << "::BackProp()" << endl;
	debug << nWeight << " = " << iCol + 1 - WeightStart << endl;
#endif // OUTPUT_TO_SCREEN

}

void Group::UpdateWeight()
{
	int iRow, nRow;
	
	iRow = WeightStart;
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		nRow = nNeuron[iLayer] * nNeuron[iLayer - 1];
		W[iLayer] += Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
			(pNetwork->dWeight.segment(iRow, nRow).data(), nNeuron[iLayer], nNeuron[iLayer - 1]);
		iRow += nRow;
		b[iLayer] += pNetwork->dWeight.segment(iRow, nNeuron[iLayer]);
		iRow += nNeuron[iLayer];
	}

#ifdef DEBUG_GROUP
	debug << "Group" << iGroup << "::UpdateWeight()" << endl;
	debug << iRow - WeightStart << " = " << nWeight << endl;
#endif // DEBUG_GROUP

}

void Group::BackupWeight()
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		W_copy[iLayer] = W[iLayer];
		b_copy[iLayer] = b[iLayer];
	}
}

void Group::RestoreWeight()
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		W[iLayer] = W_copy[iLayer];
		b[iLayer] = b_copy[iLayer];
	}
}

void Group::CalDevEnergy()
{
	for (int i = 1; i < nLayer - 1; ++i) {
		for (int j = 0; j < nNet; ++j) {
			devA[i][j] = ((W[i] * devA[i - 1][j]).colwise() + b[i]).array() /
				sqrt(((W[i] * devA[i - 1][j]).colwise() + b[i]).array().square() + 1);
		}
	}

	for (int iNet = 0; iNet < nNet; ++iNet) {
		devA[nLayer - 1][iNet] = (W[nLayer - 1] * devA[nLayer - 2][iNet]).colwise() + b[nLayer - 1];
		pNetwork->vEnergy += devA[nLayer - 1][iNet];
	}
}

void Group::OutputWeight(std::ostream & outWb)
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		Map<RowVectorXd> tW(W[iLayer].data(), W[iLayer].size());
		Map<RowVectorXd> tb(b[iLayer].data(), b[iLayer].size());
		outWb << std::setprecision(16) << tW << " ";
		outWb << std::setprecision(16) << tb << " ";
	}
}

void Group::SaveWeight(std::ostream & fout)
{
	int iLayer;

	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		Map<RowVectorXd> tW(W[iLayer].data(), W[iLayer].size());
		Map<RowVectorXd> tb(b[iLayer].data(), b[iLayer].size());
		for (size_t i = 0; i < tW.size(); ++i)
			fout << setw(25) << left << tW(i);
		fout << endl;
		for (size_t i = 0; i < tb.size(); ++i)
			fout << setw(25) << left << tb(i);
		fout << endl;
	}
}


void Group::RandWeight()
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		W[iLayer].setRandom();
		b[iLayer].setRandom();
	}
}
