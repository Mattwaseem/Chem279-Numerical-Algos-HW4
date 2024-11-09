// FockMatrix.h
#ifndef FOCKMATRIX_H
#define FOCKMATRIX_H

#include <armadillo>
#include <vector>
#include "DensityMatrix.h"
#include "OverlapMatrix.h"

class FockMatrix
{
public:
    FockMatrix(const DensityMatrix &densityMatrix, const OverlapMatrix &overlapMatrix,
               const std::vector<int> &atomicNumbersPerBasisFunction, const std::vector<double> &alphas,
               const std::vector<double> &d_total);

    arma::mat computeDiagonalElements();
    arma::mat computeOffDiagonalElements();

private:
    arma::mat fockMatrix;
    arma::mat overlapMatrix;
    const DensityMatrix &densityMatrix;
    const std::vector<int> atomicNumbersPerBasisFunction;
    const std::vector<double> alphas;
    const std::vector<double> d_total;
};

#endif // FOCKMATRIX_H
