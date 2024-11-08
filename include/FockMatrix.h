#ifndef FOCKMATRIX_H
#define FOCKMATRIX_H

#include <vector>
#include <armadillo>
#include "DensityMatrix.h"
#include "OverlapMatrix.h"
#include "GammaCalculator.h"
#include "Constants.h"

class FockMatrix
{
public:
    FockMatrix(const DensityMatrix &densityMatrix, const OverlapMatrix &overlapMatrix,
               const std::vector<int> &atomicNumbers, const std::vector<double> &alphas,
               const std::vector<double> &d_total);

    arma::mat computeDiagonalElements();
    arma::mat computeOffDiagonalElements();
    double calculateDiagonalElement(int mu);
    double calculateOffDiagonalElement(int mu, int nu);

private:
    DensityMatrix densityMatrix;
    arma::mat overlapMatrix;
    arma::mat fockMatrix;
    std::vector<int> atomicNumbers;
    std::vector<double> alphas;
    std::vector<double> d_total;
};

#endif // FOCKMATRIX_H
