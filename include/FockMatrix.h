// FockMatrix.h
#ifndef FOCKMATRIX_H
#define FOCKMATRIX_H

#include "DensityMatrix.h"
#include "OverlapMatrix.h"
#include "GammaCalculator.h"
#include "CartesianGaussian.h"
#include <armadillo>
#include <vector>

class FockMatrix
{
public:
    // Updated constructor to include basisFunctions
    FockMatrix(const DensityMatrix &densityMatrix,
               const OverlapMatrix &overlapMatrix,
               const arma::mat &H_core,
               const std::vector<int> &atomicNumbersPerBasisFunction,
               const std::vector<double> &alphas,
               const std::vector<double> &d_total,
               const std::vector<CartesianGaussian> &basisFunctions); // 7th argument

    arma::mat computeFAlpha();
    arma::mat computeFBeta();

private:
    DensityMatrix densityMatrix_;
    OverlapMatrix overlapMatrix_;
    arma::mat H_core_;
    std::vector<int> atomicNumbersPerBasisFunction_;
    std::vector<double> alphas_;
    std::vector<double> d_total_;
    std::vector<CartesianGaussian> basisFunctions_;
    GammaCalculator gammaCalculator_;
};

#endif // FOCKMATRIX_H
