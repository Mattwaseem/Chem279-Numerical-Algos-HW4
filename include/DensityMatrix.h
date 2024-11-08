#ifndef DENSITYMATRIX_H
#define DENSITYMATRIX_H

#include <armadillo>

class DensityMatrix
{
public:
    DensityMatrix(const arma::mat &alphaCoeffs, const arma::mat &betaCoeffs);
    arma::mat computeAlphaDensityMatrix() const;
    arma::mat computeBetaDensityMatrix() const;
    arma::mat computeTotalDensityMatrix() const;

private:
    arma::mat alphaCoeffs;
    arma::mat betaCoeffs;
};

#endif // DENSITYMATRIX_H
