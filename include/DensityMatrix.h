#ifndef DENSITYMATRIX_H
#define DENSITYMATRIX_H

#include <armadillo>

class DensityMatrix
{
public:
    DensityMatrix(const arma::mat &alphaCoeffs, const arma::mat &betaCoeffs, int p, int q);
    arma::mat computeAlphaDensityMatrix() const;
    arma::mat computeBetaDensityMatrix() const;
    arma::mat computeTotalDensityMatrix() const;

private:
    arma::mat alphaCoeffs;
    arma::mat betaCoeffs;
    int p; // Number of alpha electrons
    int q; // Number of beta electrons
};

#endif
