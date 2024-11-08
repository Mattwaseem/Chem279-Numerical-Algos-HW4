#include "DensityMatrix.h"

DensityMatrix::DensityMatrix(const arma::mat &alphaCoeffs, const arma::mat &betaCoeffs)
    : alphaCoeffs(alphaCoeffs), betaCoeffs(betaCoeffs) {}

arma::mat DensityMatrix::computeAlphaDensityMatrix() const
{
    return alphaCoeffs * alphaCoeffs.t();
}

arma::mat DensityMatrix::computeBetaDensityMatrix() const
{
    return betaCoeffs * betaCoeffs.t();
}

arma::mat DensityMatrix::computeTotalDensityMatrix() const
{
    return computeAlphaDensityMatrix() + computeBetaDensityMatrix();
}
