#include "DensityMatrix.h"
#include <iostream>

DensityMatrix::DensityMatrix(const arma::mat &alphaCoeffs, const arma::mat &betaCoeffs, int p, int q)
{
    std::cout << "alphaCoeffs size: " << alphaCoeffs.n_rows << " x " << alphaCoeffs.n_cols << std::endl;
    std::cout << "betaCoeffs size: " << betaCoeffs.n_rows << " x " << betaCoeffs.n_cols << std::endl;
    std::cout << "p: " << p << ", q: " << q << std::endl;

    if (p <= alphaCoeffs.n_cols)
    {
        this->alphaCoeffs = alphaCoeffs.cols(0, p - 1);
    }
    else
    {
        std::cerr << "Error: p exceeds alphaCoeffs column count!" << std::endl;
        exit(1);
    }

    if (q <= betaCoeffs.n_cols)
    {
        this->betaCoeffs = betaCoeffs.cols(0, q - 1);
    }
    else
    {
        std::cerr << "Error: q exceeds betaCoeffs column count!" << std::endl;
        exit(1);
    }

    this->p = p;
    this->q = q;
}

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
