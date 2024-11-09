// FockMatrix.cpp
#include "FockMatrix.h"
#include <armadillo>
#include <iostream>

FockMatrix::FockMatrix(const DensityMatrix &densityMatrix,
                       const OverlapMatrix &overlapMatrix,
                       const arma::mat &H_core,
                       const std::vector<int> &atomicNumbersPerBasisFunction,
                       const std::vector<double> &alphas,
                       const std::vector<double> &d_total,
                       const std::vector<CartesianGaussian> &basisFunctions)
    : densityMatrix_(densityMatrix),
      overlapMatrix_(overlapMatrix),
      H_core_(H_core),
      atomicNumbersPerBasisFunction_(atomicNumbersPerBasisFunction),
      alphas_(alphas),
      d_total_(d_total),
      basisFunctions_(basisFunctions)
{
}

arma::mat FockMatrix::computeFAlpha()
{
    arma::mat Fa = H_core_;
    // Example loop to compute Fock matrix elements
    for (size_t mu = 0; mu < basisFunctions_.size(); ++mu)
    {
        for (size_t nu = 0; nu < basisFunctions_.size(); ++nu)
        {
            // Compute Gamma_AB
            double gammaAB = gammaCalculator_.calculateGamma(mu, nu, alphas_, d_total_, basisFunctions_[mu].getCenter(), basisFunctions_[nu].getCenter());

            // Update Fa with Gamma_AB and other terms as per CNDO/2
            Fa(mu, nu) += gammaAB * densityMatrix_.getP_alpha()(mu, nu); // Placeholder
        }
    }
    return Fa;
}

arma::mat FockMatrix::computeFBeta()
{
    arma::mat Fb = H_core_;

    for (size_t mu = 0; mu < basisFunctions_.size(); ++mu)
    {
        for (size_t nu = 0; nu < basisFunctions_.size(); ++nu)
        {

            double gammaAB = gammaCalculator_.calculateGamma(mu, nu, alphas_, d_total_, basisFunctions_[mu].getCenter(), basisFunctions_[nu].getCenter());

            Fb(mu, nu) += gammaAB * densityMatrix_.getP_beta()(mu, nu);
        }
    }
    return Fb;
}
