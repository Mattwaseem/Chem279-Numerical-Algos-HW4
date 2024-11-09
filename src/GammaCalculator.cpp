#include "GammaCalculator.h"
#include <cmath>
#include <armadillo>
#include <iostream>

double GammaCalculator::Zero_zero(double distance, double sigmaA, double sigmaB)
{
    double U = M_PI * sigmaA * sigmaB * std::sqrt(M_PI);
    double V2 = 1.0 / (sigmaA + sigmaB);
    double T = V2 * std::pow(distance, 2);

    if (std::abs(distance) < 1e-8)
    {
        return U * 2.0 * std::sqrt(V2 / M_PI);
    }
    else
    {
        return U / std::abs(distance) * std::erf(std::sqrt(T));
    }
}

// Calculates Gamma_AB integral between basis functions mu and nu
double GammaCalculator::calculateGamma(size_t mu, size_t nu, const std::vector<double> &alphas, const std::vector<double> &d_total, const arma::vec &position_mu, const arma::vec &position_nu)
{
    // Compute the Euclidean distance between basis functions mu and nu
    double distance = arma::norm(position_mu - position_nu);

    // Retrieve sigma (exponents) for basis functions mu and nu
    double sigmaA = alphas[mu];
    double sigmaB = alphas[nu];

    // Compute Gamma_AB using Zero_zero function
    double Gamma_AB = Zero_zero(distance, sigmaA, sigmaB);

    return Gamma_AB;
}
