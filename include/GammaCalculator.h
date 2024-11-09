#ifndef GAMMACALCULATOR_H
#define GAMMACALCULATOR_H

#include <cmath>
#include <vector>
#include <armadillo>

class GammaCalculator
{
public:
    // Calculates Gamma_AB integral between basis functions mu and nu
    double calculateGamma(size_t mu, size_t nu, const std::vector<double> &alphas, const std::vector<double> &d_total, const arma::vec &position_mu, const arma::vec &position_nu);

private:
    // Helper function to compute the integral based on distance and sigmas
    double Zero_zero(double distance, double sigmaA, double sigmaB);
};

#endif // GAMMACALCULATOR_H
