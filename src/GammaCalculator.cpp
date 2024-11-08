#include "GammaCalculator.h"
#include <cmath>
#include <armadillo>

double GammaCalculator::Zero_zero(double Ra, double Rb, double sigmaA, double sigmaB)
{
    double U = M_PI * sigmaA * sigmaB * std::sqrt(M_PI);
    double V2 = 1.0 / (sigmaA + sigmaB);
    double T = V2 * std::pow(Ra - Rb, 2);

    if (std::abs(Ra - Rb) < 1e-8)
    {
        return U * 2.0 * std::sqrt(V2 / M_PI);
    }
    else
    {
        return U / std::abs(Ra - Rb) * std::erf(std::sqrt(T));
    }
}

double GammaCalculator::calculateGamma(double Ra, double Rb, const std::vector<double> &alphas, const std::vector<double> &d_total)
{
    double gamma = 0.0;
    for (size_t k1 = 0; k1 < alphas.size(); ++k1)
    {
        for (size_t k2 = 0; k2 < alphas.size(); ++k2)
        {
            double sigmaA = 1.0 / (alphas[k1] + alphas[k2]);
            for (size_t l1 = 0; l1 < alphas.size(); ++l1)
            {
                for (size_t l2 = 0; l2 < alphas.size(); ++l2)
                {
                    double sigmaB = 1.0 / (alphas[l1] + alphas[l2]);
                    gamma += d_total[k1] * d_total[k2] * d_total[l1] * d_total[l2] * Zero_zero(Ra, Rb, sigmaA, sigmaB);
                }
            }
        }
    }
    return gamma * 27.211;
}
