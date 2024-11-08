#ifndef GAMMACALCULATOR_H
#define GAMMACALCULATOR_H

#include <cmath>
#include <vector>

class GammaCalculator
{
public:
    double calculateGamma(double Ra, double Rb, const std::vector<double> &alphas, const std::vector<double> &d_total);

private:
    double Zero_zero(double Ra, double Rb, double sigmaA, double sigmaB);
};

#endif // GAMMACALCULATOR_H
