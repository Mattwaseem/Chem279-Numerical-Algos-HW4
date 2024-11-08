#include "Molecule.h"
#include "CartesianGaussian.h" // Assuming needed for Gaussian functions
#include <vector>

Molecule::Molecule() : numBasisFunctions(0), numElectrons(0) {}

void Molecule::computeBasisFunctions()
{
    basisFunctions.clear();
    numElectrons = 0;

    for (const auto &atom : atoms)
    {
        arma::vec center = {atom.x, atom.y, atom.z};
        arma::ivec angularMomentumS = {0, 0, 0};

        // Use pre-calculated basis set information instead of defining exponents and coefficients here
        std::vector<double> exponents;
        std::vector<double> coefficients;

        // Determine the basis set for each element
        if (atom.element == "H")
        {
            exponents = {3.42525091, 0.62391373, 0.16885540};
            coefficients = {0.15432897, 0.53532814, 0.44463454};
            numElectrons += 1;
        }
        else if (atom.element == "C")
        {
            exponents = {71.6168370, 13.0450963, 3.5305122};
            coefficients = {0.15432897, 0.53532814, 0.44463454};
            numElectrons += 4;

            // Add p orbitals for carbon
            std::vector<arma::ivec> angularMomentumP = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            for (const auto &angularMomentum : angularMomentumP)
            {
                basisFunctions.emplace_back(center, exponents, angularMomentum, coefficients);
            }
        }

        // Create s-orbital CartesianGaussian for the atom
        CartesianGaussian s_orbital(center, exponents, angularMomentumS, coefficients);
        basisFunctions.push_back(s_orbital);
    }

    numBasisFunctions = basisFunctions.size();
}
