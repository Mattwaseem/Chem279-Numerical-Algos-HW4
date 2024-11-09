// FockMatrix.cpp
#include "FockMatrix.h"
#include "GammaCalculator.h"
#include "Constants.h"
#include <iostream>

FockMatrix::FockMatrix(const DensityMatrix &densityMatrix, const OverlapMatrix &overlapMatrix,
                       const std::vector<int> &atomicNumbersPerBasisFunction, const std::vector<double> &alphas,
                       const std::vector<double> &d_total)
    : densityMatrix(densityMatrix), overlapMatrix(overlapMatrix.getMatrix()),
      atomicNumbersPerBasisFunction(atomicNumbersPerBasisFunction), alphas(alphas), d_total(d_total)
{
    fockMatrix = arma::mat(densityMatrix.computeTotalDensityMatrix().n_rows,
                           densityMatrix.computeTotalDensityMatrix().n_cols, arma::fill::zeros);
}

arma::mat FockMatrix::computeDiagonalElements()
{
    std::cout << "FockMatrix computeDiagonalElements: fockMatrix size: " << fockMatrix.n_rows << " x " << fockMatrix.n_cols << std::endl;
    std::cout << "Density matrix size: " << densityMatrix.computeTotalDensityMatrix().n_rows << " x " << densityMatrix.computeTotalDensityMatrix().n_cols << std::endl;
    std::cout << "Overlap matrix size: " << overlapMatrix.n_rows << " x " << overlapMatrix.n_cols << std::endl;

    Constants constants;
    GammaCalculator gammaCalculator;

    for (int mu = 0; mu < fockMatrix.n_rows; ++mu)
    {
        int atomicNumber_mu;
        if (mu < atomicNumbersPerBasisFunction.size())
            atomicNumber_mu = atomicNumbersPerBasisFunction[mu];
        else
        {
            std::cerr << "Error: atomicNumber_mu index out of range for mu = " << mu << std::endl;
            throw std::out_of_range("atomicNumbersPerBasisFunction index out of range");
        }

        double I_mu = constants.getIonizationPotential(atomicNumber_mu);
        double A_mu = constants.getElectronAffinity(atomicNumber_mu);
        double Z_A = constants.getValenceElectrons(atomicNumber_mu);
        double P_AA = densityMatrix.computeTotalDensityMatrix()(mu, mu);
        double P_alpha_mu_mu = densityMatrix.computeAlphaDensityMatrix()(mu, mu);

        double gamma_AA = gammaCalculator.calculateGamma(mu, mu, alphas, d_total);

        double sum_B = 0.0;
        for (int B = 0; B < fockMatrix.n_cols; ++B)
        {
            if (B != mu)
            {
                if (B >= atomicNumbersPerBasisFunction.size())
                {
                    std::cerr << "Error: atomicNumber_B index out of range for B = " << B << std::endl;
                    throw std::out_of_range("atomicNumbersPerBasisFunction index out of range");
                }

                int atomicNumber_B = atomicNumbersPerBasisFunction[B];
                double Z_B = constants.getValenceElectrons(atomicNumber_B);
                double gamma_AB = gammaCalculator.calculateGamma(mu, B, alphas, d_total);
                double P_BB = densityMatrix.computeTotalDensityMatrix()(B, B);
                sum_B += (P_BB - Z_B) * gamma_AB;
            }
        }

        fockMatrix(mu, mu) = -(I_mu + A_mu) + (P_AA - Z_A - (P_alpha_mu_mu - 0.5)) * gamma_AA + sum_B;
    }
    return fockMatrix;
}

arma::mat FockMatrix::computeOffDiagonalElements()
{
    std::cout << "FockMatrix computeOffDiagonalElements: fockMatrix size: " << fockMatrix.n_rows << " x " << fockMatrix.n_cols << std::endl;
    std::cout << "Density matrix size: " << densityMatrix.computeTotalDensityMatrix().n_rows << " x " << densityMatrix.computeTotalDensityMatrix().n_cols << std::endl;
    std::cout << "Overlap matrix size: " << overlapMatrix.n_rows << " x " << overlapMatrix.n_cols << std::endl;

    Constants constants;
    GammaCalculator gammaCalculator;

    for (int mu = 0; mu < fockMatrix.n_rows; ++mu)
    {
        if (mu >= atomicNumbersPerBasisFunction.size())
        {
            std::cerr << "Error: atomicNumber_mu index out of range for mu = " << mu << std::endl;
            throw std::out_of_range("atomicNumbersPerBasisFunction index out of range");
        }
        int atomicNumber_mu = atomicNumbersPerBasisFunction[mu];

        for (int nu = 0; nu < mu; ++nu)
        {
            if (nu >= atomicNumbersPerBasisFunction.size())
            {
                std::cerr << "Error: atomicNumber_nu index out of range for nu = " << nu << std::endl;
                throw std::out_of_range("atomicNumbersPerBasisFunction index out of range");
            }
            int atomicNumber_nu = atomicNumbersPerBasisFunction[nu];

            double beta_A = constants.getBondingParameter(atomicNumber_mu);
            double beta_B = constants.getBondingParameter(atomicNumber_nu);
            double s_mu_nu = overlapMatrix(mu, nu);
            double P_alpha_mu_nu = densityMatrix.computeAlphaDensityMatrix()(mu, nu);

            double gamma_AB = gammaCalculator.calculateGamma(mu, nu, alphas, d_total);

            fockMatrix(mu, nu) = 0.5 * (beta_A + beta_B) * s_mu_nu - P_alpha_mu_nu * gamma_AB;
            fockMatrix(nu, mu) = fockMatrix(mu, nu); // Ensure symmetry
        }
    }
    return fockMatrix;
}
