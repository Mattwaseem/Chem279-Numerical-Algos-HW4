#include "FockMatrix.h"
#include "GammaCalculator.h"
#include "Constants.h"
#include "OverlapMatrix.h"

FockMatrix::FockMatrix(const DensityMatrix &densityMatrix, const OverlapMatrix &overlapMatrix,
                       const std::vector<int> &atomicNumbers, const std::vector<double> &alphas,
                       const std::vector<double> &d_total)
    : densityMatrix(densityMatrix), overlapMatrix(overlapMatrix.getMatrix()),
      atomicNumbers(atomicNumbers), alphas(alphas), d_total(d_total)
{
    fockMatrix = arma::mat(densityMatrix.computeTotalDensityMatrix().n_rows,
                           densityMatrix.computeTotalDensityMatrix().n_cols, arma::fill::zeros);
}

arma::mat FockMatrix::computeDiagonalElements()
{
    Constants constants;
    GammaCalculator gammaCalculator;

    for (int mu = 0; mu < fockMatrix.n_rows; ++mu)
    {
        int atomicNumber_mu = atomicNumbers[mu];
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
                int atomicNumber_B = atomicNumbers[B];
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
    Constants constants;
    GammaCalculator gammaCalculator;

    for (int mu = 0; mu < fockMatrix.n_rows; ++mu)
    {
        int atomicNumber_mu = atomicNumbers[mu];
        for (int nu = 0; nu < mu; ++nu)
        {
            int atomicNumber_nu = atomicNumbers[nu];
            double beta_A = constants.getBondingParameter(atomicNumber_mu);
            double beta_B = constants.getBondingParameter(atomicNumber_nu);
            double s_mu_nu = overlapMatrix(mu, nu);
            double P_alpha_mu_nu = densityMatrix.computeAlphaDensityMatrix()(mu, nu);

            double gamma_AB = gammaCalculator.calculateGamma(mu, nu, alphas, d_total);

            fockMatrix(mu, nu) = 0.5 * (beta_A + beta_B) * s_mu_nu - P_alpha_mu_nu * gamma_AB;
            fockMatrix(nu, mu) = fockMatrix(mu, nu);
        }
    }
    return fockMatrix;
}
