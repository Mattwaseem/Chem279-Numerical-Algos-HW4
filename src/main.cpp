#include "DensityMatrix.h"
#include "FockMatrix.h"
#include "GammaCalculator.h"
#include "Constants.h"
#include "FileInputParser.h"
#include "OverlapMatrix.h"
#include <iostream>
#include <armadillo>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include "STO3GBasis.h"

namespace fs = std::filesystem;

void processMolecule(const std::string &inputFilename)
{
    FileInputParser parser;
    Molecule molecule = parser.parseMolecule(inputFilename);
    const auto &atoms = molecule.getAtoms();

    std::unordered_map<int, STO3GBasis> basisSets;
    std::vector<int> atomicNumbers;

    for (const auto &atom : atoms)
    {
        int atomicNumber = parser.getAtomicNumberFromSymbol(atom.element);
        if (basisSets.find(atomicNumber) == basisSets.end())
        {
            std::string basisFilename = parser.getBasisSetFilename(atomicNumber);
            basisSets[atomicNumber] = parser.parseBasisSet(basisFilename);
        }
        atomicNumbers.push_back(atomicNumber);
    }

    std::vector<CartesianGaussian> basisFunctions;
    for (const auto &atom : atoms)
    {
        int atomicNumber = parser.getAtomicNumberFromSymbol(atom.element);
        if (basisSets.find(atomicNumber) != basisSets.end())
        {
            const auto &basisSet = basisSets[atomicNumber];
            arma::vec position = {atom.x, atom.y, atom.z};
            arma::Col<long long> angularMomentum = {0, 0, 0};
            basisFunctions.emplace_back(position, basisSet.exponents, angularMomentum, basisSet.coefficients);
        }
    }

    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();

    arma::mat alphaCoeffs = arma::zeros(atoms.size(), atoms.size());
    arma::mat betaCoeffs = arma::zeros(atoms.size(), atoms.size());
    DensityMatrix densityMatrix(alphaCoeffs, betaCoeffs);

    arma::mat pAlpha_old = arma::zeros(atoms.size(), atoms.size());
    arma::mat pBeta_old = arma::zeros(atoms.size(), atoms.size());
    double tolerance = 1e-6;
    int maxIterations = 100;
    bool converged = false;

    arma::mat fAlpha;
    int p = 1;
    int q = 0;
    std::vector<double> alphas, d_total;

    for (int iter = 0; iter < maxIterations; ++iter)
    {
        FockMatrix fockMatrix(densityMatrix, overlapMatrix, atomicNumbers, alphas, d_total);
        fAlpha = fockMatrix.computeDiagonalElements();
        fAlpha += fockMatrix.computeOffDiagonalElements();

        arma::mat C_alpha, C_beta;
        arma::vec epsilon_alpha, epsilon_beta;
        arma::eig_sym(epsilon_alpha, C_alpha, fAlpha);

        std::cout << "fAlpha matrix size: " << fAlpha.n_rows << " x " << fAlpha.n_cols << std::endl;
        std::cout << "C_alpha matrix size: " << C_alpha.n_rows << " x " << C_alpha.n_cols << std::endl;

        p = std::min(p, static_cast<int>(C_alpha.n_cols));
        q = std::min(q, static_cast<int>(C_beta.n_cols));

        if (p > 0 && p <= C_alpha.n_cols && q > 0 && q <= C_beta.n_cols)
        {
            arma::mat pAlpha_new = C_alpha.cols(0, p - 1) * C_alpha.cols(0, p - 1).t();
            arma::mat pBeta_new = C_beta.cols(0, q - 1) * C_beta.cols(0, q - 1).t();

            if (arma::approx_equal(pAlpha_new, pAlpha_old, "absdiff", tolerance) &&
                arma::approx_equal(pBeta_new, pBeta_old, "absdiff", tolerance))
            {
                converged = true;
                break;
            }
            pAlpha_old = pAlpha_new;
            pBeta_old = pBeta_new;
        }
        else
        {
            std::cerr << "Error: Dimension mismatch in C_alpha or C_beta." << std::endl;
            return;
        }
    }

    if (!converged)
    {
        std::cerr << "SCF did not converge after " << maxIterations << " iterations." << std::endl;
    }

    arma::mat pTotal = pAlpha_old + pBeta_old;
    std::string outputDirectory = "calculated_outputs";
    if (!fs::exists(outputDirectory))
        fs::create_directory(outputDirectory);
    std::string outputFilename = outputDirectory + "/" + fs::path(inputFilename).stem().string() + "_calculated_output.txt";

    std::ofstream outputFile(outputFilename);
    if (outputFile.is_open())
    {
        outputFile << "Fock Matrix fAlpha for " << inputFilename << ":\n"
                   << fAlpha << "\n";
        outputFile << "Density Matrix (P total):\n"
                   << pTotal << std::endl;
        outputFile.close();
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file1> <input_file2> ... <input_fileN>" << std::endl;
        return 1;
    }
    for (int i = 1; i < argc; ++i)
    {
        std::string inputFilename = argv[i];
        processMolecule(inputFilename);
    }
    return 0;
}
