// main.cpp
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
    std::vector<int> atomicNumbersPerBasisFunction; // One entry per basis function
    Constants constants;

    int totalValenceElectrons = 0;
    for (const auto &atom : atoms)
    {
        int atomicNumber = parser.getAtomicNumberFromSymbol(atom.element);
        if (basisSets.find(atomicNumber) == basisSets.end())
        {
            std::string basisFilename = parser.getBasisSetFilename(atomicNumber);
            basisSets[atomicNumber] = parser.parseBasisSet(basisFilename);
        }
        totalValenceElectrons += constants.getValenceElectrons(atomicNumber);

        // Determine the number of basis functions for each atom
        if (atom.element == "H")
        {
            atomicNumbersPerBasisFunction.push_back(atomicNumber); // 1 basis function (1s)
        }
        else if (atom.element == "C" || atom.element == "N" || atom.element == "O" || atom.element == "F")
        {
            atomicNumbersPerBasisFunction.push_back(atomicNumber); // 1 basis function (s)
            atomicNumbersPerBasisFunction.push_back(atomicNumber); // 1 basis function (p_x)
            atomicNumbersPerBasisFunction.push_back(atomicNumber); // 1 basis function (p_y)
            atomicNumbersPerBasisFunction.push_back(atomicNumber); // 1 basis function (p_z)
        }
        // Add similar blocks for other elements if needed
    }

    int p = totalValenceElectrons / 2;
    int q = totalValenceElectrons - p;

    // Override p and q for specific molecules if necessary
    if (inputFilename.find("H2") != std::string::npos)
    {
        p = 1;
        q = 1;
    }
    else if (inputFilename.find("HF") != std::string::npos)
    {
        p = 4;
        q = 4;
    }
    else if (inputFilename.find("HO") != std::string::npos)
    {
        p = 4;
        q = 3;
    }

    // Compute basis functions using Molecule class
    molecule.computeBasisFunctions();
    int numBasisFunctions = molecule.getNumBasisFunctions();

    // Retrieve basis functions from Molecule
    const std::vector<CartesianGaussian> &basisFunctions = molecule.getBasisFunctions();

    // Validate the size of atomicNumbersPerBasisFunction
    if (atomicNumbersPerBasisFunction.size() != numBasisFunctions)
    {
        std::cerr << "Error: Number of atomicNumbersPerBasisFunction (" << atomicNumbersPerBasisFunction.size()
                  << ") does not match numBasisFunctions (" << numBasisFunctions << ")." << std::endl;
        exit(1);
    }

    // **Debugging: Print atomicNumbersPerBasisFunction**
    std::cout << "Atomic Numbers per Basis Function:" << std::endl;
    for (size_t i = 0; i < atomicNumbersPerBasisFunction.size(); ++i)
    {
        std::cout << "Basis Function " << i << ": Atomic Number = " << atomicNumbersPerBasisFunction[i] << std::endl;
    }

    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();

    std::cout << "Overlap matrix calculated, dimensions: " << overlapMatrix.getMatrix().n_rows << " x " << overlapMatrix.getMatrix().n_cols << std::endl;

    arma::mat alphaCoeffs = arma::zeros(numBasisFunctions, numBasisFunctions);
    arma::mat betaCoeffs = arma::zeros(numBasisFunctions, numBasisFunctions);
    DensityMatrix densityMatrix(alphaCoeffs, betaCoeffs, p, q);

    arma::mat pAlpha_old = arma::zeros(numBasisFunctions, numBasisFunctions);
    arma::mat pBeta_old = arma::zeros(numBasisFunctions, numBasisFunctions);
    double tolerance = 1e-6;
    int maxIterations = 100;
    bool converged = false;

    arma::mat fAlpha;
    std::vector<double> alphas, d_total;

    for (int iter = 0; iter < maxIterations; ++iter)
    {
        FockMatrix fockMatrix(densityMatrix, overlapMatrix, atomicNumbersPerBasisFunction, alphas, d_total);
        fAlpha = fockMatrix.computeDiagonalElements();
        fAlpha += fockMatrix.computeOffDiagonalElements();

        std::cout << "Fock matrix size: " << fAlpha.n_rows << " x " << fAlpha.n_cols << std::endl;

        if (!fAlpha.is_symmetric())
        {
            std::cerr << "Warning: fAlpha matrix is not symmetric!" << std::endl;
        }
        else
        {
            std::cout << "fAlpha matrix is symmetric." << std::endl;
        }

        arma::mat C_alpha, C_beta;
        arma::vec epsilon_alpha, epsilon_beta;

        if (fAlpha.n_rows > 0 && fAlpha.n_cols > 0)
        {
            arma::eig_sym(epsilon_alpha, C_alpha, fAlpha);
            std::cout << "Eigen decomposition successful." << std::endl;
        }
        else
        {
            std::cerr << "Error: fAlpha matrix is empty or has invalid dimensions." << std::endl;
            return;
        }

        int pEffective = std::min(p, static_cast<int>(C_alpha.n_cols));
        int qEffective = (q > 0) ? std::min(q, static_cast<int>(C_beta.n_cols)) : 0;

        std::cout << "C_alpha size: " << C_alpha.n_rows << " x " << C_alpha.n_cols << std::endl;
        std::cout << "pEffective: " << pEffective << ", qEffective: " << qEffective << std::endl;

        if (pEffective > C_alpha.n_cols || qEffective > C_beta.n_cols)
        {
            std::cerr << "Error: pEffective or qEffective exceeds available columns in C_alpha or C_beta." << std::endl;
            exit(1);
        }

        arma::mat pAlpha_new = C_alpha.cols(0, pEffective - 1) * C_alpha.cols(0, pEffective - 1).t();
        arma::mat pBeta_new = arma::zeros(C_alpha.n_rows, C_alpha.n_cols);

        if (qEffective > 0 && qEffective <= C_beta.n_cols)
        {
            pBeta_new = C_beta.cols(0, qEffective - 1) * C_beta.cols(0, qEffective - 1).t();
        }

        if (arma::approx_equal(pAlpha_new, pAlpha_old, "absdiff", tolerance) &&
            arma::approx_equal(pBeta_new, pBeta_old, "absdiff", tolerance))
        {
            converged = true;
            break;
        }
        pAlpha_old = pAlpha_new;
        pBeta_old = pBeta_new;
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
