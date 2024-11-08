#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include "CartesianGaussian.h" // Assuming needed for Gaussian functions
#include <vector>
#include <string>
#include <armadillo>

struct Atom
{
    std::string element;
    double x, y, z;
};

class Molecule
{
public:
    Molecule(); // Constructor

    // Function to compute basis functions and electrons based on the molecule's atoms
    void computeBasisFunctions();

    // Getter functions for basis functions, electrons, and atoms
    const std::vector<CartesianGaussian> &getBasisFunctions() const { return basisFunctions; }
    int getNumBasisFunctions() const { return numBasisFunctions; }
    int getNumElectrons() const { return numElectrons; }
    std::vector<Atom> &getAtoms() { return atoms; }

    // Function to add an atom (used by FileInputParser)
    void addAtom(const Atom &atom) { atoms.push_back(atom); }

private:
    std::vector<Atom> atoms;
    std::vector<CartesianGaussian> basisFunctions;
    int numBasisFunctions;
    int numElectrons;
};

#endif
