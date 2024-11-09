#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <unordered_map>
#include <vector>
#include <string>

struct AtomicParams
{
    double I_s;  // Ionization energy (s-orbital)
    double I_p;  // Ionization
    double beta; // Bonding
};

class Constants
{
public:
    Constants();

    double getIonizationPotential(int atomicNumber) const;
    double getElectronAffinity(int atomicNumber) const;
    double getValenceElectrons(int atomicNumber) const;
    double getBondingParameter(int atomicNumber) const;

private:
    std::unordered_map<int, double> electronAffinities;
    std::unordered_map<int, double> ionizationPotentials;
    std::unordered_map<int, double> bondingParameters;
    std::unordered_map<int, double> valenceElectrons;
};

#endif
