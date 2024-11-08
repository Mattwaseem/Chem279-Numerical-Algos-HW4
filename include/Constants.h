#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <unordered_map>
#include <vector>
#include <string>

struct AtomicParams
{
    double I_s;  // Ionization energy (s-orbital)
    double I_p;  // Ionization energy (p-orbital)
    double beta; // Bonding parameter
};

class Constants
{
public:
    Constants(); // Constructor to initialize the constants

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

#endif // CONSTANTS_H
