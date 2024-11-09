# Implementation of CNDO/2 Method for Molecular Calculations
## Numerical Algorithms Applied to Computational Quantum Chemistry
### Homework 4: Implementation of Complete Neglect of Differential Overlap (CNDO/2) Method for Molecular Calculations
#### GitHub Repository: Chem279-Numerical-Algos-HW4

#### Repository: Publicly Available (TA invited)

#### Username: MattWaseem
Assignment Overview
The goal of this assignment was to extend concepts from previous homework (EHT, Matrix Overlap calculations, etc.) to implement the CNDO/2 method, which is a semi-empirical Hartree-Fock calculation.

The key steps for this assignment included:

Parsing molecule input files (code recycled from previous homework)
Building the CNDO/2 Fock Matrix
Solving the Self Consistent Field (SCF) equations
Calculating the total energy
The assignment was organized into two main problems with a supplementary section:

Problem 1: Correctly Building the CNDO/2 Fock Matrix
Problem 2: Calculating the SCF equations with the open-shell CNDO/2 method
Problem 3: Applying the CNDO/2 method to more complex molecules such as N₂ and O₂


#### How to Compile and Run
##### Dependencies:

Ensure you have the Armadillo library installed on your system.
Preparation:

Make sure you have the input files under the sample_input directory.
Compilation:

Navigate to the build directory (create it if it doesn't exist):

1.
> mkdir build
> cd build
2.
> cmake ..
3.
> make
4.
> make run


##### Directory Structure
.
├── CMakeLists.txt
├── GammaPseudoCode.cpp
├── README.md
├── STO3G_info
│   ├── C_STO3G.txt
│   ├── F_STO3G.txt
│   ├── H_STO3G.txt
│   ├── N_STO3G.txt
│   └── O_STO3G.txt
├── bin
│   └── CNDO2
├── build
├── calculated_outputs
│   ├── H2_calculated_output.txt
│   ├── HF_calculated_output.txt
│   ├── HO_calculated_output.txt
│   ├── N2_calculated_output.txt
│   └── O2_calculated_output.txt
├── homework4.pdf
├── include
│   ├── CartesianGaussian.h
│   ├── Constants.h
│   ├── DensityMatrix.h
│   ├── FileInputParser.h
│   ├── FockMatrix.h
│   ├── GammaCalculator.h
│   ├── Molecule.h
│   ├── OverlapMatrix.h
│   └── STO3GBasis.h
├── sample_input
│   ├── H2.txt
│   ├── HF.txt
│   ├── HO.txt
│   ├── N2.txt
│   └── O2.txt
├── sample_output
│   ├── H2.txt
│   ├── HF.txt
│   └── HO.txt
└── src
    ├── CartesianGaussian.cpp
    ├── Constants.cpp
    ├── DensityMatrix.cpp
    ├── FileInputParser.cpp
    ├── FockMatrix.cpp
    ├── GammaCalculator.cpp
    ├── Molecule.cpp
    ├── OverlapMatrix.cpp
    ├── STO3GBasis.cpp
    └── main.cpp
#### Problem 1: Building the CNDO/2 Fock Matrix for Simple Molecules
#### Problem 2: Evaluating the Open-Shell CNDO/2 Self-Consistent Field Equations
