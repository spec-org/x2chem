#pragma once

#include <array>
#include <complex>

namespace X2Chem {

  // Inputs to x2c_hamiltonian
  struct Integrals {
    const double* S;                  ///< Overlap
    const double* T;                  ///< Kinetic
    const double* V;                  ///< Nuclear potential
    const std::array<double*,4> pVp;  ///< Spin-orbit
  };

  // Operators output by x2c_hamiltonian
  struct X2COperators {
    std::complex<double>* U;     ///< Transformation matrix
    std::complex<double>* coreH; ///< Core hamiltonian
  };
  
  // Compute X2C Core Hamiltonian
  void x2c_hamiltonian(const unsigned int, const Integrals&, X2COperators&);
  
  // Boettger 2e SOC scaling factor
  void boettger_2e_soc(double*, double*);

}
