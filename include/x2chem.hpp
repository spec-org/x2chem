#pragma once

#include <array>
#include <complex>

namespace X2Chem {

  // Inputs to x2c_hamiltonian
  struct Integrals {
    const std::complex<double>* S;                  ///< Overlap
    const std::complex<double>* T;                  ///< Kinetic
    const std::complex<double>* V;                  ///< Nuclear potential
    const std::array<std::complex<double>*,4> pVp;  ///< Spin-orbit
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
