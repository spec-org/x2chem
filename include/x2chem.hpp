#pragma once

#include <array>
#include <complex>

namespace X2Chem {

  // Physical constants
  constexpr double LIGHTSPEED = 137.035999074; // Atomic units

  // Inputs to x2c_hamiltonian
  // FIXME: Change to const
  struct Integrals {
    double* S;                  ///< Overlap
    double* T;                  ///< Kinetic
    double* V;                  ///< Nuclear potential
    std::array<double*,4> pVp;  ///< Spin-orbit
  };

  // Operators output by x2c_hamiltonian
  struct X2COperators {
    std::complex<double>* U;     ///< Transformation matrix
    std::complex<double>* coreH; ///< Core hamiltonian
  };
  
  // Compute X2C Core Hamiltonian
  void x2c_hamiltonian(const unsigned int, const Integrals&, X2COperators&, std::complex<double>*);
  
  // Boettger 2e SOC scaling factor
  void boettger_2e_soc(double*, double*);

  // Construct 4C Core Hamiltonian
  void _build_4c_core_ham(const unsigned int, double*, std::complex<double>*, 
                          std::complex<double>*, std::complex<double>*);

  // Form spin-orbit coupling matrix (W) 
  void _form_1e_soc_matrix(const unsigned int, std::complex<double>*, std::array<double*,4>, bool);


  // Auxilary functions
  void _set_submat_complex(unsigned int, unsigned int, const double*, unsigned int, 
    std::complex<double>*, unsigned int); 
  void _set_submat_complex(unsigned int, unsigned int, const std::complex<double>*, unsigned int, 
    std::complex<double>*, unsigned int);

  // Dev functions
  void _print_matrix(unsigned int, const std::complex<double>*);
  void _print_matrix(unsigned int, const double*);


}
