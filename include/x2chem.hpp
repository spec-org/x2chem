#pragma once

#include <array>
#include <complex>
#include <stdexcept>
#include <string>

namespace X2Chem {

  //
  // Input and output structures
  //

  // Physical constants
  constexpr double LIGHTSPEED = 137.035999074; // Atomic units

  // Inputs to x2c_hamiltonian
  struct Integrals {
    const double* S;                  ///< Overlap
    const double* T;                  ///< Kinetic
    const double* V;                  ///< Nuclear potential
    std::array<const double*,4> pVp;  ///< Spin-orbit PV.P, PVxP(X), PVxP(Y), PVxP(Z) 
  };

  // Operators output by x2c_hamiltonian
  struct X2COperators {
    std::complex<double>* UL;     ///< L Transformation matrix
    std::complex<double>* US;     ///< S Transformation matrix
    std::complex<double>* coreH; ///< Core hamiltonian
  };

  // Exception from bad linear algebra
  class LinAlgExcept : virtual public std::runtime_error {

    int64_t code;
    std::string name;
    std::string message;

    public:

    LinAlgExcept(std::string function_name, int64_t exit_code) :
      std::runtime_error(""), code(exit_code), name(function_name) {
      message = "Linear algebra function " + function_name
                + " exited with code " + std::to_string(exit_code);
    }

    LinAlgExcept() = delete;

    int64_t exit_code() { return code; }
    const std::string& function_name() { return name; }
    const char* what() const noexcept { return message.c_str(); }

  };


  //
  // Main API calls
  //

  // Compute X2C Core Hamiltonian
  void x2c_hamiltonian(const int64_t, const Integrals&, X2COperators&, void*);

  // Compute X2C Core Hamiltonian with internal scratch
  void x2c_hamiltonian(const int64_t, const Integrals&, X2COperators&);

  // Compute X2C Core Hamiltonian with nonorthogonal integrals
  void x2c_hamiltonian_ao(const int64_t, const Integrals&, X2COperators&, void*, void*);

  // Compute X2C Core Hamiltonian with nonorthogonal integrals and internal
  //   scratch
  void x2c_hamiltonian_ao(const int64_t, const Integrals&, X2COperators&);
  
  // Boettger 2e SOC scaling factor
  void boettger_2e_soc(int64_t, std::complex<double>*, double*, int64_t*);

  // Orthonormalize matrix and give transformation
  // Only implemented for double and std::complex<double>
  template <typename T>
  int64_t orthonormalize(int64_t, T*, double*, double);

  
}
