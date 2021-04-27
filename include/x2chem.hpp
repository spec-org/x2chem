#pragma once

#include <array>
#include <complex>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace X2Chem {

  // Physical constants
  constexpr double LIGHTSPEED = 137.035999074; // Atomic units

  // Inputs to x2c_hamiltonian
  // FIXME: Change to const
  struct Integrals {
    double* S;                  ///< Overlap
    double* T;                  ///< Kinetic
    double* V;                  ///< Nuclear potential
    std::array<double*,4> pVp;  ///< Spin-orbit PV.P, PVxP(X), PVxP(Y), PVxP(Z) 
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
  void x2c_hamiltonian(const int64_t, const Integrals&, X2COperators&, std::complex<double>*);
  
  // Boettger 2e SOC scaling factor
  void boettger_2e_soc(int64_t, std::complex<double>*, double*, int64_t*);

  template <typename T>
  int64_t orthonormalize(int64_t, T*, double*, double);


  //
  // Internal subroutines
  //

  // Construct 4C Core Hamiltonian
  void _build_4c_core_ham(const int64_t, double*, std::complex<double>*, 
                          std::complex<double>*, std::complex<double>*);

  // Form spin-orbit coupling matrix (W) 
  void _form_1e_soc_matrix(const int64_t, std::complex<double>*, std::array<double*,4>, bool);

  // Form picture change unitary matrices UL and US
  void _form_U(const int64_t, std::complex<double>*, std::complex<double>*,
               double*, double*, std::complex<double>*, std::complex<double>*, 
               std::complex<double>*, std::complex<double>*);


  //
  // Auxilary functions
  //

  namespace detail {

    // Inverse of a matrix using LU factorization (getrf + getri)
    void LUinv_square(const int64_t, std::complex<double>*, int64_t, int64_t*);

    // Set submats of larger matrix
    template <typename sourceT, typename destT>
    void set_submat(const int64_t n, const int64_t m, const sourceT* A,
      const int64_t LDA, destT* B, const int64_t LDB)
    {
      for (auto i = 0; i < m; i++)
      for (auto j = 0; j < n; j++) {
          B[i*LDB + j] = destT(A[i*LDA + j]);
      }
    }

    // Dev functions
    template <typename T>
    void print_matrix(const int64_t N, const T* matrix)
    {

      std::cout << std::scientific << std::setprecision(5);
      // Print matrix column major
      for (auto i = 0; i < N; i++) {
        std::cout << "Row " << i << ":  ";
        for (auto j = 0; j < N; j++) {
          std::cout << std::setw(8) << matrix[j*N + i] << " ";
        }
        std::cout << std::endl;
      }

    }

  }

}
