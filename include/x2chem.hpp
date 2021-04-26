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
    std::array<double*,4> pVp;  ///< Spin-orbit PV.P, PVxP(X), PVxP(Y), PVxP(Z) 
  };

  // Operators output by x2c_hamiltonian
  struct X2COperators {
    std::complex<double>* UL;     ///< L Transformation matrix
    std::complex<double>* US;     ///< S Transformation matrix
    std::complex<double>* coreH; ///< Core hamiltonian
  };
  
  //
  // Main API calls
  //

  // Compute X2C Core Hamiltonian
  void x2c_hamiltonian(const unsigned int, const Integrals&, X2COperators&, std::complex<double>*);
  
  // Boettger 2e SOC scaling factor
  void boettger_2e_soc(double*, double*);



  //
  // Internal subroutines
  //

  // Construct 4C Core Hamiltonian
  void _build_4c_core_ham(const unsigned int, double*, std::complex<double>*, 
                          std::complex<double>*, std::complex<double>*);

  // Form spin-orbit coupling matrix (W) 
  void _form_1e_soc_matrix(const unsigned int, std::complex<double>*, std::array<double*,4>, bool);

  // Form picture change unitary matrices UL and US
  void _form_U(const unsigned int, std::complex<double>*, std::complex<double>*,
               double*, double*, std::complex<double>*, std::complex<double>*, 
               std::complex<double>*, std::complex<double>*);


  //
  // Auxilary functions
  //

  namespace detail {

    // Inverse of a matrix using LU factorization (getrf + getri)
    void LUinv_square(int64_t, std::complex<double>*, int64_t, int64_t*);

    // Set submats of larger matrix
    template <typename sourceT, typename destT>
    void set_submat(unsigned int n, unsigned int m, const sourceT* A,
      unsigned int LDA, destT* B, unsigned int LDB)
    {
      for (auto i = 0; i < m; i++)
      for (auto j = 0; j < n; j++) {
          B[i*LDB + j] = destT(A[i*LDA + j]);
      }
    }

    // Dev functions
    template <typename T>
    void print_matrix(unsigned int N, const T* matrix)
    {
      // Print matrix column major
      for (auto i = 0; i < N; i++) {
        std::cout << "Row " << i << ":  ";
        for (auto j = 0; j < N; j++) {
          std::cout << matrix[j*N + i] << " ";
        }
        std::cout << std::endl;
      }

    }

  }

}
