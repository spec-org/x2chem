#pragma once

#include <iomanip>

namespace X2Chem {

  namespace detail {

    //
    // Internal subroutines
    //

    // Construct 4C Core Hamiltonian
    void build_4c_core_ham(const int64_t, double*, double*,
      std::complex<double>*, std::complex<double>*); 

    // Form spin-orbit coupling matrix (W) 
    void form_1e_soc_matrix(const int64_t, std::complex<double>*, int64_t,
      std::array<double*,4>, bool);

    // Form picture change unitary matrices UL and US
    void form_picture_change(const int64_t, std::complex<double>*,
      std::complex<double>*, double*, std::complex<double>*,
      std::complex<double>*, double*, std::complex<double>*);


    //
    // Auxilary functions
    //

    // Inverse of a matrix using LU factorization (getrf + getri)
    void LUinv_square(const int64_t, std::complex<double>*, int64_t, int64_t*);

    // Transforms into or out of a basis
    template <typename OpT, typename TransT>
    void transform(int64_t n, int64_t m, OpT* A, int64_t LDA, TransT* U,
      int64_t LDU, OpT* SCR, int64_t LDS, OpT* B, int64_t LDB, bool forward);

    // Specific transform to 2 component operators with a block diagonal transform
    template <typename OpT, typename TransT>
    void blockTransform(int64_t n, int64_t m, OpT* A, int64_t LDA, TransT* U,
      int64_t LDU, OpT* SCR, int64_t LDS, OpT* B, int64_t LDB, bool forward);

    // Set submatrix of larger matrix
    template <typename sourceT, typename destT>
    void set_submat(const int64_t n, const int64_t m, const sourceT* A,
      const int64_t LDA, destT* B, const int64_t LDB)
    {
      for (auto i = 0; i < m; i++)
      for (auto j = 0; j < n; j++) {
          B[i*LDB + j] = destT(A[i*LDA + j]);
      }
    }

    // Print a matrix
    template <typename T>
    void print_matrix(const int64_t N, const T* matrix)
    {

      std::cout << std::scientific << std::setprecision(8);
      // Print matrix column major
      for (auto i = 0; i < N; i++) {
        std::cout << "Row " << i << ":  ";
        for (auto j = 0; j < N; j++) {
          std::cout << std::setw(11) << matrix[j*N + i] << " ";
        }
        std::cout << std::endl;
      }

    }

  }

}
