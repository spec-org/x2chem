#include <iostream>
#include <lapack.hh>
#include <x2chem.hpp>

namespace X2Chem {

  // Real part of A is written to arbitrary subblock of B
  //
  //                         [. . .]  
  // real(A) -> complex(B) = [. A .]
  //                         [. . .]
  //
  void _real_submat_to_complex(unsigned int n, unsigned int m, const double* A, unsigned int LDA, 
    std::complex<double>* B, unsigned int LDB) 
  {
    for (auto i = 0; i < m; i++)
    for (auto j = 0; j < n; j++) {
        B[i*LDB + j] = std::complex<double>(A[i*LDA + j],0.0);
    }
  }

  void _build_4c_core_ham(const unsigned int nb, const Integrals& ints, std::complex<double>* core4c)
  {
  
    // Set spatial ints to subblocks of 4C
    _real_submat_to_complex(nb, nb, ints.S, nb, core4c, 4*nb);

    for (auto i = 0; i < 4*nb; i++) 
    for (auto j = 0; j < 4*nb; j++) {
      std::cout << "4C[" << i << "," << j << "]: " << core4c[i*4*nb + j] << std::endl;
    }

    

    return;
  }

  void x2c_hamiltonian(const unsigned int nb, const Integrals& ints,
    X2COperators& output, std::complex<double>* core4c)
  {
  
    //const double* S = ints.S;
    std::complex<double>* U = output.U;
    
    // Build 4C Core Hamiltonian
    _build_4c_core_ham(nb, ints, core4c);

  
    return;
  }
  
  void boettger_2e_soc(double *basis, double *coreH)
  {
    return;
  }

} // namespace X2Chem
