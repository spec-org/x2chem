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
  void _set_submat_complex(unsigned int n, unsigned int m, const double* A, unsigned int LDA, 
    std::complex<double>* B, unsigned int LDB) 
  {
    for (auto i = 0; i < m; i++)
    for (auto j = 0; j < n; j++) {
        B[i*LDB + j] = std::complex<double>(A[i*LDA + j],0.0);
    }
  }

  // Complex matrix A is written to arbitrary subblock of B
  //
  //                            [. . .]  
  // complex(A) -> complex(B) = [. A .]
  //                            [. . .]
  //
  void _set_submat_complex(unsigned int n, unsigned int m, const std::complex<double>* A, unsigned int LDA, 
    std::complex<double>* B, unsigned int LDB) 
  {
    for (auto i = 0; i < m; i++)
    for (auto j = 0; j < n; j++) {
        B[i*LDB + j] = A[i*LDA + j];
    }
  }

  void _build_4c_core_ham(const unsigned int nb, double* V, double* p, 
                          std::complex<double>* W, std::complex<double>* core4c)
  {
    // Zero out 4C Core Ham memory
    std::fill_n(core4c,16*nb*nb,std::complex<double>(0.));


    // Set proper matrices to subblocks of 4C Hamiltonian
    
    // Transformed potential V
    _set_submat_complex(nb, nb, V, nb, core4c, 4*nb);
    _set_submat_complex(nb, nb, V, nb, core4c + nb + (4*nb * nb) , 4*nb);
    
    // Off-diagonal c*p terms
    // FIXME pick better variable names
    std::complex<double>* CP11 = core4c + 8*nb*nb;
    std::complex<double>* CP12 = CP11   + 4*nb*nb + nb;
    std::complex<double>* CP21 = core4c + 2*nb;
    std::complex<double>* CP22 = CP21   + 4*nb*nb + nb;

    for (auto i = 0; i < nb; i++) {
      CP11[i + 4*nb*i] = LIGHTSPEED + p[i];
      CP12[i + 4*nb*i] = LIGHTSPEED + p[i];
      CP21[i + 4*nb*i] = LIGHTSPEED + p[i];
      CP22[i + 4*nb*i] = LIGHTSPEED + p[i];
    }

    // Spin orbit matrix W
    _set_submat_complex(2*nb, 2*nb, W, 2*nb, core4c + 2*nb + (4*nb * 2*nb) , 4*nb);

    // Print matrix column major
    for (auto i = 0; i < 4*nb; i++) {
      std::cout << "Row " << i << ":  ";
      for (auto j = 0; j < 4*nb; j++) {
        std::cout << core4c[j*4*nb + i] << " ";
      }
      std::cout << std::endl;
    }

    

    return;
  }

  void x2c_hamiltonian(const unsigned int nb, const Integrals& ints,
    X2COperators& output, std::complex<double>* core4c)
  {
  
    //const double* S = ints.S;
    std::complex<double>* U = output.U;
    
    // Orthonormalize Basis for the transformation K
    // TK = SKt
    std::complex<double>* W = new std::complex<double>[4*nb*nb];
    std::fill_n(W,4*nb*nb,std::complex<double>(0.,1.));

    double* VSR = new double[nb*nb];
    double* VSL = new double[nb*nb];
    std::complex<double>* eig = new std::complex<double>[nb];
    double* beta = new double[nb];
    double* p = new double[nb];
    std::fill_n(p,nb,0.0);

    //auto res = lapack::ggev(lapack::Job::NoVec, lapack::Job::Vec, nb, V, nb, 
    //                        V, nb, eig, beta, VSL, nb, VSR, nb);

    // Transform non-rel potential V

    // Form spin-orbit coupling matrix W

    // Build 4C Core Hamiltonian
    _build_4c_core_ham(nb, ints.V, p, W, core4c);

  
    return;
  }
  
  void boettger_2e_soc(double *basis, double *coreH)
  {
    return;
  }

} // namespace X2Chem
