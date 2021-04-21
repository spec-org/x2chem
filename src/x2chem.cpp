#include <iostream>
#include <algorithm>
#include <lapack.hh>
#include <blas.hh>
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


  void _print_matrix(unsigned int N, const std::complex<double>* matrix)
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

  void _print_matrix(unsigned int N, const double* matrix)
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


  void _build_4c_core_ham(const unsigned int nb, double* V, std::complex<double>* p, 
                          std::complex<double>* W, std::complex<double>* core4c)
  {
    // Zero out 4C Core Ham memory
    std::fill_n(core4c,16*nb*nb,std::complex<double>(0.));


    // Set proper matrices to subblocks of 4C Hamiltonian
    
    // Transformed 1C potential V 
    _set_submat_complex(nb, nb, V, nb, core4c, 4*nb);
    _set_submat_complex(nb, nb, V, nb, core4c + nb + (4*nb * nb) , 4*nb);
    
    // Off-diagonal c*p terms
    // FIXME pick better variable names
    std::complex<double>* CP11 = core4c + 8*nb*nb;
    std::complex<double>* CP12 = CP11   + 4*nb*nb + nb;
    std::complex<double>* CP21 = core4c + 2*nb;
    std::complex<double>* CP22 = CP21   + 4*nb*nb + nb;

    for (auto i = 0; i < nb; i++) {
      CP11[i + 4*nb*i] = LIGHTSPEED * p[i];
      CP12[i + 4*nb*i] = LIGHTSPEED * p[i];
      CP21[i + 4*nb*i] = LIGHTSPEED * p[i];
      CP22[i + 4*nb*i] = LIGHTSPEED * p[i];
    }

    // Spin orbit matrix W
    _set_submat_complex(2*nb, 2*nb, W, 2*nb, core4c + 2*nb + (4*nb * 2*nb) , 4*nb);

    //std::cout << "4CCH" << "\n";
    //_print_matrix(4*nb, core4c); 

    return;
  }

  //FIXME: Why do these lapack calls require int64_t specifically??
  void _LUinv_square(int64_t n, std::complex<double>* A, int64_t lda, int64_t* ipiv)
  {
    int64_t info = lapack::getrf(n, n, A, lda, ipiv);
    if( info == 0 ) {
      auto info2 = lapack::getri(n, A, lda, ipiv);
    } else {
      std::cerr << "LU factorization failed" << std::endl;
      exit(1);
    }

    
  }

    
  void _form_1e_soc_matrix(const unsigned int nb, std::complex<double>* W, 
                           std::array<double*,4> pVp, bool soc)
  {
    if( soc ) {

      // W = [ W1  W2 ]
      //     [ W3  W4 ]
      std::complex<double> *W1 = W;
      std::complex<double> *W2 = W1 + 2*nb*nb;
      std::complex<double> *W3 = W1 + nb;
      std::complex<double> *W4 = W2 + nb;

      for(auto i = 0; i < nb; i++) {
        for(auto j = 0; j < nb; j++) {

          // W1 = pV.p + i (pVxp)(Z)
          W1[j + i*2*nb] =  pVp[0][j + i*nb] + std::complex<double>(0.,1.) * pVp[3][j + i*nb];

          // W2 = (pVxp)(Y) + i (pVxp)(X)
          W2[j + i*2*nb] =  pVp[2][j + i*nb] + std::complex<double>(0.,1.) * pVp[1][j + i*nb];

          // W3 = -(pVxp)(Y) + i (pVxp)(X)
          W3[j + i*2*nb] = -pVp[2][j + i*nb] + std::complex<double>(0.,1.) * pVp[1][j + i*nb];

          // W4 = pV.p - i (pVxp)(Z)
          W4[j + i*2*nb] =  pVp[0][j + i*nb] - std::complex<double>(0.,1.) * pVp[3][j + i*nb];
        }
      }

    } else {
      // FIXME: Scalar X2C calcs should be completely real 
      // W = [ W1  0  ]
      //     [ 0   W4 ]
      std::complex<double> *W1 = W;
      std::complex<double> *W4 = W1 + 2*nb*nb + nb;

      for(auto i = 0; i < nb; i++) {
        for(auto j = 0; j < nb; j++) {

          // W1 = pV.p 
          W1[j + i*2*nb] =  pVp[0][j + i*nb];

          // W4 = pV.p
          W4[j + i*2*nb] =  pVp[0][j + i*nb];
        }
      }
    
    }
      
    return;
  } 


  void x2c_hamiltonian(const unsigned int nb, const Integrals& ints,
    X2COperators& output, std::complex<double>* core4c)
  {
  
    //const double* S = ints.S;
    std::complex<double>* U = output.U;
    
    std::complex<double>* W = new std::complex<double>[4*nb*nb];
    std::fill_n(W,4*nb*nb,std::complex<double>(0.,1.));

    double* K = new double[nb*nb];
    double* SCR1 = new double[nb*nb];
    double* SCR2 = new double[nb*nb];
    std::complex<double>* CSCR = new std::complex<double>[4*nb*nb];
    std::complex<double>* X = new std::complex<double>[4*nb*nb];
    std::complex<double>* R = new std::complex<double>[4*nb*nb];
    std::complex<double>* eig = new std::complex<double>[nb];
    double* beta = new double[4*nb];
    int64_t* ipiv = new int64_t[2*nb];
    std::complex<double>* p = new std::complex<double>[nb];
    double* S = ints.S;
    double* T = ints.T;
    double* V = ints.V;
    double* T_copy = new double[nb*nb]; 
    double* S_copy = new double[nb*nb]; 
    double* V_tilde = new double[nb*nb]; 
    std::copy_n(S, nb*nb, S_copy);
    std::copy_n(T, nb*nb, T_copy);
    //std::fill_n(p,nb,0.0);

    //std::cout << "Overlap" << "\n";
    //_print_matrix(nb, S); 
    //std::cout << "Kinetic" << "\n";
    //_print_matrix(nb, T);

    // General eigen TK = SKt
    // to solve for transformation matrix K
    // T = 1C kinetic; S = 1C overlap
    // FIXME: orthonormalize T and S first, eigvals (p) should then be real
    auto res = lapack::ggev(lapack::Job::NoVec, lapack::Job::Vec, nb, T_copy, nb, 
                            S_copy, nb, eig, beta, SCR1, nb, K, nb);



    // K^+ S K transform
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
               nb, nb, nb, 1.0, K, nb, S, nb, 0.0, SCR1, nb);

    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nb, nb, 1.0, SCR1, nb, K, nb, 0.0, SCR2, nb);

    std::cout << "SCR2" << std::endl;
    _print_matrix(nb, SCR2);

    // Scale K by 1 / sqrt(diag( K^+ S K ))
    for (auto i = 0; i < nb; i++) 
    for (auto j = 0; j < nb; j++) 
      K[j + nb*i] /= std::sqrt(SCR2[i + i*nb]);

    // Transform non-rel potential V
    // K^+ V K
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
               nb, nb, nb, 1.0, K, nb, V, nb, 0.0, SCR1, nb);

    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nb, nb, 1.0, SCR1, nb, K, nb, 0.0, V_tilde, nb);

    std::cout << "V_tilde" << std::endl;
    _print_matrix(nb, V_tilde);

    // Transform pV.p and pVxp integrals
    // K^+ w K
    for(auto mat = 0; mat < 4; mat++) {
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nb, nb, nb, 1.0, K, nb, ints.pVp[mat], nb, 0.0, SCR1, nb);

      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nb, nb, nb, 1.0, SCR1, nb, K, nb, 0.0, ints.pVp[mat], nb);
    }

    // p^-1 = 1 / sqrt(2 * t) 
    for (auto i = 0; i < nb; i++) p[i] = 1.0 / std::sqrt( 2 * eig[i] / beta[i] );

    // Momentum transform p^-1 K^+ w K p^-1 
    // where p^-1 is a diagonal matrix
    for(auto mat = 0; mat < 4; mat++) {
      for (auto i = 0; i < nb; i++) {
        for (auto j = 0; j < nb; j++) {
          // FIXME: once orthnormal transform is fixed, real cast wont be necessary
          // anymore
          ints.pVp[mat][j + i*nb] *= std::real(p[i] * p[j]);
        }
      }
    }

    // Form full spin-orbit coupling matrix W
    _form_1e_soc_matrix(nb, W, ints.pVp, true);


    // Subtract out 2mc^2 from W diagonals
    const double Wscale = 2. * LIGHTSPEED * LIGHTSPEED;
    for(auto i = 0; i < 2*nb; i++) W[i + i*2*nb] -= Wscale;

    // P^-1 -> P
    for(auto i = 0; i < nb; i++) p[i] = 1./p[i];

    // Build 4C Core Hamiltonian
    _build_4c_core_ham(nb, V_tilde, p, W, core4c);
    std::cout << "4C" << std::endl;
    _print_matrix(4*nb, core4c);

    // Diagonalize 4C Core Hamiltonian
    lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 4*nb, core4c, 4*nb, beta);

    std::cout << "4C Diag" << std::endl;
    _print_matrix(4*nb, core4c);
    for (auto i = 0; i < 4*nb; i++) std::cout << "eig " << i << ": " << beta[i] << std::endl;
    std::cout << "MOVEC" << std::endl;
    for (auto i = 0; i < 2*nb; i++) std::cout << i << ": " << core4c[i + (2*nb)*4*nb] << " ";
    std::cout << std::endl;

    // Get large and small components of the eigenvectors
    std::complex<double>* large = core4c + 8*nb*nb;
    std::complex<double>* small = large + 2*nb;

    std::cout << "MOVEC" << std::endl;
    for (auto i = 0; i < 4*nb; i++) std::cout << i << ": " << large[i] << " ";
    std::cout << std::endl;


    // Compute inverse of large -> large^-1
    _LUinv_square(2*nb, large, 4*nb, ipiv);

    //_print_matrix(4*nb, core4c);
  
    /*
    std::complex<double>* TEST = new std::complex<double>[nb*nb]; 
    TEST[0] = 2;
    TEST[1] = 1;
    TEST[2] = 6;
    TEST[3] = 8;
    std::cout << "TEST before" << std::endl;
    _print_matrix(nb, TEST);
    _LUinv_square(nb, TEST, nb, ipiv);
    std::cout << "TEST after" << std::endl;
    _print_matrix(nb, TEST);
    */


    // X = small * large^-1
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, small, 4*nb, large, 4*nb, 0.0, X, 2*nb);

    _print_matrix(2*nb, X);

    
    // Form renormalization matrix
    // R = sqrt(1 + X^+ * X)

    // R = X^+ * X
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, X, 2*nb, X, 2*nb, 0.0, R, 2*nb);

    // R = R + I
    for(auto i = 0; i < 2*nb; i++) R[i + 2*nb*i] += 1.0;

    // Diagonalize R: R -> V * r * V^+
    lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, R, 2*nb, beta);

    // CSCR -> V * r^-0.25
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < 2*nb; j++)
      CSCR[j + 2*nb*i] = R[j + 2*nb*i] * std::pow(beta[i],-0.25);

    // R = SCR1 * SCR1**H
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, CSCR, 2*nb, CSCR, 2*nb, 0.0, R, 2*nb);

    
    std::cout << "R" << std::endl;
    _print_matrix(2*nb, R);



    // Construct 2C CH
    // 2C CH = R * (V' + cp * X + X^+ * cp + X^+ * W' * X) * R

    std::fill_n(output.coreH,4*nb*nb,std::complex<double>(0.));

    // Copy V_tilde into spin diagonal blocks of 2C CH
    _set_submat_complex(nb, nb, V_tilde, nb, output.coreH, 2*nb);
    _set_submat_complex(nb, nb, V_tilde, nb, output.coreH + nb + (2*nb * nb) , 2*nb);

    // SCR1 = cp * X
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < nb; j++) {
      CSCR[j + 2*nb*i] = LIGHTSPEED * p[j] * X[j + 2*nb*i];
      CSCR[j + nb + 2*nb*i] = LIGHTSPEED * p[j] * X[j + nb + 2*nb*i];
    }

    // 2C CH += SCR1 + SCR1^+
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < 2*nb; j++) {
      output.coreH[j + i*2*nb] += CSCR[j + i*2*nb] + std::conj(CSCR[i + j*2*nb]);
    }

    // CSCR = X**H * W
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, X, 2*nb, W, 2*nb, 0.0, CSCR, 2*nb);

    // 2C CH += CSCR * X
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, CSCR, 2*nb, X, 2*nb, 1.0, output.coreH, 2*nb);

    // CSCR = 2C CH * R
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, output.coreH, 2*nb, R, 2*nb, 0.0, CSCR, 2*nb);

    // 2C CH = Y * SCR1
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, R, 2*nb, CSCR, 2*nb, 0.0, output.coreH, 2*nb);


    // Eig test
    //lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, output.coreH, 2*nb, beta);
    //std::cout << "2C eigs" << std::endl;
    //for (auto i = 0; i < 2*nb; i++) std::cout << "eig " << i << ": " << beta[i] << " ";
    //std::cout << std::endl;
 

 
    // Transform 2C CH back

    // Form U_X2C for picture change

    return;
  }
  
  void boettger_2e_soc(double *basis, double *coreH)
  {
    return;
  }

} // namespace X2Chem
