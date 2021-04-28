#include <iostream>
#include <algorithm>
#include <lapack.hh>
#include <blas.hh>
#include <x2chem.hpp>


#ifdef X2CHEM_DEBUG_PRINT
  #define x2chem_dbgout(name, size, location) \
    std::cout << "-----------------------------------------" << std::endl; \
    std::cout << name << std::endl; \
    detail::print_matrix(size, location); \
    std::cout << "-----------------------------------------" << std::endl;
#else
  #define x2chem_dbgout(name, size, location)
#endif


namespace X2Chem {

  namespace detail {

    template <typename OpT, typename TransT>
    void transform(int64_t n, OpT* A, int64_t LDA, TransT* U, int64_t LDU,
      OpT* SCR, int64_t LDS, OpT* B, int64_t LDB, bool forward)
    {

      blas::Op first = forward ? blas::Op::NoTrans : blas::Op::ConjTrans;
      blas::Op second = forward ? blas::Op::ConjTrans : blas::Op::NoTrans;

      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, first, 
                 n, n, n, 1.0, A, LDA, U, LDU, 0.0, SCR, LDS);

      blas::gemm(blas::Layout::ColMajor, second, blas::Op::NoTrans, 
                 n, n, n, 1.0, U, LDU, SCR, LDS, 0.0, B, LDB);

    }

    template void transform<double,double>(int64_t, double*, int64_t, double*,
        int64_t, double*, int64_t, double*, int64_t, bool);

    template void transform<std::complex<double>,double>(int64_t,
        std::complex<double>*, int64_t, double*, int64_t,
        std::complex<double>*, int64_t, std::complex<double>*, int64_t, bool);

    template void transform<std::complex<double>,std::complex<double>>(int64_t,
        std::complex<double>*, int64_t, std::complex<double>*, int64_t,
        std::complex<double>*, int64_t, std::complex<double>*, int64_t, bool);


    void LUinv_square(int64_t n, std::complex<double>* A, int64_t lda, int64_t* ipiv)
    {

      auto info = lapack::getrf(n, n, A, lda, ipiv);
      if(info != 0)
        throw LinAlgExcept("GETRF", info);

      info = lapack::getri(n, A, lda, ipiv);
      if(info != 0)
        throw LinAlgExcept("GETRI", info);

    }

  } // namespace detail


  void _build_4c_core_ham(const int64_t nb, double* V, double* p, 
                          std::complex<double>* W, std::complex<double>* core4c)
  {

    // Zero out 4C Core Ham memory
    std::fill_n(core4c,8*nb*nb,std::complex<double>(0.));
    for( auto i = 0; i < 2*nb; i++ )
      std::fill_n(core4c + 8*nb*nb + 4*nb*i, 2*nb, 0.);

    // Set proper matrices to subblocks of 4C Hamiltonian
    
    // Transformed 1C potential V 
    detail::set_submat(nb, nb, V, nb, core4c, 4*nb);
    detail::set_submat(nb, nb, V, nb, core4c + nb + (4*nb * nb) , 4*nb);
    
    // Off-diagonal c*p terms
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
    detail::set_submat(2*nb, 2*nb, W, 2*nb, core4c + 2*nb + (4*nb * 2*nb) , 4*nb);
  }


  /**
   * @brief Builds the 1e spin-orbit coupling (SOC) matrix from 
   *        real 1C pVp integrals, corresponding to the bottom right 
   *        block of the 4C Hamiltonian. 
   *
   * Modifies the W pointer, and assumes W is 2*nb x 2*nb 
   * in dimension and is contiguous in memory. The W matrix
   * is formed by 4 different nb x nb blocks. 
   *
   *      [ pV.p + i (pVxp)(Z)        |  (pVxp)(Y) + i (pVxp)(X) ]
   *  W = [ ------------------------  |  ----------------------- ]
   *      [ -(pVxp)(Y) + i (pVxp)(X)  |  pV.p - i (pVxp)(Z)      ]
   *
   *
   * @param[in]   nb        Number of spatial basis functions
   * @param[in]   W         Pointer to the W 1e SOC Matrix
   * @param[in]   pVp       Array of 4 pointers to each real pVp 
   *                        matrix input by the user, ordered 
   *                        pV.p, pVxp(X), pVxp(Y), pVxp(Z)
   * @param[in]   soc       Boolean where true computes the full
   *                        1e SOC matrix, false only computes the 
   *                        on diagonal blocks for scalar relativity
   *                        without SOC.
   *
   **/
  void _form_1e_soc_matrix(const int64_t nb, std::complex<double>* W, int64_t LDW,
                           std::array<double*,4> pVp, bool soc)
  {
    if( soc ) {

      // W = [ W1  W2 ]
      //     [ W3  W4 ]
      std::complex<double> *W1 = W;
      std::complex<double> *W2 = W1 + nb*LDW;
      std::complex<double> *W3 = W1 + nb;
      std::complex<double> *W4 = W2 + nb;

      for(auto i = 0; i < nb; i++) {
        for(auto j = 0; j < nb; j++) {

          // W1 = pV.p + i (pVxp)(Z)
          W1[j + i*LDW] =  pVp[0][j + i*nb] + std::complex<double>(0.,1.) * pVp[3][j + i*nb];

          // W2 = (pVxp)(Y) + i (pVxp)(X)
          W2[j + i*LDW] =  pVp[2][j + i*nb] + std::complex<double>(0.,1.) * pVp[1][j + i*nb];

          // W3 = -(pVxp)(Y) + i (pVxp)(X)
          W3[j + i*LDW] = -pVp[2][j + i*nb] + std::complex<double>(0.,1.) * pVp[1][j + i*nb];

          // W4 = pV.p - i (pVxp)(Z)
          W4[j + i*LDW] =  pVp[0][j + i*nb] - std::complex<double>(0.,1.) * pVp[3][j + i*nb];
        }
      }

    } else {
      // FIXME: Scalar X2C calcs should be completely real 
      // W = [ W1  0  ]
      //     [ 0   W4 ]
      std::complex<double> *W1 = W;
      std::complex<double> *W2 = W1 + nb*LDW;
      std::complex<double> *W3 = W1 + nb;
      std::complex<double> *W4 = W2 + nb;

      for(auto i = 0; i < nb; i++) {
        for(auto j = 0; j < nb; j++) {

          // W1 = pV.p 
          W1[j + i*LDW] = pVp[0][j + i*nb];

          // W4 = pV.p
          W4[j + i*LDW] = pVp[0][j + i*nb];

          // Zero out W2 and W3 blocks
          W2[j + i*LDW] = 0.0;
          W3[j + i*LDW] = 0.0;
        }
      }
    
    }

  } 


  void x2c_hamiltonian(const int64_t nb, const Integrals& ints,
    X2COperators& output, void* scratch)
  {
      
    std::cout.precision(11);
    int64_t nbsq = nb*nb;
    
    //
    // Aliases into the memory provided by the user
    //

    // First chunks are for things that need to last for the whole function
    double* p = reinterpret_cast<double*>(scratch);
    double* K = p + nb;
    double* V_tilde = K + nbsq;

    // Remaining memory is for the core hamiltonian / scratch space
    std::complex<double>* core4c = reinterpret_cast<std::complex<double>*>(V_tilde + nbsq);

    // These are continguous chunks of the core hamiltonian
    std::complex<double>* SCR1 = core4c;
    std::complex<double>* SCR2 = SCR1 + 4*nbsq;
    std::complex<double>* SCR3 = SCR2 + 4*nbsq;
    std::complex<double>* SCR4 = SCR3 + 4*nbsq;

    // Double alias for SCR1
    double* DSCR1 = reinterpret_cast<double*>(SCR1);

    //
    // Other aliases
    //

    // Use US for eigenvalues and transformed integrals
    double* eig = reinterpret_cast<double*>(output.US);

    // Use UL for W
    std::complex<double>* W = output.UL;

    // Integral short names
    double* S = ints.S;
    double* T = ints.T;
    double* V = ints.V;


    //
    // Transform into K basis
    //

    // Eigenproblem TK = Kt to solve for transformation matrix K
    // T = 1C kinetic
    std::copy_n(T, nbsq, K);
    lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, nb, K, nb, eig);

    x2chem_dbgout("K Matrix", nb, K);

    // p^-1 = 1 / sqrt(2 * t) 
    for (auto i = 0; i < nb; i++) p[i] = 1. / std::sqrt( 2. * eig[i] );


    // Transform non-rel potential V -> V_tilde
    // K^+ V K
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
               nb, nb, nb, 1.0, K, nb, V, nb, 0.0, DSCR1, nb);

    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nb, nb, 1.0, DSCR1, nb, K, nb, 0.0, V_tilde, nb);

    // Transform pV.p and pVxp integrals
    // K^+ w K
    for(auto mat = 0; mat < 4; mat++) {
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nb, nb, nb, 1.0, K, nb, ints.pVp[mat], nb, 0.0, DSCR1, nb);

      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nb, nb, nb, 1.0, DSCR1, nb, K, nb, 0.0, ints.pVp[mat], nb);
    }

    // Momentum transform p^-1 K^+ w K p^-1 
    // where p^-1 is a diagonal matrix
    for(auto mat = 0; mat < 4; mat++) {
      for (auto i = 0; i < nb; i++) {
        for (auto j = 0; j < nb; j++) {
          ints.pVp[mat][j + i*nb] *= p[i] * p[j];
        }
      }
    }

    //
    // Form 4c hamiltonian and diagonalize
    //

    // Form full spin-orbit coupling matrix W
    int64_t LDW = 2*nb;
    _form_1e_soc_matrix(nb, W, LDW, ints.pVp, true);
    x2chem_dbgout("W Matrix", 2*nb, W);

    // Subtract out 2mc^2 from W diagonals
    const double Wscale = 2. * LIGHTSPEED * LIGHTSPEED;
    for(auto i = 0; i < 2*nb; i++) W[i + i*LDW] -= Wscale;

    // P^-1 -> P
    for(auto i = 0; i < nb; i++) p[i] = 1./p[i];

    // Build 4C Core Hamiltonian
    _build_4c_core_ham(nb, V_tilde, p, W, core4c);
    x2chem_dbgout("4C Core", 4*nb, core4c);

    // Diagonalize 4C Core Hamiltonian
    lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 4*nb, core4c, 4*nb, eig);
    x2chem_dbgout("4C Diag", 4*nb, core4c);

    //
    // Compute X2C transformations
    //

    // Get large and small components of the eigenvectors
    // [ _ L ]
    // [ _ S ]
    std::complex<double>* large = core4c + 8*nbsq;
    std::complex<double>* small = core4c + 8*nbsq+2*nb;

    // Compute inverse of large -> large^-1
    detail::LUinv_square(2*nb, large, 4*nb, reinterpret_cast<int64_t*>(eig));

    // X = small * large^-1
    std::complex<double>* X = SCR2;
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, small, 4*nb, large, 4*nb, 0.0, X, 2*nb);

    
    // Form renormalization matrix
    // R = sqrt(1 + X^+ * X)

    // R = X^+ * X
    std::complex<double>* R = SCR3;
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
               2*nb, 2*nb, 2*nb, 1.0, X, 2*nb, X, 2*nb, 0.0, R, 2*nb);

    // R = R + I
    for(auto i = 0; i < 2*nb; i++) R[i + 2*nb*i] += 1.0;

    // Diagonalize R: R -> V * r * V^+
    lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, R, 2*nb, eig);

    // CSCR -> V * r^-0.25
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < 2*nb; j++)
      SCR1[j + i*2*nb] = R[j + i*2*nb] * std::pow(eig[i],-0.25);

    // R = SCR1 * SCR1**H
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, SCR1, 2*nb, SCR1, 2*nb, 0.0, R, 2*nb);

    x2chem_dbgout("R", 2*nb, R);



    //
    // Construct 2C CH
    //

    // 2C CH = R * (V' + cp * X + X^+ * cp + X^+ * W' * X) * R

    std::fill_n(output.coreH,4*nbsq,std::complex<double>(0.));

    // Copy V_tilde into spin diagonal blocks of 2C CH
    detail::set_submat(nb, nb, V_tilde, nb, output.coreH, 2*nb);
    detail::set_submat(nb, nb, V_tilde, nb, output.coreH + nb + (2*nbsq) , 2*nb);

    // CSCR = cp * X
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < nb; j++) {
      SCR1[j + 2*nb*i] = LIGHTSPEED * p[j] * X[j + 2*nb*i];
      SCR1[j + nb + 2*nb*i] = LIGHTSPEED * p[j] * X[j + nb + 2*nb*i];
    }

    // 2C CH += CSCR + CSCR^+
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < 2*nb; j++) {
      output.coreH[j + i*2*nb] += SCR1[j + i*2*nb] + std::conj(SCR1[i + j*2*nb]);
    }

    // CSCR = X**H * W
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, X, 2*nb, W, LDW, 0.0, SCR1, 2*nb);

    // 2C CH += CSCR * X
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, SCR1, 2*nb, X, 2*nb, 1.0, output.coreH, 2*nb);

    // CSCR = 2C CH * R
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, output.coreH, 2*nb, R, 2*nb, 0.0, SCR1, 2*nb);

    // 2C CH = R * CSCR
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 2*nb, 2*nb, 2*nb, 1.0, R, 2*nb, SCR1, 2*nb, 0.0, output.coreH, 2*nb);

    x2chem_dbgout("2C Core in K basis", 2*nb, output.coreH);

    //
    // Back Transform 2C CH
    //

    // K^-1 = K\dag

    // Transform 2C CH by block
    // K H K\dag
    // Top-left
    detail::transform(nb, output.coreH, 2*nb, K, nb, SCR1, nb, output.coreH, 2*nb, false);
    // Bottom-left
    detail::transform(nb, output.coreH + nb, 2*nb, K, nb, SCR1, nb, output.coreH + nb, 2*nb, false);
    // Top-right
    detail::transform(nb, output.coreH + 2*nb*nb, 2*nb, K, nb, SCR1, nb, output.coreH + 2*nbsq, 2*nb, false);
    // Bottom-right
    detail::transform(nb, output.coreH + 2*nb*nb + nb, 2*nb, K, nb, SCR1, nb, output.coreH + 2*nbsq + nb, 2*nb, false);

    x2chem_dbgout("2C Core in ortho basis", 2*nb, output.coreH);


    //
    // Form U_X2C for picture change
    //

    _form_U(nb, output.UL, output.US, K, X, R, p, SCR1); 

  }
  
  // Form picture change unitary matrices UL and US
  void _form_U(const int64_t nb, std::complex<double>* UL, std::complex<double>* US,
               double* K, std::complex<double>* X, std::complex<double>* R, 
               double* p, std::complex<double>* SCR) 
  {

    // Form UL = K * R * K^-1
    // Top-left
    detail::transform(nb, R, 2*nb, K, nb, SCR, nb, UL, 2*nb, false);
    // Bottom-left
    detail::transform(nb, R + nb, 2*nb, K, nb, SCR, nb, UL + nb, 2*nb, false);
    // Top-right
    detail::transform(nb, R + 2*nb*nb, 2*nb, K, nb, SCR, nb, UL + 2*nb*nb, 2*nb, false);
    // Bottom-right
    detail::transform(nb, R + 2*nb*nb + nb, 2*nb, K, nb, SCR, nb, UL + 2*nb*nb + nb, 2*nb, false);


    // Form US = K * 2cp^-1 * X * R * K^-1
    // (Overwrites R)

    // X*R -> R
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               2*nb, 2*nb, 2*nb, 1.0, X, 2*nb, R, 2*nb, 0.0, R, 2*nb);

    x2chem_dbgout("X*R", 2*nb, R);

    // R = 2 * c * p^-1 * R  
    for(auto i = 0; i < 2*nb; i++)
    for(auto j = 0; j < nb; j++) {
      R[j + 2*nb*i] = 2*LIGHTSPEED * (1./p[j]) * R[j + 2*nb*i];
      R[j + nb + 2*nb*i] = 2*LIGHTSPEED * (1./p[j]) * R[j + nb + 2*nb*i];
    }

    // Transform K * R * K^-1
    // where R = 2 * c * p^-1 * R * X
    // Top-left
    detail::transform(nb, R, 2*nb, K, nb, SCR, nb, US, 2*nb, false);
    // Bottom-left
    detail::transform(nb, R + nb, 2*nb, K, nb, SCR, nb, US + nb, 2*nb, false);
    // Top-right
    detail::transform(nb, R + 2*nb*nb, 2*nb, K, nb, SCR, nb, US + 2*nb*nb, 2*nb, false);
    // Bottom-right
    detail::transform(nb, R + 2*nb*nb + nb, 2*nb, K, nb, SCR, nb, US + 2*nb*nb + nb, 2*nb, false);

  }

  /**
   * @brief Scales the X2C core Hamiltonian using
   *        the Boettger 2-electron spin-orbit coupling
   *        formula. 
   *
   * Modifies the X2C core Hamiltonian in place with computed
   * Boettger values to approximate two-electron spin-orbit
   * coupling. Assumes coreH is 2*nb x 2*nb in dimension and
   * is contiguous in memory.
   *
   * @param[in]   nb        Number of spatial basis functions
   * @param[in]   coreH     Pointer to the X2C Core Hamiltonian
   * @param[in]   nucList   Pointer to an nb-length array of atomic nuclear charge
   *                        values associated with each spatial basis function.
   * @param[in]   angList   Pointer to an nb-length array of angular momentum
   *                        values for each spatial basis function.
   *
   **/
  void boettger_2e_soc(int64_t nb, std::complex<double>* coreH, double* nucList, int64_t* angList)
  {

    std::complex<double> factor;
    std::complex<double> on_diag;

    // Blocks of the X2C Core Hamiltonian
    std::complex<double>* HAA = coreH;
    std::complex<double>* HAB = coreH + nb;
    std::complex<double>* HBA = coreH + nb * 2*nb;
    std::complex<double>* HBB = coreH + nb * 2*nb + nb;

    // Q factors indexed by angular momentum s, p, d...
    std::array<double,6> Ql={0.,2.,10.,28.,60.,110.};
 
    // Modify each block of the X2C Core Hamiltonian
    for( auto i = 0; i < nb; i++ ) {
      for( auto j = 0; j < nb; j++ ) {

        // Compute factor f
        factor = -1.0 * std::sqrt( Ql[angList[j]] * Ql[angList[i]] / nucList[j] / nucList[i]);

        // Scale elements of each block
        on_diag = (factor / 2.0) * (HAA[j + i*2*nb] - HBB[j + i*2*nb]);
        HAA[j + i*2*nb] += on_diag;
        HBB[j + i*2*nb] -= on_diag;
        HAB[j + i*2*nb] += factor * HAB[j + i*2*nb];
        HBA[j + i*2*nb] += factor * HBA[j + i*2*nb];

      }
    }   

    return;
  }

  /**
   * @brief Orthonormalizes a matrix and gives the transformation matrix
   *
   * Orthonormalizes basisMat and gives the transformation matrix.
   * Assumes square matrix that is contiguous in memory.
   * Truncates vectors with less than screen singular value,
   * and returns the number of basis vectors left after truncation.
   *
   * @param[in]   N         Dimension of basis
   * @param[in]   basisMat  Basis matrix for the orthonormalization (N x N)
   * @param[in]   scratch   Scratch space (N)
   * @param[in]   screen    Minimum singular value for keeping orthogonal
   *                        vectors
   *
   * @returns     Number of significant basis functions remaining
   **/
  template <typename T>
  int64_t orthonormalize(int64_t N, T* basisMat, double* scratch, double screen)
  {

    // Do the SVD
    int64_t code = lapack::gesvd(lapack::Job::OverwriteVec, lapack::Job::NoVec,
                                 N, N, basisMat, N, scratch, nullptr, N, nullptr, N);
    if( code != 0 )
      throw LinAlgExcept("GESVD", code);

    // Find number of singular values above tolerance
    int64_t nSigVec = std::distance(
        scratch,
        std::find_if(scratch, scratch+N,
          [&](double x){ return x < screen; }
        )
    );

    // Scale columns of basisMat by 1/sqrt(s)
    for(auto i = 0; i < nSigVec; i++) {
      blas::scal(N, T(1.)/std::sqrt(scratch[i]), basisMat + i*N, 1);
    }

    return nSigVec;

  }

  template int64_t orthonormalize<double>(int64_t, double*, double*, double);
  template int64_t orthonormalize<std::complex<double>>(int64_t, std::complex<double>*, double*, double);

} // namespace X2Chem
