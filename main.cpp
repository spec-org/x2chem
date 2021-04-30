#include <algorithm>
#include <iostream>
#include <stdio.h>
#include "x2chem.hpp"
#include <lapack.hh>

int main() {

    // Allocate test matrices
    int64_t nb = 4;
    double *S    = new double[nb*nb];
    double *T    = new double[nb*nb];
    double *V    = new double[nb*nb];
    std::array<double*,4> pVp;
    pVp[0] = new double[nb*nb];
    pVp[1] = new double[nb*nb];
    pVp[2] = new double[nb*nb];
    pVp[3] = new double[nb*nb];
    std::complex<double> *UL        = new std::complex<double>[2*nb*2*nb];
    std::complex<double> *US        = new std::complex<double>[2*nb*2*nb];
    std::complex<double> *coreX2C  = new std::complex<double>[2*nb*2*nb];
    std::complex<double> *core4C   = new std::complex<double>[4*nb*4*nb+2*nb*nb+nb];

    double* SCR1 = new double[nb*nb];
    double* SCR2 = new double[nb*nb];


    for (int64_t i = 0; i < 4*nb; i++)
    for (int64_t j = 0; j < 4*nb; j++) {
      core4C[i*4*nb + j] = std::complex<double>(0.);
    }

    // Load UH91+ custom basis 4 primitive (CQ)
    std::fill_n(S,nb*nb,0.0); 
    std::fill_n(T,nb*nb,0.0); 
    std::fill_n(V,nb*nb,0.0); 
    std::fill_n(pVp[0],nb*nb,0.0); 
    std::fill_n(pVp[1],nb*nb,0.0); 
    std::fill_n(pVp[2],nb*nb,0.0); 
    std::fill_n(pVp[3],nb*nb,0.0); 

    S[0] = 1.0;
    S[3] = -5.90215761e-01;
    S[5] = 1.0;
    S[10] = 1.0;
    S[12] = -5.90215761e-01;
    S[15] = 1.0;

    T[0] = 9.35870595e-01;
    T[3] = -3.68625145e-01;
    T[5] = 0.5;
    T[10] = 0.5;
    T[12] = -3.68625145e-01;
    T[15] = 0.5;

    V[0] = -5.49638952e+01;
    V[3] = 2.99331686e+01;
    V[5] = -4.41585567e+01;
    V[10] = -4.41585567e+01;
    V[12] = 2.99331686e+01;
    V[15] = -4.43078077e+01;

    pVp[0][0] = -1.00133396e+02;
    pVp[0][3] = 4.37558149e+01;
    pVp[0][5] = -5.29624772e+01;
    pVp[0][10] = -5.29624772e+01;
    pVp[0][12] = 4.37558149e+01;
    pVp[0][15] = -5.29481798e+01;


    pVp[1][2] = 1.43713517e+01;
    pVp[1][8] = -1.43713517e+01;
    pVp[1][11] = 1.75481407e+01;
    pVp[1][14] = -1.75481407e+01;

    pVp[2][1] = -1.43713517e+01;
    pVp[2][4] = 1.43713517e+01;
    pVp[2][7] = -1.75481407e+01;
    pVp[2][13] = 1.75481407e+01;

    pVp[3][6] = 1.76078411e+01;
    pVp[3][9] = -1.76078411e+01;


    // Load HeH+ custom basis 2 primitive (CQ)
    
    /*
    S[0] = 1.0;
    S[1] = 0.421524100;
    S[2] = 0.421524100;
    S[3] = 1.0;

    T[0] = 9.35870595e-01;
    T[1] = 2.41700233e-01;
    T[2] = 2.41700233e-01;
    T[3] = 1.73838450;

    V[0] = -2.65182646;
    V[1] = -1.53418322;
    V[2] = -1.53418322;
    V[3] = -4.14903710;
 
    //pVdotp transformed  
    //pVp[0][0] = -2.98297859;
    //pVp[0][1] = 1.88891345e-01;
    //pVp[0][2] = 1.88891345e-01;
    //pVp[0][3] = -2.20165241;

    pVp[0][0] = -4.04003775;
    pVp[0][1] = -1.06076334;
    pVp[0][2] = -1.06076334;
    pVp[0][3] = -1.04220128e+01;

    std::fill_n(pVp[1],nb*nb,0.0); 
    std::fill_n(pVp[2],nb*nb,0.0); 
    std::fill_n(pVp[3],nb*nb,0.0); 
    */


    // Load HeH+ sto-3g (pyscf)
    /*
    S[0] = 1.0;
    S[1] = 0.56088979;
    S[2] = 0.56088979;
    S[3] = 1.0;

    T[0] = 0.76003188;
    T[1] = 0.22156963;
    T[2] = 0.22156963;
    T[3] = 1.41176318;

    V[0] = -2.53536895;
    V[1] = -1.73230852;
    V[2] = -1.73230852;
    V[3] = -4.03736535;
    */

    std::copy_n(S, nb*nb, SCR1);
    auto nbu = X2Chem::orthonormalize(nb, SCR1, SCR2, 1e-12);
    X2Chem::detail::print_matrix(nb, SCR1);

    X2Chem::detail::transform(nb, nbu, T, nb, SCR1, nb, SCR2, nb, T, nbu, true);
    X2Chem::detail::transform(nb, nbu, V, nb, SCR1, nb, SCR2, nb, V, nbu, true);
    X2Chem::detail::transform(nb, nbu, pVp[0], nb, SCR1, nb, SCR2, nb, pVp[0], nbu, true);
    X2Chem::detail::transform(nb, nbu, pVp[1], nb, SCR1, nb, SCR2, nb, pVp[1], nbu, true);
    X2Chem::detail::transform(nb, nbu, pVp[2], nb, SCR1, nb, SCR2, nb, pVp[2], nbu, true);
    X2Chem::detail::transform(nb, nbu, pVp[3], nb, SCR1, nb, SCR2, nb, pVp[3], nbu, true);

    // initialize structs
    std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
    X2Chem::Integrals integrals{S,T,V,pVp_c};
    X2Chem::X2COperators output{UL, US, coreX2C};

    // X2C routine
    X2Chem::x2c_hamiltonian(nbu, integrals, output, core4C);

    double* transform = new double[nb*nb];
    std::complex<double>* eig = new std::complex<double>[2*nb];
    std::complex<double>* extraMat = new std::complex<double>[2*nb*2*nb];
    std::complex<double>* SCR3 = new std::complex<double>[2*nb*2*nb];
    std::complex<double>* SCR4 = new std::complex<double>[2*nb*2*nb];
    std::complex<double>* ULS = new std::complex<double>[2*nb*2*nb];
    std::complex<double>* USS = new std::complex<double>[2*nb*2*nb];

    std::copy_n(SCR1, nb*nb, transform);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, S, nb, SCR1, nb, 0.0, transform, nb);

    // X2Chem::detail::set_submat(nb, nb, SCR1, nb, transform, 2*nb);
    // X2Chem::detail::set_submat(nb, nb, SCR1, nb, transform+2*nb*nb+nb, 2*nb);
    // X2Chem::detail::print_matrix(2*nb, transform);
    // X2Chem::detail::transform(2*nb, extraMat, 2*nb, transform, 2*nb, SCR3, 2*nb, extraMat, 2*nb, false);
    // Top-left
    //
    std::cout << "core Ham" << std::endl;
    X2Chem::detail::print_matrix(2*nbu, coreX2C);

    std::cout << "Transform" << std::endl;
    X2Chem::detail::print_matrix(2*nbu, transform);


    X2Chem::detail::transform(nbu, nb, coreX2C, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
    X2Chem::detail::transform(nbu, nb, coreX2C + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + nb, 2*nb, false);
    X2Chem::detail::transform(nbu, nb, coreX2C + 2*nbu*nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb, 2*nb, false);
    X2Chem::detail::transform(nbu, nb, coreX2C + 2*nbu*nbu + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb + nb, 2*nb, false);

    std::cout << std::setprecision(14) << std::endl;
    for(auto i = 0; i < 4*nb*nb; i++) {
      if( std::abs(extraMat[i]) > 1e-13 ) {
        std::cout << "expectCore[" << i << "] = std::complex<double>(";
        double real = extraMat[i].real();
        double imag = extraMat[i].imag();
        std::cout << (std::abs(real) > 1e-14 ? real : 0.) << ",";
        std::cout << (std::abs(imag) > 1e-14 ? imag : 0.) << ");" << std::endl;
      }
    }

    lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, extraMat, 2*nb, transform);

    std::cout << std::setprecision(10) << std::endl;
    for( auto i = 0; i < 2*nb; i++ ) {
      std::cout << "Back eig " << i << ": " << transform[i] << std::endl;
    }

    std::copy_n(SCR1, nb*nb, transform);
    std::cout << "UL" << "\n";
    X2Chem::detail::print_matrix(2*nb, UL);

    X2Chem::detail::blockTransform(nbu, nb, UL, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
    // X2Chem::detail::transform(nbu, nb, UL, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
    // X2Chem::detail::transform(nbu, nb, UL + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + nb, 2*nb, false);
    // X2Chem::detail::transform(nbu, nb, UL + 2*nbu*nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb, 2*nb, false);
    // X2Chem::detail::transform(nbu, nb, UL + 2*nbu*nbu + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb + nb, 2*nb, false);

    // multiply S on right
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat, 2*nb, S, nb, 0.0, ULS, 2*nb);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat + nb, 2*nb, S, nb, 0.0, ULS + nb, 2*nb);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat + 2*nb*nb, 2*nb, S, nb, 0.0, ULS + 2*nb*nb, 2*nb);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat + 2*nb*nb +nb, 2*nb, S, nb, 0.0, ULS + 2*nb*nb +nb, 2*nb);

    std::cout << "ULS backtransform" << "\n";
    X2Chem::detail::print_matrix(2*nb, ULS);
    //lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, extraMat, 2*nb, transform);
    lapack::geev(lapack::Job::Vec, lapack::Job::NoVec, 2*nb, ULS, 2*nb, eig, SCR4, 2*nb, SCR3, 2*nb);
    std::cout << std::setprecision(10) << std::endl;
    for( auto i = 0; i < 2*nb; i++ ) {
      std::cout << "UL eig " << i << ": " << eig[i] << std::endl;
    }


    std::cout << "US" << "\n";
    X2Chem::detail::print_matrix(2*nb, US);

    X2Chem::detail::transform(nbu, nb, US, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
    X2Chem::detail::transform(nbu, nb, US + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + nb, 2*nb, false);
    X2Chem::detail::transform(nbu, nb, US + 2*nbu*nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb, 2*nb, false);
    X2Chem::detail::transform(nbu, nb, US + 2*nbu*nbu + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb + nb, 2*nb, false);

    // multiply S on right
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat, 2*nb, S, nb, 0.0, USS, 2*nb);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat + nb, 2*nb, S, nb, 0.0, USS + nb, 2*nb);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat + 2*nb*nb, 2*nb, S, nb, 0.0, USS + 2*nb*nb, 2*nb);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
               nb, nbu, nb, 1.0, extraMat + 2*nb*nb +nb, 2*nb, S, nb, 0.0, USS + 2*nb*nb +nb, 2*nb);

    std::cout << "USS backtransform" << "\n";
    X2Chem::detail::print_matrix(2*nb, USS);

    return 0;
}



