#include <stdio.h>
#include <iostream>
#include "x2chem.hpp"

int main() {

    // Allocate test matrices
    size_t nb = 4;
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
    std::complex<double> *core4C   = new std::complex<double>[4*nb*4*nb];


    for (auto i = 0; i < 4*nb; i++)
    for (auto j = 0; j < 4*nb; j++) {
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

    // initialize structs
    X2Chem::Integrals integrals{S,T,V,pVp};
    X2Chem::X2COperators output{UL, US, coreX2C};

    // X2C routine
    X2Chem::x2c_hamiltonian(nb, integrals, output, core4C);


    return 0;
}



