#include <x2chem.hpp>

void x2c_hamiltonian(int NB, double *S, double *T, double *V, double *pVp, double *pVxp, double *U, double *coreH) {
    for (auto i = 0; i < NB; i++) 
    for (auto j = 0; j < NB; j++) {
      U[i + j*NB] += 1.0;
    }
    return;
}

void boettger_2e_soc(double *basis, double *coreH) {
    return;
}
