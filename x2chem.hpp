// Compute X2C Core Hamiltonian
void x2c_hamiltonian(int NB, double *S, double *T, double *V, double *pVp, 
                     double *pVxp, double *U, double *coreH);

// Boettger 2e SOC scaling factor
void boettger_2e_soc(double *basis, double *coreH);
