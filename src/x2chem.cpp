#include <x2chem.hpp>

namespace X2Chem {

  void x2c_hamiltonian(const unsigned int nb, const Integrals& ints,
    X2COperators& output)
  {
  
    std::complex<double>* U = output.U;
    
    for (auto i = 0; i < nb; i++) 
    for (auto j = 0; j < nb; j++) {
      U[i + j*nb] += 1.0;
    }
  
    return;
  }
  
  void boettger_2e_soc(double *basis, double *coreH)
  {
    return;
  }

} // namespace X2Chem
