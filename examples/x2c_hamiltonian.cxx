#include <iostream>
#include <vector>

#include <x2chem.hpp>

// This is _only_ included for matrix printing
#include <x2chem/detail.hpp> 


void setup_ints(std::vector<double>& S, std::vector<double>& T,
                std::vector<double>& V, std::array<double*,4> pVp);

int main() {

  int64_t nb = 4;
  int64_t nbsq = nb*nb;

  // Form and allocate integral objects
  std::vector<double> S(nbsq, 0.);
  std::vector<double> T(nbsq, 0.);
  std::vector<double> V(nbsq, 0.);

  std::vector<double> pvdp(nbsq, 0.);
  std::vector<double> pvxp_z(nbsq, 0.);
  std::vector<double> pvxp_y(nbsq, 0.);
  std::vector<double> pvxp_x(nbsq, 0.);

  std::array<double*,4> pVp =
    {pvdp.data(), pvxp_x.data(), pvxp_y.data(), pvxp_z.data()};

  // Form and allocate output operator objects
  std::vector<std::complex<double>> UL(4*nbsq, 0.);
  std::vector<std::complex<double>> US(4*nbsq, 0.);
  std::vector<std::complex<double>> coreX2C(4*nbsq, 0.);

  // Populate integrals
  setup_ints(S, T, V, pVp);

  // Create structs for x2c_hamiltonian
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  X2Chem::Integrals ints{S.data(), T.data(), V.data(), pVp_c};
  X2Chem::X2COperators out{UL.data(), US.data(), coreX2C.data()};

  // Form the X2C core Hamiltonian and picture change matrices
  X2Chem::x2c_hamiltonian_ao(nb, ints, out);

  // Print output operators
  std::cout << "------- X2C core Hamiltonian -------" << std::endl;
  X2Chem::detail::print_matrix(2*nb, out.coreH);
  std::cout << "------------------------------------\n" << std::endl;

  std::cout << "------- Large Picture Change Transform -------" << std::endl;
  X2Chem::detail::print_matrix(2*nb, out.UL);
  std::cout << "----------------------------------------------\n" << std::endl;

  std::cout << "------- Small Picture Change Transform -------" << std::endl;
  X2Chem::detail::print_matrix(2*nb, out.US);
  std::cout << "----------------------------------------------\n" << std::endl;

}


// Assumes size of basis is 4
void setup_ints(std::vector<double>& S, std::vector<double>& T,
                std::vector<double>& V, std::array<double*,4> pVp)
{
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

  std::cout << "------- Overlap Integral -------" << std::endl;
  X2Chem::detail::print_matrix(4, S.data());
  std::cout << "--------------------------------\n" << std::endl;

  std::cout << "------- Kinetic Integral -------" << std::endl;
  X2Chem::detail::print_matrix(4, T.data());
  std::cout << "--------------------------------\n" << std::endl;

  std::cout << "------- Potential Integral -------" << std::endl;
  X2Chem::detail::print_matrix(4, V.data());
  std::cout << "----------------------------------\n" << std::endl;

  std::cout << "------- pV dot p Integral -------" << std::endl;
  X2Chem::detail::print_matrix(4, pVp[0]);
  std::cout << "---------------------------------\n" << std::endl;

  std::cout << "------- pV cross p Integral (X) -------" << std::endl;
  X2Chem::detail::print_matrix(4, pVp[1]);
  std::cout << "---------------------------------------\n" << std::endl;

  std::cout << "------- pV cross p Integral (Y) -------" << std::endl;
  X2Chem::detail::print_matrix(4, pVp[2]);
  std::cout << "---------------------------------------\n" << std::endl;

  std::cout << "------- pV cross p Integral (Z) -------" << std::endl;
  X2Chem::detail::print_matrix(4, pVp[3]);
  std::cout << "---------------------------------------\n" << std::endl;

}

