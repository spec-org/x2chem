#include <algorithm>

#include <gtest/gtest.h>

#include <lapack.hh>

#include <x2chem.hpp>


using namespace X2Chem;

TEST( X2C_Hamiltonian, UH_91Plus ) {

  int64_t nb = 4;
  int64_t nbsq = nb*nb;

  //
  // Allocate memory
  //

  std::vector<double> S(nbsq, 0.);
  std::vector<double> T(nbsq, 0.);
  std::vector<double> V(nbsq, 0.);

  std::vector<double> pvdp(nbsq, 0.);
  std::vector<double> pvxp_z(nbsq, 0.);
  std::vector<double> pvxp_y(nbsq, 0.);
  std::vector<double> pvxp_x(nbsq, 0.);

  std::array<double*,4> pVp =
    {pvdp.data(), pvxp_z.data(), pvxp_y.data(), pvxp_x.data()};


  std::vector<std::complex<double>> UL(4*nbsq, 0.);
  std::vector<std::complex<double>> US(4*nbsq, 0.);
  std::vector<std::complex<double>> coreX2C(4*nbsq, 0.);

  std::vector<std::complex<double>> SCR(16*nbsq, 0.);

  //
  // Insert integral elements
  //

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

  //
  // Orthonormalize prior to X2C Hamiltonian call
  //

  // Get transformation - because we don't need this after the hamiltonian call
  //  we can store it in scratch
  double* SCR1 = reinterpret_cast<double*>(SCR.data());
  double* SCR2 = SCR1 + nbsq;

  std::copy_n(S.data(), nbsq, SCR1);
  int64_t nbuse = orthonormalize(nb, SCR1, SCR2, 1e-12);

  // Transform all integrals in place
  detail::transform(nb, T.data(), nb, SCR1, nb, SCR2, nb, T.data(), nb, true);
  detail::transform(nb, V.data(), nb, SCR1, nb, SCR2, nb, V.data(), nb, true);
  detail::transform(nb, pVp[0], nb, SCR1, nb, SCR2, nb, pVp[0], nb, true);
  detail::transform(nb, pVp[1], nb, SCR1, nb, SCR2, nb, pVp[1], nb, true);
  detail::transform(nb, pVp[2], nb, SCR1, nb, SCR2, nb, pVp[2], nb, true);
  detail::transform(nb, pVp[3], nb, SCR1, nb, SCR2, nb, pVp[3], nb, true);

  //
  // Do main work
  //

  // Create structs
  Integrals ints{S.data(), T.data(), V.data(), pVp};
  X2COperators out{UL.data(), US.data(), coreX2C.data()};

  x2c_hamiltonian(nb, ints, out, SCR.data());

  //
  // Expected values
  //

  // X2C Core elements
  std::vector<std::complex<double>> expectCore(4*nbsq, 0.);

  expectCore[0] = std::complex<double>(-4.935340376492e+01, 0.000000000000e+00);
  expectCore[3] = std::complex<double>(6.330215103947e+00, 0.000000000000e+00);
  expectCore[5] = std::complex<double>(-2.380561465990e-04, 0.000000000000e+00);
  expectCore[6] = std::complex<double>(0.000000000000e+00, -2.380561465987e-04);
  expectCore[9] = std::complex<double>(-4.365868030995e+01, 0.000000000000e+00);
  expectCore[10] = std::complex<double>(0.000000000000e+00, 2.341816771659e-04);
  expectCore[12] = std::complex<double>(2.380561465989e-04, 0.000000000000e+00);
  expectCore[15] = std::complex<double>(4.666481801232e-05, 0.000000000000e+00);
  expectCore[17] = std::complex<double>(0.000000000000e+00, -2.341816771659e-04);
  expectCore[18] = std::complex<double>(-4.365868030995e+01, 0.000000000000e+00);
  expectCore[20] = std::complex<double>(0.000000000000e+00, 2.380561465969e-04);
  expectCore[23] = std::complex<double>(0.000000000000e+00, 4.666481801378e-05);
  expectCore[24] = std::complex<double>(6.330215103947e+00, 0.000000000000e+00);
  expectCore[27] = std::complex<double>(-4.722814937611e+01, 0.000000000000e+00);
  expectCore[29] = std::complex<double>(-4.666481801667e-05, 0.000000000000e+00);
  expectCore[30] = std::complex<double>(0.000000000000e+00, -4.666481801691e-05);
  expectCore[33] = std::complex<double>(2.380561465989e-04, 0.000000000000e+00);
  expectCore[34] = std::complex<double>(0.000000000000e+00, -2.380561465969e-04);
  expectCore[36] = std::complex<double>(-4.935340376492e+01, 0.000000000000e+00);
  expectCore[39] = std::complex<double>(6.330215103947e+00, 0.000000000000e+00);
  expectCore[40] = std::complex<double>(-2.380561465990e-04, 0.000000000000e+00);
  expectCore[43] = std::complex<double>(-4.666481801667e-05, 0.000000000000e+00);
  expectCore[45] = std::complex<double>(-4.365868030995e+01, 0.000000000000e+00);
  expectCore[46] = std::complex<double>(0.000000000000e+00, -2.341816771289e-04);
  expectCore[48] = std::complex<double>(0.000000000000e+00, 2.380561465987e-04);
  expectCore[51] = std::complex<double>(0.000000000000e+00, 4.666481801691e-05);
  expectCore[53] = std::complex<double>(0.000000000000e+00, 2.341816771289e-04);
  expectCore[54] = std::complex<double>(-4.365868030995e+01, 0.000000000000e+00);
  expectCore[57] = std::complex<double>(4.666481801232e-05, 0.000000000000e+00);
  expectCore[58] = std::complex<double>(0.000000000000e+00, -4.666481801378e-05);
  expectCore[60] = std::complex<double>(6.330215103947e+00, 0.000000000000e+00);
  expectCore[63] = std::complex<double>(-4.722814937611e+01, 0.000000000000e+00);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(coreX2C[i].real(), expectCore[i].real(), 1e-8);
    EXPECT_NEAR(coreX2C[i].imag(), expectCore[i].imag(), 1e-8);
  }

  // X2C Core eigenvalues
  std::vector<double> expectedEig(2*nb, 0.);

  expectedEig[0] = -5.4709561492e+01;
  expectedEig[1] = -5.4709561492e+01;
  expectedEig[2] = -4.3658914528e+01;
  expectedEig[3] = -4.3658914528e+01;
  expectedEig[4] = -4.3658446128e+01;
  expectedEig[5] = -4.3658446128e+01;
  expectedEig[6] = -4.1871991613e+01;
  expectedEig[7] = -4.1871991613e+01;


  double* eig = new double[2*nbsq];
  std::complex<double>* moreSCR = new std::complex<double>[4*nbsq];

  std::copy_n(coreX2C.data(), 4*nbsq, moreSCR);
  lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, moreSCR, 2*nb, eig);

  for(auto i = 0; i < 2*nb; i++) {
    EXPECT_NEAR(eig[i], expectedEig[i], 1e-8);
  }

  delete[] eig;
  delete[] moreSCR;

}
