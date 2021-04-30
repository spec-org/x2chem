#include <algorithm>
#include <gtest/gtest.h>
#include <lapack.hh>
#include <x2chem.hpp>


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
}

std::vector<std::complex<double>> getExpectedOrtho(int64_t size) {
  std::vector<std::complex<double>> expectCore(size, 0.);

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

  return expectCore;

}

std::vector<std::complex<double>> getExpectedAO(int64_t size) {
  std::vector<std::complex<double>> expectCore(size, 0.);

  expectCore[0] = std::complex<double>(-5.40279937735942e+01,0.00000000000000e+00);
  expectCore[3] = std::complex<double>(2.95646046372509e+01,0.00000000000000e+00);
  expectCore[5] = std::complex<double>(1.91149026587235e-04,0.00000000000000e+00);
  expectCore[6] = std::complex<double>(0.00000000000000e+00,1.91149026586882e-04);
  expectCore[9] = std::complex<double>(-4.36586803099474e+01,0.00000000000000e+00);
  expectCore[10] = std::complex<double>(2.13162820728030e-14,2.34181677165934e-04);
  expectCore[12] = std::complex<double>(-1.91149026589130e-04,0.00000000000000e+00);
  expectCore[15] = std::complex<double>(2.33394696335650e-04,-1.04504163654950e-14);
  expectCore[17] = std::complex<double>(1.77635683940025e-14,-2.34181677165934e-04);
  expectCore[18] = std::complex<double>(-4.36586803099474e+01,0.00000000000000e+00);
  expectCore[20] = std::complex<double>(0.00000000000000e+00,-1.91149026586712e-04);
  expectCore[23] = std::complex<double>(0.00000000000000e+00,2.33394696334558e-04);
  expectCore[24] = std::complex<double>(2.95646046372509e+01,0.00000000000000e+00);
  expectCore[27] = std::complex<double>(-4.38079180038424e+01,0.00000000000000e+00);
  expectCore[29] = std::complex<double>(-2.33394696337697e-04,0.00000000000000e+00);
  expectCore[30] = std::complex<double>(0.00000000000000e+00,-2.33394696337559e-04);
  expectCore[33] = std::complex<double>(-1.91149026589130e-04,0.00000000000000e+00);
  expectCore[34] = std::complex<double>(0.00000000000000e+00,1.91149026586712e-04);
  expectCore[36] = std::complex<double>(-5.40279937735941e+01,0.00000000000000e+00);
  expectCore[39] = std::complex<double>(2.95646046372509e+01,0.00000000000000e+00);
  expectCore[40] = std::complex<double>(1.91149026587235e-04,0.00000000000000e+00);
  expectCore[43] = std::complex<double>(-2.33394696337697e-04,0.00000000000000e+00);
  expectCore[45] = std::complex<double>(-4.36586803099474e+01,0.00000000000000e+00);
  expectCore[46] = std::complex<double>(0.00000000000000e+00,-2.34181677128945e-04);
  expectCore[48] = std::complex<double>(0.00000000000000e+00,-1.91149026586882e-04);
  expectCore[51] = std::complex<double>(0.00000000000000e+00,2.33394696337559e-04);
  expectCore[53] = std::complex<double>(-1.42108547152020e-14,2.34181677128945e-04);
  expectCore[54] = std::complex<double>(-4.36586803099473e+01,0.00000000000000e+00);
  expectCore[57] = std::complex<double>(2.33394696335650e-04,1.04504131798162e-14);
  expectCore[58] = std::complex<double>(0.00000000000000e+00,-2.33394696334558e-04);
  expectCore[60] = std::complex<double>(2.95646046372509e+01,0.00000000000000e+00);
  expectCore[63] = std::complex<double>(-4.38079180038425e+01,0.00000000000000e+00);

  return expectCore;
}

std::vector<double> getEigValsOrtho(int64_t size) {

  std::vector<double> expectedEig(size, 0.);

  expectedEig[0] = -5.4709561492e+01;
  expectedEig[1] = -5.4709561492e+01;
  expectedEig[2] = -4.3658914528e+01;
  expectedEig[3] = -4.3658914528e+01;
  expectedEig[4] = -4.3658446128e+01;
  expectedEig[5] = -4.3658446128e+01;
  expectedEig[6] = -4.1871991613e+01;
  expectedEig[7] = -4.1871991613e+01;

  return expectedEig;
}

std::vector<double> getEigValsAO(int64_t size) {
  std::vector<double> expectedEig(size, 0.);

  expectedEig[0] = -7.8920927989e+01;
  expectedEig[1] = -7.8920927989e+01;
  expectedEig[2] = -4.3658914487e+01;
  expectedEig[3] = -4.3658914487e+01;
  expectedEig[4] = -4.3658446128e+01;
  expectedEig[5] = -4.3658446128e+01;
  expectedEig[6] = -1.8914983793e+01;
  expectedEig[7] = -1.8914983793e+01;

  return expectedEig;
}

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
    {pvdp.data(), pvxp_x.data(), pvxp_y.data(), pvxp_z.data()};


  std::vector<std::complex<double>> UL(4*nbsq, 0.);
  std::vector<std::complex<double>> US(4*nbsq, 0.);
  std::vector<std::complex<double>> coreX2C(4*nbsq, 0.);

  std::vector<std::complex<double>> SCR(18*nbsq+nb, 0.);

  //
  // Insert integral elements
  //
  setup_ints(S, T, V, pVp);

  //
  // Orthonormalize prior to X2C Hamiltonian call
  //

  // Get transformation - because we don't need this after the hamiltonian call
  //  we can store it in scratch
  double* SCR1 = reinterpret_cast<double*>(SCR.data());
  double* SCR2 = SCR1 + nbsq;

  std::copy_n(S.data(), nbsq, SCR1);
  orthonormalize(nb, SCR1, SCR2, 1e-12);

  // Transform all integrals in place
  detail::transform(nb, nb, T.data(), nb, SCR1, nb, SCR2, nb, T.data(), nb, true);
  detail::transform(nb, nb, V.data(), nb, SCR1, nb, SCR2, nb, V.data(), nb, true);
  detail::transform(nb, nb, pVp[0], nb, SCR1, nb, SCR2, nb, pVp[0], nb, true);
  detail::transform(nb, nb, pVp[1], nb, SCR1, nb, SCR2, nb, pVp[1], nb, true);
  detail::transform(nb, nb, pVp[2], nb, SCR1, nb, SCR2, nb, pVp[2], nb, true);
  detail::transform(nb, nb, pVp[3], nb, SCR1, nb, SCR2, nb, pVp[3], nb, true);

  //
  // Do main work
  //

  // Create structs
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  Integrals ints{S.data(), T.data(), V.data(), pVp_c};
  X2COperators out{UL.data(), US.data(), coreX2C.data()};

  x2c_hamiltonian(nb, ints, out, SCR.data());

  //
  // Expected values
  //
  std::vector<std::complex<double>> expectCore = getExpectedOrtho(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(coreX2C[i].real(), expectCore[i].real(), 1e-8);
    EXPECT_NEAR(coreX2C[i].imag(), expectCore[i].imag(), 1e-8);
  }

  // X2C Core eigenvalues
  std::vector<double> expectedEig = getEigValsOrtho(2*nb);

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


TEST( X2C_Hamiltonian, UH_91Plus_Mem ) {
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
    {pvdp.data(), pvxp_x.data(), pvxp_y.data(), pvxp_z.data()};


  std::vector<std::complex<double>> UL(4*nbsq, 0.);
  std::vector<std::complex<double>> US(4*nbsq, 0.);
  std::vector<std::complex<double>> coreX2C(4*nbsq, 0.);

  //
  // Insert integral elements
  //
  setup_ints(S, T, V, pVp);

  //
  // Orthonormalize prior to X2C Hamiltonian call
  //

  // Get transformation - because we don't need this after the hamiltonian call
  //  we can store it in scratch
  double* SCR1 = new double[nbsq];
  double* SCR2 = new double[nbsq];

  std::copy_n(S.data(), nbsq, SCR1);
  orthonormalize(nb, SCR1, SCR2, 1e-12);

  // Transform all integrals in place
  detail::transform(nb, nb, T.data(), nb, SCR1, nb, SCR2, nb, T.data(), nb, true);
  detail::transform(nb, nb, V.data(), nb, SCR1, nb, SCR2, nb, V.data(), nb, true);
  detail::transform(nb, nb, pVp[0], nb, SCR1, nb, SCR2, nb, pVp[0], nb, true);
  detail::transform(nb, nb, pVp[1], nb, SCR1, nb, SCR2, nb, pVp[1], nb, true);
  detail::transform(nb, nb, pVp[2], nb, SCR1, nb, SCR2, nb, pVp[2], nb, true);
  detail::transform(nb, nb, pVp[3], nb, SCR1, nb, SCR2, nb, pVp[3], nb, true);

  //
  // Do main work
  //

  // Create structs
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  Integrals ints{S.data(), T.data(), V.data(), pVp_c};
  X2COperators out{UL.data(), US.data(), coreX2C.data()};

  x2c_hamiltonian(nb, ints, out);

  //
  // Expected values
  //
  std::vector<std::complex<double>> expectCore = getExpectedOrtho(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(coreX2C[i].real(), expectCore[i].real(), 1e-8);
    EXPECT_NEAR(coreX2C[i].imag(), expectCore[i].imag(), 1e-8);
  }

  // X2C Core eigenvalues
  std::vector<double> expectedEig = getEigValsOrtho(2*nb);

  double* eig = new double[2*nbsq];
  std::complex<double>* moreSCR = new std::complex<double>[4*nbsq];

  std::copy_n(coreX2C.data(), 4*nbsq, moreSCR);
  lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, moreSCR, 2*nb, eig);

  for(auto i = 0; i < 2*nb; i++) {
    EXPECT_NEAR(eig[i], expectedEig[i], 1e-8);
  }

  delete[] eig;
  delete[] moreSCR;
  delete[] SCR1, SCR2;

}


TEST( X2C_Hamiltonian, UH_91Plus_AO ) {
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
    {pvdp.data(), pvxp_x.data(), pvxp_y.data(), pvxp_z.data()};


  std::vector<std::complex<double>> UL(4*nbsq, 0.);
  std::vector<std::complex<double>> US(4*nbsq, 0.);
  std::vector<std::complex<double>> coreX2C(4*nbsq, 0.);

  std::vector<std::complex<double>> SCR(18*nbsq+nb, 0.);
  std::vector<double> SCR2(7*nbsq, 0.);

  //
  // Insert integral elements
  //
  setup_ints(S, T, V, pVp);

  // Create structs
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  Integrals ints{S.data(), T.data(), V.data(), pVp_c};
  X2COperators out{UL.data(), US.data(), coreX2C.data()};

  x2c_hamiltonian_ao(nb, ints, out, SCR2.data(), SCR.data());

  //
  // Expected values
  //

  std::vector<std::complex<double>> expectCore = getExpectedAO(4*nbsq);
  detail::print_matrix(2*nb, expectCore.data());

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(coreX2C[i].real(), expectCore[i].real(), 1e-8);
    EXPECT_NEAR(coreX2C[i].imag(), expectCore[i].imag(), 1e-8);
  }

  // X2C Core eigenvalues
  std::vector<double> expectedEig = getEigValsAO(2*nb);

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


TEST( X2C_Hamiltonian, UH_91Plus_AO_Mem ) {
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
    {pvdp.data(), pvxp_x.data(), pvxp_y.data(), pvxp_z.data()};


  std::vector<std::complex<double>> UL(4*nbsq, 0.);
  std::vector<std::complex<double>> US(4*nbsq, 0.);
  std::vector<std::complex<double>> coreX2C(4*nbsq, 0.);

  //
  // Insert integral elements
  //
  setup_ints(S, T, V, pVp);

  // Create structs
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  Integrals ints{S.data(), T.data(), V.data(), pVp_c};
  X2COperators out{UL.data(), US.data(), coreX2C.data()};

  x2c_hamiltonian_ao(nb, ints, out);

  //
  // Expected values
  //

  std::vector<std::complex<double>> expectCore = getExpectedAO(4*nbsq);
  detail::print_matrix(2*nb, expectCore.data());

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(coreX2C[i].real(), expectCore[i].real(), 1e-8);
    EXPECT_NEAR(coreX2C[i].imag(), expectCore[i].imag(), 1e-8);
  }

  // X2C Core eigenvalues
  std::vector<double> expectedEig = getEigValsAO(2*nb);

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
