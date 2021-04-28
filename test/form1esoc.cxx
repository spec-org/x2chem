#include <gtest/gtest.h>

#include <x2chem.hpp>

using namespace X2Chem;


TEST( Form1eSOC, dense_W ) {

  // Boettger test using the following CQ input

  int64_t nb = 2;
  std::array<double*,4> pVp;
  pVp[0] = new double[nb*nb];
  pVp[1] = new double[nb*nb];
  pVp[2] = new double[nb*nb];
  pVp[3] = new double[nb*nb];
  std::complex<double>* computed_W = new std::complex<double>[4*nb*nb];
  std::complex<double>* expected   = new std::complex<double>[4*nb*nb];

  pVp[0][0] = 1.0;
  pVp[0][1] = 1.0;
  pVp[0][2] = 1.0;
  pVp[0][3] = 1.0;

  pVp[1][0] = 2.0;
  pVp[1][1] = 2.0;
  pVp[1][2] = 2.0;
  pVp[1][3] = 2.0;

  pVp[2][0] = 3.0;
  pVp[2][1] = 3.0;
  pVp[2][2] = 3.0;
  pVp[2][3] = 3.0;

  pVp[3][0] = 4.0;
  pVp[3][1] = 4.0;
  pVp[3][2] = 4.0;
  pVp[3][3] = 4.0;

  // Expected 
  expected[0] = std::complex<double>(1.0,4.0);
  expected[1] = std::complex<double>(1.0,4.0);
  expected[2] = std::complex<double>(-3.0,2.0);
  expected[3] = std::complex<double>(-3.0,2.0);

  expected[4] = std::complex<double>(1.0,4.0);
  expected[5] = std::complex<double>(1.0,4.0);
  expected[6] = std::complex<double>(-3.0,2.0);
  expected[7] = std::complex<double>(-3.0,2.0);

  expected[8] = std::complex<double>(3.0,2.0);
  expected[9] = std::complex<double>(3.0,2.0);
  expected[10] = std::complex<double>(1.0,-4.0);
  expected[11] = std::complex<double>(1.0,-4.0);

  expected[12] = std::complex<double>(3.0,2.0);
  expected[13] = std::complex<double>(3.0,2.0);
  expected[14] = std::complex<double>(1.0,-4.0);
  expected[15] = std::complex<double>(1.0,-4.0);

  // Form test W
  int64_t LDW = 2*nb;
  _form_1e_soc_matrix(nb, computed_W, LDW, pVp, true);

  // Compare SOC matrices
  for( auto i = 0; i < 2*nb; i++ ) {
    for( auto j = 0; j < 2*nb; j++ ) {
      double difference = std::abs(expected[j + i*2*nb] - computed_W[j + i*2*nb]);
      std::cout << "Computed W[" << j << "," << i << "]: " << computed_W[j + i*2*nb] <<"\n";
      std::cout << "Expected W[" << j << "," << i << "]: " << expected[j + i*2*nb] <<"\n";
      EXPECT_LT(difference, 1e-12);
    }
  }

  free(pVp[0]);
  free(pVp[1]);
  free(pVp[2]);
  free(pVp[3]);
  free(computed_W);
  free(expected);
}

TEST( Form1eSOC, dense_W_noSOC ) {

  // Boettger test using the following CQ input

  int64_t nb = 2;
  std::array<double*,4> pVp;
  pVp[0] = new double[nb*nb];
  pVp[1] = new double[nb*nb];
  pVp[2] = new double[nb*nb];
  pVp[3] = new double[nb*nb];
  std::complex<double>* computed_W = new std::complex<double>[4*nb*nb];
  std::complex<double>* expected   = new std::complex<double>[4*nb*nb];

  pVp[0][0] = 1.0;
  pVp[0][1] = 1.0;
  pVp[0][2] = 1.0;
  pVp[0][3] = 1.0;

  pVp[1][0] = 2.0;
  pVp[1][1] = 2.0;
  pVp[1][2] = 2.0;
  pVp[1][3] = 2.0;

  pVp[2][0] = 3.0;
  pVp[2][1] = 3.0;
  pVp[2][2] = 3.0;
  pVp[2][3] = 3.0;

  pVp[3][0] = 4.0;
  pVp[3][1] = 4.0;
  pVp[3][2] = 4.0;
  pVp[3][3] = 4.0;

  // Expected 
  expected[0] = std::complex<double>(1.0);
  expected[1] = std::complex<double>(1.0);
  expected[2] = std::complex<double>(0.0);
  expected[3] = std::complex<double>(0.0);

  expected[4] = std::complex<double>(1.0);
  expected[5] = std::complex<double>(1.0);
  expected[6] = std::complex<double>(0.0);
  expected[7] = std::complex<double>(0.0);

  expected[8] = std::complex<double>(0.0);
  expected[9] = std::complex<double>(0.0);
  expected[10] = std::complex<double>(1.0);
  expected[11] = std::complex<double>(1.0);

  expected[12] = std::complex<double>(0.0);
  expected[13] = std::complex<double>(0.0);
  expected[14] = std::complex<double>(1.0);
  expected[15] = std::complex<double>(1.0);

  // Form test W
  int64_t LDW = 2*nb;
  _form_1e_soc_matrix(nb, computed_W, LDW, pVp, false);

  // Compare SOC matrices
  for( auto i = 0; i < 2*nb; i++ ) {
    for( auto j = 0; j < 2*nb; j++ ) {
      double difference = std::abs(expected[j + i*2*nb] - computed_W[j + i*2*nb]);
      std::cout << "Computed W[" << j << "," << i << "]: " << computed_W[j + i*2*nb] <<"\n";
      std::cout << "Expected W[" << j << "," << i << "]: " << expected[j + i*2*nb] <<"\n";
      EXPECT_LT(difference, 1e-12);
    }
  }

  free(pVp[0]);
  free(pVp[1]);
  free(pVp[2]);
  free(pVp[3]);
  free(computed_W);
  free(expected);
}


