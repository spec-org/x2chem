#include <gtest/gtest.h>
#include <x2chem.hpp>

using namespace X2Chem;


TEST( Boettger, UH_91Plus ) {

  // Boettger test using the following CQ input
  // [Molecule]
  // charge = 91
  // mult = 1
  // geom:
  //  H   0.0  0.0  0.0
  //  U   0.0  0.0  0.9
  //
  // [BASIS]
  // definebasis = on
  // basisdef:
  //  ****
  //  H     0
  //  S    1   1.00
  //        0.6239137298     1.0
  //  ****
  //  U     0
  //  P    1   1.00
  //        0.2     1.0
  //  ****

  // Define test parameters
  int64_t nb = 4;
  double*  nucList = new double[nb];
  int64_t* angList = new int64_t[nb];
  std::complex<double>* core_test = new std::complex<double>[4*nb*nb];
  std::complex<double>* expected  = new std::complex<double>[4*nb*nb];

  // basis function order: S P P P
  nucList[0] = 1.0;
  nucList[1] = 92.0;
  nucList[2] = 92.0;
  nucList[3] = 92.0;

  angList[0] = 0;
  angList[1] = 1;
  angList[2] = 1;
  angList[3] = 1;

  // Raw X2C Core Hamiltonian Output from CQ
  std::fill_n(core_test,4*nb*nb,0.0);
  core_test[0] = -5.40279938e+01;
  core_test[3] =  2.95646047e+01;
  core_test[5] =  1.91149027e-04;
  core_test[6] =  std::complex<double>(0.,1.91149027e-04);

  core_test[9] =  -4.36586803e+01;
  core_test[10] = std::complex<double>(0.,2.34181677e-04);
  core_test[12] = -1.91149027e-04;
  core_test[15] =  2.33394696e-04;

  core_test[17] = std::complex<double>(0.,-2.34181677e-04);
  core_test[18] =  -4.36586803e+01;
  core_test[20] = std::complex<double>(0.,-1.91149027e-04);
  core_test[23] = std::complex<double>(0.,2.33394696e-04);

  core_test[24] =  2.95646047e+01;
  core_test[27] = -4.38079180e+01;
  core_test[29] = -2.33394696e-04;
  core_test[30] = std::complex<double>(0.,-2.33394696e-04);


  core_test[33] =  -1.91149027e-04;
  core_test[34] = std::complex<double>(0.,1.91149027e-04);
  core_test[36] =  -5.40279938e+01;
  core_test[39] =   2.95646047e+01;

  core_test[40] =   1.91149027e-04;
  core_test[43] =  -2.33394696e-04;
  core_test[45] =  -4.36586803e+01;
  core_test[46] = std::complex<double>(0.,-2.34181677e-04);

  core_test[48] = std::complex<double>(0.,-1.91149027e-04);
  core_test[51] = std::complex<double>(0.,2.33394696e-04);
  core_test[53] = std::complex<double>(0.,2.34181677e-04);
  core_test[54] =  -4.36586803e+01;

  core_test[57] =  2.33394696e-04;
  core_test[60] =  2.95646047e+01;
  core_test[63] =  -4.38079180e+01;
  core_test[58] = std::complex<double>(0.,-2.33394696e-04);

  // Expected output after Boettger from CQ
  std::fill_n(expected,4*nb*nb,0.0);
  expected[0] = -5.40279938e+01;
  expected[3] =  2.95646047e+01;
  expected[5] =  1.91149027e-04;
  expected[6] =  std::complex<double>(0.,1.91149027e-04);

  expected[9] =  -4.36586803e+01;
  expected[10] = std::complex<double>(0.,2.29090771e-04);
  expected[12] = -1.91149027e-04;
  expected[15] =  2.28320898e-04;

  expected[17] = std::complex<double>(0.,-2.29090771e-04);
  expected[18] =  -4.36586803e+01;
  expected[20] = std::complex<double>(0.,-1.91149027e-04);
  expected[23] = std::complex<double>(0.,2.28320898e-04);

  expected[24] =  2.95646047e+01;
  expected[27] = -4.38079180e+01;
  expected[29] = -2.28320898e-04;
  expected[30] = std::complex<double>(0.,-2.28320898e-04);

  expected[33] =  -1.91149027e-04;
  expected[34] = std::complex<double>(0.,1.91149027e-04);
  expected[36] =  -5.40279938e+01;
  expected[39] =   2.95646047e+01;

  expected[40] =   1.91149027e-04;
  expected[43] =  -2.28320898e-04;
  expected[45] =  -4.36586803e+01;
  expected[46] = std::complex<double>(0.,-2.29090771e-04);

  expected[48] = std::complex<double>(0.,-1.91149027e-04);
  expected[51] = std::complex<double>(0.,2.28320898e-04);
  expected[53] = std::complex<double>(0.,2.29090771e-04);
  expected[54] =  -4.36586803e+01;

  expected[57] =  2.28320898e-04;
  expected[60] =  2.95646047e+01;
  expected[63] =  -4.38079180e+01;
  expected[58] = std::complex<double>(0.,-2.28320898e-04);


  // Run boettger
  boettger_2e_soc(nb, core_test, nucList, angList);

  // Compare Hamiltonian matrices
  for( auto i = 0; i < 2*nb; i++ ) {
    for( auto j = 0; j < 2*nb; j++ ) {
      double difference = std::abs(expected[j + i*2*nb] - core_test[j + i*2*nb]);
      std::cout << "Computed CoreH[" << j << "," << i << "]: " << core_test[j + i*2*nb] <<"\n";
      std::cout << "Expected CoreH[" << j << "," << i << "]: " << expected[j + i*2*nb] <<"\n";
      EXPECT_LT(difference, 1e-12);
    }
  }

  delete[] nucList;
  delete[] angList;
  delete[] core_test;
  delete[] expected;
}

