#include <algorithm>
#include <gtest/gtest.h>
#include <blas.hh>
#include <x2chem.hpp>


using namespace X2Chem;

template <typename T>
void check_3x3_ortho(int64_t nVec, T* in, T* expected) {

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Transformation matrix" << std::endl;
  detail::print_matrix(3, in);

  std::cout << "\nExpected matrix" << std::endl;
  detail::print_matrix(3, in);

  EXPECT_EQ(nVec, 3);

  for( auto i = 0; i < 9; i++ ) {
    // We don't care if there is a sign difference
    double difference = std::min(std::abs(in[i] - expected[i]),
                                 std::abs(in[i] + expected[i]));
    EXPECT_LT(difference, 1e-12);
  }

  std::cout << "-----------------------------------\n" << std::endl;

}

TEST( Orthonormal, I_double ) {

  double in[9];
  double scr[3];
  double expected[9];

  std::fill_n(in, 9, 0.);
  std::fill_n(expected, 9, 0.);

  in[0] = 1.;
  in[4] = 1.;
  in[8] = 1.;

  expected[0] = 1.;
  expected[4] = 1.;
  expected[8] = 1.;

  auto nVec = orthonormalize(3, in, scr, 0.);
  check_3x3_ortho(nVec, in, expected);

}

TEST( Orthonormal, I_dcomplex ) {

  std::complex<double> in[9];
  double scr[3];
  std::complex<double> expected[9];

  std::fill_n(in, 9, std::complex<double>(0.));
  std::fill_n(expected, 9, std::complex<double>(0.));

  in[0] = 1.;
  in[4] = 1.;
  in[8] = 1.;

  expected[0] = 1.;
  expected[4] = 1.;
  expected[8] = 1.;

  auto nVec = orthonormalize(3, in, scr, 0.);
  check_3x3_ortho(nVec, in, expected);

}

TEST( Orthonormal, Scale_double ) {
  
  double in[9];
  double scr[3];
  double expected[9];

  std::fill_n(in, 9, 0.);
  std::fill_n(expected, 9, 0.);

  // Regular scaling
  in[0] = 7.;
  in[4] = 3.;
  in[8] = 2.;

  expected[0] = 1./std::sqrt(7.);
  expected[4] = 1./std::sqrt(3.);
  expected[8] = 1./std::sqrt(2.);

  auto nVec = orthonormalize(3, in, scr, 0.);
  check_3x3_ortho(nVec, in, expected);

  // Reordered scaling
  std::fill_n(in, 9, 0.);
  std::fill_n(expected, 9, 0.);
  in[0] = 2.;
  in[4] = 3.;
  in[8] = 7.;

  expected[2] = 1./std::sqrt(7.);
  expected[4] = 1./std::sqrt(3.);
  expected[6] = 1./std::sqrt(2.);

  nVec = orthonormalize(3, in, scr, 0.);
  check_3x3_ortho(nVec, in, expected);

}

TEST( Orthonormal, Scale_dcomplex ) {

  std::complex<double> in[9];
  double scr[3];
  std::complex<double> expected[9];

  std::fill_n(in, 9, 0.);
  std::fill_n(expected, 9, 0.);

  // Regular scaling
  in[0] = 7.;
  in[4] = std::complex<double>(0.,3.);
  in[8] = 2.;

  expected[0] = 1./std::sqrt(7.);
  expected[4] = std::complex<double>(0.,1./std::sqrt(3.));
  expected[8] = 1./std::sqrt(2.);

  auto nVec = orthonormalize(3, in, scr, 0.);
  check_3x3_ortho(nVec, in, expected);

  // Reordered scaling
  std::fill_n(in, 9, 0.);
  std::fill_n(expected, 9, 0.);
  in[0] = std::complex<double>(0.,2.);
  in[4] = std::complex<double>(0.,3.);
  in[8] = 7.;

  expected[2] = 1./std::sqrt(7.);
  expected[4] = std::complex<double>(0.,1./std::sqrt(3.));
  expected[6] = std::complex<double>(0.,1./std::sqrt(2.));

  nVec = orthonormalize(3, in, scr, 0.);
  check_3x3_ortho(nVec, in, expected);

}

TEST( Orthonormal, Throws ) {

  double in[9];
  double scr[3];

  std::fill_n(in, 9, std::nan(""));
  EXPECT_THROW(orthonormalize(3, in, scr, 0.), LinAlgExcept);

  try {
    orthonormalize(3, in, scr, 0.);
  }
  catch(LinAlgExcept& e) {
    std::cout << "Code: " << e.exit_code() << std::endl;
    std::cout << "Name: " << e.function_name() << std::endl;
    std::cout << "What: " << e.what() << std::endl;
  }

}

TEST( Orthonormal, Trim ) {

  double in[9];
  double scr[3];

  std::fill_n(in, 9, 0.);
  in[0] = 1.;
  in[4] = 1.;
  in[5] = 1.;
  in[7] = 1.;
  in[8] = 1.;

  auto nVec = orthonormalize(3, in, scr, 1e-12);

  std::cout << "nVec: " << nVec << std::endl;
  detail::print_matrix(3, in);

  EXPECT_EQ(nVec, 2);

}

TEST( Orthonormal, Orthonormalize ) {

  double in[25];
  double scr[5];
  double copy[25];

  std::fill_n(in, 25, 0.);
  
  // set diagonal to 1,2,3,4,5
  for( auto i = 0; i < 5; i++ )
    in[i*6] = double(i+1);

  // Add some off diagonal "noise"
  in[1] = 0.5;
  in[5] = 0.5;

  in[4] = 0.75;
  in[20] = 0.75;

  // Make vectors 2 and 3 linearly dependent
  in[18] = 1./3.;
  in[13] = 1.;
  in[17] = 1.;

  in[2] = 1.5;
  in[10] = 1.5;
  in[7] = 0.75;
  in[11] = 0.75;

  in[3] = 0.5;
  in[15] = 0.5;
  in[8] = 0.25;
  in[16] = 0.25;

  // Create copy (since orthonormalize destroys it)
  std::copy_n(in, 25, copy);

  std::cout << "Original matrix:" << std::endl;
  detail::print_matrix(5, in);

  auto nVec = orthonormalize(5, in, scr, 1e-12);

  std::cout << "\nnVec: " << nVec << std::endl;
  std::cout << "Transformation matrix:" << std::endl;
  detail::print_matrix(5, in);

  // Check that linearly dependent vectors are trimmed
  EXPECT_EQ(nVec, 4);

  double out[16];
  double tmp[20];

  // Verify that U\dag S U is identity
  blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
             4, 5, 5, 1.0, in, 5, copy, 5, 0., tmp, 4);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
             4, 4, 5, 1.0, tmp, 4, in, 5, 0., out, 4);

  for( auto i = 0; i < 16; i++ ) {
    EXPECT_NEAR(out[i], i % 5 == 0 ? 1. : 0., 1e-12);
  }

}
