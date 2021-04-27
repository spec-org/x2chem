#include <gtest/gtest.h>

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

  std::fill_n(in, 9, 0.);
  in[0] = 1.;
  in[4] = 1.;
  in[7] = 1.;
  in[8] = 1e-16;

  auto nVec = orthonormalize(3, in, scr, 1e-12);
  std::cout << "nVec: " << nVec;

  detail::print_matrix(3, in);
}

TEST( Orthonormal, Trim ) {


}
