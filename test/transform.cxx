#include <algorithm>

#include <gtest/gtest.h>

#include <x2chem.hpp>
#include <x2chem/detail.hpp>


using namespace X2Chem;

TEST( Transform, Identity ) {
  int64_t nb = 3;
  int64_t nbsq = nb*nb;

  std::vector<double> A_storage(nbsq, 0.);
  std::vector<double> B_storage(nbsq, 0.);
  std::vector<double> U_storage(nbsq, 0.);
  std::vector<double> SCR_storage(nbsq, 0.);

  double* A = A_storage.data();
  double* B = B_storage.data();
  double* U = U_storage.data();
  double* SCR = SCR_storage.data();

  for(auto i = 0; i < nbsq; i++)
    A[i] = i;

  U[0] = 1.;
  U[4] = 1.;
  U[8] = 1.;

  detail::transform(nb, nb, A, nb, U, nb, SCR, nb, B, nb, true);

  for(auto i = 0; i < nbsq; i++)
    EXPECT_NEAR(A[i], B[i], 1e-14);
}

TEST( Transform, Scale ) {
  int64_t nb = 3;
  int64_t nbsq = nb*nb;

  std::vector<double> A_storage(nbsq, 0.);
  std::vector<double> B_storage(nbsq, 0.);
  std::vector<double> U_storage(nbsq, 0.);
  std::vector<double> SCR_storage(nbsq, 0.);

  double* A = A_storage.data();
  double* B = B_storage.data();
  double* U = U_storage.data();
  double* SCR = SCR_storage.data();

  for(auto i = 0; i < nbsq; i++)
    A[i] = i;

  U[0] = 7.;
  U[4] = 7.;
  U[8] = 7.;

  detail::transform(nb, nb, A, nb, U, nb, SCR, nb, B, nb, true);

  for(auto i = 0; i < nbsq; i++)
    EXPECT_NEAR(A[i], B[i]/49., 1e-14);
}

TEST( Transform, RectanglarScale ) {
  int64_t nb = 3;
  int64_t nbsq = nb*nb;
  int64_t nbu = 2;
  int64_t nbusq = nbu*nbu;
  

  std::vector<double> A_storage(nbsq, 0.);
  std::vector<double> B_storage(nbsq, 0.);
  std::vector<double> U_storage(nbu*nb, 0.);
  std::vector<double> SCR_storage(nbu*nb, 0.);

  double* A = A_storage.data();
  double* B = B_storage.data();
  double* U = U_storage.data();
  double* SCR = SCR_storage.data();

  for(auto i = 0; i < nbsq; i++)
    A[i] = i;

  std::copy_n(A, nbsq, B);

  U[0] = 7.;
  U[4] = 7.;

  detail::transform(nb, nbu, A, nb, U, nb, SCR, nb, B, nb, true);

  detail::print_matrix(nb, A);
  detail::print_matrix(nb, B);

  // Check that the transformed subblock is correct
  for(auto i = 0; i < nbu; i++)
    for(auto j = 0; j < nbu; j++)
      EXPECT_NEAR(A[j+i*nb], B[j+i*nb]/49., 1e-14);

  // Check that the remaining matrix is unchanged
  for(auto i = 0; i < nb; i++) {
    // Bottom row
    EXPECT_NEAR(A[(i+1)*nb-1], B[(i+1)*nb-1], 1e-14);
    // Left-most column
    EXPECT_NEAR(A[(nb-1)*nb+i], B[(nb-1)*nb+i], 1e-14);
  }
}

TEST( Transform, Orthonormalize ) {
  int64_t nb = 3;
  int64_t nbsq = nb*nb;

  std::vector<double> A_storage(nbsq, 0.);
  std::vector<double> B_storage(nbsq, 0.);
  std::vector<double> U_storage(nbsq, 0.);
  std::vector<double> SCR_storage(nbsq, 0.);
  
  double* A = A_storage.data();
  double* B = B_storage.data();
  double* U = U_storage.data();
  double* SCR = SCR_storage.data();

  A[0] = 1.;
  A[1] = 0.5;
  A[2] = 0.75;
  A[3] = 0.5;
  A[4] = 2.;
  A[5] = 0.3;
  A[6] = 0.75;
  A[7] = 0.3;
  A[8] = 1.;
  
  U[0] = -0.3031622603828735;
  U[1] = -0.4980067202729349;
  U[2] = -0.2580907179231789;
  U[3] = 0.43297016027207674;
  U[4] = -0.541229326159526;
  U[5] = 0.5357636659020664;
  U[6] = -1.4945860225929999;
  U[7] = 0.1863277161474962;
  U[8] = 1.396058041260623;

  detail::transform(nb, nb, A, nb, U, nb, SCR, nb, B, nb, true);

  detail::print_matrix(nb, B);
  
  for(auto i = 0; i < nbsq; i++) {
    EXPECT_NEAR(B[i], i%(nb+1) == 0 ? 1. : 0., 1e-14);
  }
}

TEST( Transform, Unorthonormalize ) {
  int64_t nb = 4;
  int64_t nbu = 3;
  int64_t nbsq = nb*nb;
  int64_t nbusq = nbu*nbu;
  int64_t nbbu = nb*nbu;

  std::vector<double> A_storage(nbusq, 0.);
  std::vector<double> B_storage(nbsq, 0.);
  std::vector<double> U_storage(nbbu, 0.);
  std::vector<double> SCR_storage(nbbu, 0.);

  double* A = A_storage.data();
  double* B = B_storage.data();
  double* U = U_storage.data();
  double* SCR = SCR_storage.data();

  // A is 3x3 identity
  std::fill_n(A, 9, 0.);
  A[0] = 1.;
  A[4] = 1.;
  A[8] = 1.;

  // U is U.I.UT = original mat
  U[0] = -0.5114592331735963;
  U[1] = -1.3899267933208561;
  U[2] = -1.3899267933208561;
  U[3] = -0.37296819622002875;
  U[4] = 0.782505774870358;
  U[5] = -0.26015072937258954;
  U[6] = -0.26015072937258954;
  U[7] = 0.8659213783317278;
  U[8] = 0.3550974022659386;
  U[9] = -0.02061812832708612;
  U[10] = -0.02061812832708612;
  U[11] = -0.33327899897302093;

  detail::transform(nbu, nb, A, nbu, U, nb, SCR, nbu, B, nb, false);

  std::vector<double> expected = {
    1., 0.5, 0.5, 0.75,
    0.5, 2., 2., 0.3,
    0.5, 2., 2., 0.3,
    0.75, 0.3, 0.3, 1.};

  for(auto i = 0; i < nbsq; i++)
    EXPECT_NEAR(B[i], expected[i], 1e-14);

}
