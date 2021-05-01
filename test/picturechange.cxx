#include <algorithm>
#include <gtest/gtest.h>
#include <lapack.hh>
#include <x2chem.hpp>
#include <x2chem/detail.hpp>


using namespace X2Chem;

TEST( Picture_Change, UH_91Plus ) {

  int64_t nb = 4;
  int64_t nbsq = nb*nb;

  //
  // Allocate memory
  //

  // Allocate required matrices
  double* S = new double[nbsq];
  std::complex<double>* UL = new std::complex<double>[4*nb*nb];
  std::complex<double>* US = new std::complex<double>[4*nb*nb];
  std::complex<double>* X = new std::complex<double>[4*nb*nb];
  std::complex<double>* R = new std::complex<double>[4*nb*nb];
  double* p = new double[nb];
  double* K = new double[nbsq];
  
  // Allocate Scratch
  double* SCR1 = new double[nbsq];
  double* SCR2 = new double[nbsq];
  std::complex<double>* CSCR1 = new std::complex<double>[4*nb*nb];
  std::complex<double>* CSCR2 = new std::complex<double>[4*nb*nb];

  // Overwite required matrices with zeros
  std::fill_n(S,nb*nb,0.0); 
  std::fill_n(K,nb*nb,0.0); 
  std::fill_n(p,nb,0.0); 
  std::fill_n(X,4*nb*nb,std::complex<double>(0.0)); 
  std::fill_n(R,4*nb*nb,std::complex<double>(0.0)); 


  //
  // Data
  //

  // Overlap matrix
  S[0] = 1.0;
  S[3] = -5.90215761e-01;
  S[5] = 1.0;
  S[10] = 1.0;
  S[12] = -5.90215761e-01;
  S[15] = 1.0;


  // Intermediate Matrices
  K[0] = -8.05896581568e-01;
  K[3] = -5.92056331625e-01;
  K[5] = -3.56253103175e-01;
  K[6] = 9.34389493990e-01;
  K[9] = -9.34389493990e-01;
  K[10] = -3.56253103175e-01;
  K[12] = -5.92056331625e-01;
  K[15] = 8.05896581568e-01;

  X[0] = std::complex<double>(3.59252231497e-03,-9.21851146041e-25);
  X[3] = std::complex<double>(-3.62517127147e-08,-3.30564911701e-23);
  X[5] = std::complex<double>(-5.70579015790e-07,1.49652882486e-06);
  X[6] = std::complex<double>(-1.49652882486e-06,-5.70579015790e-07);
  X[9] = std::complex<double>(3.64777425954e-03,-8.22367347830e-17);
  X[10] = std::complex<double>(1.12757025938e-17,1.70890626348e-06);
  X[12] = std::complex<double>(5.79362276429e-07,-1.51956577881e-06);
  X[15] = std::complex<double>(1.85230523985e-07,-4.85827222118e-07);
  X[17] = std::complex<double>(3.93565388612e-17,-1.70890626341e-06);
  X[18] = std::complex<double>(3.64777425954e-03,6.25584653524e-17);
  X[20] = std::complex<double>(1.51956577885e-06,5.79362276343e-07);
  X[23] = std::complex<double>(4.85827222080e-07,1.85230523856e-07);
  X[24] = std::complex<double>(-3.78161437513e-07,-2.11986924011e-22);
  X[27] = std::complex<double>(5.28983263983e-03,-1.41434896858e-19);
  X[29] = std::complex<double>(-2.68762245270e-07,7.04916297220e-07);
  X[30] = std::complex<double>(-7.04916297218e-07,-2.68762245271e-07);
  X[33] = std::complex<double>(5.70579015793e-07,1.49652882489e-06);
  X[34] = std::complex<double>(1.49652882489e-06,-5.70579015814e-07);
  X[36] = std::complex<double>(3.59252231497e-03,5.67037736210e-17);
  X[39] = std::complex<double>(-3.62517127469e-08,-4.35849273339e-17);
  X[40] = std::complex<double>(-5.79362276303e-07,-1.51956577886e-06);
  X[43] = std::complex<double>(-1.85230523878e-07,-4.85827222109e-07);
  X[45] = std::complex<double>(3.64777425954e-03,-8.89045781438e-18);
  X[46] = std::complex<double>(1.66967134563e-17,-1.70890626351e-06);
  X[48] = std::complex<double>(-1.51956577884e-06,5.79362276305e-07);
  X[51] = std::complex<double>(-4.85827222090e-07,1.85230523883e-07);
  X[53] = std::complex<double>(9.10729824888e-18,1.70890626341e-06);
  X[54] = std::complex<double>(3.64777425954e-03,1.73472347598e-17);
  X[57] = std::complex<double>(2.68762245317e-07,7.04916297199e-07);
  X[58] = std::complex<double>(7.04916297285e-07,-2.68762245291e-07);
  X[60] = std::complex<double>(-3.78161437574e-07,4.66206934169e-17);
  X[63] = std::complex<double>(5.28983263983e-03,3.46944695195e-17);

  R[0] = std::complex<double>(9.99993546952e-01,0.00000000000e+00);
  R[3] = std::complex<double>(7.73927484993e-10,-1.74716412121e-24);
  R[5] = std::complex<double>(2.08082590640e-09,-5.45764190767e-09);
  R[6] = std::complex<double>(5.45764191148e-09,2.08082591099e-09);
  R[9] = std::complex<double>(9.99993346935e-01,0.00000000000e+00);
  R[10] = std::complex<double>(0.00000000000e+00,-6.23212209616e-09);
  R[12] = std::complex<double>(-2.08082612078e-09,5.45764201106e-09);
  R[15] = std::complex<double>(-9.79741797240e-10,2.56969096819e-09);
  R[17] = std::complex<double>(0.00000000000e+00,6.23212209616e-09);
  R[18] = std::complex<double>(9.99993346935e-01,0.00000000000e+00);
  R[20] = std::complex<double>(-5.45764196423e-09,-2.08082556394e-09);
  R[23] = std::complex<double>(-2.56969089895e-09,-9.79741742515e-10);
  R[24] = std::complex<double>(7.73927484993e-10,1.74716412121e-24);
  R[27] = std::complex<double>(9.99986009128e-01,0.00000000000e+00);
  R[29] = std::complex<double>(9.79741711645e-10,-2.56969096367e-09);
  R[30] = std::complex<double>(2.56969085793e-09,9.79741710639e-10);
  R[33] = std::complex<double>(-2.08082612078e-09,-5.45764201106e-09);
  R[34] = std::complex<double>(-5.45764196423e-09,2.08082556394e-09);
  R[36] = std::complex<double>(9.99993546952e-01,0.00000000000e+00);
  R[39] = std::complex<double>(7.73927504175e-10,-3.25802752832e-17);
  R[40] = std::complex<double>(2.08082590640e-09,5.45764190767e-09);
  R[43] = std::complex<double>(9.79741711645e-10,2.56969096367e-09);
  R[45] = std::complex<double>(9.99993346935e-01,0.00000000000e+00);
  R[46] = std::complex<double>(1.66533453694e-16,6.23212198514e-09);
  R[48] = std::complex<double>(5.45764191148e-09,-2.08082591099e-09);
  R[51] = std::complex<double>(2.56969085793e-09,-9.79741710639e-10);
  R[53] = std::complex<double>(1.66533453694e-16,-6.23212198514e-09);
  R[54] = std::complex<double>(9.99993346935e-01,0.00000000000e+00);
  R[57] = std::complex<double>(-9.79741797240e-10,-2.56969096819e-09);
  R[58] = std::complex<double>(-2.56969089895e-09,9.79741742515e-10);
  R[60] = std::complex<double>(7.73927504175e-10,3.25802752832e-17);
  R[63] = std::complex<double>(9.99986009128e-01,0.00000000000e+00);

  p[0] = 9.84825628954e-01;
  p[1] = 1.00000000000e+00;
  p[2] = 1.00000000000e+00;
  p[3] = 1.44966370909e+00;


  // Orthonormalize overlap to get the ortho matrix
  std::copy_n(S, nbsq, SCR1);
  orthonormalize(nb, SCR1, SCR2, 1e-12);
 
  // Rename SCR1 to ortho
  double* ortho = SCR1;


  //
  // Do main work
  //

  // Run picture change routine
  detail::form_picture_change(nb, UL, US, K, X, R, p, CSCR1);

  // Back transform UL to AO basis
  detail::transform(nb, nb, UL, 2*nb, ortho, nb, CSCR2, nb, CSCR1, 2*nb, false);
  detail::transform(nb, nb, UL + nb, 2*nb, ortho, nb, CSCR2, nb, CSCR1 + nb, 2*nb, false);
  detail::transform(nb, nb, UL + 2*nb*nb, 2*nb, ortho, nb, CSCR2, nb, CSCR1 + 2*nb*nb, 2*nb, false);
  detail::transform(nb, nb, UL + 2*nb*nb + nb, 2*nb, ortho, nb, CSCR2, nb, CSCR1 + 2*nb*nb + nb, 2*nb, false);

  // Multiply 1C overlap on right of UL for each block
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1, 2*nb, S, nb, 0.0, UL, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1 + nb, 2*nb, S, nb, 0.0, UL + nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1 + 2*nb*nb, 2*nb, S, nb, 0.0, UL + 2*nb*nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1 + 2*nb*nb +nb, 2*nb, S, nb, 0.0, UL + 2*nb*nb +nb, 2*nb);

  // Back transform US to AO basis
  detail::transform(nb, nb, US, 2*nb, ortho, nb, CSCR2, nb, CSCR1, 2*nb, false);
  detail::transform(nb, nb, US + nb, 2*nb, ortho, nb, CSCR2, nb, CSCR1 + nb, 2*nb, false);
  detail::transform(nb, nb, US + 2*nb*nb, 2*nb, ortho, nb, CSCR2, nb, CSCR1 + 2*nb*nb, 2*nb, false);
  detail::transform(nb, nb, US + 2*nb*nb + nb, 2*nb, ortho, nb, CSCR2, nb, CSCR1 + 2*nb*nb + nb, 2*nb, false);

  // Multiply 1C overlap on right of US for each block
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1, 2*nb, S, nb, 0.0, US, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1 + nb, 2*nb, S, nb, 0.0, US + nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1 + 2*nb*nb, 2*nb, S, nb, 0.0, US + 2*nb*nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, CSCR1 + 2*nb*nb +nb, 2*nb, S, nb, 0.0, US + 2*nb*nb +nb, 2*nb);

  //
  // Expected values
  //

  // UL and US elements
  std::complex<double>* expectUL = new std::complex<double>[4*nb*nb];
  std::complex<double>* expectUS = new std::complex<double>[4*nb*nb];
  
  std::fill_n(expectUL,4*nb*nb,std::complex<double>(0.0)); 
  std::fill_n(expectUS,4*nb*nb,std::complex<double>(0.0)); 

  expectUL[0] = std::complex<double>(9.99985323e-01,1.16959699e-16);
  expectUL[3] = std::complex<double>(-3.75686529e-06,1.98164310e-16);
  expectUL[5] = std::complex<double>(-5.08706717e-09,-1.90715497e-16);
  expectUL[6] = std::complex<double>(-5.45736219e-17,-5.08706712e-09);
  expectUL[9] = std::complex<double>(9.99993347e-01,0.00000000e+00);
  expectUL[10] = std::complex<double>(-5.55111512e-17,-6.23212232e-09);
  expectUL[12] = std::complex<double>(2.18075291e-09,2.48116579e-17);
  expectUL[15] = std::complex<double>(-4.92415551e-09,1.13322914e-17);
  expectUL[17] = std::complex<double>(-5.55111512e-17,6.23212232e-09);
  expectUL[18] = std::complex<double>(9.99993347e-01,0.00000000e+00);
  expectUL[20] = std::complex<double>(2.48120742e-17,2.18075285e-09);
  expectUL[23] = std::complex<double>(1.13324748e-17,-4.92415553e-09);
  expectUL[24] = std::complex<double>(1.50204187e-06,-1.98164310e-16);
  expectUL[27] = std::complex<double>(9.99994233e-01,-1.16959699e-16);
  expectUL[29] = std::complex<double>(6.21127043e-09,9.76957673e-17);
  expectUL[30] = std::complex<double>(2.98751170e-17,6.21127023e-09);
  expectUL[33] = std::complex<double>(5.08706710e-09,-1.81231609e-17);
  expectUL[34] = std::complex<double>(1.81234690e-17,-5.08706705e-09);
  expectUL[36] = std::complex<double>(9.99985323e-01,-1.52963820e-27);
  expectUL[39] = std::complex<double>(-3.75686529e-06,-2.59165936e-27);
  expectUL[40] = std::complex<double>(-2.18075285e-09,2.04181483e-16);
  expectUL[43] = std::complex<double>(4.92415573e-09,2.28153618e-17);
  expectUL[45] = std::complex<double>(9.99993347e-01,0.00000000e+00);
  expectUL[46] = std::complex<double>(-2.49800181e-16,6.23212266e-09);
  expectUL[48] = std::complex<double>(-5.66885909e-17,2.18075295e-09);
  expectUL[51] = std::complex<double>(-3.58338282e-18,-4.92415547e-09);
  expectUL[53] = std::complex<double>(-2.49800181e-16,-6.23212266e-09);
  expectUL[54] = std::complex<double>(9.99993347e-01,0.00000000e+00);
  expectUL[57] = std::complex<double>(-6.21127025e-09,3.31194013e-18);
  expectUL[58] = std::complex<double>(-3.31200243e-18,6.21127023e-09);
  expectUL[60] = std::complex<double>(1.50204187e-06,2.59165936e-27);
  expectUL[63] = std::complex<double>(9.99994233e-01,1.52963820e-27);

  expectUS[0] = std::complex<double>(1.00011931e+00,-1.10392934e-14);
  expectUS[3] = std::complex<double>(2.52772346e-04,2.89430985e-15);
  expectUS[5] = std::complex<double>(3.82373043e-04,5.81306680e-16);
  expectUS[6] = std::complex<double>(1.31871662e-14,3.82373043e-04);
  expectUS[9] = std::complex<double>(9.99746129e-01,1.33287112e-14);
  expectUS[10] = std::complex<double>(-1.34139058e-14,4.68354010e-04);
  expectUS[12] = std::complex<double>(-2.86731500e-05,-1.84490503e-17);
  expectUS[15] = std::complex<double>(4.45622281e-04,-8.85585078e-18);
  expectUS[17] = std::complex<double>(3.33383672e-15,-4.68354010e-04);
  expectUS[18] = std::complex<double>(9.99746129e-01,3.42669213e-15);
  expectUS[20] = std::complex<double>(-1.84431714e-17,-2.86731500e-05);
  expectUS[23] = std::complex<double>(-8.85317235e-18,4.45622281e-04);
  expectUS[24] = std::complex<double>(-5.54639035e-05,-1.16875911e-14);
  expectUS[27] = std::complex<double>(9.99731629e-01,-1.54022001e-15);
  expectUS[29] = std::complex<double>(-4.66798108e-04,-1.23854426e-14);
  expectUS[30] = std::complex<double>(-1.76337170e-15,-4.66798108e-04);
  expectUS[33] = std::complex<double>(-3.82373043e-04,-1.56517037e-17);
  expectUS[34] = std::complex<double>(1.54163930e-17,3.82373043e-04);
  expectUS[36] = std::complex<double>(1.00011931e+00,-5.94492278e-22);
  expectUS[39] = std::complex<double>(2.52772346e-04,-2.71384782e-22);
  expectUS[40] = std::complex<double>(2.86731500e-05,-1.81146638e-14);
  expectUS[43] = std::complex<double>(-4.45622281e-04,-3.66613131e-15);
  expectUS[45] = std::complex<double>(9.99746129e-01,-2.92096468e-14);
  expectUS[46] = std::complex<double>(-2.02775459e-14,-4.68354010e-04);
  expectUS[48] = std::complex<double>(-1.84304366e-14,-2.86731500e-05);
  expectUS[51] = std::complex<double>(-4.28812213e-15,4.45622281e-04);
  expectUS[53] = std::complex<double>(-4.76640035e-15,4.68354010e-04);
  expectUS[54] = std::complex<double>(9.99746129e-01,-4.01139381e-15);
  expectUS[57] = std::complex<double>(4.66798108e-04,-8.20267022e-19);
  expectUS[58] = std::complex<double>(1.77330553e-18,-4.66798108e-04);
  expectUS[60] = std::complex<double>(-5.54639035e-05,1.11081893e-22);
  expectUS[63] = std::complex<double>(9.99731629e-01,5.07067549e-23);



  // Check values for UL and US
  
  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(UL[i].real(), expectUL[i].real(), 1e-9);
    EXPECT_NEAR(UL[i].imag(), expectUL[i].imag(), 1e-9);
  }

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(US[i].real(), expectUS[i].real(), 1e-9);
    EXPECT_NEAR(US[i].imag(), expectUS[i].imag(), 1e-9);
  }

  // Free memory
  delete[] S;
  delete[] UL;
  delete[] US;
  delete[] X; 
  delete[] R; 
  delete[] p; 
  delete[] K; 
  delete[] SCR1; 
  delete[] SCR2; 
  delete[] CSCR1; 
  delete[] CSCR2; 
  delete[] expectUL;
  delete[] expectUS;
  
}
