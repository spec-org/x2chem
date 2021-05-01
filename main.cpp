#include <algorithm>
#include <iostream>
#include <stdio.h>
#include "x2chem.hpp"
#include <lapack.hh>


#define x2chem_dbgout(name, size, location) \
  std::cout << "-----------------------------------------" << std::endl; \
  std::cout << name << std::endl; \
  X2Chem::detail::print_matrix(size, location); \
  std::cout << "-----------------------------------------" << std::endl;

int main() {

  // Allocate test matrices
  int64_t nb = 7;
  double *S    = new double[nb*nb];
  double *T    = new double[nb*nb];
  double *V    = new double[nb*nb];
  std::array<double*,4> pVp;
  pVp[0] = new double[nb*nb];
  pVp[1] = new double[nb*nb];
  pVp[2] = new double[nb*nb];
  pVp[3] = new double[nb*nb];
  std::complex<double> *UL        = new std::complex<double>[2*nb*2*nb];
  std::complex<double> *US        = new std::complex<double>[2*nb*2*nb];
  std::complex<double> *coreX2C  = new std::complex<double>[2*nb*2*nb];
  std::complex<double> *core4C   = new std::complex<double>[4*nb*4*nb+2*nb*nb+nb];

  double* SCR1 = new double[nb*nb];
  double* SCR2 = new double[nb*nb];


  for (int64_t i = 0; i < 4*nb; i++)
  for (int64_t j = 0; j < 4*nb; j++) {
    core4C[i*4*nb + j] = std::complex<double>(0.);
  }

  // Load UH91+ custom basis 4 primitive (CQ)
  std::fill_n(S,nb*nb,0.0); 
  std::fill_n(T,nb*nb,0.0); 
  std::fill_n(V,nb*nb,0.0); 
  std::fill_n(pVp[0],nb*nb,0.0); 
  std::fill_n(pVp[1],nb*nb,0.0); 
  std::fill_n(pVp[2],nb*nb,0.0); 
  std::fill_n(pVp[3],nb*nb,0.0); 

  S[0] = 1.00000000000000e+00;
  S[1] = 9.99999597260658e-01;
  S[2] = 9.99999996415623e-01;
  S[3] = 9.99999432622440e-01;
  S[6] = -5.90215760754225e-01;
  S[7] = 9.99999597260658e-01;
  S[8] = 1.00000000000000e+00;
  S[9] = 9.99999517687640e-01;
  S[10] = 9.99998073840605e-01;
  S[13] = -5.90431348498698e-01;
  S[14] = 9.99999996415623e-01;
  S[15] = 9.99999517687640e-01;
  S[16] = 1.00000000000000e+00;
  S[17] = 9.99999519231039e-01;
  S[20] = -5.90195399822279e-01;
  S[21] = 9.99999432622440e-01;
  S[22] = 9.99998073840605e-01;
  S[23] = 9.99999519231039e-01;
  S[24] = 1.00000000000000e+00;
  S[27] = -5.89959310190970e-01;
  S[32] = 1.00000000000000e+00;
  S[40] = 1.00000000000000e+00;
  S[42] = -5.90215760754225e-01;
  S[43] = -5.90431348498698e-01;
  S[44] = -5.90195399822279e-01;
  S[45] = -5.89959310190970e-01;
  S[48] = 1.00000000000000e+00;

  T[0] = 9.35870594700001e-01;
  T[1] = 9.35184418533109e-01;
  T[2] = 9.35935289522261e-01;
  T[3] = 9.36684057290625e-01;
  T[6] = -3.68625144860845e-01;
  T[7] = 9.35184418533109e-01;
  T[8] = 9.34500000000000e-01;
  T[9] = 9.35248947474192e-01;
  T[10] = 9.35995793273283e-01;
  T[13] = -3.68656425937958e-01;
  T[14] = 9.35935289522261e-01;
  T[15] = 9.35248947474192e-01;
  T[16] = 9.36000000000000e-01;
  T[17] = 9.36748949159580e-01;
  T[20] = -3.68622170858745e-01;
  T[21] = 9.36684057290625e-01;
  T[22] = 9.35995793273283e-01;
  T[23] = 9.36748949159580e-01;
  T[24] = 9.37500000000001e-01;
  T[27] = -3.68587439655044e-01;
  T[32] = 5.00000000000000e-01;
  T[40] = 5.00000000000000e-01;
  T[42] = -3.68625144860845e-01;
  T[43] = -3.68656425937958e-01;
  T[44] = -3.68622170858745e-01;
  T[45] = -3.68587439655044e-01;
  T[48] = 5.00000000000000e-01;
  
  V[0] = -5.49638951773469e+01;
  V[1] = -5.49622605142595e+01;
  V[2] = -5.49640470402716e+01;
  V[3] = -5.49657762953119e+01;
  V[6] = 2.99331686357291e+01;
  V[7] = -5.49622605142595e+01;
  V[8] = -5.49606664840632e+01;
  V[9] = -5.49624085435321e+01;
  V[10] = -5.49640933950872e+01;
  V[13] = 2.99383351458594e+01;
  V[14] = -5.49640470402716e+01;
  V[15] = -5.49624085435321e+01;
  V[16] = -5.49641992648960e+01;
  V[17] = -5.49659327093469e+01;
  V[20] = 2.99326796394455e+01;
  V[21] = -5.49657762953119e+01;
  V[22] = -5.49640933950872e+01;
  V[23] = -5.49659327093469e+01;
  V[24] = -5.49677146779820e+01;
  V[27] = 2.99269964300843e+01;
  V[32] = -4.41585566940449e+01;
  V[40] = -4.41585566940449e+01;
  V[42] = 2.99331686357291e+01;
  V[43] = 2.99383351458594e+01;
  V[44] = 2.99326796394455e+01;
  V[45] = 2.99269964300843e+01;
  V[48] = -4.43078077293299e+01;

  pVp[0][0] = -1.00133395770494e+02;
  pVp[0][1] = -1.00052789627644e+02;
  pVp[0][2] = -1.00140996221105e+02;
  pVp[0][3] = -1.00228970844626e+02;
  pVp[0][6] = 4.37558148766396e+01;
  pVp[0][7] = -1.00052789627644e+02;
  pVp[0][8] = -9.99723657385543e+01;
  pVp[0][9] = -1.00060372876946e+02;
  pVp[0][10] = -1.00148148190726e+02;
  pVp[0][13] = 4.37527331577816e+01;
  pVp[0][14] = -1.00140996221105e+02;
  pVp[0][15] = -1.00060372876946e+02;
  pVp[0][16] = -1.00148598295196e+02;
  pVp[0][17] = -1.00236591729833e+02;
  pVp[0][20] = 4.37561020733976e+01;
  pVp[0][21] = -1.00228970844626e+02;
  pVp[0][22] = -1.00148148190726e+02;
  pVp[0][23] = -1.00236591729833e+02;
  pVp[0][24] = -1.00324803127441e+02;
  pVp[0][27] = 4.37593838551730e+01;
  pVp[0][32] = -5.29624771935904e+01;
  pVp[0][40] = -5.29624771935904e+01;
  pVp[0][42] = 4.37558148766396e+01;
  pVp[0][43] = 4.37527331577816e+01;
  pVp[0][44] = 4.37561020733976e+01;
  pVp[0][45] = 4.37593838551730e+01;
  pVp[0][48] = -5.29481797958743e+01;

  pVp[1][5] = 1.43713517368872e+01;
  pVp[1][12] = 1.43706526234754e+01;
  pVp[1][19] = 1.43714163288945e+01;
  pVp[1][26] = 1.43721472800142e+01;
  pVp[1][35] = -1.43713517368872e+01;
  pVp[1][36] = -1.43706526234754e+01;
  pVp[1][37] = -1.43714163288945e+01;
  pVp[1][38] = -1.43721472800142e+01;
  pVp[1][41] = 1.75481406677038e+01;
  pVp[1][47] = -1.75481406677038e+01;

  pVp[2][4] = -1.43713517368872e+01;
  pVp[2][11] = -1.43706526234754e+01;
  pVp[2][18] = -1.43714163288945e+01;
  pVp[2][25] = -1.43721472800142e+01;
  pVp[2][28] = 1.43713517368872e+01;
  pVp[2][29] = 1.43706526234754e+01;
  pVp[2][30] = 1.43714163288945e+01;
  pVp[2][31] = 1.43721472800142e+01;
  pVp[2][34] = -1.75481406677038e+01;
  pVp[2][46] = 1.75481406677038e+01;

  pVp[3][33] = 1.76078410818178e+01;
  pVp[3][39] = -1.76078410818178e+01;




  // Load HeH+ custom basis 2 primitive (CQ)
  
  /*
  S[0] = 1.0;
  S[1] = 0.421524100;
  S[2] = 0.421524100;
  S[3] = 1.0;

  T[0] = 9.35870595e-01;
  T[1] = 2.41700233e-01;
  T[2] = 2.41700233e-01;
  T[3] = 1.73838450;

  V[0] = -2.65182646;
  V[1] = -1.53418322;
  V[2] = -1.53418322;
  V[3] = -4.14903710;

  //pVdotp transformed  
  //pVp[0][0] = -2.98297859;
  //pVp[0][1] = 1.88891345e-01;
  //pVp[0][2] = 1.88891345e-01;
  //pVp[0][3] = -2.20165241;

  pVp[0][0] = -4.04003775;
  pVp[0][1] = -1.06076334;
  pVp[0][2] = -1.06076334;
  pVp[0][3] = -1.04220128e+01;

  std::fill_n(pVp[1],nb*nb,0.0); 
  std::fill_n(pVp[2],nb*nb,0.0); 
  std::fill_n(pVp[3],nb*nb,0.0); 
  */


  // Load HeH+ sto-3g (pyscf)
  /*
  S[0] = 1.0;
  S[1] = 0.56088979;
  S[2] = 0.56088979;
  S[3] = 1.0;

  T[0] = 0.76003188;
  T[1] = 0.22156963;
  T[2] = 0.22156963;
  T[3] = 1.41176318;

  V[0] = -2.53536895;
  V[1] = -1.73230852;
  V[2] = -1.73230852;
  V[3] = -4.03736535;
  */

  std::copy_n(S, nb*nb, SCR1);
  auto nbu = X2Chem::orthonormalize(nb, SCR1, SCR2, 1e-14);
  std::cout << "NBU: " << nbu << std::endl;
  X2Chem::detail::print_matrix(nb, SCR1);

  X2Chem::detail::transform(nb, nbu, T, nb, SCR1, nb, SCR2, nb, T, nbu, true);
  X2Chem::detail::transform(nb, nbu, V, nb, SCR1, nb, SCR2, nb, V, nbu, true);
  X2Chem::detail::transform(nb, nbu, pVp[0], nb, SCR1, nb, SCR2, nb, pVp[0], nbu, true);
  X2Chem::detail::transform(nb, nbu, pVp[1], nb, SCR1, nb, SCR2, nb, pVp[1], nbu, true);
  X2Chem::detail::transform(nb, nbu, pVp[2], nb, SCR1, nb, SCR2, nb, pVp[2], nbu, true);
  X2Chem::detail::transform(nb, nbu, pVp[3], nb, SCR1, nb, SCR2, nb, pVp[3], nbu, true);

  x2chem_dbgout("T_t", nbu, T);
  x2chem_dbgout("V_t", nbu, V);
  x2chem_dbgout("pvdp_t", nbu, pVp[0]);
  x2chem_dbgout("pvxpx_t", nbu, pVp[1]);
  x2chem_dbgout("pvxpy_t", nbu, pVp[2]);
  x2chem_dbgout("pvxpz_t", nbu, pVp[3]);

  // initialize structs
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  X2Chem::Integrals integrals{S,T,V,pVp_c};
  X2Chem::X2COperators output{UL, US, coreX2C};

  // X2C routine
  X2Chem::x2c_hamiltonian(nbu, integrals, output, core4C);

  double* transform = new double[nb*nb];
  std::complex<double>* eig = new std::complex<double>[2*nb];
  std::complex<double>* extraMat = new std::complex<double>[2*nb*2*nb];
  std::complex<double>* SCR3 = new std::complex<double>[2*nb*2*nb];
  std::complex<double>* SCR4 = new std::complex<double>[2*nb*2*nb];
  std::complex<double>* ULS = new std::complex<double>[2*nb*2*nb];
  std::complex<double>* USS = new std::complex<double>[2*nb*2*nb];

  std::copy_n(SCR1, nb*nb, transform);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nbu, nb, 1.0, S, nb, SCR1, nb, 0.0, transform, nb);

  // X2Chem::detail::set_submat(nb, nb, SCR1, nb, transform, 2*nb);
  // X2Chem::detail::set_submat(nb, nb, SCR1, nb, transform+2*nb*nb+nb, 2*nb);
  // X2Chem::detail::print_matrix(2*nb, transform);
  // X2Chem::detail::transform(2*nb, extraMat, 2*nb, transform, 2*nb, SCR3, 2*nb, extraMat, 2*nb, false);
  // Top-left
  //
  std::cout << "core Ham" << std::endl;
  X2Chem::detail::print_matrix(2*nbu, coreX2C);

  std::cout << "Transform" << std::endl;
  X2Chem::detail::print_matrix(nb, transform);


  X2Chem::detail::blockTransform(nbu, nb, coreX2C, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, coreX2C, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, coreX2C + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + nb, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, coreX2C + 2*nbu*nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, coreX2C + 2*nbu*nbu + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb + nb, 2*nb, false);

  std::cout << std::setprecision(14) << std::endl;
  for(auto i = 0; i < 4*nb*nb; i++) {
    if( std::abs(extraMat[i]) > 1e-13 ) {
      std::cout << "expectCore[" << i << "] = std::complex<double>(";
      double real = extraMat[i].real();
      double imag = extraMat[i].imag();
      std::cout << (std::abs(real) > 1e-14 ? real : 0.) << ",";
      std::cout << (std::abs(imag) > 1e-14 ? imag : 0.) << ");" << std::endl;
    }
  }

  std::cout << "core Ham" << std::endl;
  X2Chem::detail::print_matrix(2*nbu, extraMat);

  lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, extraMat, 2*nb, transform);

  std::cout << std::setprecision(10) << std::endl;
  for( auto i = 0; i < 2*nb; i++ ) {
    std::cout << "Back eig " << i << ": " << transform[i] << std::endl;
  }

  // std::copy_n(SCR1, nb*nb, transform);
  std::cout << "UL" << "\n";
  X2Chem::detail::print_matrix(2*nb, UL);


  std::copy_n(SCR1, nb*nb, transform);
  X2Chem::detail::blockTransform(nbu, nb, UL, 2*nbu, transform, nb, SCR3, nbu, extraMat, 2*nb, false);
  // X2Chem::detail::transform(nbu, nb, UL, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
  // X2Chem::detail::transform(nbu, nb, UL + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + nb, 2*nb, false);
  // X2Chem::detail::transform(nbu, nb, UL + 2*nbu*nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb, 2*nb, false);
  // X2Chem::detail::transform(nbu, nb, UL + 2*nbu*nbu + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb + nb, 2*nb, false);

  std::cout << "VULV" << std::endl;
  X2Chem::detail::print_matrix(2*nb, extraMat);
  // multiply U
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
  //           nb, nbu, nbu, 1.0, SCR1, nb, UL, 2*nbu, 0.0, extraMat, 2*nb);
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
  //           nb, nbu, nbu, 1.0, SCR1, nb, UL+nbu, 2*nbu, 0.0, extraMat + nb, 2*nb);
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
  //           nb, nbu, nbu, 1.0, SCR1, nb, UL+2*nbu*nbu, 2*nbu, 0.0, extraMat + 2*nb*nbu, 2*nb);
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
  //           nb, nbu, nbu, 1.0, SCR1, nb, UL+2*nbu*nbu+nbu, 2*nbu, 0.0, extraMat + 2*nb*nbu + nb, 2*nb);
  //// multiply U\dagS on right
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
  //           nb, nb, nbu, 1.0, extraMat, 2*nb, transform, nb, 0.0, ULS, 2*nb);
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
  //           nb, nb, nbu, 1.0, extraMat + nb, 2*nb, transform, nb, 0.0, ULS + nb, 2*nb);
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
  //           nb, nb, nbu, 1.0, extraMat + 2*nb*nbu, 2*nb, transform, nb, 0.0, ULS + 2*nb*nb, 2*nb);
  //blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
  //           nb, nb, nbu, 1.0, extraMat + 2*nb*nbu +nb, 2*nb, transform, nb, 0.0, ULS + 2*nb*nb + nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat, 2*nb, S, nb, 0.0, ULS, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat + nb, 2*nb, S, nb, 0.0, ULS + nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat + 2*nb*nb, 2*nb, S, nb, 0.0, ULS + 2*nb*nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat + 2*nb*nb+nb, 2*nb, S, nb, 0.0, ULS + 2*nb*nb + nb, 2*nb);

  std::cout << "ULS backtransform" << "\n";
  X2Chem::detail::print_matrix(2*nb, ULS);

  std::cout << std::setprecision(14) << std::endl;
  for(auto i = 0; i < 4*nb*nb; i++) {
    if( std::abs(ULS[i]) > 1e-13 ) {
      std::cout << "expectUL[" << i << "] = std::complex<double>(";
      double real = ULS[i].real();
      double imag = ULS[i].imag();
      std::cout << (std::abs(real) > 1e-14 ? real : 0.) << ",";
      std::cout << (std::abs(imag) > 1e-14 ? imag : 0.) << ");" << std::endl;
    }
  }

  //lapack::heev(lapack::Job::Vec, lapack::Uplo::Lower, 2*nb, extraMat, 2*nb, transform);
  lapack::geev(lapack::Job::Vec, lapack::Job::NoVec, 2*nb, ULS, 2*nb, eig, SCR4, 2*nb, SCR3, 2*nb);
  std::cout << std::setprecision(10) << std::endl;
  for( auto i = 0; i < 2*nb; i++ ) {
    std::cout << "UL eig " << i << ": " << eig[i] << std::endl;
  }


  std::cout << "US" << "\n";
  X2Chem::detail::print_matrix(2*nb, US);

  X2Chem::detail::blockTransform(nbu, nb, US, 2*nbu, transform, nb, SCR3, nbu, extraMat, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, US, 2*nbu, transform, nb, SCR3, nb, extraMat, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, US + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + nb, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, US + 2*nbu*nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb, 2*nb, false);
  //X2Chem::detail::transform(nbu, nb, US + 2*nbu*nbu + nbu, 2*nbu, transform, nb, SCR3, nb, extraMat + 2*nb*nb + nb, 2*nb, false);

  // multiply S on right
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat, 2*nb, S, nb, 0.0, USS, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat + nb, 2*nb, S, nb, 0.0, USS + nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat + 2*nb*nb, 2*nb, S, nb, 0.0, USS + 2*nb*nb, 2*nb);
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
             nb, nb, nb, 1.0, extraMat + 2*nb*nb +nb, 2*nb, S, nb, 0.0, USS + 2*nb*nb +nb, 2*nb);

  std::cout << "USS backtransform" << "\n";
  X2Chem::detail::print_matrix(2*nb, USS);

  std::cout << std::setprecision(14) << std::endl;
  for(auto i = 0; i < 4*nb*nb; i++) {
    if( std::abs(USS[i]) > 1e-13 ) {
      std::cout << "expectUS[" << i << "] = std::complex<double>(";
      double real = USS[i].real();
      double imag = USS[i].imag();
      std::cout << (std::abs(real) > 1e-14 ? real : 0.) << ",";
      std::cout << (std::abs(imag) > 1e-14 ? imag : 0.) << ");" << std::endl;
    }
  }


  return 0;
}



