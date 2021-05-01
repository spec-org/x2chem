#include <algorithm>
#include <gtest/gtest.h>
#include <lapack.hh>
#include <x2chem.hpp>
#include <x2chem/detail.hpp>


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

void setup_linearDependentInts(std::vector<double>& S, std::vector<double>& T,
                std::vector<double>& V, std::array<double*,4> pVp)
{

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

std::vector<std::complex<double>> getULOrtho(int64_t size) {
  std::vector<std::complex<double>> expectUL(size, 0.);

  expectUL[0] = std::complex<double>(9.99990905451666e-01,0.00000000000000e+00);
  expectUL[3] = std::complex<double>(3.59633681223936e-06,0.00000000000000e+00);
  expectUL[5] = std::complex<double>(6.33536307388156e-09,0.00000000000000e+00);
  expectUL[6] = std::complex<double>(0.00000000000000e+00,6.33536307106944e-09);
  expectUL[9] = std::complex<double>(9.99993346934954e-01,0.00000000000000e+00);
  expectUL[10] = std::complex<double>(0.00000000000000e+00,-6.23212279698970e-09);
  expectUL[12] = std::complex<double>(-6.33536306967103e-09,0.00000000000000e+00);
  expectUL[15] = std::complex<double>(-1.24180151080622e-09,0.00000000000000e+00);
  expectUL[17] = std::complex<double>(0.00000000000000e+00,6.23212279698970e-09);
  expectUL[18] = std::complex<double>(9.99993346934954e-01,0.00000000000000e+00);
  expectUL[20] = std::complex<double>(0.00000000000000e+00,-6.33536304929818e-09);
  expectUL[23] = std::complex<double>(0.00000000000000e+00,-1.24180152600221e-09);
  expectUL[24] = std::complex<double>(3.59633681229488e-06,0.00000000000000e+00);
  expectUL[27] = std::complex<double>(9.99988650628242e-01,0.00000000000000e+00);
  expectUL[29] = std::complex<double>(1.24180155943968e-09,0.00000000000000e+00);
  expectUL[30] = std::complex<double>(0.00000000000000e+00,1.24180156187250e-09);
  expectUL[33] = std::complex<double>(-6.33536306967103e-09,0.00000000000000e+00);
  expectUL[34] = std::complex<double>(0.00000000000000e+00,6.33536304929818e-09);
  expectUL[36] = std::complex<double>(9.99990905451666e-01,0.00000000000000e+00);
  expectUL[39] = std::complex<double>(3.59633681301652e-06,0.00000000000000e+00);
  expectUL[40] = std::complex<double>(6.33536307388155e-09,0.00000000000000e+00);
  expectUL[43] = std::complex<double>(1.24180155943968e-09,0.00000000000000e+00);
  expectUL[45] = std::complex<double>(9.99993346934953e-01,0.00000000000000e+00);
  expectUL[46] = std::complex<double>(0.00000000000000e+00,6.23212237371717e-09);
  expectUL[48] = std::complex<double>(0.00000000000000e+00,-6.33536307106944e-09);
  expectUL[51] = std::complex<double>(0.00000000000000e+00,-1.24180156187250e-09);
  expectUL[53] = std::complex<double>(0.00000000000000e+00,-6.23212237371717e-09);
  expectUL[54] = std::complex<double>(9.99993346934953e-01,0.00000000000000e+00);
  expectUL[57] = std::complex<double>(-1.24180151080622e-09,0.00000000000000e+00);
  expectUL[58] = std::complex<double>(0.00000000000000e+00,1.24180152600221e-09);
  expectUL[60] = std::complex<double>(3.59633681301652e-06,0.00000000000000e+00);
  expectUL[63] = std::complex<double>(9.99988650628242e-01,0.00000000000000e+00);

  return expectUL;
}

std::vector<std::complex<double>> getUSOrtho(int64_t size) {
  std::vector<std::complex<double>> expectUS(size, 0.);

  expectUS[0] = std::complex<double>(9.99826814771920e-01,0.00000000000000e+00);
  expectUS[3] = std::complex<double>(-1.76635065149588e-04,0.00000000000000e+00);
  expectUS[5] = std::complex<double>(-4.76159229310070e-04,0.00000000000000e+00);
  expectUS[6] = std::complex<double>(0.00000000000000e+00,-4.76159229309741e-04);
  expectUS[9] = std::complex<double>(9.99746128678913e-01,2.17013364743612e-14);
  expectUS[10] = std::complex<double>(0.00000000000000e+00,4.68354010751905e-04);
  expectUS[12] = std::complex<double>(4.22923649128596e-04,1.92895206400246e-14);
  expectUS[15] = std::complex<double>(1.88732070392437e-04,0.00000000000000e+00);
  expectUS[17] = std::complex<double>(-3.78030939884866e-14,-4.68354010772046e-04);
  expectUS[18] = std::complex<double>(9.99746128678898e-01,3.12503657933644e-14);
  expectUS[20] = std::complex<double>(0.00000000000000e+00,4.22923649121554e-04);
  expectUS[23] = std::complex<double>(0.00000000000000e+00,1.88732070424737e-04);
  expectUS[24] = std::complex<double>(-7.82497870212007e-05,0.00000000000000e+00);
  expectUS[27] = std::complex<double>(1.00002412321525e+00,0.00000000000000e+00);
  expectUS[29] = std::complex<double>(-9.32564309057209e-05,0.00000000000000e+00);
  expectUS[30] = std::complex<double>(0.00000000000000e+00,-9.32564309061679e-05);
  expectUS[33] = std::complex<double>(4.76159229295499e-04,2.52801523204643e-14);
  expectUS[34] = std::complex<double>(0.00000000000000e+00,-4.76159229295761e-04);
  expectUS[36] = std::complex<double>(9.99826814771909e-01,0.00000000000000e+00);
  expectUS[39] = std::complex<double>(-1.76635065151032e-04,-1.52566472698272e-14);
  expectUS[40] = std::complex<double>(-4.22923649129771e-04,0.00000000000000e+00);
  expectUS[43] = std::complex<double>(-1.88732070414317e-04,0.00000000000000e+00);
  expectUS[45] = std::complex<double>(9.99746128678922e-01,-2.67815012788136e-14);
  expectUS[46] = std::complex<double>(-2.73669975570101e-14,-4.68354010787333e-04);
  expectUS[48] = std::complex<double>(0.00000000000000e+00,4.22923649129039e-04);
  expectUS[51] = std::complex<double>(0.00000000000000e+00,1.88732070416019e-04);
  expectUS[53] = std::complex<double>(2.81441536742477e-14,4.68354010787214e-04);
  expectUS[54] = std::complex<double>(9.99746128678923e-01,-2.64672994892273e-14);
  expectUS[57] = std::complex<double>(9.32564308804507e-05,0.00000000000000e+00);
  expectUS[58] = std::complex<double>(-3.09077308018046e-14,-9.32564308930139e-05);
  expectUS[60] = std::complex<double>(-7.82497870457921e-05,-1.50147656827899e-14);
  expectUS[63] = std::complex<double>(1.00002412321525e+00,0.00000000000000e+00);

  return expectUS;
}

std::vector<std::complex<double>> getULAO(int64_t size) {
  std::vector<std::complex<double>> expectUL(size, 0.);

  expectUL[0] = std::complex<double>(9.99985322968168e-01,0.00000000000000e+00);
  expectUL[3] = std::complex<double>(-3.75686529641150e-06,0.00000000000000e+00);
  expectUL[5] = std::complex<double>(-5.08706700925181e-09,0.00000000000000e+00);
  expectUL[6] = std::complex<double>(0.00000000000000e+00,-5.08706700564307e-09);
  expectUL[9] = std::complex<double>(9.99993346934954e-01,0.00000000000000e+00);
  expectUL[10] = std::complex<double>(0.00000000000000e+00,-6.23212279698970e-09);
  expectUL[12] = std::complex<double>(2.18075283882525e-09,0.00000000000000e+00);
  expectUL[15] = std::complex<double>(-4.92415550503408e-09,0.00000000000000e+00);
  expectUL[17] = std::complex<double>(0.00000000000000e+00,6.23212279698970e-09);
  expectUL[18] = std::complex<double>(9.99993346934954e-01,0.00000000000000e+00);
  expectUL[20] = std::complex<double>(0.00000000000000e+00,2.18075281061591e-09);
  expectUL[23] = std::complex<double>(0.00000000000000e+00,-4.92415551039591e-09);
  expectUL[24] = std::complex<double>(1.50204187276959e-06,0.00000000000000e+00);
  expectUL[27] = std::complex<double>(9.99994233111741e-01,0.00000000000000e+00);
  expectUL[29] = std::complex<double>(6.21127022712266e-09,0.00000000000000e+00);
  expectUL[30] = std::complex<double>(0.00000000000000e+00,6.21127022571634e-09);
  expectUL[33] = std::complex<double>(5.08706702751128e-09,0.00000000000000e+00);
  expectUL[34] = std::complex<double>(0.00000000000000e+00,-5.08706700246658e-09);
  expectUL[36] = std::complex<double>(9.99985322968167e-01,0.00000000000000e+00);
  expectUL[39] = std::complex<double>(-3.75686529685559e-06,0.00000000000000e+00);
  expectUL[40] = std::complex<double>(-2.18075278746541e-09,0.00000000000000e+00);
  expectUL[43] = std::complex<double>(4.92415556111588e-09,0.00000000000000e+00);
  expectUL[45] = std::complex<double>(9.99993346934953e-01,0.00000000000000e+00);
  expectUL[46] = std::complex<double>(0.00000000000000e+00,6.23212237371717e-09);
  expectUL[48] = std::complex<double>(0.00000000000000e+00,2.18075278320127e-09);
  expectUL[51] = std::complex<double>(0.00000000000000e+00,-4.92415556222633e-09);
  expectUL[53] = std::complex<double>(0.00000000000000e+00,-6.23212237371717e-09);
  expectUL[54] = std::complex<double>(9.99993346934953e-01,0.00000000000000e+00);
  expectUL[57] = std::complex<double>(-6.21127020135423e-09,0.00000000000000e+00);
  expectUL[58] = std::complex<double>(0.00000000000000e+00,6.21127019006647e-09);
  expectUL[60] = std::complex<double>(1.50204187332470e-06,0.00000000000000e+00);
  expectUL[63] = std::complex<double>(9.99994233111741e-01,0.00000000000000e+00);

  return expectUL;
}

std::vector<std::complex<double>> getUSAO(int64_t size) {
  std::vector<std::complex<double>> expectUS(size, 0.);

  expectUS[0] = std::complex<double>(1.00011930926894e+00,0.00000000000000e+00);
  expectUS[3] = std::complex<double>(2.52772345584407e-04,0.00000000000000e+00);
  expectUS[5] = std::complex<double>(3.82373041861601e-04,0.00000000000000e+00);
  expectUS[6] = std::complex<double>(0.00000000000000e+00,3.82373041861105e-04);
  expectUS[9] = std::complex<double>(9.99746128678913e-01,2.17013364743612e-14);
  expectUS[10] = std::complex<double>(0.00000000000000e+00,4.68354010751905e-04);
  expectUS[12] = std::complex<double>(-2.86731487259412e-05,-1.09536159847631e-14);
  expectUS[15] = std::complex<double>(4.45622283195567e-04,1.06789683892014e-14);
  expectUS[17] = std::complex<double>(-3.78030939884866e-14,-4.68354010772046e-04);
  expectUS[18] = std::complex<double>(9.99746128678898e-01,3.12503657933644e-14);
  expectUS[20] = std::complex<double>(1.31604573211798e-14,-2.86731486863145e-05);
  expectUS[23] = std::complex<double>(0.00000000000000e+00,4.45622283227297e-04);
  expectUS[24] = std::complex<double>(-5.54639022590786e-05,0.00000000000000e+00);
  expectUS[27] = std::complex<double>(9.99731628718228e-01,0.00000000000000e+00);
  expectUS[29] = std::complex<double>(-4.66798108744805e-04,0.00000000000000e+00);
  expectUS[30] = std::complex<double>(0.00000000000000e+00,-4.66798108744714e-04);
  expectUS[33] = std::complex<double>(-3.82373041860047e-04,-2.50253116813995e-14);
  expectUS[34] = std::complex<double>(-1.43616032868030e-14,3.82373041854593e-04);
  expectUS[36] = std::complex<double>(1.00011930926894e+00,1.65525987502224e-14);
  expectUS[39] = std::complex<double>(2.52772345585850e-04,1.47034842449987e-14);
  expectUS[40] = std::complex<double>(2.86731487024326e-05,0.00000000000000e+00);
  expectUS[43] = std::complex<double>(-4.45622283220394e-04,0.00000000000000e+00);
  expectUS[45] = std::complex<double>(9.99746128678922e-01,-2.67815012788136e-14);
  expectUS[46] = std::complex<double>(-2.73669975570101e-14,-4.68354010787333e-04);
  expectUS[48] = std::complex<double>(0.00000000000000e+00,-2.86731487001410e-05);
  expectUS[51] = std::complex<double>(0.00000000000000e+00,4.45622283221865e-04);
  expectUS[53] = std::complex<double>(2.81441536742477e-14,4.68354010787214e-04);
  expectUS[54] = std::complex<double>(9.99746128678923e-01,-2.64672994892273e-14);
  expectUS[57] = std::complex<double>(4.66798108720373e-04,2.00587163328820e-14);
  expectUS[58] = std::complex<double>(-1.36191720313990e-14,-4.66798108726293e-04);
  expectUS[60] = std::complex<double>(-5.54639022488645e-05,0.00000000000000e+00);
  expectUS[63] = std::complex<double>(9.99731628718215e-01,-2.11238874234042e-14);

  return expectUS;
}

std::vector<std::complex<double>> getExpectedAOLinDep(int64_t size) {
  std::vector<std::complex<double>> expectCore(size, 0.);

  expectCore[0] = std::complex<double>(-5.40279901804441e+01,0.00000000000000e+00);
  expectCore[1] = std::complex<double>(-5.40270415792123e+01,0.00000000000000e+00);
  expectCore[2] = std::complex<double>(-5.40280773598087e+01,0.00000000000000e+00);
  expectCore[3] = std::complex<double>(-5.40290579719522e+01,0.00000000000000e+00);
  expectCore[6] = std::complex<double>(2.95645627976302e+01,0.00000000000000e+00);
  expectCore[11] = std::complex<double>(1.91146102773418e-04,0.00000000000000e+00);
  expectCore[12] = std::complex<double>(0.00000000000000e+00,1.91146102776360e-04);
  expectCore[14] = std::complex<double>(-5.40270415792123e+01,0.00000000000000e+00);
  expectCore[15] = std::complex<double>(-5.40261318532464e+01,0.00000000000000e+00);
  expectCore[16] = std::complex<double>(-5.40271250908074e+01,0.00000000000000e+00);
  expectCore[17] = std::complex<double>(-5.40280632214250e+01,0.00000000000000e+00);
  expectCore[20] = std::complex<double>(2.95696980410912e+01,0.00000000000000e+00);
  expectCore[25] = std::complex<double>(1.91136816336149e-04,0.00000000000000e+00);
  expectCore[26] = std::complex<double>(0.00000000000000e+00,1.91136816339096e-04);
  expectCore[28] = std::complex<double>(-5.40280773598087e+01,0.00000000000000e+00);
  expectCore[29] = std::complex<double>(-5.40271250908074e+01,0.00000000000000e+00);
  expectCore[30] = std::complex<double>(-5.40281648852170e+01,0.00000000000000e+00);
  expectCore[31] = std::complex<double>(-5.40291495053782e+01,0.00000000000000e+00);
  expectCore[34] = std::complex<double>(2.95640767774153e+01,0.00000000000000e+00);
  expectCore[39] = std::complex<double>(1.91146960820344e-04,0.00000000000000e+00);
  expectCore[40] = std::complex<double>(0.00000000000000e+00,1.91146960823285e-04);
  expectCore[42] = std::complex<double>(-5.40290579719522e+01,0.00000000000000e+00);
  expectCore[43] = std::complex<double>(-5.40280632214250e+01,0.00000000000000e+00);
  expectCore[44] = std::complex<double>(-5.40291495053782e+01,0.00000000000000e+00);
  expectCore[45] = std::complex<double>(-5.40301805480547e+01,0.00000000000000e+00);
  expectCore[48] = std::complex<double>(2.95584282826761e+01,0.00000000000000e+00);
  expectCore[53] = std::complex<double>(1.91156669675008e-04,0.00000000000000e+00);
  expectCore[54] = std::complex<double>(0.00000000000000e+00,1.91156669677943e-04);
  expectCore[60] = std::complex<double>(-4.36586802997236e+01,0.00000000000000e+00);
  expectCore[61] = std::complex<double>(-3.46389583683049e-14,2.34177408265159e-04);
  expectCore[63] = std::complex<double>(-1.91146102791568e-04,-1.78879334248712e-14);
  expectCore[64] = std::complex<double>(-1.91136816354322e-04,-1.78858415835939e-14);
  expectCore[65] = std::complex<double>(-1.91146960838491e-04,-1.78881296278397e-14);
  expectCore[66] = std::complex<double>(-1.91156669693129e-04,-1.78903875155872e-14);
  expectCore[69] = std::complex<double>(2.33413073566227e-04,2.78052254314883e-14);
  expectCore[74] = std::complex<double>(-3.55271367880050e-14,-2.34177408265159e-04);
  expectCore[75] = std::complex<double>(-4.36586802997238e+01,0.00000000000000e+00);
  expectCore[77] = std::complex<double>(-3.14888808569579e-14,-1.91146102789623e-04);
  expectCore[78] = std::complex<double>(-3.14974780229584e-14,-1.91136816352372e-04);
  expectCore[79] = std::complex<double>(-3.14880688116237e-14,-1.91146960836547e-04);
  expectCore[80] = std::complex<double>(-3.14786518040782e-14,-1.91156669691191e-04);
  expectCore[83] = std::complex<double>(3.89264088950723e-14,2.33413073569753e-04);
  expectCore[84] = std::complex<double>(2.95645627976302e+01,0.00000000000000e+00);
  expectCore[85] = std::complex<double>(2.95696980410912e+01,0.00000000000000e+00);
  expectCore[86] = std::complex<double>(2.95640767774153e+01,0.00000000000000e+00);
  expectCore[87] = std::complex<double>(2.95584282826761e+01,0.00000000000000e+00);
  expectCore[90] = std::complex<double>(-4.38078931100226e+01,0.00000000000000e+00);
  expectCore[95] = std::complex<double>(-2.33413073539535e-04,0.00000000000000e+00);
  expectCore[96] = std::complex<double>(0.00000000000000e+00,-2.33413073545220e-04);
  expectCore[102] = std::complex<double>(-1.91146102791568e-04,1.78879375863756e-14);
  expectCore[103] = std::complex<double>(-3.14888743768380e-14,1.91146102789624e-04);
  expectCore[105] = std::complex<double>(-5.40279901804442e+01,0.00000000000000e+00);
  expectCore[106] = std::complex<double>(-5.40270415792123e+01,0.00000000000000e+00);
  expectCore[107] = std::complex<double>(-5.40280773598087e+01,0.00000000000000e+00);
  expectCore[108] = std::complex<double>(-5.40290579719523e+01,0.00000000000000e+00);
  expectCore[111] = std::complex<double>(2.95645627976302e+01,0.00000000000000e+00);
  expectCore[116] = std::complex<double>(-1.91136816354314e-04,1.78858457440334e-14);
  expectCore[117] = std::complex<double>(-3.14974715474461e-14,1.91136816352364e-04);
  expectCore[119] = std::complex<double>(-5.40270415792123e+01,0.00000000000000e+00);
  expectCore[120] = std::complex<double>(-5.40261318532464e+01,0.00000000000000e+00);
  expectCore[121] = std::complex<double>(-5.40271250908075e+01,0.00000000000000e+00);
  expectCore[122] = std::complex<double>(-5.40280632214250e+01,0.00000000000000e+00);
  expectCore[125] = std::complex<double>(2.95696980410913e+01,0.00000000000000e+00);
  expectCore[130] = std::complex<double>(-1.91146960838492e-04,1.78881337894435e-14);
  expectCore[131] = std::complex<double>(-3.14880623310701e-14,1.91146960836548e-04);
  expectCore[133] = std::complex<double>(-5.40280773598087e+01,0.00000000000000e+00);
  expectCore[134] = std::complex<double>(-5.40271250908075e+01,0.00000000000000e+00);
  expectCore[135] = std::complex<double>(-5.40281648852171e+01,0.00000000000000e+00);
  expectCore[136] = std::complex<double>(-5.40291495053783e+01,0.00000000000000e+00);
  expectCore[139] = std::complex<double>(2.95640767774153e+01,0.00000000000000e+00);
  expectCore[144] = std::complex<double>(-1.91156669693140e-04,1.78903916783277e-14);
  expectCore[145] = std::complex<double>(-3.14786453185140e-14,1.91156669691201e-04);
  expectCore[147] = std::complex<double>(-5.40290579719523e+01,0.00000000000000e+00);
  expectCore[148] = std::complex<double>(-5.40280632214250e+01,0.00000000000000e+00);
  expectCore[149] = std::complex<double>(-5.40291495053783e+01,0.00000000000000e+00);
  expectCore[150] = std::complex<double>(-5.40301805480548e+01,0.00000000000000e+00);
  expectCore[153] = std::complex<double>(2.95584282826761e+01,0.00000000000000e+00);
  expectCore[154] = std::complex<double>(1.91146102773418e-04,0.00000000000000e+00);
  expectCore[155] = std::complex<double>(1.91136816336157e-04,0.00000000000000e+00);
  expectCore[156] = std::complex<double>(1.91146960820343e-04,0.00000000000000e+00);
  expectCore[157] = std::complex<double>(1.91156669674997e-04,0.00000000000000e+00);
  expectCore[160] = std::complex<double>(-2.33413073539527e-04,0.00000000000000e+00);
  expectCore[165] = std::complex<double>(-4.36586802997236e+01,0.00000000000000e+00);
  expectCore[166] = std::complex<double>(-2.39808173319034e-14,-2.34177408272667e-04);
  expectCore[168] = std::complex<double>(0.00000000000000e+00,-1.91146102776360e-04);
  expectCore[169] = std::complex<double>(0.00000000000000e+00,-1.91136816339105e-04);
  expectCore[170] = std::complex<double>(0.00000000000000e+00,-1.91146960823284e-04);
  expectCore[171] = std::complex<double>(0.00000000000000e+00,-1.91156669677932e-04);
  expectCore[174] = std::complex<double>(0.00000000000000e+00,2.33413073545212e-04);
  expectCore[179] = std::complex<double>(-2.39808173319034e-14,2.34177408272667e-04);
  expectCore[180] = std::complex<double>(-4.36586802997237e+01,0.00000000000000e+00);
  expectCore[186] = std::complex<double>(2.33413073566235e-04,-2.78052334998654e-14);
  expectCore[187] = std::complex<double>(3.89264039168798e-14,-2.33413073569761e-04);
  expectCore[189] = std::complex<double>(2.95645627976302e+01,0.00000000000000e+00);
  expectCore[190] = std::complex<double>(2.95696980410913e+01,0.00000000000000e+00);
  expectCore[191] = std::complex<double>(2.95640767774153e+01,0.00000000000000e+00);
  expectCore[192] = std::complex<double>(2.95584282826761e+01,0.00000000000000e+00);
  expectCore[195] = std::complex<double>(-4.38078931100226e+01,0.00000000000000e+00);

  return expectCore;
}

std::vector<std::complex<double>> getULAOLinDep(int64_t size) {
  std::vector<std::complex<double>> expectUL(size, 0.);

  expectUL[0] = std::complex<double>(4.95403838387574e-01,1.72359968462704e-13);
  expectUL[1] = std::complex<double>(3.05803449555242e-02,-1.81723695463658e-13);
  expectUL[2] = std::complex<double>(4.97184372457923e-01,1.75942062542370e-13);
  expectUL[3] = std::complex<double>(-2.32003980290756e-02,-1.66577918710591e-13);
  expectUL[6] = std::complex<double>(-3.09158276756705e-09,0.00000000000000e+00);
  expectUL[7] = std::complex<double>(1.28035992972699e-12,-1.04350824591134e-12);
  expectUL[8] = std::complex<double>(-1.34812044689135e-12,1.09827966431408e-12);
  expectUL[9] = std::complex<double>(1.30681838939236e-12,-1.06503417054234e-12);
  expectUL[10] = std::complex<double>(-1.23906434832351e-12,1.01026824301364e-12);
  expectUL[11] = std::complex<double>(-5.08693540878997e-09,0.00000000000000e+00);
  expectUL[12] = std::complex<double>(0.00000000000000e+00,-5.08693544015407e-09);
  expectUL[14] = std::complex<double>(1.28780759514484e-02,1.72581567363360e-13);
  expectUL[15] = std::complex<double>(1.01573863100202e+00,-1.81957063357879e-13);
  expectUL[16] = std::complex<double>(-3.52674083151214e-02,1.76168244254412e-13);
  expectUL[17] = std::complex<double>(6.58184160784003e-03,-1.66792331414001e-13);
  expectUL[20] = std::complex<double>(-3.09652992136478e-09,0.00000000000000e+00);
  expectUL[21] = std::complex<double>(1.28144991826270e-12,-1.04470406117591e-12);
  expectUL[22] = std::complex<double>(-1.34925609979872e-12,1.09953062641103e-12);
  expectUL[23] = std::complex<double>(1.30792989757561e-12,-1.06625401615812e-12);
  expectUL[24] = std::complex<double>(-1.24013015829325e-12,1.01143294047854e-12);
  expectUL[25] = std::complex<double>(-5.08668842186820e-09,0.00000000000000e+00);
  expectUL[26] = std::complex<double>(0.00000000000000e+00,-5.08668845324297e-09);
  expectUL[28] = std::complex<double>(4.98341112670460e-01,1.72339056772012e-13);
  expectUL[29] = std::complex<double>(-1.75103533292713e-02,-1.81701673147732e-13);
  expectUL[30] = std::complex<double>(5.03933297874028e-01,1.75920718382684e-13);
  expectUL[31] = std::complex<double>(1.51638546612958e-02,-1.66557685177609e-13);
  expectUL[34] = std::complex<double>(-3.11086778559400e-09,0.00000000000000e+00);
  expectUL[35] = std::complex<double>(1.28025708050389e-12,-1.04339543899561e-12);
  expectUL[36] = std::complex<double>(-1.34801328812937e-12,1.09816165490104e-12);
  expectUL[37] = std::complex<double>(1.30671350953033e-12,-1.06491909682676e-12);
  expectUL[38] = std::complex<double>(-1.23896378124536e-12,1.01015837189166e-12);
  expectUL[39] = std::complex<double>(-5.08695822943571e-09,0.00000000000000e+00);
  expectUL[40] = std::complex<double>(0.00000000000000e+00,-5.08695826079880e-09);
  expectUL[42] = std::complex<double>(-1.00038251475780e-02,1.72096793687153e-13);
  expectUL[43] = std::complex<double>(-3.80562665668549e-03,-1.81446543192833e-13);
  expectUL[44] = std::complex<double>(2.87381255802757e-02,1.75673445102864e-13);
  expectUL[45] = std::complex<double>(9.84935118491194e-01,-1.66323278785094e-13);
  expectUL[48] = std::complex<double>(-3.08331982168397e-09,0.00000000000000e+00);
  expectUL[49] = std::complex<double>(1.27906569782247e-12,-1.04208906301978e-12);
  expectUL[50] = std::complex<double>(-1.34677197491335e-12,1.09679502631832e-12);
  expectUL[51] = std::complex<double>(1.30549860379845e-12,-1.06358646818748e-12);
  expectUL[52] = std::complex<double>(-1.23779884287600e-12,1.00888599732619e-12);
  expectUL[53] = std::complex<double>(-5.08721644505975e-09,0.00000000000000e+00);
  expectUL[54] = std::complex<double>(0.00000000000000e+00,-5.08721647641117e-09);
  expectUL[56] = std::complex<double>(1.91844161094351e-10,-1.42584476468741e-11);
  expectUL[57] = std::complex<double>(-2.02087542219061e-10,1.50130550701988e-11);
  expectUL[58] = std::complex<double>(1.95816201374716e-10,-1.45530998822994e-11);
  expectUL[59] = std::complex<double>(-1.85572748205515e-10,1.37984907534094e-11);
  expectUL[60] = std::complex<double>(9.99993346934782e-01,0.00000000000000e+00);
  expectUL[61] = std::complex<double>(0.00000000000000e+00,-6.23195186037414e-09);
  expectUL[63] = std::complex<double>(2.10326952639744e-03,3.08809968976107e-11);
  expectUL[64] = std::complex<double>(-2.21487662188125e-03,-3.25537141771557e-11);
  expectUL[65] = std::complex<double>(2.14675844280671e-03,3.15223673633889e-11);
  expectUL[66] = std::complex<double>(-2.03515165744504e-03,-2.98496898542689e-11);
  expectUL[69] = std::complex<double>(-5.68950445581233e-09,0.00000000000000e+00);
  expectUL[70] = std::complex<double>(-3.24846809710427e-11,2.40041423573131e-12);
  expectUL[71] = std::complex<double>(3.41007875685963e-11,-2.52746926096450e-12);
  expectUL[72] = std::complex<double>(-3.31473487545405e-11,2.45002042370163e-12);
  expectUL[73] = std::complex<double>(3.15312692931559e-11,-2.32296516228677e-12);
  expectUL[74] = std::complex<double>(0.00000000000000e+00,6.23195186037414e-09);
  expectUL[75] = std::complex<double>(9.99993346934785e-01,0.00000000000000e+00);
  expectUL[77] = std::complex<double>(-3.79805592105628e-11,2.10326959998787e-03);
  expectUL[78] = std::complex<double>(4.00042472320616e-11,-2.21487669939401e-03);
  expectUL[79] = std::complex<double>(-3.87665656228151e-11,2.14675851792021e-03);
  expectUL[80] = std::complex<double>(3.67430305403554e-11,-2.03515172863636e-03);
  expectUL[83] = std::complex<double>(0.00000000000000e+00,-5.68950455910091e-09);
  expectUL[84] = std::complex<double>(8.53108909788716e-01,-4.63174803914542e-13);
  expectUL[85] = std::complex<double>(-8.98546356333100e-01,4.88036292972044e-13);
  expectUL[86] = std::complex<double>(8.70823637469584e-01,-4.72775559975489e-13);
  expectUL[87] = std::complex<double>(-8.25345606303017e-01,4.47913532438073e-13);
  expectUL[90] = std::complex<double>(9.99993789603317e-01,0.00000000000000e+00);
  expectUL[91] = std::complex<double>(-3.02683894535147e-12,2.86733098653561e-12);
  expectUL[92] = std::complex<double>(3.17310890559773e-12,-3.00903486799116e-12);
  expectUL[93] = std::complex<double>(-3.08822468030112e-12,2.92574352621000e-12);
  expectUL[94] = std::complex<double>(2.94192454831813e-12,-2.78404585710540e-12);
  expectUL[95] = std::complex<double>(6.21199627589057e-09,0.00000000000000e+00);
  expectUL[96] = std::complex<double>(0.00000000000000e+00,6.21199633663242e-09);
  expectUL[98] = std::complex<double>(-8.10357067150789e-12,-2.49173682530823e-12);
  expectUL[99] = std::complex<double>(8.54137979186578e-12,2.62659171375909e-12);
  expectUL[100] = std::complex<double>(-8.27178063673228e-12,-2.54347874075896e-12);
  expectUL[101] = std::complex<double>(7.83395591878996e-12,2.40861849565788e-12);
  expectUL[102] = std::complex<double>(5.08693560245097e-09,0.00000000000000e+00);
  expectUL[103] = std::complex<double>(0.00000000000000e+00,-5.08693558274887e-09);
  expectUL[105] = std::complex<double>(4.95464873583842e-01,-4.30925342936662e-11);
  expectUL[106] = std::complex<double>(3.05193097592564e-02,4.53064733161328e-11);
  expectUL[107] = std::complex<double>(4.97245407656010e-01,-4.39774629141780e-11);
  expectUL[108] = std::complex<double>(-2.32614332235244e-02,4.17635710804895e-11);
  expectUL[111] = std::complex<double>(-3.09522074637414e-09,0.00000000000000e+00);
  expectUL[112] = std::complex<double>(-8.10858337453462e-12,-2.49627539069674e-12);
  expectUL[113] = std::complex<double>(8.54666308128764e-12,2.63137518700637e-12);
  expectUL[114] = std::complex<double>(-8.27689737148339e-12,-2.54811149069615e-12);
  expectUL[115] = std::complex<double>(7.83880205761083e-12,2.41300632903898e-12);
  expectUL[116] = std::complex<double>(5.08668861514931e-09,0.00000000000000e+00);
  expectUL[117] = std::complex<double>(0.00000000000000e+00,-5.08668859590316e-09);
  expectUL[119] = std::complex<double>(1.29391111477162e-02,-4.29532823548250e-11);
  expectUL[120] = std::complex<double>(1.01567759580576e+00,4.51598768608249e-11);
  expectUL[121] = std::complex<double>(-3.52063731170347e-02,-4.38353354331269e-11);
  expectUL[122] = std::complex<double>(6.52080641157227e-03,4.16287880151789e-11);
  expectUL[125] = std::complex<double>(-3.10016790017187e-09,0.00000000000000e+00);
  expectUL[126] = std::complex<double>(-8.10309739927798e-12,-2.49130851681916e-12);
  expectUL[127] = std::complex<double>(8.54088097232738e-12,2.62614029305910e-12);
  expectUL[128] = std::complex<double>(-8.27129754237284e-12,-2.54304154390356e-12);
  expectUL[129] = std::complex<double>(7.83349837265734e-12,2.40820441181805e-12);
  expectUL[130] = std::complex<double>(5.08695842309162e-09,0.00000000000000e+00);
  expectUL[131] = std::complex<double>(0.00000000000000e+00,-5.08695840338743e-09);
  expectUL[133] = std::complex<double>(4.98402147866727e-01,-4.31056693930127e-11);
  expectUL[134] = std::complex<double>(-1.75713885255391e-02,4.53203012078748e-11);
  expectUL[135] = std::complex<double>(5.03994333070295e-01,-4.39908692533988e-11);
  expectUL[136] = std::complex<double>(1.51028194650280e-02,4.17762846368083e-11);
  expectUL[139] = std::complex<double>(-3.11450576440109e-09,0.00000000000000e+00);
  expectUL[140] = std::complex<double>(-8.09761153259977e-12,-2.48634632494266e-12);
  expectUL[141] = std::complex<double>(8.53509897763852e-12,2.62091033384926e-12);
  expectUL[142] = std::complex<double>(-8.26569782414757e-12,-2.53797637636354e-12);
  expectUL[143] = std::complex<double>(7.82819479289499e-12,2.40340702116709e-12);
  expectUL[144] = std::complex<double>(5.08721663822309e-09,0.00000000000000e+00);
  expectUL[145] = std::complex<double>(0.00000000000000e+00,-5.08721661849451e-09);
  expectUL[147] = std::complex<double>(-9.94278995131026e-03,-4.32577683074414e-11);
  expectUL[148] = std::complex<double>(-3.86666185295326e-03,4.54804225240264e-11);
  expectUL[149] = std::complex<double>(2.87991607765434e-02,-4.41461090243998e-11);
  expectUL[150] = std::complex<double>(9.84874083296745e-01,4.19235021158490e-11);
  expectUL[153] = std::complex<double>(-3.08695780049106e-09,0.00000000000000e+00);
  expectUL[154] = std::complex<double>(-2.10326950933934e-03,1.60568056028156e-11);
  expectUL[155] = std::complex<double>(2.21487660399848e-03,-1.69003904923199e-11);
  expectUL[156] = std::complex<double>(-2.14675842540265e-03,1.63881018711455e-11);
  expectUL[157] = std::complex<double>(2.03515164086568e-03,-1.55445252312553e-11);
  expectUL[160] = std::complex<double>(5.68950420631519e-09,0.00000000000000e+00);
  expectUL[161] = std::complex<double>(2.01940228932063e-10,3.39618373528690e-11);
  expectUL[162] = std::complex<double>(-2.12641915271647e-10,-3.58023869510588e-11);
  expectUL[163] = std::complex<double>(2.06114534208972e-10,3.46672785240637e-11);
  expectUL[164] = std::complex<double>(-1.95412990295599e-10,-3.28266513551320e-11);
  expectUL[165] = std::complex<double>(9.99993346934783e-01,0.00000000000000e+00);
  expectUL[166] = std::complex<double>(0.00000000000000e+00,6.23195194580928e-09);
  expectUL[168] = std::complex<double>(-8.23750195349594e-13,2.10326957571497e-03);
  expectUL[169] = std::complex<double>(8.69601586825995e-13,-2.21487667383947e-03);
  expectUL[170] = std::complex<double>(-8.40961966172108e-13,2.14675849314596e-03);
  expectUL[171] = std::complex<double>(7.95108484899033e-13,-2.03515170514370e-03);
  expectUL[174] = std::complex<double>(0.00000000000000e+00,-5.68950425578133e-09);
  expectUL[175] = std::complex<double>(-2.47951182122632e-11,4.34453379966003e-11);
  expectUL[176] = std::complex<double>(2.60518066372342e-11,-4.57672193481438e-11);
  expectUL[177] = std::complex<double>(-2.53028614552387e-11,4.43450338182889e-11);
  expectUL[178] = std::complex<double>(2.40461627075377e-11,-4.20231015825245e-11);
  expectUL[179] = std::complex<double>(0.00000000000000e+00,-6.23195194580928e-09);
  expectUL[180] = std::complex<double>(9.99993346934784e-01,0.00000000000000e+00);
  expectUL[182] = std::complex<double>(2.91749727708479e-12,4.14347724117908e-12);
  expectUL[183] = std::complex<double>(-3.07562400295459e-12,-4.36724322985783e-12);
  expectUL[184] = std::complex<double>(2.97809949837698e-12,4.22947782008611e-12);
  expectUL[185] = std::complex<double>(-2.81996590266408e-12,-4.00570296351434e-12);
  expectUL[186] = std::complex<double>(-6.21199658480199e-09,0.00000000000000e+00);
  expectUL[187] = std::complex<double>(0.00000000000000e+00,6.21199661523219e-09);
  expectUL[189] = std::complex<double>(8.53108909719595e-01,-4.34888673738973e-11);
  expectUL[190] = std::complex<double>(-8.98515838685853e-01,4.58465930885105e-11);
  expectUL[191] = std::complex<double>(8.70762602244213e-01,-4.43922674096305e-11);
  expectUL[192] = std::complex<double>(-8.25315088659409e-01,4.20345696981940e-11);
  expectUL[195] = std::complex<double>(9.99993789603317e-01,0.00000000000000e+00);

  return expectUL;
}

std::vector<std::complex<double>> getUSAOLinDep(int64_t size) {
  std::vector<std::complex<double>> expectUS(size, 0.);

  expectUS[0] = std::complex<double>(-1.82565084655889e+01,-2.94711454335115e-10);
  expectUS[1] = std::complex<double>(1.98660173448516e+01,3.10621135930229e-10);
  expectUS[2] = std::complex<double>(-1.86496590963943e+01,-3.00827873360676e-10);
  expectUS[3] = std::complex<double>(1.80405872737319e+01,2.84917690091023e-10);
  expectUS[6] = std::complex<double>(4.14657751177661e-04,0.00000000000000e+00);
  expectUS[7] = std::complex<double>(8.59755824016740e-12,1.58842622352774e-12);
  expectUS[8] = std::complex<double>(-9.04891910949334e-12,-1.67355679102301e-12);
  expectUS[9] = std::complex<double>(8.77492210872133e-12,1.62134071069098e-12);
  expectUS[10] = std::complex<double>(-8.32355484043569e-12,-1.53620162810569e-12);
  expectUS[11] = std::complex<double>(3.82373445508800e-04,0.00000000000000e+00);
  expectUS[12] = std::complex<double>(0.00000000000000e+00,3.82373445536124e-04);
  expectUS[14] = std::complex<double>(-1.87505338084138e+01,-2.94391748917018e-10);
  expectUS[15] = std::complex<double>(2.08634013093069e+01,3.10284164139345e-10);
  expectUS[16] = std::complex<double>(-1.91939535511156e+01,-3.00501532171330e-10);
  expectUS[17] = std::complex<double>(1.80814249801788e+01,2.84608615865408e-10);
  expectUS[20] = std::complex<double>(4.14912215713659e-04,0.00000000000000e+00);
  expectUS[21] = std::complex<double>(8.61078456225753e-12,1.58609348608151e-12);
  expectUS[22] = std::complex<double>(-9.06282884954254e-12,-1.67110637514427e-12);
  expectUS[23] = std::complex<double>(8.78842036974535e-12,1.61896025086231e-12);
  expectUS[24] = std::complex<double>(-8.33636962570629e-12,-1.53393883761175e-12);
  expectUS[25] = std::complex<double>(3.82354839726443e-04,0.00000000000000e+00);
  expectUS[26] = std::complex<double>(0.00000000000000e+00,3.82354839753767e-04);
  expectUS[28] = std::complex<double>(-1.82524983343865e+01,-2.94741611311575e-10);
  expectUS[29] = std::complex<double>(1.98168130180966e+01,3.10652921622093e-10);
  expectUS[30] = std::complex<double>(-1.86418683193242e+01,-3.00858656273973e-10);
  expectUS[31] = std::complex<double>(1.80778894181421e+01,2.84946844258742e-10);
  expectUS[34] = std::complex<double>(4.14633691344424e-04,0.00000000000000e+00);
  expectUS[35] = std::complex<double>(8.59630851722967e-12,1.58864640825516e-12);
  expectUS[36] = std::complex<double>(-9.04760480658721e-12,-1.67378808396149e-12);
  expectUS[37] = std::complex<double>(8.77364668913068e-12,1.62156540002786e-12);
  expectUS[38] = std::complex<double>(-8.32234400396596e-12,-1.53641520982253e-12);
  expectUS[39] = std::complex<double>(3.82375164711435e-04,0.00000000000000e+00);
  expectUS[40] = std::complex<double>(0.00000000000000e+00,3.82375164738759e-04);
  expectUS[42] = std::complex<double>(-1.87480370135472e+01,-2.95090822932839e-10);
  expectUS[43] = std::complex<double>(1.98169451793619e+01,3.11020993200026e-10);
  expectUS[44] = std::complex<double>(-1.91038843241549e+01,-3.01215116084304e-10);
  expectUS[45] = std::complex<double>(1.90353089374148e+01,2.85284443518416e-10);
  expectUS[48] = std::complex<double>(4.14354742477974e-04,0.00000000000000e+00);
  expectUS[49] = std::complex<double>(8.58181011441341e-12,1.59119792970487e-12);
  expectUS[50] = std::complex<double>(-9.03235721129383e-12,-1.67646834149036e-12);
  expectUS[51] = std::complex<double>(8.75885018852935e-12,1.62416912097194e-12);
  expectUS[52] = std::complex<double>(-8.30829676188096e-12,-1.53889020437888e-12);
  expectUS[53] = std::complex<double>(3.82394618155034e-04,0.00000000000000e+00);
  expectUS[54] = std::complex<double>(0.00000000000000e+00,3.82394618210780e-04);
  expectUS[56] = std::complex<double>(-1.92110183926785e-09,1.71248533048747e-09);
  expectUS[57] = std::complex<double>(2.02442394177252e-09,-1.80359606890455e-09);
  expectUS[58] = std::complex<double>(-1.96093984998991e-09,1.74791421763438e-09);
  expectUS[59] = std::complex<double>(1.85761540605593e-09,-1.65680246870428e-09);
  expectUS[60] = std::complex<double>(9.99746137214690e-01,-1.69527749043635e-14);
  expectUS[61] = std::complex<double>(0.00000000000000e+00,4.68345474582400e-04);
  expectUS[63] = std::complex<double>(-2.57884907045832e+01,-8.08169777659495e-10);
  expectUS[64] = std::complex<double>(2.71733394018661e+01,8.50272671068658e-10);
  expectUS[65] = std::complex<double>(-2.63230870981119e+01,-8.24815678203182e-10);
  expectUS[66] = std::complex<double>(2.49382902085979e+01,7.82687943638499e-10);
  expectUS[69] = std::complex<double>(4.80923329377406e-04,-2.39532618641014e-14);
  expectUS[70] = std::complex<double>(5.42526314063400e-09,-5.11174542810685e-10);
  expectUS[71] = std::complex<double>(-5.71049280871369e-09,5.38300496111903e-10);
  expectUS[72] = std::complex<double>(5.53721799873804e-09,-5.21744093113514e-10);
  expectUS[73] = std::complex<double>(-5.25199197935845e-09,4.94618490870881e-10);
  expectUS[74] = std::complex<double>(-2.50355292052973e-14,-4.68345474565394e-04);
  expectUS[75] = std::complex<double>(9.99746137214751e-01,0.00000000000000e+00);
  expectUS[77] = std::complex<double>(1.73018293524625e-09,-2.57884907273977e+01);
  expectUS[78] = std::complex<double>(-1.83586112399710e-09,2.71733394258952e+01);
  expectUS[79] = std::complex<double>(1.76711774785976e-09,-2.63230871213985e+01);
  expectUS[80] = std::complex<double>(-1.66146029760919e-09,2.49382902306698e+01);
  expectUS[83] = std::complex<double>(-2.11342561318926e-14,4.80923329431474e-04);
  expectUS[84] = std::complex<double>(-3.19210165719633e+01,1.58772133123609e-10);
  expectUS[85] = std::complex<double>(3.36044823188786e+01,-1.67337678145864e-10);
  expectUS[86] = std::complex<double>(-3.25801977330993e+01,1.62066808848908e-10);
  expectUS[87] = std::complex<double>(3.08966273639126e+01,-1.53501003763101e-10);
  expectUS[90] = std::complex<double>(9.99701131538387e-01,0.00000000000000e+00);
  expectUS[91] = std::complex<double>(-2.08152883967750e-11,4.17724290695224e-12);
  expectUS[92] = std::complex<double>(2.18905361880492e-11,-4.39016178687162e-12);
  expectUS[93] = std::complex<double>(-2.12432331344250e-11,4.26288265309616e-12);
  expectUS[94] = std::complex<double>(2.01679242596827e-11,-4.04998172487447e-12);
  expectUS[95] = std::complex<double>(-4.66889306715609e-04,0.00000000000000e+00);
  expectUS[96] = std::complex<double>(0.00000000000000e+00,-4.66889306727960e-04);
  expectUS[98] = std::complex<double>(-8.12964494469848e-09,-1.17281436150124e-09);
  expectUS[99] = std::complex<double>(8.56470144194595e-09,1.23510876574342e-09);
  expectUS[100] = std::complex<double>(-8.29804727388031e-09,-1.19706949533496e-09);
  expectUS[101] = std::complex<double>(7.86298919812372e-09,1.13477546843686e-09);
  expectUS[102] = std::complex<double>(-3.82373445551561e-04,1.47283047822977e-14);
  expectUS[103] = std::complex<double>(2.69059073103515e-14,3.82373445529456e-04);
  expectUS[105] = std::complex<double>(-1.82563863963733e+01,4.11326593028287e-09);
  expectUS[106] = std::complex<double>(1.98660173460084e+01,-4.33287794205588e-09);
  expectUS[107] = std::complex<double>(-1.86496590975130e+01,4.19842820339689e-09);
  expectUS[108] = std::complex<double>(1.80406483099450e+01,-3.97881003659324e-09);
  expectUS[111] = std::complex<double>(4.14657754814085e-04,0.00000000000000e+00);
  expectUS[112] = std::complex<double>(-8.13406668192462e-09,-1.17334651947230e-09);
  expectUS[113] = std::complex<double>(8.56935996680688e-09,1.23566925025488e-09);
  expectUS[114] = std::complex<double>(-8.30256061898788e-09,-1.19761266405959e-09);
  expectUS[115] = std::complex<double>(7.86726575404894e-09,1.13529031030722e-09);
  expectUS[116] = std::complex<double>(-3.82354839769209e-04,1.47512167452681e-14);
  expectUS[117] = std::complex<double>(2.69309615496147e-14,3.82354839718675e-04);
  expectUS[119] = std::complex<double>(-1.87504117392000e+01,4.11427847307245e-09);
  expectUS[120] = std::complex<double>(2.08634013104638e+01,-4.33394272656651e-09);
  expectUS[121] = std::complex<double>(-1.91939535522342e+01,4.19946155717443e-09);
  expectUS[122] = std::complex<double>(1.80814860163919e+01,-3.97979115654027e-09);
  expectUS[125] = std::complex<double>(4.14912219350083e-04,0.00000000000000e+00);
  expectUS[126] = std::complex<double>(-8.12922756545683e-09,-1.17276408140628e-09);
  expectUS[127] = std::complex<double>(8.56426171176463e-09,1.23505580929955e-09);
  expectUS[128] = std::complex<double>(-8.29762124753170e-09,-1.19701817491253e-09);
  expectUS[129] = std::complex<double>(7.86258552282600e-09,1.13472682437507e-09);
  expectUS[130] = std::complex<double>(-3.82375164725773e-04,1.47261412711283e-14);
  expectUS[131] = std::complex<double>(2.69035407093643e-14,3.82375164703670e-04);
  expectUS[133] = std::complex<double>(-1.82523762651708e+01,4.11316990250813e-09);
  expectUS[134] = std::complex<double>(1.98168130192535e+01,-4.33277695898704e-09);
  expectUS[135] = std::complex<double>(-1.86418683204429e+01,4.19833020155474e-09);
  expectUS[136] = std::complex<double>(1.80779504543552e+01,-3.97871698983246e-09);
  expectUS[139] = std::complex<double>(4.14633694980848e-04,0.00000000000000e+00);
  expectUS[140] = std::complex<double>(-8.12439081723748e-09,-1.17218080328805e-09);
  expectUS[141] = std::complex<double>(8.55916595352391e-09,1.23444148426134e-09);
  expectUS[142] = std::complex<double>(-8.29268429361391e-09,-1.19642282842722e-09);
  expectUS[143] = std::complex<double>(7.85790758058837e-09,1.13416252513542e-09);
  expectUS[144] = std::complex<double>(-3.82394618254632e-04,1.47010594113994e-14);
  expectUS[145] = std::complex<double>(2.68760943645165e-14,3.82394618232537e-04);
  expectUS[147] = std::complex<double>(-1.87479149443316e+01,4.11205136845048e-09);
  expectUS[148] = std::complex<double>(1.98169451805188e+01,-4.33160070162567e-09);
  expectUS[149] = std::complex<double>(-1.91038843252718e+01,4.19718867722424e-09);
  expectUS[150] = std::complex<double>(1.90353699736261e+01,-3.97763317939987e-09);
  expectUS[153] = std::complex<double>(4.14354746114398e-04,0.00000000000000e+00);
  expectUS[154] = std::complex<double>(2.57884907169934e+01,4.14149491290891e-10);
  expectUS[155] = std::complex<double>(-2.71733394149397e+01,-4.36178752433148e-10);
  expectUS[156] = std::complex<double>(2.63230871107791e+01,4.22717237171377e-10);
  expectUS[157] = std::complex<double>(-2.49382902206017e+01,-4.00687575703760e-10);
  expectUS[160] = std::complex<double>(-4.80923329424853e-04,0.00000000000000e+00);
  expectUS[161] = std::complex<double>(1.20950246001308e-08,-5.68787773925303e-09);
  expectUS[162] = std::complex<double>(-1.27509207906029e-08,5.99046952271961e-09);
  expectUS[163] = std::complex<double>(1.23462878483846e-08,-5.80554903623443e-09);
  expectUS[164] = std::complex<double>(-1.16904686170741e-08,5.50297290026826e-09);
  expectUS[165] = std::complex<double>(9.99746137214699e-01,0.00000000000000e+00);
  expectUS[166] = std::complex<double>(0.00000000000000e+00,-4.68345474535757e-04);
  expectUS[168] = std::complex<double>(-9.06316941316847e-10,-2.57884907170662e+01);
  expectUS[169] = std::complex<double>(9.54503499267870e-10,2.71733394150165e+01);
  expectUS[170] = std::complex<double>(-9.25064593940314e-10,-2.63230871108534e+01);
  expectUS[171] = std::complex<double>(8.76877449416729e-10,2.49382902206720e+01);
  expectUS[174] = std::complex<double>(0.00000000000000e+00,4.80923329441206e-04);
  expectUS[175] = std::complex<double>(4.91559576741757e-09,1.18690115721054e-08);
  expectUS[176] = std::complex<double>(-5.17657446673207e-09,-1.25128209219887e-08);
  expectUS[177] = std::complex<double>(5.01724579996139e-09,1.21155937696797e-08);
  expectUS[178] = std::complex<double>(-4.75628115198183e-09,-1.14718605903389e-08);
  expectUS[179] = std::complex<double>(0.00000000000000e+00,4.68345474528796e-04);
  expectUS[180] = std::complex<double>(9.99746137214693e-01,0.00000000000000e+00);
  expectUS[182] = std::complex<double>(5.56926390759873e-09,-1.45395859968984e-11);
  expectUS[183] = std::complex<double>(-5.86712190845098e-09,1.53620524352336e-11);
  expectUS[184] = std::complex<double>(5.68461390372731e-09,-1.48444749319978e-11);
  expectUS[185] = std::complex<double>(-5.38675448161316e-09,1.40222180661254e-11);
  expectUS[186] = std::complex<double>(4.66889306754114e-04,-3.37127812376875e-14);
  expectUS[187] = std::complex<double>(-3.31466718976541e-14,-4.66889306702452e-04);
  expectUS[189] = std::complex<double>(-3.19211081228386e+01,-2.42560483040041e-09);
  expectUS[190] = std::complex<double>(3.36044517993396e+01,2.55285878982795e-09);
  expectUS[191] = std::complex<double>(-3.25801977312040e+01,-2.47563645858887e-09);
  expectUS[192] = std::complex<double>(3.08966273621227e+01,2.34838701797425e-09);
  expectUS[195] = std::complex<double>(9.99701131534752e-01,2.46764745394750e-14);

  return expectUS;
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

  // X2C Large picture change
  std::vector<std::complex<double>> expectUL = getULOrtho(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(UL[i].real(), expectUL[i].real(), 1e-8);
    EXPECT_NEAR(UL[i].imag(), expectUL[i].imag(), 1e-8);
  }

  // X2C small picture change
  std::vector<std::complex<double>> expectUS = getUSOrtho(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(US[i].real(), expectUS[i].real(), 1e-8);
    EXPECT_NEAR(US[i].imag(), expectUS[i].imag(), 1e-8);
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

  // X2C Large picture change
  std::vector<std::complex<double>> expectUL = getULOrtho(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(UL[i].real(), expectUL[i].real(), 1e-8);
    EXPECT_NEAR(UL[i].imag(), expectUL[i].imag(), 1e-8);
  }

  // X2C small picture change
  std::vector<std::complex<double>> expectUS = getUSOrtho(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(US[i].real(), expectUS[i].real(), 1e-8);
    EXPECT_NEAR(US[i].imag(), expectUS[i].imag(), 1e-8);
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

  // X2C Core
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

  // X2C Large picture change
  std::vector<std::complex<double>> expectUL = getULAO(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(UL[i].real(), expectUL[i].real(), 1e-8);
    EXPECT_NEAR(UL[i].imag(), expectUL[i].imag(), 1e-8);
  }

  // X2C small picture change
  std::vector<std::complex<double>> expectUS = getUSAO(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(US[i].real(), expectUS[i].real(), 1e-8);
    EXPECT_NEAR(US[i].imag(), expectUS[i].imag(), 1e-8);
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

  // X2C Core matrix
  std::vector<std::complex<double>> expectCore = getExpectedAO(4*nbsq);

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

  // X2C Large picture change
  std::vector<std::complex<double>> expectUL = getULAO(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(UL[i].real(), expectUL[i].real(), 1e-8);
    EXPECT_NEAR(UL[i].imag(), expectUL[i].imag(), 1e-8);
  }

  // X2C small picture change
  std::vector<std::complex<double>> expectUS = getUSAO(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(US[i].real(), expectUS[i].real(), 1e-8);
    EXPECT_NEAR(US[i].imag(), expectUS[i].imag(), 1e-8);
  }

  delete[] eig;
  delete[] moreSCR;

}

TEST( X2C_Hamiltonian, UH_91Plus_AO_LinDep ) {
  int64_t nb = 7;
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
  setup_linearDependentInts(S, T, V, pVp);

  // Create structs
  std::array<const double*, 4> pVp_c = {pVp[0], pVp[1], pVp[2], pVp[3]};
  Integrals ints{S.data(), T.data(), V.data(), pVp_c};
  X2COperators out{UL.data(), US.data(), coreX2C.data()};

  x2c_hamiltonian_ao(nb, ints, out, SCR2.data(), SCR.data());

  //
  // Expected values
  //

  // X2C Core
  std::vector<std::complex<double>> expectCore = getExpectedAOLinDep(4*nbsq);
  detail::print_matrix(2*nb, expectCore.data());

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(coreX2C[i].real(), expectCore[i].real(), 1e-8);
    EXPECT_NEAR(coreX2C[i].imag(), expectCore[i].imag(), 1e-8);
  }

  // X2C Large picture change
  std::vector<std::complex<double>> expectUL = getULAOLinDep(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    EXPECT_NEAR(UL[i].real(), expectUL[i].real(), 1e-8);
    EXPECT_NEAR(UL[i].imag(), expectUL[i].imag(), 1e-8);
  }

  // X2C small picture change
  std::vector<std::complex<double>> expectUS = getUSAOLinDep(4*nbsq);

  for(auto i = 0; i < 4*nbsq; i++) {
    std::cout << "i: " << i << std::endl;
    EXPECT_NEAR(US[i].real(), expectUS[i].real(), 1e-8);
    EXPECT_NEAR(US[i].imag(), expectUS[i].imag(), 1e-8);
  }

}
