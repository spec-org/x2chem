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
