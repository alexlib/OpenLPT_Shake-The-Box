//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// chaccum.cpp
//
// Code generation for function 'chaccum'
//

// Include files
#include "chaccum.h"
#include "colon.h"
#include "conv2AXPYValidCMP.h"
#include "imfilter.h"
#include "minOrMax.h"
#include "multithresh.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Declarations
static float rt_hypotf_snf(float u0, float u1);

// Function Definitions
static float rt_hypotf_snf(float u0, float u1)
{
  float a;
  float y;
  a = std::abs(u0);
  y = std::abs(u1);
  if (a < y) {
    a /= y;
    y *= std::sqrt(a * a + 1.0F);
  } else if (a > y) {
    y /= a;
    y = a * std::sqrt(y * y + 1.0F);
  } else if (!std::isnan(y)) {
    y = a * 1.41421354F;
  }
  return y;
}

namespace coder {
void chaccum(const ::coder::array<unsigned char, 2U> &varargin_1,
             const double varargin_2_data[], const int varargin_2_size[2],
             ::coder::array<creal_T, 2U> &accumMatrix,
             ::coder::array<float, 2U> &gradientImg)
{
  array<creal_T, 2U> r1;
  array<double, 2U> b_a;
  array<double, 2U> r;
  array<double, 2U> radiusRange;
  array<double, 2U> w;
  array<double, 2U> w0;
  array<double, 1U> Ex_chunk;
  array<float, 2U> A;
  array<float, 2U> Gx;
  array<float, 2U> Gy;
  array<float, 2U> a;
  array<float, 1U> b_gradientImg;
  array<int, 1U> Ex;
  array<int, 1U> Ey;
  array<int, 1U> idxE;
  array<int, 1U> idxE_chunk;
  array<int, 1U> xckeep;
  array<int, 1U> yckeep;
  array<bool, 2U> inside;
  array<bool, 1U> rows_to_keep;
  double pad[2];
  double radiusRangeIn_data[2];
  int csz_idx_0;
  int csz_idx_1;
  int i;
  int i1;
  int i2;
  int loop_ub;
  int radiusRangeIn_size_idx_1;
  bool exitg1;
  bool y;
  if (varargin_2_size[1] == 2) {
    if (varargin_2_data[0] == varargin_2_data[1]) {
      radiusRangeIn_size_idx_1 = 1;
      radiusRangeIn_data[0] = varargin_2_data[0];
    } else {
      radiusRangeIn_size_idx_1 = 2;
      for (i = 0; i < 2; i++) {
        radiusRangeIn_data[i] = varargin_2_data[i];
      }
    }
  } else {
    radiusRangeIn_size_idx_1 = 1;
    radiusRangeIn_data[0] = varargin_2_data[0];
  }
  rows_to_keep.set_size(varargin_1.size(0) * varargin_1.size(1));
  csz_idx_1 = varargin_1.size(0) * varargin_1.size(1);
  for (i = 0; i < csz_idx_1; i++) {
    rows_to_keep[i] = (varargin_1[i] == varargin_1[0]);
  }
  y = true;
  csz_idx_0 = 1;
  exitg1 = false;
  while ((!exitg1) && (csz_idx_0 <= rows_to_keep.size(0))) {
    if (!rows_to_keep[csz_idx_0 - 1]) {
      y = false;
      exitg1 = true;
    } else {
      csz_idx_0++;
    }
  }
  if (y) {
    accumMatrix.set_size(varargin_1.size(0), varargin_1.size(1));
    csz_idx_1 = varargin_1.size(0) * varargin_1.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      accumMatrix[i].re = 0.0;
      accumMatrix[i].im = 0.0;
    }
    gradientImg.set_size(varargin_1.size(0), varargin_1.size(1));
    csz_idx_1 = varargin_1.size(0) * varargin_1.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      gradientImg[i] = 0.0F;
    }
  } else {
    double xcStep;
    float Gmax;
    float edgeThresh;
    int M;
    int idx;
    int nx;
    A.set_size(varargin_1.size(0), varargin_1.size(1));
    csz_idx_1 = varargin_1.size(0) * varargin_1.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      A[i] = static_cast<float>(varargin_1[i]) / 255.0F;
    }
    pad[0] = 1.0;
    pad[1] = 1.0;
    if ((A.size(0) == 0) || (A.size(1) == 0)) {
      Gx.set_size(A.size(0), A.size(1));
      csz_idx_1 = A.size(0) * A.size(1);
      for (i = 0; i < csz_idx_1; i++) {
        Gx[i] = A[i];
      }
    } else {
      padImage_outSize(A, pad, a);
      b_a.set_size(a.size(0), a.size(1));
      csz_idx_1 = a.size(0) * a.size(1);
      for (i = 0; i < csz_idx_1; i++) {
        b_a[i] = a[i];
      }
      internal::conv2AXPYValidCMP(b_a, w);
      Gx.set_size(w.size(0), w.size(1));
      csz_idx_1 = w.size(0) * w.size(1);
      for (i = 0; i < csz_idx_1; i++) {
        Gx[i] = static_cast<float>(w[i]);
      }
    }
    pad[0] = 1.0;
    pad[1] = 1.0;
    if ((A.size(0) == 0) || (A.size(1) == 0)) {
      Gy.set_size(A.size(0), A.size(1));
      csz_idx_1 = A.size(0) * A.size(1);
      for (i = 0; i < csz_idx_1; i++) {
        Gy[i] = A[i];
      }
    } else {
      padImage_outSize(A, pad, a);
      b_a.set_size(a.size(0), a.size(1));
      csz_idx_1 = a.size(0) * a.size(1);
      for (i = 0; i < csz_idx_1; i++) {
        b_a[i] = a[i];
      }
      internal::b_conv2AXPYValidCMP(b_a, w);
      Gy.set_size(w.size(0), w.size(1));
      csz_idx_1 = w.size(0) * w.size(1);
      for (i = 0; i < csz_idx_1; i++) {
        Gy[i] = static_cast<float>(w[i]);
      }
    }
    if (Gx.size(0) <= Gy.size(0)) {
      csz_idx_0 = Gx.size(0);
    } else {
      csz_idx_0 = Gy.size(0);
    }
    if (Gx.size(1) <= Gy.size(1)) {
      csz_idx_1 = Gx.size(1);
    } else {
      csz_idx_1 = Gy.size(1);
    }
    gradientImg.set_size(csz_idx_0, csz_idx_1);
    nx = csz_idx_0 * csz_idx_1;
    for (csz_idx_0 = 0; csz_idx_0 < nx; csz_idx_0++) {
      gradientImg[csz_idx_0] = rt_hypotf_snf(Gx[csz_idx_0], Gy[csz_idx_0]);
    }
    csz_idx_0 = gradientImg.size(0) * gradientImg.size(1);
    b_gradientImg = gradientImg.reshape(csz_idx_0);
    Gmax = internal::maximum(b_gradientImg);
    a.set_size(gradientImg.size(0), gradientImg.size(1));
    csz_idx_1 = gradientImg.size(0) * gradientImg.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      a[i] = gradientImg[i] / Gmax;
    }
    edgeThresh = multithresh(a);
    Gmax *= edgeThresh;
    inside.set_size(gradientImg.size(0), gradientImg.size(1));
    csz_idx_1 = gradientImg.size(0) * gradientImg.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      inside[i] = (gradientImg[i] > Gmax);
    }
    nx = inside.size(0) * inside.size(1);
    if (nx == 0) {
      xckeep.set_size(0);
      yckeep.set_size(0);
    } else {
      idx = 0;
      xckeep.set_size(nx);
      yckeep.set_size(nx);
      csz_idx_0 = 1;
      csz_idx_1 = 1;
      exitg1 = false;
      while ((!exitg1) && (csz_idx_1 <= inside.size(1))) {
        bool guard1{false};
        guard1 = false;
        if (inside[(csz_idx_0 + inside.size(0) * (csz_idx_1 - 1)) - 1]) {
          idx++;
          xckeep[idx - 1] = csz_idx_0;
          yckeep[idx - 1] = csz_idx_1;
          if (idx >= nx) {
            exitg1 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
        if (guard1) {
          csz_idx_0++;
          if (csz_idx_0 > inside.size(0)) {
            csz_idx_0 = 1;
            csz_idx_1++;
          }
        }
      }
      if (nx == 1) {
        if (idx == 0) {
          xckeep.set_size(0);
          yckeep.set_size(0);
        }
      } else {
        if (1 > idx) {
          i = 0;
        } else {
          i = idx;
        }
        xckeep.set_size(i);
        if (1 > idx) {
          idx = 0;
        }
        yckeep.set_size(idx);
      }
    }
    Ey.set_size(xckeep.size(0));
    csz_idx_1 = xckeep.size(0);
    for (i = 0; i < csz_idx_1; i++) {
      Ey[i] = xckeep[i];
    }
    Ex.set_size(yckeep.size(0));
    csz_idx_1 = yckeep.size(0);
    for (i = 0; i < csz_idx_1; i++) {
      Ex[i] = yckeep[i];
    }
    csz_idx_0 = gradientImg.size(0);
    idxE.set_size(Ey.size(0));
    csz_idx_1 = Ey.size(0);
    for (i = 0; i < csz_idx_1; i++) {
      idxE[i] = Ey[i] + csz_idx_0 * (Ex[i] - 1);
    }
    if (radiusRangeIn_size_idx_1 > 1) {
      if (std::isnan(radiusRangeIn_data[0]) ||
          std::isnan(radiusRangeIn_data[1])) {
        radiusRange.set_size(1, 1);
        radiusRange[0] = rtNaN;
      } else if (radiusRangeIn_data[1] < radiusRangeIn_data[0]) {
        radiusRange.set_size(1, 0);
      } else if ((std::isinf(radiusRangeIn_data[0]) ||
                  std::isinf(radiusRangeIn_data[1])) &&
                 (radiusRangeIn_data[0] == radiusRangeIn_data[1])) {
        radiusRange.set_size(1, 1);
        radiusRange[0] = rtNaN;
      } else {
        eml_float_colon(radiusRangeIn_data[0], radiusRangeIn_data[1],
                        radiusRange);
      }
    } else {
      radiusRange.set_size(1, 1);
      radiusRange[0] = radiusRangeIn_data[0];
    }
    w0.set_size(1, radiusRange.size(1));
    csz_idx_1 = radiusRange.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      w0[i] = 1.0 / (6.2831853071795862 * radiusRange[i]);
    }
    xcStep = std::floor(1.0E+6 / static_cast<double>(radiusRange.size(1)));
    radiusRangeIn_size_idx_1 = A.size(1);
    M = A.size(0);
    accumMatrix.set_size(A.size(0), A.size(1));
    csz_idx_1 = A.size(0) * A.size(1);
    for (i = 0; i < csz_idx_1; i++) {
      accumMatrix[i].re = 0.0;
      accumMatrix[i].im = 0.0;
    }
    i = static_cast<int>((static_cast<double>(Ex.size(0)) + (xcStep - 1.0)) /
                         xcStep);
    if (0 <= i - 1) {
      i1 = radiusRange.size(1);
      i2 = radiusRange.size(1);
      pad[0] = static_cast<unsigned int>(A.size(0));
      pad[1] = static_cast<unsigned int>(A.size(1));
      loop_ub = A.size(0) * A.size(1);
    }
    for (int b_i{0}; b_i < i; b_i++) {
      double c_i;
      c_i = static_cast<double>(b_i) * xcStep + 1.0;
      csz_idx_0 = static_cast<int>(
          std::fmin((c_i + xcStep) - 1.0, static_cast<double>(Ex.size(0))));
      if (std::isnan(c_i)) {
        csz_idx_1 = 1;
      } else if (csz_idx_0 < c_i) {
        csz_idx_1 = 0;
      } else if (c_i == c_i) {
        csz_idx_1 = static_cast<int>(static_cast<double>(csz_idx_0) - c_i) + 1;
      } else {
        b_eml_float_colon(c_i, static_cast<double>(csz_idx_0), r);
        csz_idx_1 = r.size(1);
      }
      Ex_chunk.set_size(csz_idx_1);
      xckeep.set_size(csz_idx_1);
      idxE_chunk.set_size(csz_idx_1);
      for (idx = 0; idx < csz_idx_1; idx++) {
        csz_idx_0 = (static_cast<int>(c_i) + idx) - 1;
        Ex_chunk[idx] = Ex[csz_idx_0];
        xckeep[idx] = Ey[csz_idx_0];
        idxE_chunk[idx] = idxE[csz_idx_0];
      }
      a.set_size(idxE_chunk.size(0), radiusRange.size(1));
      A.set_size(idxE_chunk.size(0), radiusRange.size(1));
      w.set_size(idxE_chunk.size(0), radiusRange.size(1));
      inside.set_size(idxE_chunk.size(0), radiusRange.size(1));
      rows_to_keep.set_size(idxE_chunk.size(0));
      csz_idx_1 = idxE_chunk.size(0);
      for (idx = 0; idx < csz_idx_1; idx++) {
        rows_to_keep[idx] = false;
      }
      for (nx = 0; nx < i1; nx++) {
        idx = idxE_chunk.size(0);
        for (csz_idx_1 = 0; csz_idx_1 < idx; csz_idx_1++) {
          c_i = radiusRange[nx];
          Gmax = static_cast<float>(Ex_chunk[csz_idx_1]) +
                 static_cast<float>(-c_i) *
                     (Gx[idxE_chunk[csz_idx_1] - 1] /
                      gradientImg[idxE_chunk[csz_idx_1] - 1]);
          if (Gmax > 0.0F) {
            edgeThresh = Gmax + 0.5F;
          } else if (Gmax < 0.0F) {
            edgeThresh = Gmax - 0.5F;
          } else {
            edgeThresh = 0.0F;
          }
          a[csz_idx_1 + a.size(0) * nx] = edgeThresh;
          Gmax = static_cast<float>(xckeep[csz_idx_1]) +
                 static_cast<float>(-c_i) *
                     (Gy[idxE_chunk[csz_idx_1] - 1] /
                      gradientImg[idxE_chunk[csz_idx_1] - 1]);
          if (Gmax > 0.0F) {
            Gmax += 0.5F;
          } else if (Gmax < 0.0F) {
            Gmax -= 0.5F;
          } else {
            Gmax = 0.0F;
          }
          A[csz_idx_1 + A.size(0) * nx] = Gmax;
          w[csz_idx_1 + w.size(0) * nx] = w0[nx];
          y = ((edgeThresh >= 1.0F) &&
               (static_cast<double>(edgeThresh) <= radiusRangeIn_size_idx_1) &&
               (Gmax >= 1.0F) && (static_cast<double>(Gmax) < M));
          inside[csz_idx_1 + inside.size(0) * nx] = y;
          if (y) {
            rows_to_keep[csz_idx_1] = true;
          }
        }
      }
      xckeep.set_size(a.size(0) * a.size(1));
      yckeep.set_size(A.size(0) * A.size(1));
      Ex_chunk.set_size(w.size(0) * w.size(1));
      csz_idx_0 = -1;
      for (nx = 0; nx < i2; nx++) {
        idx = idxE_chunk.size(0);
        for (csz_idx_1 = 0; csz_idx_1 < idx; csz_idx_1++) {
          if (rows_to_keep[csz_idx_1] &&
              inside[csz_idx_1 + inside.size(0) * nx]) {
            csz_idx_0++;
            xckeep[csz_idx_0] = static_cast<int>(a[csz_idx_1 + a.size(0) * nx]);
            yckeep[csz_idx_0] = static_cast<int>(A[csz_idx_1 + A.size(0) * nx]);
            Ex_chunk[csz_idx_0] = w[csz_idx_1 + w.size(0) * nx];
          }
        }
      }
      r1.set_size(static_cast<int>(pad[0]), static_cast<int>(pad[1]));
      for (idx = 0; idx < loop_ub; idx++) {
        r1[idx].re = 0.0;
        r1[idx].im = 0.0;
      }
      for (idx = 0; idx <= csz_idx_0; idx++) {
        r1[(yckeep[idx] + r1.size(0) * (xckeep[idx] - 1)) - 1].re =
            r1[(yckeep[idx] + r1.size(0) * (xckeep[idx] - 1)) - 1].re +
            Ex_chunk[idx];
      }
      csz_idx_1 = accumMatrix.size(0) * accumMatrix.size(1);
      for (idx = 0; idx < csz_idx_1; idx++) {
        accumMatrix[idx].re = accumMatrix[idx].re + r1[idx].re;
        accumMatrix[idx].im = accumMatrix[idx].im + r1[idx].im;
      }
    }
  }
}

} // namespace coder

// End of code generation (chaccum.cpp)
