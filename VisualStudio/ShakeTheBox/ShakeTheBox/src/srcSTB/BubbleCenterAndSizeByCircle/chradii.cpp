//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// chradii.cpp
//
// Code generation for function 'chradii'
//

// Include files
#include "chradii.h"
#include "BubbleCenterAndSizeByCircle_rtwutil.h"
#include "colon.h"
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
void chradii(const ::coder::array<double, 2U> &varargin_1,
             const ::coder::array<float, 2U> &varargin_2,
             const double varargin_3_data[],
             ::coder::array<double, 1U> &r_estimated)
{
  array<creal_T, 1U> h;
  array<double, 2U> b_y;
  array<double, 2U> r;
  array<double, 2U> xx;
  array<double, 2U> y;
  array<double, 2U> yy;
  array<double, 1U> b_r;
  array<double, 1U> bins;
  array<float, 2U> b_varargin_2;
  array<float, 1U> gradientImg;
  array<bool, 1U> c_r;
  array<bool, 1U> x;
  creal_T ex;
  double b_varargin_3;
  double varargin_3;
  int i;
  int loop_ub;
  int ny;
  r_estimated.set_size(varargin_1.size(0));
  loop_ub = varargin_1.size(0);
  for (i = 0; i < loop_ub; i++) {
    r_estimated[i] = 0.0;
  }
  i = varargin_1.size(0);
  if (0 <= varargin_1.size(0) - 1) {
    varargin_3 = varargin_3_data[0];
    b_varargin_3 = varargin_3_data[1];
  }
  for (int k{0}; k < i; k++) {
    double bottom;
    double br;
    double d;
    double left;
    double right;
    double top;
    int b_i;
    int b_loop_ub;
    int i1;
    int i2;
    int j;
    int nx;
    bool exitg1;
    br = varargin_1[k];
    left = std::fmax(std::floor(br - varargin_3_data[1]), 1.0);
    right = std::fmin(std::ceil(br + varargin_3_data[1]),
                      static_cast<double>(varargin_2.size(1)));
    d = varargin_1[k + varargin_1.size(0)];
    top = std::fmax(std::floor(d - varargin_3_data[1]), 1.0);
    bottom = std::fmin(std::ceil(d + varargin_3_data[1]),
                       static_cast<double>(varargin_2.size(0)));
    if (top > bottom) {
      i1 = -1;
      j = -1;
    } else {
      i1 = static_cast<int>(top) - 2;
      j = static_cast<int>(bottom) - 1;
    }
    if (left > right) {
      i2 = -1;
      ny = -1;
    } else {
      i2 = static_cast<int>(left) - 2;
      ny = static_cast<int>(right) - 1;
    }
    right = (br - left) + 1.0;
    left = (d - top) + 1.0;
    loop_ub = ny - i2;
    if (loop_ub < 1) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, loop_ub);
      b_loop_ub = loop_ub - 1;
      for (ny = 0; ny <= b_loop_ub; ny++) {
        y[ny] = static_cast<double>(ny) + 1.0;
      }
    }
    b_loop_ub = j - i1;
    if (b_loop_ub < 1) {
      b_y.set_size(1, 0);
    } else {
      b_y.set_size(1, b_loop_ub);
      nx = b_loop_ub - 1;
      for (j = 0; j <= nx; j++) {
        b_y[j] = static_cast<double>(j) + 1.0;
      }
    }
    nx = y.size(1);
    ny = b_y.size(1);
    xx.set_size(b_y.size(1), y.size(1));
    yy.set_size(b_y.size(1), y.size(1));
    if ((y.size(1) != 0) && (b_y.size(1) != 0)) {
      for (j = 0; j < nx; j++) {
        for (b_i = 0; b_i < ny; b_i++) {
          xx[b_i + xx.size(0) * j] = y[j];
          yy[b_i + yy.size(0) * j] = b_y[b_i];
        }
      }
    }
    nx = xx.size(0) * xx.size(1);
    for (j = 0; j < nx; j++) {
      xx[j] = xx[j] - right;
    }
    nx = yy.size(0) * yy.size(1);
    for (j = 0; j < nx; j++) {
      yy[j] = yy[j] - left;
    }
    if (xx.size(0) <= yy.size(0)) {
      j = xx.size(0);
    } else {
      j = yy.size(0);
    }
    if (xx.size(1) <= yy.size(1)) {
      ny = xx.size(1);
    } else {
      ny = yy.size(1);
    }
    r.set_size(j, ny);
    nx = j * ny;
    for (b_i = 0; b_i < nx; b_i++) {
      r[b_i] = rt_hypotd_snf(xx[b_i], yy[b_i]);
    }
    nx = r.size(0) * r.size(1);
    for (b_i = 0; b_i < nx; b_i++) {
      r[b_i] = std::round(r[b_i]);
    }
    b_r.set_size(r.size(0) * r.size(1));
    nx = r.size(0) * r.size(1);
    for (j = 0; j < nx; j++) {
      b_r[j] = r[j];
    }
    b_varargin_2.set_size(b_loop_ub, loop_ub);
    for (j = 0; j < loop_ub; j++) {
      for (ny = 0; ny < b_loop_ub; ny++) {
        b_varargin_2[ny + b_varargin_2.size(0) * j] =
            varargin_2[((i1 + ny) + varargin_2.size(0) * ((i2 + j) + 1)) + 1];
      }
    }
    ny = b_loop_ub * loop_ub;
    gradientImg.set_size(ny);
    for (i1 = 0; i1 < ny; i1++) {
      gradientImg[i1] = b_varargin_2[i1];
    }
    x.set_size(r.size(0) * r.size(1));
    loop_ub = r.size(0) * r.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      x[i1] = (r[i1] >= varargin_3);
    }
    c_r.set_size(r.size(0) * r.size(1));
    loop_ub = r.size(0) * r.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_r[i1] = (r[i1] <= b_varargin_3);
    }
    ny = x.size(0) - 1;
    j = 0;
    for (b_i = 0; b_i <= ny; b_i++) {
      if (x[b_i] && c_r[b_i]) {
        j++;
      }
    }
    nx = 0;
    for (b_i = 0; b_i <= ny; b_i++) {
      if (x[b_i] && c_r[b_i]) {
        gradientImg[nx] = gradientImg[b_i];
        nx++;
      }
    }
    gradientImg.set_size(j);
    ny = x.size(0) - 1;
    j = 0;
    for (b_i = 0; b_i <= ny; b_i++) {
      if (x[b_i] && c_r[b_i]) {
        j++;
      }
    }
    nx = 0;
    for (b_i = 0; b_i <= ny; b_i++) {
      if (x[b_i] && c_r[b_i]) {
        b_r[nx] = b_r[b_i];
        nx++;
      }
    }
    b_r.set_size(j);
    if (j <= 2) {
      if (j == 1) {
        right = b_r[0];
      } else {
        right = b_r[j - 1];
        if ((!(b_r[0] > right)) &&
            ((!std::isnan(b_r[0])) || std::isnan(right))) {
          right = b_r[0];
        }
      }
    } else {
      if (!std::isnan(b_r[0])) {
        nx = 1;
      } else {
        nx = 0;
        b_i = 2;
        exitg1 = false;
        while ((!exitg1) && (b_i <= j)) {
          if (!std::isnan(b_r[b_i - 1])) {
            nx = b_i;
            exitg1 = true;
          } else {
            b_i++;
          }
        }
      }
      if (nx == 0) {
        right = b_r[0];
      } else {
        right = b_r[nx - 1];
        i1 = nx + 1;
        for (b_i = i1; b_i <= j; b_i++) {
          br = b_r[b_i - 1];
          if (right > br) {
            right = br;
          }
        }
      }
    }
    left = internal::maximum(b_r);
    if (std::isnan(right) || std::isnan(left)) {
      y.set_size(1, 1);
      y[0] = rtNaN;
    } else if (left < right) {
      y.set_size(1, 0);
    } else if ((std::isinf(right) || std::isinf(left)) && (right == left)) {
      y.set_size(1, 1);
      y[0] = rtNaN;
    } else if (std::floor(right) == right) {
      loop_ub = static_cast<int>(std::floor(left - right));
      y.set_size(1, loop_ub + 1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        y[i1] = right + static_cast<double>(i1);
      }
    } else {
      b_eml_float_colon(right, left, y);
    }
    bins.set_size(y.size(1));
    loop_ub = y.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      bins[i1] = y[i1];
    }
    for (i1 = 0; i1 < j; i1++) {
      b_r[i1] = (b_r[i1] - bins[0]) + 1.0;
    }
    ny = b_r.size(0);
    if (b_r.size(0) <= 2) {
      if (b_r.size(0) == 1) {
        right = b_r[0];
      } else if ((b_r[0] < b_r[b_r.size(0) - 1]) ||
                 (std::isnan(b_r[0]) && (!std::isnan(b_r[b_r.size(0) - 1])))) {
        right = b_r[b_r.size(0) - 1];
      } else {
        right = b_r[0];
      }
    } else {
      if (!std::isnan(b_r[0])) {
        nx = 1;
      } else {
        nx = 0;
        b_i = 2;
        exitg1 = false;
        while ((!exitg1) && (b_i <= ny)) {
          if (!std::isnan(b_r[b_i - 1])) {
            nx = b_i;
            exitg1 = true;
          } else {
            b_i++;
          }
        }
      }
      if (nx == 0) {
        right = b_r[0];
      } else {
        right = b_r[nx - 1];
        i1 = nx + 1;
        for (b_i = i1; b_i <= ny; b_i++) {
          br = b_r[b_i - 1];
          if (right < br) {
            right = br;
          }
        }
      }
    }
    ny = static_cast<int>(right);
    h.set_size(ny);
    for (i1 = 0; i1 < ny; i1++) {
      h[i1].re = 0.0;
      h[i1].im = 0.0;
    }
    i1 = b_r.size(0);
    for (nx = 0; nx < i1; nx++) {
      ny = static_cast<int>(b_r[nx]) - 1;
      h[ny].re = static_cast<float>(h[ny].re) + gradientImg[nx];
      h[ny].im = static_cast<float>(h[ny].im);
    }
    i1 = h.size(0);
    for (nx = 0; nx < i1; nx++) {
      left = h[nx].re;
      bottom = h[nx].im;
      br = 6.2831853071795862 * bins[nx];
      if (bottom == 0.0) {
        right = left / br;
        left = 0.0;
      } else if (left == 0.0) {
        right = 0.0;
        left = bottom / br;
      } else {
        right = left / br;
        left = bottom / br;
      }
      h[nx].re = right;
      h[nx].im = left;
    }
    internal::maximum(h, &ex, &ny);
    r_estimated[k] = bins[ny - 1];
  }
}

} // namespace coder

// End of code generation (chradii.cpp)
