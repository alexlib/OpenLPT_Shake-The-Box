//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// minOrMax.cpp
//
// Code generation for function 'minOrMax'
//

// Include files
#include "minOrMax.h"
#include "BubbleCenterAndSizeByCircle_rtwutil.h"
#include "relop.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "rt_defines.h"
#include <cmath>

// Function Declarations
static double rt_atan2d_snf(double u0, double u1);

// Function Definitions
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int b_u0;
    int b_u1;
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }
    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }
    y = std::atan2(static_cast<double>(b_u0), static_cast<double>(b_u1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }
  return y;
}

namespace coder {
namespace internal {
float maximum(const ::coder::array<float, 1U> &x)
{
  float ex;
  int last;
  last = x.size(0);
  if (x.size(0) <= 2) {
    if (x.size(0) == 1) {
      ex = x[0];
    } else if ((x[0] < x[x.size(0) - 1]) ||
               (std::isnan(x[0]) && (!std::isnan(x[x.size(0) - 1])))) {
      ex = x[x.size(0) - 1];
    } else {
      ex = x[0];
    }
  } else {
    int idx;
    int k;
    if (!std::isnan(x[0])) {
      idx = 1;
    } else {
      bool exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x[0];
    } else {
      ex = x[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        float f;
        f = x[k - 1];
        if (ex < f) {
          ex = f;
        }
      }
    }
  }
  return ex;
}

double maximum(const double x_data[], const int x_size[2])
{
  double ex;
  int last;
  last = x_size[1];
  if (x_size[1] <= 2) {
    if (x_size[1] == 1) {
      ex = x_data[0];
    } else {
      ex = x_data[x_size[1] - 1];
      if ((!(x_data[0] < ex)) && ((!std::isnan(x_data[0])) || std::isnan(ex))) {
        ex = x_data[0];
      }
    }
  } else {
    int idx;
    int k;
    if (!std::isnan(x_data[0])) {
      idx = 1;
    } else {
      bool exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x_data[0];
    } else {
      ex = x_data[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        double d;
        d = x_data[k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }
  return ex;
}

double maximum(const ::coder::array<double, 1U> &x)
{
  double ex;
  int last;
  last = x.size(0);
  if (x.size(0) <= 2) {
    if (x.size(0) == 1) {
      ex = x[0];
    } else if ((x[0] < x[x.size(0) - 1]) ||
               (std::isnan(x[0]) && (!std::isnan(x[x.size(0) - 1])))) {
      ex = x[x.size(0) - 1];
    } else {
      ex = x[0];
    }
  } else {
    int idx;
    int k;
    if (!std::isnan(x[0])) {
      idx = 1;
    } else {
      bool exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x[0];
    } else {
      ex = x[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        double d;
        d = x[k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }
  return ex;
}

void maximum(const ::coder::array<creal_T, 1U> &x, creal_T *ex, int *idx)
{
  int istop;
  istop = x.size(0);
  *idx = 1;
  *ex = x[0];
  for (int k{2}; k <= istop; k++) {
    creal_T dc;
    bool SCALEA;
    dc = x[k - 1];
    if (std::isnan(dc.re) || std::isnan(x[k - 1].im)) {
      SCALEA = false;
    } else if (std::isnan(ex->re) || std::isnan(ex->im)) {
      SCALEA = true;
    } else {
      double b_x;
      double bi;
      bool SCALEB;
      if ((std::abs(ex->re) > 8.9884656743115785E+307) ||
          (std::abs(ex->im) > 8.9884656743115785E+307)) {
        SCALEA = true;
      } else {
        SCALEA = false;
      }
      if ((std::abs(x[k - 1].re) > 8.9884656743115785E+307) ||
          (std::abs(x[k - 1].im) > 8.9884656743115785E+307)) {
        SCALEB = true;
      } else {
        SCALEB = false;
      }
      if (SCALEA || SCALEB) {
        b_x = rt_hypotd_snf(ex->re / 2.0, ex->im / 2.0);
        bi = rt_hypotd_snf(x[k - 1].re / 2.0, x[k - 1].im / 2.0);
      } else {
        b_x = rt_hypotd_snf(ex->re, ex->im);
        bi = rt_hypotd_snf(x[k - 1].re, x[k - 1].im);
      }
      if (iseq(b_x, bi)) {
        b_x = rt_atan2d_snf(ex->im, ex->re);
        bi = rt_atan2d_snf(x[k - 1].im, x[k - 1].re);
        if (iseq(b_x, bi)) {
          double br;
          br = x[k - 1].re;
          bi = x[k - 1].im;
          if (ex->re != br) {
            if (b_x >= 0.0) {
              b_x = br;
              bi = ex->re;
            } else {
              b_x = ex->re;
              bi = br;
            }
          } else if (ex->re < 0.0) {
            b_x = bi;
            bi = ex->im;
          } else {
            b_x = ex->im;
          }
          if (iseq(b_x, bi)) {
            b_x = 0.0;
            bi = 0.0;
          }
        }
      }
      SCALEA = (b_x < bi);
    }
    if (SCALEA) {
      *ex = dc;
      *idx = k;
    }
  }
}

} // namespace internal
} // namespace coder

// End of code generation (minOrMax.cpp)
