//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// relop.cpp
//
// Code generation for function 'relop'
//

// Include files
#include "relop.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <math.h>

// Function Definitions
namespace coder {
namespace internal {
bool iseq(double x, double y)
{
  double absx;
  int exponent;
  bool p;
  absx = std::abs(y / 2.0);
  if ((!std::isinf(absx)) && (!std::isnan(absx))) {
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = std::ldexp(1.0, exponent - 53);
    }
  } else {
    absx = rtNaN;
  }
  if ((std::abs(y - x) < absx) ||
      (std::isinf(x) && std::isinf(y) && ((x > 0.0) == (y > 0.0)))) {
    p = true;
  } else {
    p = false;
  }
  return p;
}

} // namespace internal
} // namespace coder

// End of code generation (relop.cpp)
