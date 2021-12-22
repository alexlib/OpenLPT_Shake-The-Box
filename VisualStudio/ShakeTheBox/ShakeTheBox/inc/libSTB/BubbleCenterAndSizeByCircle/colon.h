//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// colon.h
//
// Code generation for function 'colon'
//

#ifndef COLON_H
#define COLON_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void b_eml_float_colon(double a, double b, ::coder::array<double, 2U> &y);

void eml_float_colon(double a, double b, ::coder::array<double, 2U> &y);

} // namespace coder

#endif
// End of code generation (colon.h)
