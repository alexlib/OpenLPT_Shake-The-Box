//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// imhmax.h
//
// Code generation for function 'imhmax'
//

#ifndef IMHMAX_H
#define IMHMAX_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void imhmax(const ::coder::array<double, 2U> &b_I, double H,
            ::coder::array<double, 2U> &J);

}

#endif
// End of code generation (imhmax.h)
