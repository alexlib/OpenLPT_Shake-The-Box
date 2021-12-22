//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// chradii.h
//
// Code generation for function 'chradii'
//

#ifndef CHRADII_H
#define CHRADII_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void chradii(const ::coder::array<double, 2U> &varargin_1,
             const ::coder::array<float, 2U> &varargin_2,
             const double varargin_3_data[],
             ::coder::array<double, 1U> &r_estimated);

}

#endif
// End of code generation (chradii.h)
