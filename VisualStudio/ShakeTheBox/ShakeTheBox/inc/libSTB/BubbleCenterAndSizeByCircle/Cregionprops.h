//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  regionprops.h
//
//  Code generation for function 'regionprops'
//


#ifndef CREGIONPROPS_H
#define CREGIONPROPS_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "BubbleCenterAndSizeByCircle_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void Cregionprops(const coder::array<bool, 2U> &varargin_1, const coder::
  array<double, 2U> &varargin_2, coder::array<struct_TT, 1U> &outstats);

#endif

// End of code generation (regionprops.h)
