//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  CircleIdentifier.cpp
//
//  Code generation for function 'CircleIdentifier'
//


// Include files
#include "CircleIdentifier.h"
#include "BubbleCenterAndSizeByCircle_data.h"
#include "BubbleCenterAndSizeByCircle_internal_types.h"
#include "BubbleCenterAndSizeByCircle_rtwutil.h"
#include "NeighborhoodProcessor.h"
#include "chaccum.h"
#include "chradii.h"
#include "imhmax.h"
#include "medfilt2.h"
#include "regionprops.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>
#include <math.h>

#include "Frame.h"
#include "NumDataIO.h"

// Function Definitions
void CircleIdentifier::BubbleCenterAndSizeByCircle(double rmin, double rmax, double sense)
{
    coder::array<double, 2U> centers;
    coder::array<double, 2U> radii;
    coder::images::internal::coder::NeighborhoodProcessor np;
    coder::array<struct_T, 1U> s;
    coder::array<creal_T, 2U> accumMatrix;
    coder::array<double, 2U> Hd;
    coder::array<double, 2U> accumMatrixRe;
    coder::array<double, 2U> b_accumMatrixRe;
    coder::array<double, 2U> b_centers;
    coder::array<double, 2U> b_metric;
    coder::array<double, 2U> c_x;
    coder::array<double, 2U> metric;
    coder::array<double, 1U> d_x;
    coder::array<double, 1U> idx2Keep;
    coder::array<float, 2U> gradientImg;
    coder::array<int, 2U> iidx;
    coder::array<int, 1U> ii;
    coder::array<bool, 2U> b_x;
    coder::array<bool, 2U> bwpre;
    coder::array<bool, 2U> expl_temp_bw;
    coder::array<bool, 1U> x;
    b_struct_T expl_temp;
    double radiusRange_data[2];
    int b_np[2];
    int radiusRange_size[2];
    int exponent;
    int nrows;
    bool c_np[9];
    bool continuePropagation;
    bool exitg1;
    //      s = regionprops(img,'Centroid','MajorAxisLength','MinorAxisLength');
    if (rmin == rmax) {
        radiusRange_size[0] = 1;
        radiusRange_size[1] = 1;
        radiusRange_data[0] = rmin;
    }
    else {
        radiusRange_size[0] = 1;
        radiusRange_size[1] = 2;
        radiusRange_data[0] = rmin;
        radiusRange_data[1] = rmax;
    }
    centers.set_size(0, 0);
    radii.set_size(0, 0);
    coder::chaccum(img, radiusRange_data, radiusRange_size, accumMatrix,
        gradientImg);
    continuePropagation = false;
    nrows = 0;
    exitg1 = false;
    while ((!exitg1) &&
        (nrows + 1 <= accumMatrix.size(0) * accumMatrix.size(1))) {
        if (((accumMatrix[nrows].re == 0.0) && (accumMatrix[nrows].im == 0.0)) ||
            (std::isnan(accumMatrix[nrows].re) ||
                std::isnan(accumMatrix[nrows].im))) {
            nrows++;
        }
        else {
            continuePropagation = true;
            exitg1 = true;
        }
    }
    if (continuePropagation) {
        int i;
        int i1;
        int idx;
        int loop_ub;
        int nx;
        centers.set_size(0, 0);
        metric.set_size(0, 0);
        nx = accumMatrix.size(0) * accumMatrix.size(1);
        accumMatrixRe.set_size(accumMatrix.size(0), accumMatrix.size(1));
        for (nrows = 0; nrows < nx; nrows++) {
            accumMatrixRe[nrows] =
                rt_hypotd_snf(accumMatrix[nrows].re, accumMatrix[nrows].im);
        }
        x.set_size(accumMatrixRe.size(0) * accumMatrixRe.size(1));
        loop_ub = accumMatrixRe.size(0) * accumMatrixRe.size(1);
        for (i = 0; i < loop_ub; i++) {
            x[i] = (accumMatrixRe[i] == accumMatrixRe[0]);
        }
        continuePropagation = true;
        nrows = 1;
        exitg1 = false;
        while ((!exitg1) && (nrows <= x.size(0))) {
            if (!x[nrows - 1]) {
                continuePropagation = false;
                exitg1 = true;
            }
            else {
                nrows++;
            }
        }
        if (!continuePropagation) {
            double absx;
            if (static_cast<unsigned int>(accumMatrixRe.size(0)) >
                static_cast<unsigned int>(accumMatrixRe.size(1))) {
                i = accumMatrixRe.size(1);
            }
            else {
                i = accumMatrixRe.size(0);
            }
            if (i > 5) {
                if ((accumMatrixRe.size(0) == 0) || (accumMatrixRe.size(1) == 0)) {
                    Hd.set_size(accumMatrixRe.size(0), accumMatrixRe.size(1));
                    loop_ub = accumMatrixRe.size(0) * accumMatrixRe.size(1);
                    for (i = 0; i < loop_ub; i++) {
                        Hd[i] = accumMatrixRe[i];
                    }
                }
                else {
                    b_accumMatrixRe.set_size(accumMatrixRe.size(0),
                        accumMatrixRe.size(1));
                    loop_ub = accumMatrixRe.size(0) * accumMatrixRe.size(1) - 1;
                    for (i = 0; i <= loop_ub; i++) {
                        b_accumMatrixRe[i] = accumMatrixRe[i];
                    }
                    coder::images::internal::coder::optimized::medfilt2(b_accumMatrixRe,
                        Hd);
                }
            }
            else {
                Hd.set_size(accumMatrixRe.size(0), accumMatrixRe.size(1));
                loop_ub = accumMatrixRe.size(0) * accumMatrixRe.size(1);
                for (i = 0; i < loop_ub; i++) {
                    Hd[i] = accumMatrixRe[i];
                }
            }
            absx = std::abs(1.0 - sense);
            if ((!std::isinf(absx)) && (!std::isnan(absx))) {
                if (absx <= 2.2250738585072014E-308) {
                    absx = 4.94065645841247E-324;
                }
                else {
                    frexp(absx, &exponent);
                    absx = std::ldexp(1.0, exponent - 53);
                }
            }
            else {
                absx = rtNaN;
            }
            b_accumMatrixRe.set_size(Hd.size(0), Hd.size(1));
            loop_ub = Hd.size(0) * Hd.size(1) - 1;
            for (i = 0; i <= loop_ub; i++) {
                b_accumMatrixRe[i] = Hd[i];
            }
            coder::imhmax(b_accumMatrixRe, std::fmax((1.0 - sense) - absx, 0.0), Hd);
            np.ImageSize[0] = Hd.size(0);
            np.ImageSize[1] = Hd.size(1);
            np.Padding = 1.0;
            np.ProcessBorder = true;
            np.NeighborhoodCenter = 1.0;
            np.PadValue = 0.0;
            for (exponent = 0; exponent < 9; exponent++) {
                np.Neighborhood[exponent] = true;
                np.ImageNeighborLinearOffsets[exponent] = 0;
                np.NeighborLinearIndices[exponent] = 0;
            }
            std::memset(&np.NeighborSubscriptOffsets[0], 0, 18U * sizeof(int));
            expl_temp_bw.set_size(Hd.size(0), Hd.size(1));
            loop_ub = Hd.size(0) * Hd.size(1);
            for (i = 0; i < loop_ub; i++) {
                expl_temp_bw[i] = true;
            }
            continuePropagation = true;
            while (continuePropagation) {
                bool p;
                bwpre.set_size(expl_temp_bw.size(0), expl_temp_bw.size(1));
                loop_ub = expl_temp_bw.size(0) * expl_temp_bw.size(1);
                for (i = 0; i < loop_ub; i++) {
                    bwpre[i] = expl_temp_bw[i];
                }
                b_np[0] = np.ImageSize[0];
                b_np[1] = np.ImageSize[1];
                for (int i2{ 0 }; i2 < 9; i2++) {
                    c_np[i2] = np.Neighborhood[i2];
                }
                coder::images::internal::coder::NeighborhoodProcessor::
                    computeParameters(b_np, c_np, np.ImageNeighborLinearOffsets,
                        np.NeighborLinearIndices,
                        np.NeighborSubscriptOffsets, np.InteriorStart,
                        np.InteriorEnd);
                expl_temp.bw = expl_temp_bw;
                np.process2D(Hd, expl_temp_bw, &expl_temp);
                p = false;
                if ((bwpre.size(0) == expl_temp_bw.size(0)) &&
                    (bwpre.size(1) == expl_temp_bw.size(1))) {
                    p = true;
                }
                if (p && ((bwpre.size(0) != 0) && (bwpre.size(1) != 0)) &&
                    ((expl_temp_bw.size(0) != 0) && (expl_temp_bw.size(1) != 0))) {
                    nrows = 0;
                    exitg1 = false;
                    while ((!exitg1) &&
                        (nrows <= expl_temp_bw.size(0) * expl_temp_bw.size(1) - 1)) {
                        if (bwpre[nrows] != expl_temp_bw[nrows]) {
                            p = false;
                            exitg1 = true;
                        }
                        else {
                            nrows++;
                        }
                    }
                }
                continuePropagation = p;
                continuePropagation = !continuePropagation;
            }
            coder::regionprops(expl_temp_bw, accumMatrixRe, s);
            if (s.size(0) != 0) {
                centers.set_size(s.size(0), 2);
                i = s.size(0);
                for (idx = 0; idx < i; idx++) {
                    centers[idx] = s[idx].WeightedCentroid[0];
                    centers[idx + centers.size(0)] = s[idx].WeightedCentroid[1];
                }
                i = centers.size(0) - 1;
                i1 = static_cast<int>(
                    ((-1.0 - static_cast<double>(centers.size(0))) + 1.0) / -1.0);
                for (idx = 0; idx < i1; idx++) {
                    nx = i - idx;
                    if (std::isnan(centers[nx]) ||
                        std::isnan(centers[nx + centers.size(0)])) {
                        c_x.set_size(centers.size(0), 2);
                        loop_ub = centers.size(0) * 2;
                        for (exponent = 0; exponent < loop_ub; exponent++) {
                            c_x[exponent] = centers[exponent];
                        }
                        nrows = centers.size(0) - 1;
                        for (loop_ub = 0; loop_ub < 2; loop_ub++) {
                            for (exponent = nx + 1; exponent <= nrows; exponent++) {
                                c_x[(exponent + c_x.size(0) * loop_ub) - 1] =
                                    c_x[exponent + c_x.size(0) * loop_ub];
                            }
                        }
                        if (1 > centers.size(0) - 1) {
                            loop_ub = 0;
                        }
                        else {
                            loop_ub = centers.size(0) - 1;
                        }
                        for (exponent = 0; exponent < 2; exponent++) {
                            for (nx = 0; nx < loop_ub; nx++) {
                                c_x[nx + loop_ub * exponent] = c_x[nx + c_x.size(0) * exponent];
                            }
                        }
                        c_x.set_size(loop_ub, 2);
                        centers.set_size(loop_ub, 2);
                        nrows = loop_ub << 1;
                        for (exponent = 0; exponent < nrows; exponent++) {
                            centers[exponent] = c_x[exponent];
                        }
                    }
                }
                if (centers.size(0) != 0) {
                    loop_ub = centers.size(0);
                    idx2Keep.set_size(centers.size(0));
                    for (i = 0; i < loop_ub; i++) {
                        idx2Keep[i] = centers[i + centers.size(0)];
                    }
                    i = centers.size(0) - 1;
                    for (nrows = 0; nrows <= i; nrows++) {
                        idx2Keep[nrows] = std::round(idx2Keep[nrows]);
                    }
                    loop_ub = centers.size(0);
                    d_x.set_size(centers.size(0));
                    for (i = 0; i < loop_ub; i++) {
                        d_x[i] = centers[i];
                    }
                    i = centers.size(0) - 1;
                    for (nrows = 0; nrows <= i; nrows++) {
                        d_x[nrows] = std::round(d_x[nrows]);
                    }
                    nrows = Hd.size(0);
                    loop_ub = idx2Keep.size(0);
                    for (i = 0; i < loop_ub; i++) {
                        idx2Keep[i] = static_cast<int>(idx2Keep[i]) +
                            nrows * (static_cast<int>(d_x[i]) - 1);
                    }
                    b_metric.set_size(idx2Keep.size(0), 1);
                    loop_ub = idx2Keep.size(0);
                    for (i = 0; i < loop_ub; i++) {
                        b_metric[i] = Hd[static_cast<int>(idx2Keep[i]) - 1];
                    }
                    coder::internal::sort(b_metric, iidx);
                    metric.set_size(b_metric.size(0), 1);
                    loop_ub = b_metric.size(0);
                    for (i = 0; i < loop_ub; i++) {
                        metric[i] = b_metric[i];
                    }
                    b_centers.set_size(iidx.size(0), 2);
                    loop_ub = iidx.size(0);
                    for (i = 0; i < 2; i++) {
                        for (i1 = 0; i1 < loop_ub; i1++) {
                            b_centers[i1 + b_centers.size(0) * i] =
                                centers[(iidx[i1] + centers.size(0) * i) - 1];
                        }
                    }
                    centers.set_size(b_centers.size(0), 2);
                    loop_ub = b_centers.size(0) * 2;
                    for (i = 0; i < loop_ub; i++) {
                        centers[i] = b_centers[i];
                    }
                }
            }
        }
        if ((centers.size(0) != 0) && (centers.size(1) != 0)) {
            b_x.set_size(metric.size(0), metric.size(1));
            loop_ub = metric.size(0) * metric.size(1);
            for (i = 0; i < loop_ub; i++) {
                b_x[i] = (metric[i] >= 1.0 - sense);
            }
            nx = b_x.size(0) * b_x.size(1);
            idx = 0;
            ii.set_size(nx);
            nrows = 0;
            exitg1 = false;
            while ((!exitg1) && (nrows <= nx - 1)) {
                if (b_x[nrows]) {
                    idx++;
                    ii[idx - 1] = nrows + 1;
                    if (idx >= nx) {
                        exitg1 = true;
                    }
                    else {
                        nrows++;
                    }
                }
                else {
                    nrows++;
                }
            }
            if (nx == 1) {
                if (idx == 0) {
                    ii.set_size(0);
                }
            }
            else {
                if (1 > idx) {
                    idx = 0;
                }
                ii.set_size(idx);
            }
            idx2Keep.set_size(ii.size(0));
            loop_ub = ii.size(0);
            for (i = 0; i < loop_ub; i++) {
                idx2Keep[i] = ii[i];
            }
            nrows = centers.size(1) - 1;
            b_centers.set_size(idx2Keep.size(0), centers.size(1));
            for (i = 0; i <= nrows; i++) {
                loop_ub = idx2Keep.size(0);
                for (i1 = 0; i1 < loop_ub; i1++) {
                    b_centers[i1 + b_centers.size(0) * i] =
                        centers[(static_cast<int>(idx2Keep[i1]) + centers.size(0) * i) -
                        1];
                }
            }
            centers.set_size(b_centers.size(0), b_centers.size(1));
            loop_ub = b_centers.size(0) * b_centers.size(1);
            for (i = 0; i < loop_ub; i++) {
                centers[i] = b_centers[i];
            }
            if (idx2Keep.size(0) == 0) {
                centers.set_size(0, 0);
            }
            else if (radiusRange_size[1] == 1) {
                b_metric.set_size(idx2Keep.size(0), 1);
                nrows = idx2Keep.size(0);
                for (nx = 0; nx < nrows; nx++) {
                    b_metric[nx] = radiusRange_data[0];
                }
                radii.set_size(b_metric.size(0), 1);
                loop_ub = b_metric.size(0);
                for (i = 0; i < loop_ub; i++) {
                    radii[i] = b_metric[i];
                }
            }
            else {
                coder::chradii(centers, gradientImg, radiusRange_data, idx2Keep);
                radii.set_size(idx2Keep.size(0), 1);
                loop_ub = idx2Keep.size(0);
                for (i = 0; i < loop_ub; i++) {
                    radii[i] = idx2Keep[i];
                }
            }
        }
    }
    int n_center = *centers.size();
    for (int i = 0; i < n_center; i++) {
        Position bubble(centers[i] - 2, centers[n_center + i] - 2, 0, radii[i]);
        center.push_back(bubble);
        radius.push_back(radii[i]);
    }
}
Frame CircleIdentifier::CreateFrame() {
    return Frame(center);
}

CircleIdentifier::~CircleIdentifier()
{
  omp_destroy_nest_lock(&emlrtNestLockGlobal);
}

CircleIdentifier::CircleIdentifier(int**& p, int rows, int cols) : img_p(p), rows(rows), cols(cols)
{
  //rt_InitInfAndNaN();
  omp_init_nest_lock(&emlrtNestLockGlobal);

  img.set_size(rows, cols);
  // Loop over the array to initialize each element.
  for (int idx0 = 0; idx0 < img.size(0); idx0++) {
      for (int idx1 = 0; idx1 < img.size(1); idx1++) {
          // Set the value of the array element.
          // Change this value to the value that the application requires.
          //if (img_p[idx0][idx1] > threshold)
          //    img[idx0 + img.size(0) * idx1] = true;
          //else
          //    img[idx0 + img.size(0) * idx1] = false;
          img[idx0 + img.size(0) * idx1] = p[idx0][idx1];
      }
  }
}

void CircleIdentifier::SaveCenter(string file_path) {
    NumDataIO<double> data_io;
    file_path.erase(file_path.end() - 3, file_path.end()); // erase "tif"
    data_io.SetFilePath(file_path + "txt");

    //	double particle_2Dcenter[x.size()][2];
    int n_center = center.size();
    double* particle_2Dcenter = new double[n_center * 2];
    //	for (int i = 0; i < x.size(); ++i)
    //		particle_2Dcenter[i] = new double[2];

    for (int i = 0; i < n_center; i++) {
        particle_2Dcenter[2 * i] = center[i].X();
        particle_2Dcenter[2 * i + 1] = center[i].Y();
    }
    data_io.SetTotalNumber(n_center * 2);
    data_io.WriteData(particle_2Dcenter);

    delete[] particle_2Dcenter;
}

void CircleIdentifier::SaveRadius(string file_path) {
    NumDataIO<double> data_io;
    file_path.erase(file_path.end() - 4, file_path.end()); // erase "tif"
    data_io.SetFilePath(file_path + "_r.txt");

    //	double particle_2Dcenter[x.size()][2];
    int n_center = radius.size();
    double* particle_2Dcenter = new double[n_center];
    //	for (int i = 0; i < x.size(); ++i)
    //		particle_2Dcenter[i] = new double[2];

    for (int i = 0; i < n_center; i++) {
        particle_2Dcenter[i] = radius[i];
    }
    data_io.SetTotalNumber(n_center);
    data_io.WriteData(particle_2Dcenter);

    delete[] particle_2Dcenter;
}

// End of code generation (CircleIdentifier.cpp)
