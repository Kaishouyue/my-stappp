/*****************************************************************************/
/*                              gaussIntegral                                */
/*              funtion for Q4 element stiffmatrix                           */
/*****************************************************************************/

#pragma once
#include "MatrixMultiplication.h"

const double GaussPoints1[1][2] = {
    {0.0, 0.0}
};

const double GaussWeights1[1] = {4.0};

const double GaussPoints2[4][2] = {
    {-0.5773502691896257, -0.5773502691896257},
    {0.5773502691896257, -0.5773502691896257},
    {0.5773502691896257, 0.5773502691896257},
    {-0.5773502691896257, 0.5773502691896257}
};
const double GaussWeights2[4] = {1.0, 1.0, 1.0, 1.0};

myMatrix GaussInt4Q4(const myMatrix& D, const double ksi, const double eta, const double weight, const double x[4], const double y[4]);