//
//  GlobalVarb.h
//  gcta_code
//
//  Created by Zhihong Zhu on 31/08/2014.
//  Copyright (c) 2014 Queensland Brain Institute. All rights reserved.
//

#ifndef _GLOBALVARB_H
#define _GLOBALVARB_H

#endif

#ifndef EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <string>

using namespace Eigen;
using namespace std;

#ifdef SINGLE_PRECISION
typedef DiagonalMatrix<float, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;
typedef SparseMatrix<float> eigenSparseMat;
typedef DynamicSparseMatrix<float> eigenDynSparseMat;
#else
typedef DiagonalMatrix<double, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;
typedef SparseMatrix<double> eigenSparseMat;
typedef DynamicSparseMatrix<double> eigenDynSparseMat;
#endif

// for missing value
struct missValue
{
    int num;
    string pos;
};