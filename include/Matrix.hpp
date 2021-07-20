#ifndef GCTA2_MATRIX_H
#define GCTA2_MATRIX_H

#include "cpu.h"
#include <Eigen/Eigen>
#include <iostream>
#include <Logger.h>

static_assert(std::numeric_limits<double>::is_iec559, "Not a supported compiler");

// two step to inverse the matrix

enum INVmethod{
    INV_LLT = 1,
    INV_LU = 2,
    INV_QR = 3,
    INV_FQR = 4,
    INV_SVD = 5,
    INV_ERR = 100
};

template<typename MatrixType>
bool _LLT(MatrixType &A, double &logdet){
    //MatrixType::Scalar * vi = A.data();
    auto * vi = A.data();
    Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, 1> diag = A.diagonal();
    //auto diag = A.diagonal();

    int info, cols = (int)A.cols();
    char uplo = 'L';
    LOGGER.ts("LLT");
#if GCTA_CPU_x86
    dpotrf(&uplo, &cols, vi, &cols, &info);
#else
    dpotrf_(&uplo, &cols, vi, &cols, &info);
#endif    
    //LOGGER << "  LLT time: " << LOGGER.tp("LLT") << std::endl;
    if(info == 0){
        logdet = A.diagonal().array().square().log().sum();
        //LOGGER.ts("LLT_INV");
#if GCTA_CPU_x86        
        dpotri(&uplo, &cols, vi, &cols, &info);
#else
        dpotri_(&uplo, &cols, vi, &cols, &info);
#endif
        //LOGGER << "  LLT inverse time: " << LOGGER.tp("LLT_INV") << std::endl;
        if(info == 0){
            A.template triangularView<Eigen::Upper>() = A.transpose();
            return true;
        }
    }

    A.template triangularView<Eigen::Lower>() = A.transpose();
    A.diagonal() = diag;
    return false;
}


template<typename MatrixType>
bool SquareMatrixInverse(MatrixType &A, double &logdet, int &rank, INVmethod &method){
    int n = A.rows();
    bool ret = false;
    switch(method){
        case INV_LLT:
            {
                if(_LLT(A, logdet)){
                    method = INV_LLT;
                    ret = true;
                    break;
                }
            }
        case INV_LU:
            {
                Eigen::PartialPivLU<MatrixType> lu(A);
                double det = lu.determinant();
                //std::cout << "LU det: " << std::scientific << det << std::endl;
                if(det >= 1e-10 || det <= -1e-10){
                    logdet = lu.matrixLU().diagonal().array().abs().log().sum();
                    A = lu.inverse();
                    method = INV_LU;
                    ret = true;
                    break;
                }
            }
        case INV_QR:
            {
                Eigen::HouseholderQR<MatrixType> qr(A);
                double det = qr.absDeterminant();
                //std::cout << "QR det: " << std::scientific << det << std::endl;
                if(det > 1e-16){
                    logdet = qr.logAbsDeterminant();
                    A = qr.solve(MatrixType::Identity(n, n));
                    method = INV_QR;
                    ret = true;
                    break;
                }
            }
        case INV_FQR:
            // not necessary 
            // Eigen::HouseholderQR<MatrixType> qr(A);
            {
                Eigen::ColPivHouseholderQR<MatrixType> qr(A);
                if(qr.isInvertible()){
                    logdet = qr.logAbsDeterminant();
                    A = qr.inverse();
                    method = INV_QR;
                    ret = true;
                    break;
                }else{
                    rank = qr.rank();
                    // it will be extreme slow or not accurate
                    if(n > 50000 || 1.0 * rank / n < 0.99){
                        ret = false;
                        break;
                    }
                }
            }
        case INV_SVD:
            //Eigen::BDCSVD<MatrixType> svd(A, Eigen::ComputeThinU|Eigen::ComputeThinV);
            ;
        default:
            rank = 0;
            ret = false;
            method = INV_ERR;
    }
    return ret;
}

template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif //header

