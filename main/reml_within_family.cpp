/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Within-family REML analysis
 *
 * 2012 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"
#include <Eigen/SparseCholesky>

void gcta::detect_family()
{
    LOGGER<<"Accelarating for family based GRM ..."<<endl;
    int i=0, j=0, k=0, l=0, prev_pnt=0;
    double d_buf1=0.0, d_buf2=0.0;
    _fam_brk_pnt.clear();

    _Asp.resize(_r_indx.size());

    int pos;
    for(pos = 0; pos < _r_indx.size() - 1; pos++){
        _Asp[pos] = _A[pos].sparseView();
    }

    pos=_r_indx[_r_indx.size()-1];
    (_Asp[pos]).resize(_n, _n);
    _Asp[pos].setIdentity();

    _A.clear();
}

bool gcta::calcu_Vi_within_family(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter)
{
    int i=0, j=0, k=0;
    double logdet_buf=0.0;
    string errmsg="\n  the V (variance-covariance) matrix is not invertible.";
    
    if(_r_indx.size()==1){
        Vi = eigenMatrix::Zero(_n, _n);
        Vi.diagonal()=eigenVector::Constant(_n, 1.0/prev_varcmp[0]);
        logdet=_n*log(prev_varcmp[0]);
    }
    else{
        eigenSparseMat Vit(_n, _n);
        for(i=0; i<_r_indx.size(); i++){
            eigenSparseMat vit_temp = (_Asp[_r_indx[i]])*prev_varcmp[i];
            Vit = Vit + vit_temp;
        }
        int prev_pnt=0, subn=0;
        //logdet=0.0;

        Eigen::SimplicialLDLT<eigenSparseMat> solver;
        solver.compute(Vit);

        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the matrix within familiy");
        }
        
        double logdet = solver.vectorD().array().square().log().sum();

        eigenSparseMat vit_r = solver.solve(_Asp[_r_indx[_r_indx.size() - 1]]);

        Vi = eigenMatrix(vit_r);

    }
    
    return true;
}
