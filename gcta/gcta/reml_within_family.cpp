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

void gcta::detect_family()
{
    cout<<"Detecting sub-matrix for each family from the GRM ..."<<endl;
    int i=0, j=0, k=0, l=0, prev_pnt=0;
    double d_buf1=0.0, d_buf2=0.0;
    _fam_brk_pnt.clear();
    for(i=0; i<_n-1; i++){
        d_buf1=_A[0].row(i).tail(_n-i-1).sum();
        d_buf2=_A[0].col(i).tail(_n-i-1).sum();
        if(CommFunc::FloatEqual(d_buf1, 0.0) && CommFunc::FloatEqual(d_buf2, 0.0)) _fam_brk_pnt.push_back(i);
    }
    
    _Adn.resize(_r_indx.size());
    for(i=0; i<_r_indx.size(); i++) (_Adn[i]).setZero(_n, _n);
    
    int pos=0;
    for(l=0; l<_r_indx.size()-1; l++){
        pos=_r_indx[l];
        prev_pnt=0;
        for(k=0; k<_fam_brk_pnt.size(); k++){
            for(j=prev_pnt; j<=_fam_brk_pnt[k]; j++){
                for(i=prev_pnt; i<=_fam_brk_pnt[k]; i++) _Adn[pos](i,j)=(_A[pos])(i,j);
            }
            prev_pnt=_fam_brk_pnt[k]+1;
        }
        for(j=prev_pnt; j<_n; j++){
            for(i=prev_pnt; i<_n; i++) _Adn[pos](i,j)=(_A[pos])(i,j);
        }
    }
    pos=_r_indx[_r_indx.size()-1];
    for(i=0; i<_n; i++){
        _Adn[pos](i,i)=1.0;
    }
    cout<<"There are "<<_fam_brk_pnt.size()+1<<" sub-matrices detected."<<endl;
    
    // release momery
    _A.clear();
}

bool gcta::calcu_Vi_within_family(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter)
{
    int i=0, j=0, k=0;
    double logdet_buf=0.0;
    string errmsg="\nError: the V (variance-covariance) matrix is not invertible.";
    
    Vi=eigenMatrix::Zero(_n, _n);
    if(_r_indx.size()==1){
        Vi.diagonal()=eigenVector::Constant(_n, 1.0/prev_varcmp[0]);
        logdet=_n*log(prev_varcmp[0]);
    }
    else{
        for(i=0; i<_r_indx.size(); i++) Vi+=_Adn[_r_indx[i]]*prev_varcmp[i];
        int prev_pnt=0, subn=0;
        logdet=0.0;
        for(i=0; i<_fam_brk_pnt.size()+1; i++){
            if(i==_fam_brk_pnt.size()) subn=_n-prev_pnt;
            else subn=_fam_brk_pnt[i]-prev_pnt+1;
            eigenMatrix subVi=Vi.block(prev_pnt, prev_pnt, subn, subn);
            stringstream errmsg;
            errmsg<<"Error: the sub-matrix of V for the "<<i+1<<"-th family is not invertible.";
            if(!comput_inverse_logdet_LDLT(subVi, logdet_buf)) throw(errmsg.str());
            logdet+=logdet_buf;
            //logdet+=comput_inverse_logdet_LU(subVi, errmsg.str());
            Vi.block(prev_pnt, prev_pnt, subn, subn)=subVi;
            
            if(i<_fam_brk_pnt.size()) prev_pnt=_fam_brk_pnt[i]+1;
        }
    }
    
    return true;
}