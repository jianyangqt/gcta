/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for mixed linera model association analysis
 *
 * 2013 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::mlma(string grm_file, string phen_file, string qcovar_file, string covar_file, int mphen, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, bool no_constrain, bool within_family, bool inbred, bool no_adj_covar)
{
    _within_family=within_family;
    _reml_max_iter=MaxIter;
    unsigned long i=0, j=0;
    bool grm_flag=(!grm_file.empty());
    bool qcovar_flag=(!qcovar_file.empty());
    bool covar_flag=(!covar_file.empty());
    if(!qcovar_flag && !covar_flag) no_adj_covar=false;
    
    // Read data
    int qcovar_num=0, covar_num=0;
    vector<string> phen_ID, qcovar_ID, covar_ID, grm_id;
    vector< vector<string> > phen_buf, qcovar, covar; // save individuals by column
    
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if(qcovar_flag){
        qcovar_num=read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if(covar_flag){
        covar_num=read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    if(grm_flag){
        read_grm(grm_file, grm_id);
        update_id_map_kp(grm_id, _id_map, _keep);
    }
    else{
        make_grm_mkl(false, false, inbred, true, 0, true);
        for(i=0; i<_keep.size(); i++) grm_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
    }
    
    vector<string> uni_id;
	map<string, int> uni_id_map;
    map<string, int>::iterator iter;
	for(i=0; i<_keep.size(); i++){
	    uni_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
	    uni_id_map.insert(pair<string,int>(_fid[_keep[i]]+":"+_pid[_keep[i]], i));
	}
    _n=_keep.size();
    if(_n<1) throw("Error: no individual is in common in the input files.");
    cout<<_n<<" individuals are in common in these files."<<endl;
    
    // construct model terms
    _y.setZero(_n);
    for(i=0; i<phen_ID.size(); i++){
        iter=uni_id_map.find(phen_ID[i]);
        if(iter==uni_id_map.end()) continue;
        _y[iter->second]=atof(phen_buf[i][mphen-1].c_str());
    }
    
    _r_indx.clear();
    vector<int> kp;
    for(i=0; i<2; i++) _r_indx.push_back(i);
    _A.resize(_r_indx.size());
    StrFunc::match(uni_id, grm_id, kp);
    (_A[0])=eigenMatrix::Zero(_n, _n);
    if(grm_flag){
        #pragma omp parallel for private(j)
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=_grm(kp[i],kp[j]);
        }
        _grm.resize(0,0);
    }
    else{
        #pragma omp parallel for private(j)
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=_grm_mkl[kp[i]*_n+kp[j]];
        }
        delete[] _grm_mkl;
    }
    _A[_r_indx.size()-1]=eigenMatrix::Identity(_n, _n);
    
    // construct X matrix
    vector<eigenMatrix> E_float;
    eigenMatrix qE_float;
    construct_X(_n, uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);
    
    // names of variance component
    _var_name.push_back("V(G)");
    _hsq_name.push_back("V(G)/Vp");
    _var_name.push_back("V(e)");
    
    // within family
    if(_within_family) detect_family();
    
    // run REML algorithm
    cout<<"\nPerforming MLM association analyses (including the candidate SNP) ..."<<endl;
    unsigned long k=0, n=_keep.size(), m=_include.size();
	reml(false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true);
    _P.resize(0,0);
    _A.clear();
    float *y=new float[n];
    eigenVector y_buf=_y;
    if(!no_adj_covar) y_buf=_y.array()-(_X*_b).array(); // adjust phenotype for covariates
    for(i=0; i<n; i++) y[i]=y_buf[i];
    
    if(grm_flag){
        cout<<endl;
        _geno_mkl=new float[n*m];
        make_XMat_mkl(_geno_mkl, false);
        #pragma omp parallel for private(j, k)
        for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                k=i*m+j;
                if(_geno_mkl[k]<1e5) _geno_mkl[k]-=_mu[_include[j]];
                else _geno_mkl[k]=0.0;
            }
        }
    }
    
    eigenVector beta, se, pval;
    if(no_adj_covar) mlma_calcu_stat_covar(y, _geno_mkl, n, m, beta, se, pval);
    else mlma_calcu_stat(y, _geno_mkl, n, m, beta, se, pval);
    delete[] y;
    delete[] _geno_mkl;
    
    string filename=_out+".mlma";
    cout<<"\nSaving the results of the mixed linear model association analyses of "<<m<<" SNPs to ["+filename+"] ..."<<endl;
    ofstream ofile(filename.c_str());
    if(!ofile) throw("Can not open the file ["+filename+"] to write.");
    ofile<<"Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\tp"<<endl;
	for(i=0; i<m; i++){
        j=_include[i];
        ofile<<_chr[j]<<"\t"<<_snp_name[j]<<"\t"<<_bp[j]<<"\t"<<_ref_A[j]<<"\t"<<_other_A[j]<<"\t";
        if(pval[i]>1.5) ofile<<"NA\tNA\tNA\tNA"<<endl;
        else ofile<<0.5*_mu[j]<<"\t"<<beta[i]<<"\t"<<se[i]<<"\t"<<pval[i]<<endl;
    }
    ofile.close();
}

void gcta::mlma_calcu_stat(float *y, float *geno_mkl, unsigned long n, unsigned long m, eigenVector &beta, eigenVector &se, eigenVector &pval)
{
    unsigned long i=0, j=0;
    double Xt_Vi_X=0.0, chisq=0.0;
    float *X=new float[n];
    float *Vi_X=new float[n];
    float *Vi=new float[n*n];
    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        for(j=0; j<n; j++) Vi[i*n+j]=_Vi(i,j);
    }
    _Vi.resize(0,0);
    
    beta.resize(m);
    se=eigenVector::Zero(m);
    pval=eigenVector::Constant(m,2);
    cout<<"\nRunning association tests for "<<m<<" SNPs ..."<<endl;
    for(i=0; i<m; i++){
        for(j=0; j<n; j++) X[j]=geno_mkl[j*m+i];
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, Vi, n, X, 1, 0.0, Vi_X, 1);
        Xt_Vi_X=cblas_sdot(n, X, 1, Vi_X, 1);
        se[i]=1.0/Xt_Vi_X;
        beta[i]=se[i]*cblas_sdot(n, y, 1, Vi_X, 1);
        if(se[i]>1.0e-30){
            se[i]=sqrt(se[i]);
            chisq=beta[i]/se[i];
            pval[i]=StatFunc::pchisq(chisq*chisq, 1);
        }
    }
    delete[] X;
    delete[] Vi_X;
    delete[] Vi;
}

void gcta::mlma_calcu_stat_covar(float *y, float *geno_mkl, unsigned long n, unsigned long m, eigenVector &beta, eigenVector &se, eigenVector &pval)
{
    unsigned long i=0, j=0, col_num=_X_c+1;
    double chisq=0.0, d_buf=0.0;
    float *Vi=new float[n*n];
    float *X=new float[n*col_num];
    float *Vi_X=new float[n*col_num];
    float *Xt_Vi_X=new float[col_num*col_num];
    float *Xt_Vi_y=new float[col_num];
    float *b_vec=new float[col_num];
    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        for(j=0; j<n; j++) Vi[i*n+j]=_Vi(i,j);
    }
    _Vi.resize(0,0);
    for(i=0; i<n; i++){
        for(j=0; j<_X_c; j++) X[i*col_num+j]=_X(i,j);
        X[i*col_num+_X_c]=0.0;
    }

    beta.resize(m);
    se=eigenVector::Zero(m);
    pval=eigenVector::Constant(m,2);
    cout<<"\nRunning association tests for "<<m<<" SNPs ..."<<endl;
    for(i=0; i<m; i++){
        for(j=0; j<n; j++) X[j*col_num+_X_c]=geno_mkl[j*m+i];
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, col_num, n, 1.0, Vi, n, X, col_num, 0.0, Vi_X, col_num);
        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, col_num, col_num, n, 1.0, X, col_num, Vi_X, col_num, 0.0, Xt_Vi_X, col_num);
        if(!comput_inverse_logdet_LU_mkl_array(col_num, Xt_Vi_X, d_buf)) throw("Error: Xt_Vi_X is not invertable.");
        cblas_sgemv(CblasRowMajor, CblasTrans, n, col_num, 1.0, Vi_X, col_num, y, 1, 0.0, Xt_Vi_y, 1);
        cblas_sgemv(CblasRowMajor, CblasNoTrans, col_num, col_num, 1.0, Xt_Vi_X, col_num, Xt_Vi_y, 1, 0.0, b_vec, 1);
        se[i]=Xt_Vi_X[_X_c*col_num+_X_c];
        beta[i]=b_vec[_X_c];
        if(se[i]>1.0e-30){
            se[i]=sqrt(se[i]);
            chisq=beta[i]/se[i];
            pval[i]=StatFunc::pchisq(chisq*chisq, 1);
        }
    }
    delete[] Vi;
    delete[] X;
    delete[] Vi_X;
    delete[] Xt_Vi_X;
    delete[] Xt_Vi_y;
    delete[] b_vec;
}

void gcta::mlma_loco(string phen_file, string qcovar_file, string covar_file, int mphen, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, bool no_constrain, bool inbred, bool no_adj_covar)
{
    unsigned long i=0, j=0, k=0, c1=0, c2=0, n=0;
    _reml_max_iter=MaxIter;
    bool qcovar_flag=(!qcovar_file.empty());
    bool covar_flag=(!covar_file.empty());
    if(!qcovar_flag && !covar_flag) no_adj_covar=false;
    
    // Read data
    int qcovar_num=0, covar_num=0;
    vector<string> phen_ID, qcovar_ID, covar_ID, grm_id;
    vector< vector<string> > phen_buf, qcovar, covar; // save individuals by column
    
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if(qcovar_flag){
        qcovar_num=read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if(covar_flag){
        covar_num=read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    n=_keep.size();
    _n=_keep.size();
    if(_n<1) throw("Error: no individual is in common in the input files.");
    cout<<_n<<" individuals are in common in these files."<<endl;
    
    vector<int> chrs, vi_buf(_chr);
    stable_sort(vi_buf.begin(), vi_buf.end());
	vi_buf.erase(unique(vi_buf.begin(), vi_buf.end()), vi_buf.end());
    if(vi_buf.size()<2) throw("Error: There is only one chromosome. The MLM leave-on-chromosome-out (LOCO) analysis requires at least two chromosome.");
    for(i=0; i<vi_buf.size(); i++){
        if(vi_buf[i]<=_autosome_num) chrs.push_back(vi_buf[i]);
    }
    vector<int> include_o(_include);
    map<string, int> snp_name_map_o(_snp_name_map);
    vector<float> m_chrs_f(chrs.size());
    vector<float *> grm_chrs(chrs.size());
    vector<float *> geno_chrs(chrs.size());
    vector< vector<int> > icld_chrs(chrs.size());
    cout<<endl;
    if(_mu.empty()) calcu_mu();
    cout<<"\nCalculating the genetic relationship matrix for each of the "<<chrs.size()<<" chromosomes ... "<<endl;
    for(c1=0; c1<chrs.size(); c1++){
        cout<<"Chr "<<chrs[c1]<<":"<<endl;
        extract_chr(chrs[c1], chrs[c1]);
        make_grm_mkl(false, false, inbred, true, 0, true);
        
        m_chrs_f[c1]=(float)_include.size();
        icld_chrs[c1]=_include;
        _include=include_o;
        _snp_name_map=snp_name_map_o;
        
        geno_chrs[c1]=_geno_mkl;
        _geno_mkl=NULL;
        grm_chrs[c1]=_grm_mkl;
        _grm_mkl=NULL;
    }
    for(i=0; i<_keep.size(); i++) grm_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
    
    vector<string> uni_id;
	map<string, int> uni_id_map;
    map<string, int>::iterator iter;
	for(i=0; i<_keep.size(); i++){
	    uni_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
	    uni_id_map.insert(pair<string,int>(_fid[_keep[i]]+":"+_pid[_keep[i]], i));
	}
    
    // construct model terms
    _y.setZero(_n);
    for(i=0; i<phen_ID.size(); i++){
        iter=uni_id_map.find(phen_ID[i]);
        if(iter==uni_id_map.end()) continue;
        _y[iter->second]=atof(phen_buf[i][mphen-1].c_str());
    }
    
    // construct X matrix
    vector<eigenMatrix> E_float;
    eigenMatrix qE_float;
    construct_X(_n, uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);
    
    // names of variance component
    _var_name.push_back("V(G)");
    _hsq_name.push_back("V(G)/Vp");
    _var_name.push_back("V(e)");
    
    // MLM association
    cout<<"\nPerforming MLM association analyses (leave-one-chromosome-out) ..."<<endl;
    
    vector<int> kp;
    StrFunc::match(uni_id, grm_id, kp);
    _r_indx.resize(2);
    for(i=0; i<2; i++) _r_indx[i]=i;
    _A.resize(_r_indx.size());
    _A[1]=eigenMatrix::Identity(_n, _n);
    
    eigenVector y_buf=_y;
    float *y=new float[_n];
    vector<eigenVector> beta(chrs.size()), se(chrs.size()), pval(chrs.size());
    for(c1=0; c1<chrs.size(); c1++){
        cout<<"\n-----------------------------------\n#Chr "<<chrs[c1]<<":"<<endl;
        extract_chr(chrs[c1], chrs[c1]);
        
        _A[0]=eigenMatrix::Zero(_n, _n);
        double d_buf=0;
        for(c2=0; c2<chrs.size(); c2++){
            if(chrs[c1]==chrs[c2]) continue;
            #pragma omp parallel for private(j)
            for(i=0; i<_n; i++){
                for(j=0; j<=i; j++){
                    (_A[0])(i,j)+=(grm_chrs[c2])[kp[i]*_n+kp[j]]*m_chrs_f[c2];
                }
            }
            d_buf+=m_chrs_f[c2];
        }
        
        #pragma omp parallel for private(j)
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++){
                (_A[0])(i,j)/=d_buf;
                (_A[0])(j,i)=(_A[0])(i,j);
            }
        }
        
        // run REML algorithm
        reml(false, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true);
        if(!no_adj_covar) y_buf=_y.array()-(_X*_b).array(); // adjust phenotype for covariates
        for(i=0; i<_n; i++) y[i]=y_buf[i];
        reml_priors.clear();
        reml_priors_var=_varcmp;
        _P.resize(0,0);
        _A[0].resize(0,0);
        
        mlma_calcu_stat(y, (geno_chrs[c1]), n, _include.size(), beta[c1], se[c1], pval[c1]);
        
        _include=include_o;
        _snp_name_map=snp_name_map_o;
        cout<<"-----------------------------------"<<endl;
    }
    
    delete[] y;
    for(c1=0; c1<chrs.size(); c1++){
        delete[] (grm_chrs[c1]);
        delete[] (geno_chrs[c1]);
    }
    
    string filename=_out+".loco.mlma";
    cout<<"\nSaving the results of the mixed linear model association analyses of "<<_include.size()<<" SNPs to ["+filename+"] ..."<<endl;
    ofstream ofile(filename.c_str());
    if(!ofile) throw("Can not open the file ["+filename+"] to write.");
    ofile<<"Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\tp"<<endl;
    for(c1=0; c1<chrs.size(); c1++){
        for(i=0; i<icld_chrs[c1].size(); i++){
            j=icld_chrs[c1][i];
            ofile<<_chr[j]<<"\t"<<_snp_name[j]<<"\t"<<_bp[j]<<"\t"<<_ref_A[j]<<"\t"<<_other_A[j]<<"\t";
            if(pval[c1][i]>1.5) ofile<<"NA\tNA\tNA\tNA"<<endl;
            else ofile<<0.5*_mu[j]<<"\t"<<beta[c1][i]<<"\t"<<se[c1][i]<<"\t"<<pval[c1][i]<<endl;
        }
    }
    ofile.close();
}

void gcta::grm_minus_grm(float *grm, float *sub_grm)
{
    int i=0, j=0, k=0, n=_n;
    
    #pragma omp parallel for private(j,k)
    for(i=0; i<n; i++){
		for(j=0; j<=i; j++){
            k=i*n+j;
            sub_grm[k]=grm[k]-sub_grm[k];
		}
	}

}


