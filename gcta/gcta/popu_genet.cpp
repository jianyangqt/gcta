//
//  popu_genet.cpp
//  gcta
//
//  Created by Jian Yang on 14/03/13.
//
//

#include "gcta.h"

// paa: proportion of ancestral alleles
void gcta::paa(string aa_file)
{
    check_autosome();
    if(_mu.empty()) calcu_mu();
	//rm_high_ld();
    
    int i=0, j=0, k=0;
    
    // read ancestral alleles from a file
    ifstream i_aa(aa_file.c_str());
    if(!i_aa) throw("Error: can not open the file ["+aa_file+"] to read.");
    string cbuf=".";
    string str_buf;
	cout<<"Reading ancestral alleles of the SNPs from ["+aa_file+"]."<<endl;
    map<string, int>::iterator iter, End=_snp_name_map.end();
    vector<string> aa(_snp_num);
    for(i=0; i<_snp_num; i++) aa[i]=".";
    int icount=0;
	while(i_aa){
		i_aa>>str_buf;
		if(i_aa.eof()) break;
		iter=_snp_name_map.find(str_buf);
		i_aa>>cbuf;
		if(iter!=End && cbuf!="."){
		    aa[iter->second]=cbuf;
		    icount++;
		}
	}
    i_aa.close();
	cout<<"Ancestral alleles for "<<icount<<" SNPs are included from ["+aa_file+"]."<<endl;
    
	cout<<"Calculating proportion of ancestral alleles ..."<<endl;
	double x=0.0;
	vector<double> hom_aa_rare(_keep.size()), hom_aa_comm(_keep.size()), hom_da_rare(_keep.size()), hom_da_comm(_keep.size()), het_aa_rare(_keep.size()), het_aa_comm(_keep.size()), nomiss(_keep.size());
	for(i=0; i<_keep.size(); i++){
 		for(k=0; k<_include.size(); k++){
 		    if(aa[_include[k]]==".") continue;
            if(!_snp_1[_include[k]][_keep[i]] || _snp_2[_include[k]][_keep[i]]){
                x=_snp_1[_include[k]][_keep[i]]+_snp_2[_include[k]][_keep[i]];
                if(x<0.1){
                    if(_ref_A[_include[k]]==aa[_include[k]]){
                        if(_mu[_include[k]]>1.0) hom_da_rare[i]+=1.0;
                        else hom_da_comm[i]+=1.0;
                    }
                    else{
                        if(_mu[_include[k]]>1.0) hom_aa_rare[i]+=1.0;
                        else hom_aa_comm[i]+=1.0;
                    }
                }
                else if(x>1.9){
                    if(_ref_A[_include[k]]==aa[_include[k]]){
                        if(_mu[_include[k]]>1.0) hom_aa_comm[i]+=1.0;
                        else hom_aa_rare[i]+=1.0;
                    }
                    else{
                        if(_mu[_include[k]]>1.0) hom_da_comm[i]+=1.0;
                        else hom_da_rare[i]+=1.0;
                    }
                }
                else{
                    if(_ref_A[_include[k]]==aa[_include[k]]){
                        if(_mu[_include[k]]>1.0) het_aa_comm[i]+=1.0;
                        else het_aa_rare[i]+=1.0;
                    }
                    else{
                        if(_mu[_include[k]]>1.0) het_aa_rare[i]+=1.0;
                        else het_aa_comm[i]+=1.0;
                    }
                }
                nomiss[i]+=1.0;
            }
        }
        hom_aa_rare[i]/=nomiss[i];
        hom_aa_comm[i]/=nomiss[i];
        hom_da_rare[i]/=nomiss[i];
        hom_da_comm[i]/=nomiss[i];
        het_aa_rare[i]/=nomiss[i];
        het_aa_comm[i]/=nomiss[i];
        cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}
    
    // Save matrix A in binary file
	string paa_file=_out+".paa";
	ofstream o_paa(paa_file.c_str());
	if(!o_paa) throw("Error: can not open the file ["+paa_file+"] to write.");
	o_paa<<"FID\tIID\tNOMISS\tHOM_AA_RARE\tHOM_AA_COMM\tHOM_DA_RARE\tHOM_DA_COMM\tHET_AA_RARE\tHET_AA_COMM"<<endl;
	for(i=0; i<_keep.size(); i++) o_paa<<_fid[i]<<"\t"<<_pid[i]<<"\t"<<nomiss[i]<<"\t"<<hom_aa_rare[i]<<"\t"<<hom_aa_comm[i]<<"\t"<<hom_da_rare[i]<<"\t"<<hom_da_comm[i]<<"\t"<<het_aa_rare[i]<<"\t"<<het_aa_comm[i]<<endl;
	o_paa.close();
	cout<<"Proportion of ancestral alleles has been saved in file ["+paa_file+"]."<<endl;
}

// inbreeding coefficient
void gcta::ibc(bool ibc_all)
{
    check_autosome();
    if(_mu.empty()) calcu_mu();
	//rm_high_ld();
    
    int i=0, j=0, k=0;
    
	// Calcuate A matrix
	cout<<"Calculating the inbreeding coefficients ..."<<endl;
	vector<double> h(_include.size()), h_i(_include.size()), w(_include.size()), p(_include.size()), p_q(_include.size()), q_p(_include.size()); // variance of each SNP, 2pq
	#pragma omp parallel for
    for(j=0; j<_include.size(); j++){
	    p[j]=0.5*_mu[_include[j]];
        h[j]=2.0*p[j]*(1.0-p[j]);
        p_q[j]=p[j]/(1.0-p[j]);
        if(fabs(p_q[j])<1.0e-50) q_p[j]=0.0;
        else q_p[j]=1.0/p_q[j];
        w[j]=h[j]/(1.0-h[j]);
        if(fabs(h[j])<1.0e-50) h_i[j]=0.0;
        else h_i[j]=1.0/h[j];
	}
	vector<double> Fhat1, Fhat1_w, Fhat2, Fhat2_w, Fhat3, Fhat4, Fhat5, Fhat6, Fhat7, rare_hom, comm_hom, nomiss;
    Fhat1.resize(_keep.size());
    Fhat2.resize(_keep.size());
    Fhat3.resize(_keep.size());
    nomiss.resize(_keep.size());
    if(ibc_all){
        Fhat4.resize(_keep.size());
        Fhat5.resize(_keep.size());
        Fhat6.resize(_keep.size());
        Fhat7.resize(_keep.size());
        Fhat1_w.resize(_keep.size());
        Fhat2_w.resize(_keep.size());
        rare_hom.resize(_keep.size());
        comm_hom.resize(_keep.size());
    }
    #pragma omp parallel for private(k)
    for(i=0; i<_keep.size(); i++){
        double x=0.0, sum_w=0.0, sum_h=0.0, Fhat_buf=0.0;
		for(k=0; k<_include.size(); k++){
            if(!_snp_1[_include[k]][_keep[i]] || _snp_2[_include[k]][_keep[i]]){
                x=_snp_1[_include[k]][_keep[i]]+_snp_2[_include[k]][_keep[i]];
                if(_allele2[_include[k]]==_ref_A[_include[k]]) x=2.0-x;
                Fhat_buf=(x-_mu[_include[k]])*(x-_mu[_include[k]]);
                if(ibc_all) Fhat4[i]+=Fhat_buf;
                Fhat_buf*=h_i[k];
                Fhat1[i]+=Fhat_buf;
                if(ibc_all) Fhat1_w[i]+=Fhat_buf*w[k];
                Fhat_buf=h[k]-x*(2.0-x);
                if(ibc_all) Fhat6[i]+=Fhat_buf;
                Fhat_buf*=h_i[k];
                Fhat2[i]+=Fhat_buf;
                if(ibc_all) Fhat2_w[i]+=Fhat_buf*w[k];
                Fhat_buf=(x*(x-1.0-_mu[_include[k]])+_mu[_include[k]]*p[k]);
                Fhat3[i]+=Fhat_buf*h_i[k];
                if(ibc_all){
                    Fhat5[i]+=Fhat_buf;
                    sum_w+=w[k];
                    sum_h+=h[k];
                    if(x<0.1) Fhat7[i]+=p_q[k];
                    else if(x>1.9) Fhat7[i]+=q_p[k];
                    // Count the number of rare and common homozygotes
                    if(x<0.1){
                        if(p[k]>0.5) rare_hom[i]+=1.0;
                        else comm_hom[i]+=1.0;
                    }
                    else if(x>1.9){
                        if(p[k]<0.5) rare_hom[i]+=1.0;
                        else comm_hom[i]+=1.0;
                    }
                }
                nomiss[i]+=1.0;
            }
        }
        Fhat1[i]/=nomiss[i];
        Fhat2[i]/=nomiss[i];
        Fhat3[i]/=nomiss[i];
        if(ibc_all){
            Fhat1_w[i]/=sum_w;
            Fhat2_w[i]/=sum_w;
            Fhat4[i]/=sum_h;
            Fhat5[i]/=sum_h;
            Fhat6[i]/=sum_h;
            rare_hom[i]/=nomiss[i];
            comm_hom[i]/=nomiss[i];
            Fhat7[i]/=nomiss[i];
        }
	}
    
    // Save matrix A in binary file
	string ibc_file=_out+".ibc";
	ofstream o_ibc(ibc_file.c_str());
	if(!o_ibc) throw("Error: can not open the file ["+ibc_file+"] to write.");
	if(ibc_all){
        o_ibc<<"FID\tIID\tNOMISS\tP_RARE_HOM\tP_COMM_HOM\tFhat1\tFhat1_w\tFhat2\tFhat2_w\tFhat3\tFhat4\tFhat5\tFhat6\tFhat7"<<endl;
        for(i=0; i<_keep.size(); i++) o_ibc<<_fid[_keep[i]]<<"\t"<<_pid[_keep[i]]<<"\t"<<nomiss[i]<<"\t"<<rare_hom[i]<<"\t"<<comm_hom[i]<<"\t"<<Fhat1[i]-1.0<<"\t"<<Fhat1_w[i]-1.0<<"\t"<<Fhat2[i]<<"\t"<<Fhat2_w[i]<<"\t"<<Fhat3[i]<<"\t"<<Fhat4[i]-1.0<<"\t"<<Fhat5[i]<<"\t"<<Fhat6[i]<<"\t"<<Fhat7[i]<<endl;
	}
	else{
        o_ibc<<"FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3"<<endl;
        for(i=0; i<_keep.size(); i++) o_ibc<<_fid[_keep[i]]<<"\t"<<_pid[_keep[i]]<<"\t"<<nomiss[i]<<"\t"<<Fhat1[i]-1.0<<"\t"<<Fhat2[i]<<"\t"<<Fhat3[i]<<endl;
	}
	o_ibc.close();
	cout<<"Inbreeding coefficients have been saved in the file ["+ibc_file+"]."<<endl;
}

void gcta::Fst(string filename)
{
    vector<string> subpopu, subpopu_name;
    read_subpopu(filename, subpopu, subpopu_name);
   
    eigenMatrix S;
    coeff_mat(subpopu, S, "Error: too many subpopulations specified in the file ["+filename+"].", "Error: at least two sub-populations are required.");
    int r=subpopu_name.size();
    eigenVector x(_keep.size()), tmp, n_sub(r), xt_S, p_sub(r), h_sub(r);//Fst(_include.size()),  f_all(_include.size()), chisq(_include.size()), pval=eigenVector::Constant(_include.size(), 1.0);
    if(_mu.empty()) calcu_mu();
    
    int i = 0, j = 0;
    for(i=0; i<r; i++) n_sub[i]=S.col(i).sum();
    double n_bar = n_sub.mean();
    double n_bar_r = n_bar * r;
    double n_c = (n_bar_r - n_sub.squaredNorm() / n_bar_r) / (r - 1.0);
    
    string outfile=_out+".fst";
    cout<<"\nSaving the Fst test results"<<_include.size()<<" SNPs to ["+outfile+"] ..."<<endl;
    ofstream ofile(outfile.c_str());
    if(!ofile) throw("Can not open the file ["+outfile+"] to write.");
    ofile<<"Chr\tSNP\tbp\trefA\t";
    for(i=0; i<r; i++) ofile<<"freq_"<<subpopu_name[i]<<"(n="<<n_sub[i]<<")\t";
    ofile<<"Fst"<<endl;
    double p_bar=0.0, s_sq=0.0, h_bar=0.0, a=0.0, b=0.0, c=0.0, Fst=0.0, d_buf=0.0;
    for(i=0; i<_include.size(); i++){
        ofile<<_chr[_include[i]]<<"\t"<<_snp_name[_include[i]]<<"\t"<<_bp[_include[i]]<<"\t"<<_ref_A[_include[i]]<<"\t";
        makex_eigenVector(i, x, false, false);
        for(j=0; j<r; j++){
            p_sub[j]=x.dot(S.col(j))*0.5/n_sub[j];
            ofile<<p_sub[j]<<"\t";
        }
        p_bar = p_sub.dot(n_sub) / n_bar_r;
        tmp = p_sub.array() - p_bar;
        s_sq = (tmp.array() * tmp.array() * n_sub.array()).sum() / ((r - 1.0) * n_bar);
        h_sub = 2 * p_sub.array() * (1.0 - p_sub.array());
        h_bar = h_sub.dot(n_sub) / n_bar_r;
        d_buf = p_bar * (1.0 - p_bar) - (r - 1.0) * s_sq / r;
        a = (n_bar / n_c) * (s_sq - (1.0 / (n_bar - 1.0)) * (d_buf - 0.25 * h_bar));
        b = (n_bar / (n_bar - 1.0)) * (d_buf - (2 * n_bar - 1.0) * h_bar / (4 * n_bar));
        c = 0.5 * h_bar;

        Fst= a / (a + b + c);

        /*double p_mean = p_sub.mean();
        tmp = p_sub.array() - p_mean;
        double Fst_fixed = (tmp.array() * tmp.array() * n_sub.array()).sum();
        Fst_fixed = Fst_fixed / (p_mean * (1.0 - p_mean));
        Fst_fixed = Fst_fixed / ((r - 1.0) * n_bar);*/

        ofile<<Fst<<endl;
    }
    ofile.close();
}

void gcta::read_subpopu(string filename, vector<string> &subpopu, vector<string> &subpopu_name)
{
    cout<<"Reading sub-population information from ["+filename+"]."<<endl;
    ifstream ifstream_subpopu(filename.c_str());
    if(!ifstream_subpopu) throw("Error: can not open the file ["+filename+"] to read.");
    
    vector<string> ID;
    vector< vector<string> > subpopu_buf;
    read_fac(ifstream_subpopu, ID, subpopu_buf);
    update_id_map_kp(ID, _id_map, _keep);
    cout<<"Sub-population information for "<<_keep.size()<<" individuals matched to the genotype data."<<endl;
    
    int i=0, j=0;
    map<string, int> uni_id_map;
    map<string, int>::iterator iter;
    for(i=0; i<_keep.size(); i++) uni_id_map.insert(pair<string,int>(_fid[_keep[i]]+":"+_pid[_keep[i]], i));
    subpopu.clear();
    subpopu.resize(_keep.size());
    for(i=0; i<ID.size(); i++){
        iter=uni_id_map.find(ID[i]);
        if(iter!=uni_id_map.end()) subpopu[iter->second]=subpopu_buf[i][0];
    }
    
    subpopu_name.clear();
    subpopu_name=subpopu;
    stable_sort(subpopu_name.begin(), subpopu_name.end());
    subpopu_name.erase(unique(subpopu_name.begin(), subpopu_name.end()), subpopu_name.end());
 }

 void gcta::std_XMat_subpopu(string subpopu_file, MatrixXf &X, eigenVector &sd_SNP, bool grm_xchr_flag, bool miss_with_mu, double make_grm_scl)
{
    vector<string> subpopu, subpopu_name;
    read_subpopu(subpopu_file, subpopu, subpopu_name);
   
    eigenMatrix S;
    coeff_mat(subpopu, S, "Error: too many subpopulations specified in the file ["+subpopu_file+"].", "Error: at least two sub-populations are required.");
    int popu_num = subpopu_name.size();

    if (_mu.empty()) calcu_mu();

    unsigned long i = 0, j = 0, n = _keep.size(), m = _include.size();
    sd_SNP.resize(m);
    if (_dosage_flag) {
        for (j = 0; j < m; j++)  sd_SNP[j] = (X.col(j) - VectorXf::Constant(n, _mu[_include[j]])).squaredNorm() / (n - 1.0);
    } 
    else {
        for (j = 0; j < m; j++) sd_SNP[j] = _mu[_include[j]]*(1.0 - 0.5 * _mu[_include[j]]);
    }
    if (make_grm_scl) {
        for (j = 0; j < m; j++) {
            if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
            else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
        }
    }

    int popu = 0;
    vector<eigenVector> sd_SNP_sub(popu_num), mu_SNP_sub(popu_num);
    for(popu  = 0; popu < popu_num; popu++) {
        (mu_SNP_sub[popu]).resize(m);
        (sd_SNP_sub[popu]).resize(m);
        #pragma omp parallel for private(i)
        for (j = 0; j < m; j++){
            double m_sub = 0.0;
            (mu_SNP_sub[popu])(j) = 0.0;
            for(i = 0; i < n; i++){
                if(X(i,j) < 1e5 && S(i,popu)>0.0){
                    m_sub += 1.0;
                    (mu_SNP_sub[popu])(j) += X(i,j);
                }
            }
            (mu_SNP_sub[popu])(j) /= m_sub;
            if(make_grm_scl){
                if(_dosage_flag){
                    (sd_SNP_sub[popu])(j) = 0.0;
                    for(i = 0; i < n; i++){
                        if(X(i,j) < 1e5 && S(i,popu)>0.0) (sd_SNP_sub[popu])(j) += (X(i,j) - (mu_SNP_sub[popu])(j)) * (X(i,j) - (mu_SNP_sub[popu])(j));
                    }
                    (sd_SNP_sub[popu])(j) /= (m_sub - 1.0);                    
                }
                else (sd_SNP_sub[popu])(j) = (mu_SNP_sub[popu])(j)*(1.0 - 0.5 * (mu_SNP_sub[popu])(j));
                if ((sd_SNP_sub[popu])(j) < 1.0e-50) (sd_SNP_sub[popu])(j) = 0.0;
                else (sd_SNP_sub[popu])(j) = 1.0 / sqrt((sd_SNP_sub[popu])(j));
            }
        }
    }

    for(popu  = 0; popu < popu_num; popu++) {
        #pragma omp parallel for private(j)
        for(i = 0; i < n; i++){
            if(S(i,popu)>0.0){
                for (j = 0; j < m; j++){
                    if(X(i,j) < 1e5){
                        X(i,j) -=  (mu_SNP_sub[popu])(j);
                        if(make_grm_scl) X(i,j) *= powf((sd_SNP_sub[popu])(j), make_grm_scl);
                    }
                    else if (miss_with_mu) X(i,j) = 0.0;
                }
            }
        }
    }

    if (!grm_xchr_flag) return;

    // for the X-chromosome
    check_sex();
    double f_buf = sqrt(0.5);

    #pragma omp parallel for private(j)
    for (i = 0; i < n; i++) {
        if (_sex[_keep[i]] == 1) {
            for (j = 0; j < m; j++) {
                if (X(i,j) < 1e5) X(i,j) *= f_buf;
                else if (miss_with_mu) X(i,j) = 0.0;
            }
        }
    }
}

