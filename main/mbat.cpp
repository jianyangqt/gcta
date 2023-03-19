/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementation of gene-based association test (GBAT) in GCTA
 *
 * 2013 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues> 
#include <iostream>

//When compile under windows, M_PI is defined in corecrt_math_defines.h
#if defined _WIN32 || defined _WIN64
#include <corecrt_math_defines.h>
#endif


using namespace Eigen;
using namespace std;



void gcta::mbat_calcu_lambda(vector<int> &snp_indx, MatrixXf &rval, VectorXd &eigenval, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = snp_indx.size();

    MatrixXf X;
    make_XMat_subset(X, snp_indx, false);
    vector<int> rm_ID1;
    double R_cutoff = sbat_ld_cutoff;
    int qi = 0; //alternate index

    VectorXd sumsq_x(m);
    for (j = 0; j < m; j++) sumsq_x[j] = X.col(j).dot(X.col(j));

    MatrixXf C = X.transpose() * X;
    X.resize(0,0);
    #pragma omp parallel for private(j)
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            double d_buf = sqrt(sumsq_x[i] * sumsq_x[j]);
            if(d_buf>0.0) C(i,j) /= d_buf;
            else C(i,j) = 0.0;
        }
    }
    
    rval = C;
    if (sbat_ld_cutoff < 1) rm_cor_sbat(C, R_cutoff, m, rm_ID1);
        //Create new index
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sub_indx.push_back(i);
            else {
                if (rm_ID1[qi] == i) qi++; //Skip removed snp
                else sub_indx.push_back(i);
            }
        }
        snp_count = sub_indx.size();
        if (sub_indx.size() < C.size()) { //Build new matrix
            MatrixXf D(sub_indx.size(),sub_indx.size());
            for (i = 0 ; i < sub_indx.size() ; i++) {
               for (j = 0 ; j < sub_indx.size() ; j++) {
                   D(i,j) = C(sub_indx[i],sub_indx[j]);
               }
            }
            C = D; 
        }
    
    SelfAdjointEigenSolver<MatrixXf> saes(C);
    eigenval = saes.eigenvalues().cast<double>();
}


void gcta::mbat_ACATO(double &mbat_svd_pvalue,double &fastbat_pvalue, double &P_mBATcombo){
    double cctStat =0;
    Vector2d combinedPvalue(mbat_svd_pvalue ,fastbat_pvalue);
    Vector2i isSmall(0,0);
    for (int k =0;k <2; k++)
    {
        if(combinedPvalue[k] ==1 ) combinedPvalue[k] = 1 - 1e-16;
        if(combinedPvalue[k] < 1e-16 ) isSmall[k] = 1 ;
    }

    if(isSmall.sum() == 0){for(int k =0;k <2; k++) cctStat +=   tan((0.5-combinedPvalue[k])* M_PI );}
    else if(isSmall.sum() == 1){
        for(int k =0;k <2; k++) {
            if (isSmall[k] ==0) cctStat =   (1 / combinedPvalue[1-k])/ M_PI +  tan((0.5 - combinedPvalue[k] ) * M_PI);
            else cctStat =   (1 / combinedPvalue[k])/ M_PI +  tan((0.5 - combinedPvalue[1-k] ) * M_PI);
        }
    }else {for(int k =0;k <2; k++) cctStat +=   (1 / combinedPvalue[k])/ M_PI;}
    cctStat = cctStat / 2;

    if(cctStat > 1e+15) P_mBATcombo = (1/cctStat) /M_PI;
    else P_mBATcombo = 1- (((1.0L / M_PI) * atan((cctStat - 0L)/1.0L) + 0.5L));
}

void gcta::svdDecomposition( MatrixXf &X,double &prop, int &eigenvalueNum, VectorXd &eigenvalueUsed,MatrixXd &U_prop){
    VectorXd cumsumNonNeg; // cumulative sums of non-negative values
    double sumNonNeg;

    SelfAdjointEigenSolver<MatrixXf> eigensolver(X);
    VectorXd eigenVal = eigensolver.eigenvalues().cast<double>();
    MatrixXd eigenVec = eigensolver.eigenvectors().cast<double>();
    int revIdx = eigenVal.size();
    cumsumNonNeg.resize(revIdx);
    cumsumNonNeg.setZero(); 
    revIdx = revIdx -1;
    if(eigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = eigenVal(revIdx);
    sumNonNeg = eigenVal(revIdx);
    revIdx = revIdx -1;
    
    while( eigenVal(revIdx) > 1e-10 ){
        sumNonNeg = sumNonNeg + eigenVal(revIdx);
        cumsumNonNeg(revIdx) = eigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
        if(revIdx < 0) break; 
    }
    cumsumNonNeg = cumsumNonNeg * 1/sumNonNeg;
    bool haveValue = false;
    for (revIdx = revIdx + 1; revIdx < eigenVal.size(); revIdx ++ ){
        if(prop >= cumsumNonNeg(revIdx) ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    if(!haveValue) revIdx = eigenVal.size() - 1;
    // block matrix;

    U_prop = eigenVec.rightCols(eigenVal.size() - revIdx);
    eigenvalueUsed = eigenVal.tail(eigenVal.size() - revIdx);
    eigenvalueNum = eigenVal.size() - revIdx;
}

void gcta::mbat_gene(string mbat_sAssoc_file, string mbat_gAnno_file, int mbat_wind, double mbat_svd_gamma, double sbat_ld_cutoff, bool mbat_write_snpset, bool GC, double GC_val,bool mbat_print_all_p)
{
    int i = 0, j = 0;
    vector<int> snp_chr;
    vector<string> gene_name;
    vector<int> gene_chr, gene_bp1, gene_bp2;
    /////////////////////////////////////////
    // Step 1. Read snp and gene files
    /////////////////////////////////////////
    //Step 1.1  Read SNP association results based cojo ma format
    // GC:  If this option is specified, p-values will be adjusted by the genomic control method.
    // By default, the genomic inflation factor will be calculated from the summary-level statistics 
    // of all the SNPs unless you specify a value, e.g. --cojo-gc 1.05.
    cout << " StatFunc::pchisq(1e100, 1):" << to_string( StatFunc::pchisq(1e100, 1)) << " "<< StatFunc::pchisq(1e100, 1) << endl;
    init_massoc(mbat_sAssoc_file, GC, GC_val);
    int snp_num = _include.size();
    // re-calculate chi-square 
    vector<double> snp_chisq(snp_num);
    for (i = 0; i < snp_num; i++) snp_chisq[i] = StatFunc::qchisq(_pval[i], 1);

    map<int, string> chr_begin_snp, chr_end_snp;
    chr_begin_snp.insert(pair<int, string>(_chr[_include[0]], _snp_name[_include[0]]  ));
    for (i = 1; i < snp_num; i++) {
        if (_chr[_include[i]] != _chr[ _include[i - 1]]) {
            chr_begin_snp.insert(pair<int, string>(_chr[_include[i]], _snp_name[_include[i]] ));
            chr_end_snp.insert(pair<int, string>(_chr[_include[i-1]], _snp_name[_include[i -1]] ));
        }
    }
    chr_end_snp.insert(pair<int, string>(_chr[_include[snp_num -1]], _snp_name[_include[snp_num - 1]]));
    //Step 1.2  Read gene file
    sbat_read_geneAnno(mbat_gAnno_file, gene_name, gene_chr, gene_bp1, gene_bp2);

    /////////////////////////////////////////
    // Step 2. Map snps to genes
    /////////////////////////////////////////
    LOGGER << "Mapping the physical positions of genes to SNP data (gene boundaries: " << mbat_wind / 1000 << "Kb away from UTRs) ..." << endl;

    vector<int> gene_mapped_idx;
    int gene_num = gene_name.size();
    vector<string> gene2snp_1(gene_num), gene2snp_2(gene_num);
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chr_iter;
    vector<locus_bp> snp_vec;
    for (i = 0; i < snp_num; i++){
        snp_vec.push_back(locus_bp(_snp_name[_include[i]], _chr[_include[i]], _bp[_include[i]]));
    }
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        // find lowest snp_name in the gene
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp1[i] - mbat_wind));
        if (iter != snp_vec.end()) gene2snp_1[i] = iter->locus_name;
        else gene2snp_1[i] = "NA";
    }
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] == "NA") {
            gene2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp2[i] + mbat_wind));
        if (iter != snp_vec.end()){
            if (iter->bp ==  gene_bp2[i] + mbat_wind){
                gene2snp_2[i] = iter->locus_name;
            }else {
                if(iter!=snp_vec.begin()){
                    iter--;
                    gene2snp_2[i] = iter->locus_name;
                }
                else gene2snp_2[i] = "NA";
            }
        }
        else {
            chr_iter = chr_end_snp.find(gene_chr[i]);
            if (chr_iter == chr_end_snp.end()) gene2snp_2[i] = "NA";
            else gene2snp_2[i] = chr_iter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] != "NA" && gene2snp_2[i] != "NA") 
        {
            mapped++;
            gene_mapped_idx.push_back(i);
        }
    }
    if (mapped < 1) LOGGER.e(0, "no gene can be mapped to the SNP data. Please check the input data regarding chromosome and bp.");
    else LOGGER << mapped << " genes have been mapped to SNP data." << endl;

    /////////////////////////////////////////
    // Step 3. Run gene-based test.
    /////////////////////////////////////////
    if (_mu.empty()) calcu_mu();
    LOGGER << "\nRunning mBAT-combo analysis for " << mapped << " gene(s) ..." << endl;

    LOGGER << "For computing fastBAT test statistics, SNPs are pruned with LD rsq cutoff of = " << sbat_ld_cutoff * sbat_ld_cutoff 
        << " (set by --fastBAT-ld-cutoff 0.9 as default)." << endl;
    LOGGER << "For computing mBAT test statistics, top principal components are selected to explain at least "<< mbat_svd_gamma * 100
        << "% of variance in LD (set by --mBAT-svd-gamma 0.9 as default)." << endl;
    //int mapped = gene_mapped_idx.size();
    vector<string> min_snp_name(mapped);
    vector<double> min_snp_pval(mapped);
    vector<int> snp_num_in_gene(mapped);
    vector<int> snp_num_in_gene_mBAT(mapped);
    map<string, int>::iterator iter1, iter2;
    string rgoodsnpfile = _out + ".gene.snpset.mbat";
    ofstream rogoodsnp;
    if (mbat_write_snpset) {
        rogoodsnp.open(rgoodsnpfile.c_str());
        // rogoodsnp << "gene\tsnp" << endl;
    }
    ////////////////////////////////////////////////////
    VectorXd P_mBATcombo, P_mbat_svd_prop, fastbat_gene_pvalue, Chisq_mBAT,chisq_o;
    P_mBATcombo.resize(mapped);
    P_mbat_svd_prop.resize(mapped);
    fastbat_gene_pvalue.resize(mapped);
    Chisq_mBAT.resize(mapped);
    chisq_o.resize(mapped);
    map<int,int> include_in_gene; 
    vector<int> eigenvalueNum_mBAT(mapped);
    int gene_analyzed=0;

    /// save output results
    string filename = _out + ".gene.assoc.mbat";
    LOGGER << "Saving the results of the mBAT analysis to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Gene\tChr\tStart\tEnd\tNo.SNPs\tSNP_start\tSNP_end\tTopSNP\tTopSNP_Pvalue\tNo.Eigenvalues\tChisq_mBAT\tP_mBATcombo";
    if(mbat_print_all_p){
        ofile << "\tP_mBAT\tChisq_fastBAT\tP_fastBAT" << endl;
    } else {
        ofile << endl;
    }

    // Dectecting gene with snp > 10,000
    for (j = 0; j < mapped; j++) {
        include_in_gene.clear();
        int gene_ori_idx = gene_mapped_idx[j];
        vector<int> snp_indx;
        iter1 = _snp_name_map.find(gene2snp_1[gene_ori_idx]);
        iter2 = _snp_name_map.find(gene2snp_2[gene_ori_idx]);
        bool skip = false;
        if (iter1 == _snp_name_map.end() || iter2 == _snp_name_map.end() || iter1->second >= iter2->second) skip = true;
        // 
        std::vector<int>::iterator iter_include;
        int idx_include_in_gene;
        iter_include = std::find(_include.begin(), _include.end(), iter1->second);
        if(iter_include != _include.cend()){
           idx_include_in_gene = std::distance(_include.begin(),iter_include); 
        }
        for(int k = idx_include_in_gene; k < snp_num; k ++){
                if(_include[k] > iter2->second){break;}
                snp_indx.push_back(_include[k]);
                include_in_gene.insert(pair<int,int>(_include[k],k));
        }

        snp_num_in_gene[j] = snp_indx.size();  // assume snp_name_index is ordered and continuous
        snp_num_in_gene_mBAT[j] = snp_num_in_gene[j];
        if(!skip && snp_num_in_gene[j] > 10000){
            LOGGER << "The following genes have > 10000 SNPs which may take at least a few minutes for each test. You may consider removing these genes from the input list if it takes too much time." << endl;
            // for (int k = 0; k < snp_indx.size(); k++) {
            LOGGER << gene_name[gene_ori_idx] << endl; // "\t"<< _snp_name[snp_indx[k]]<< endl;
            // } 
        } 
    }


    // Do mBAT analysis for each gene
    for (j = 0; j < mapped; j++) {
        include_in_gene.clear();
        int gene_ori_idx = gene_mapped_idx[j];
        vector<int> snp_indx;
        iter1 = _snp_name_map.find(gene2snp_1[gene_ori_idx]);
        iter2 = _snp_name_map.find(gene2snp_2[gene_ori_idx]);
        bool skip = false;
        if (iter1 == _snp_name_map.end() || iter2 == _snp_name_map.end() || iter1->second >= iter2->second) skip = true;
        // 
        std::vector<int>::iterator iter_include;
        int idx_include_in_gene;
        iter_include = std::find(_include.begin(), _include.end(), iter1->second);
        if(iter_include != _include.cend()){
           idx_include_in_gene = std::distance(_include.begin(),iter_include);
        }
        for(int k = idx_include_in_gene; k < snp_num; k ++){
                if(_include[k] > iter2->second){break;}
                snp_indx.push_back(_include[k]);
                include_in_gene.insert(pair<int,int>(_include[k],k));
        }

        snp_num_in_gene[j] = snp_indx.size();  // assume snp_name_index is ordered and continuous
        snp_num_in_gene_mBAT[j] = snp_num_in_gene[j];
        // if(!skip && snp_num_in_gene[j] > 10000){
        //    LOGGER<<"Warning: Too many SNPs in the gene region [" << gene_name[gene_ori_idx] << "]. Maximum limit is 10000. This gene is ignored in the analysis."<<endl;
           // skip = true;  
        // } 
        if(skip){
            P_mBATcombo[j] = 2.0;
            snp_num_in_gene[j] = 0;
            continue;
        }
     
        if (mbat_write_snpset) {
            rogoodsnp << gene_name[gene_ori_idx] << endl;
            for (int k = 0; k < snp_indx.size(); k++) {
                // rogoodsnp << gene_name[gene_ori_idx] << "\t"<< _snp_name[snp_indx[k]] << endl;
                rogoodsnp << _snp_name[snp_indx[k]] << endl;
            }
            rogoodsnp << "END" << endl << endl;
        }
        // Step 3.1 run mBAT analysis
        // calculate LD matrix 
        // Step 3.2 run fastBAT
        // actually both mBAT and fatBAT use LD matrix, it may be possible to use one LD matrix to speed up later.
        int snp_count;
        vector<int> snp_include_gene;
        if (_mu.empty()) calcu_mu();
        chisq_o[j] = 0;
        for(int k = 0; k < snp_indx.size(); k++) {
              chisq_o[j] += snp_chisq[include_in_gene[snp_indx[k]]];
             snp_include_gene.push_back(include_in_gene[snp_indx[k]]);
        }
        if(snp_num_in_gene[j] == 1) {
            // consider fastbat
            fastbat_gene_pvalue[j] = StatFunc::pchisq(chisq_o[j], 1.0);
            
            // mbat
            float U_prop1;
            float eigenvalueUsed1;
            float zscore_in_set1;
            float Uprop_z1;
            float lambda_prop_diag_inv1;
            // svdDecomposition(rval, mbat_svd_gamma,eigenvalueNum_mBAT[i],eigenvalueUsed, U_prop);
            U_prop1 =1.0;
            eigenvalueNum_mBAT[i] = 1;
            eigenvalueUsed1 = 0.9999997;
            // zscore_in_set.resize(snp_indx.size());
            // // min_snp_pval[j]=2;
            // // min_snp_name[j]="na";
            // for(int k = 0; k < snp_indx.size(); k++){
            //     zscore_in_set[k] = beta[snp_indx[k]]/beta_se[snp_indx[k]];
            // }
            int k = 0;
            zscore_in_set1 = _beta[include_in_gene[snp_indx[k]]]/_beta_se[include_in_gene[snp_indx[k]]];

            Uprop_z1 = U_prop1 * zscore_in_set1; 
            lambda_prop_diag_inv1 = 1 / eigenvalueUsed1;        
            Chisq_mBAT[i] = Uprop_z1 * lambda_prop_diag_inv1 * Uprop_z1;
            P_mbat_svd_prop[i] = StatFunc::pchisq(Chisq_mBAT[i], eigenvalueNum_mBAT[i]);

        } else {
            snp_count=snp_num_in_gene[j];
            VectorXd eigenval;
            vector<int> sub_indx;
            MatrixXf rval(snp_indx.size(), snp_indx.size());;
            mbat_calcu_lambda(snp_include_gene, rval, eigenval, snp_count, sbat_ld_cutoff, sub_indx);
            ///////////////////////////////////////////////////////////////////////////////
            // Step 3.1 run mBAT analysis
            MatrixXd U_prop;
            VectorXd eigenvalueUsed;
            svdDecomposition(rval, mbat_svd_gamma,eigenvalueNum_mBAT[j],eigenvalueUsed, U_prop);
            VectorXd zscore_in_gene;
            zscore_in_gene.resize(snp_indx.size());
            min_snp_pval[j]=2;
            min_snp_name[j]="na";

            for(int k = 0; k < snp_indx.size(); k++){
                if (min_snp_pval[j] > _pval[include_in_gene[snp_indx[k]]]) { 
                    min_snp_pval[j] = _pval[include_in_gene[snp_indx[k]]];
                    min_snp_name[j] = _snp_name[snp_indx[k]]; //keep minimum value - regardless of whether SNP removed by LD pruning
                }
                zscore_in_gene[k] = _beta[include_in_gene[snp_indx[k]]]/_beta_se[include_in_gene[snp_indx[k]]];
            }

            MatrixXd Uprop_z = U_prop.transpose() * zscore_in_gene; 
            MatrixXd lambda_prop_diag_inv = eigenvalueUsed.asDiagonal().inverse();        
            Chisq_mBAT[j] = (Uprop_z.transpose() * lambda_prop_diag_inv * Uprop_z)(0);
            P_mbat_svd_prop[j] = StatFunc::pchisq(Chisq_mBAT[j], eigenvalueNum_mBAT[j]);
            ///////////////////////////////////////////////////////////////////////////////
            //Step 3.2 recalculate chisq value from low correlation snp subset
            if (sbat_ld_cutoff < 1) {
                chisq_o[j] = 0;
                for (int k = 0; k < sub_indx.size(); k++){
                    chisq_o[j] += snp_chisq[snp_include_gene[ sub_indx[k]] ];
                }

            } 
            snp_num_in_gene[j] = snp_count;
            if (snp_count==1 && chisq_o[j] ==0){
                fastbat_gene_pvalue[j] = 1;
            }else {
                fastbat_gene_pvalue[j] = StatFunc::pchisqsum(chisq_o[j], eigenval);
            }
        }
        // Step 3.3 Combine fastbat and mbat p value to obtain mbat-combo p
        mbat_ACATO(P_mbat_svd_prop[j],fastbat_gene_pvalue[j],P_mBATcombo[j]);
        
        /////////////////////////////////////////
        // Step 4. output
        /////////////////////////////////////////
        if (P_mBATcombo[j] > 1.5) continue;
        gene_analyzed++;
        ofile << gene_name[gene_ori_idx] << "\t" << gene_chr[gene_ori_idx] << "\t" << gene_bp1[gene_ori_idx] << "\t" << gene_bp2[gene_ori_idx] << "\t";
        ofile << snp_num_in_gene_mBAT[j] << "\t" << gene2snp_1[gene_ori_idx] << "\t" << gene2snp_2[gene_ori_idx]  << "\t";
        ofile << min_snp_name[j] << "\t" << min_snp_pval[j] << "\t" << eigenvalueNum_mBAT[j] << "\t" <<  Chisq_mBAT[j] << "\t" << P_mBATcombo[j];
        
        if(mbat_print_all_p){ofile << "\t" << P_mbat_svd_prop[j] << "\t" << chisq_o[j]  << "\t" << fastbat_gene_pvalue[j];}
        ofile << endl;
       // LOGGER << "mBAT-combo analysis for " << gene_ori_idx+1 << "/"<< gene_num << "-th gene [" << gene_name[gene_ori_idx] << "] has been completed." << endl;
    }
    ofile.close();
    rogoodsnp.close();
}


void gcta::mbat(string mbat_sAssoc_file, string snpset_file, double mbat_svd_gamma, double sbat_ld_cutoff, bool mbat_write_snpset, bool GC, double GC_val,bool mbat_print_all_p)
{
    /////////////////////////////////////////
    // Step 1. read snp set file and ma file.
    /////////////////////////////////////////
    int i = 0, j = 0;
    int snp_count;

    // read SNP set file
    vector<string> set_name;
    vector< vector<string> > snpset;
    sbat_read_snpset(snpset_file, set_name, snpset);
    int set_num = set_name.size();
    //cout << "snpset size: " << snpset[0].size();

    // read SNP association results
    vector<string> snp_name;
    vector<int> snp_chr, snp_bp;
    vector<double> snp_pval;
    vector<double> beta;
    vector<double> beta_se;
    init_massoc(mbat_sAssoc_file, GC, GC_val);
    int snp_num = _include.size();
    // re-calculate chi-square 
    vector<double> snp_chisq(snp_num);
    snp_pval.resize(snp_num);
    snp_name.resize(snp_num);
    snp_chr.resize(snp_num);
    snp_bp.resize(snp_num);
    beta.resize(snp_num);
    beta_se.resize(snp_num);
    for (i = 0; i < snp_num; i++) {
        snp_name[i] = _snp_name[_include[i]];
        snp_chr[i]  = _chr[_include[i]];
        snp_bp[i] =  _bp[_include[i]];
        snp_pval[i] = _pval[i];
        beta[i]  = _beta[i];
        beta_se[i] = _beta_se[i];
        snp_chisq[i] = StatFunc::qchisq(_pval[i], 1);
    }
    /////////////////////////////////////////
    // Step 2. Map snp to gene cis-region.
    /////////////////////////////////////////
    // we don't have to do this like mbat_gene,

    /////////////////////////////////////////
    // Step 3. Run gene-based test.
    /////////////////////////////////////////
    if (_mu.empty()) calcu_mu();
    LOGGER << "\nRunning mBAT-combo analysis for " << set_num << " gene(s) ..." << endl;

    LOGGER << "For computing fastBAT test statistics, SNPs are pruned with LD rsq cutoff of = " << sbat_ld_cutoff * sbat_ld_cutoff 
        << " (set by --fastBAT-ld-cutoff 0.9 as default)." << endl;
    LOGGER << "For computing mBAT test statistics, top principal components are selected to explain at least "<< mbat_svd_gamma * 100
        << "% of variance in LD (set by --mBAT-svd-gamma 0.9 as default)." << endl;
    
    vector<string> min_snp_name(set_num);
    vector<double> min_snp_pval(set_num);
    vector<int> snp_num_in_set(set_num); // = snp_num_in_set
    vector<int> snp_num_in_set_mBAT(set_num); // = snp_num_in_set_mBAT;
    // map<string, int>::iterator iter1, iter2;
    map<string, int>::iterator iter;
    map<string, int> snp_name_map;
    string rgoodsnpfile = _out + ".set.snpset.mbat";
    ofstream rogoodsnp;
    if (mbat_write_snpset) {
        rogoodsnp.open(rgoodsnpfile.c_str());
        rogoodsnp << "set\tsnp" << endl;
    }
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(pair<string,int>(snp_name[i], i));
    ////////////////////////////////////////////////////
    VectorXd P_mBATcombo, P_mbat_svd_prop, fastbat_set_pvalue, Chisq_mBAT,chisq_o;
    P_mBATcombo.resize(set_num);
    P_mbat_svd_prop.resize(set_num);
    fastbat_set_pvalue.resize(set_num);
    Chisq_mBAT.resize(set_num);
    chisq_o.resize(set_num);
    map<int,int> include_in_gene; 
    vector<int> eigenvalueNum_mBAT(set_num);
    int gene_analyzed=0;

    /// save output results
    string filename = _out + ".set.assoc.mbat";
    LOGGER << "Saving the results of the mBAT analysis to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Set\tNo.SNPs\tTopSNP\tTopSNP_Pvalue\tNo.Eigenvalues\tChisq_mBAT\tP_mBATcombo";
    if(mbat_print_all_p){
        ofile << "\tP_mBAT\tChisq_fastBAT\tP_fastBAT" << endl;
    } else {
        ofile << endl;
    }


    // Dectecting gene with snp > 10,000
    for (i = 0; i < set_num; i++) {
        bool skip = false;
        if(snpset[i].size() < 1) skip = true;
        vector<int> snp_indx;
        for(j = 0; j < snpset[i].size(); j++){
            iter = snp_name_map.find(snpset[i][j]);
            if(iter!=snp_name_map.end()) snp_indx.push_back(iter->second);
        }
        snp_num_in_set[i] = snp_indx.size();
        snp_num_in_set_mBAT[i] = snp_indx.size();
        if(!skip && snp_num_in_set[i] > 10000){
            LOGGER << "The following sets have > 10000 SNPs which may take at least a few minutes for each test. You may consider removing these genes from the input list if it takes too much time." << endl;
            // for (int k = 0; k < snp_indx.size(); k++) {
            LOGGER << set_name[i]  << endl; // "\t" << snp_name[snp_indx[k]] << endl;
            // }
            // skip = true;  
        } 
    }


    // Do mBAT analysis for each set
    for (i = 0; i < set_num; i++) {
        bool skip = false;
        if(snpset[i].size() < 1) skip = true;
        vector<int> snp_indx;
        for(j = 0; j < snpset[i].size(); j++){
            iter = snp_name_map.find(snpset[i][j]);
            if(iter!=snp_name_map.end()) snp_indx.push_back(iter->second);
        }
        snp_num_in_set[i] = snp_indx.size();
        snp_num_in_set_mBAT[i] = snp_indx.size();
        // if(!skip && snp_num_in_set[i] > 20000){
        //     LOGGER<<"Warning: Too many SNPs in the set ["<<set_name[i]<<"]. Maximum limit is 20000. This gene is ignored in the analysis."<<endl;
        //     skip = true;  
        // } 
        if(skip){
            P_mBATcombo[i] = 2.0;
            snp_num_in_set[i] = 0;
            continue;
        }
     
        // Step 3.1 run mBAT analysis
        // calculate LD matrix 
        // Step 3.2 run fastBAT
        // actually both mBAT and fatBAT use LD matrix, it may be possible to use one LD matrix to speed up later.
        chisq_o[i] = 0;
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = 0; j < snp_indx.size(); j++) {
            if (min_snp_pval[i] > snp_pval[snp_indx[j]]) {
                min_snp_pval[i] = snp_pval[snp_indx[j]];
                min_snp_name[i] = snp_name[snp_indx[j]];
            }
            chisq_o[i] += snp_chisq[snp_indx[j]];
        }

        if(snp_num_in_set[i] == 1) {
            // consider fastbat
            fastbat_set_pvalue[i] = StatFunc::pchisq(chisq_o[i], 1.0);
            // mbat
            float U_prop1;
            float eigenvalueUsed1;
            float zscore_in_set1;
            float Uprop_z1;
            float lambda_prop_diag_inv1;
            // svdDecomposition(rval, mbat_svd_gamma,eigenvalueNum_mBAT[i],eigenvalueUsed, U_prop);
            U_prop1 =1.0;
            eigenvalueNum_mBAT[i] = 1;
            eigenvalueUsed1 = 0.9999997;
            // zscore_in_set.resize(snp_indx.size());
            // // min_snp_pval[j]=2;
            // // min_snp_name[j]="na";
            // for(int k = 0; k < snp_indx.size(); k++){
            //     zscore_in_set[k] = beta[snp_indx[k]]/beta_se[snp_indx[k]];
            // }
            int k = 0;
            zscore_in_set1 = beta[snp_indx[k]]/beta_se[snp_indx[k]];

            Uprop_z1 = U_prop1 * zscore_in_set1; 
            lambda_prop_diag_inv1 = 1 / eigenvalueUsed1;        
            Chisq_mBAT[i] = Uprop_z1 * lambda_prop_diag_inv1 * Uprop_z1;
            P_mbat_svd_prop[i] = StatFunc::pchisq(Chisq_mBAT[i], eigenvalueNum_mBAT[i]);
        } else {
            snp_count=snp_num_in_set[i];
            VectorXd eigenval;
            vector<int> sub_indx;
            MatrixXf rval(snp_indx.size(), snp_indx.size());;
            // extract eigenval for fastbat
            mbat_calcu_lambda(snp_indx, rval, eigenval, snp_count, sbat_ld_cutoff, sub_indx);
            ///////////////////////////////////////////////////////////////////////////////
            // Step 3.1 run mBAT analysis
            MatrixXd U_prop;
            VectorXd eigenvalueUsed;
            VectorXd zscore_in_set;
            MatrixXd Uprop_z;
            MatrixXd lambda_prop_diag_inv;
            svdDecomposition(rval, mbat_svd_gamma,eigenvalueNum_mBAT[i],eigenvalueUsed, U_prop);
            zscore_in_set.resize(snp_indx.size());
            // min_snp_pval[j]=2;
            // min_snp_name[j]="na";

            for(int k = 0; k < snp_indx.size(); k++){
                zscore_in_set[k] = beta[snp_indx[k]]/beta_se[snp_indx[k]];
            }

            Uprop_z = U_prop.transpose() * zscore_in_set; 
            lambda_prop_diag_inv = eigenvalueUsed.asDiagonal().inverse();        
            Chisq_mBAT[i] = (Uprop_z.transpose() * lambda_prop_diag_inv * Uprop_z)(0);
            P_mbat_svd_prop[i] = StatFunc::pchisq(Chisq_mBAT[i], eigenvalueNum_mBAT[i]);
            ///////////////////////////////////////////////////////////////////////////////
            //Step 3.2 recalculate chisq value from low correlation snp subset
            if (sbat_ld_cutoff < 1) {
                chisq_o[i] = 0;
                for (int k = 0; k < sub_indx.size(); k++){
                    chisq_o[i] += snp_chisq[snp_indx[sub_indx[k]] ];
                }

            } 
            snp_num_in_set[i] = snp_count;
            if (snp_count==1 && chisq_o[i] ==0){
                fastbat_set_pvalue[i] = 1;
            }else {
                fastbat_set_pvalue[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);
            }
        
            if (mbat_write_snpset) {
                rogoodsnp << set_name[i] << endl;
                for (int k = 0; k < snp_indx.size(); k++) {
                    rogoodsnp << snp_name[snp_indx[k]] << endl;
                }
                rogoodsnp << "END" << endl << endl;
            }
        }
        // Step 3.3 Combine fastbat and mbat p value to obtain mbat-combo p
        mbat_ACATO(P_mbat_svd_prop[i],fastbat_set_pvalue[i],P_mBATcombo[i]);
        
        /////////////////////////////////////////
        // Step 4. output
        /////////////////////////////////////////
        if (P_mBATcombo[i] > 1.5) continue;
        gene_analyzed++;
        ofile << set_name[i] << "\t"  << snp_num_in_set_mBAT[i] << "\t";
        ofile << min_snp_name[i] << "\t" << min_snp_pval[i] << "\t" << eigenvalueNum_mBAT[i] << "\t" <<  Chisq_mBAT[i] << "\t" << P_mBATcombo[i];
        
        if(mbat_print_all_p){ofile << "\t" << P_mbat_svd_prop[i] << "\t" << chisq_o[i]  << "\t" << fastbat_set_pvalue[i];}
        ofile << endl;
       // LOGGER << "mBAT-combo analysis for " << gene_ori_idx+1 << "/"<< gene_num << "-th gene [" << gene_name[gene_ori_idx] << "] has been completed." << endl;
    }
    ofile.close();
    rogoodsnp.close();
}

