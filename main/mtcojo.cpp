#include "gcta.h"
#include "Logger.h"
#include <limits>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>

double abnormal = -9999999999;

template <class T>
T o_digital_number(T m) {
    T eps = 1e9;
    return(abs(m) < eps ?  m : (T)nan(""));
}

double quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob) {
    if (prob < 0 or prob > 1) LOGGER.e(0, "Requested quantile probability is invalid");
    if (vals.size() == 0) return std::numeric_limits<double>::quiet_NaN();
    double index = prob * (vals.size()-1);
    unsigned below = std::floor(index), above = std::ceil(index);
    if (below == above) return vals[above];
    return (above - index) * vals[below] + (index - below) * vals[above];
}

double quantile(const std::vector<double> &vals, double prob) {
    return quantile(Eigen::Map<const Eigen::VectorXd>(vals.data(), vals.size()), prob);
}

void update_mtcojo_snp_rm(vector<string> adjsnps, map<string,int> &snp_id_map, vector<int> &include_snp_indx) {
    
    int i=0, nsnpbuf=adjsnps.size();
    std::map<string,int>::iterator iter;
    
    for(i=0; i<nsnpbuf; i++) snp_id_map.erase(adjsnps[i]);

    include_snp_indx.clear();
    for(iter=snp_id_map.begin(), i=0; iter!=snp_id_map.end(); iter++, i++) include_snp_indx.push_back(iter->second);
    stable_sort(include_snp_indx.begin(), include_snp_indx.end());
}

vector<string> read_snp_metafile(string metafile) {
    ifstream meta_snp(metafile.c_str());
    if (!meta_snp)
         LOGGER.e(0, "Can not open the file [" + metafile + "] to read.");
    
    string strbuf="";
    vector<string> snplist;
    int line_number=0;
    // Read the summary data
    while(std::getline(meta_snp, strbuf)) {
        line_number++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8) {
            LOGGER.e(0, "The summary data [" + metafile + "]: the number of elements is not equal to 8 in line " + to_string(line_number) + ".");
        }
        if(line_number==1) continue;
        snplist.push_back(line_elements[0]);
    }
    meta_snp.close();
    return(snplist);
}

void read_multi_metafile(string metafile, map<string, int> id_map,
                         vector<string> &snp_a1, vector<string> &snp_a2,
                         vector<double> &snp_freq, vector<double> &snp_b,
                         vector<double> &snp_se, vector<double> &snp_pval,
                         vector<double> &snp_n) {
    
    ifstream meta_raw(metafile.c_str());
    if (!meta_raw)
        LOGGER.e(0, "Can not open the file [" + metafile + "] to read.");
    string strbuf="", valbuf="";
    int line_number=0, snp_indx=0;
    map<string, int>::iterator iter;
    // Read the summary data
    while(std::getline(meta_raw, strbuf)) {
        line_number++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8) {
            LOGGER.e(0, "The summary data [" + metafile + "]: the number of elements is not equal to 8 in line " + to_string(line_number) + ".");
        }
        // Read the summary data
        if(line_number==1) continue;
        iter = id_map.find(line_elements[0]);
        if(iter == id_map.end()) continue;
        snp_indx = iter->second;
        snp_a1[snp_indx] = line_elements[1]; snp_a2[snp_indx] = line_elements[2];
        StrFunc::to_upper(snp_a1[snp_indx]); StrFunc::to_upper(snp_a2[snp_indx]);
        valbuf =line_elements[3];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_freq[snp_indx] = atof(valbuf.c_str());
        else snp_freq[snp_indx] = abnormal;
        valbuf = line_elements[4];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_b[snp_indx] = atof(valbuf.c_str());
        else snp_b[snp_indx] = abnormal;
        valbuf = line_elements[5];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_se[snp_indx] = atof(valbuf.c_str());
        else snp_se[snp_indx] = abnormal;
        valbuf = line_elements[6];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_pval[snp_indx] = atof(valbuf.c_str());
        else snp_pval[snp_indx] = abnormal;
        valbuf = line_elements[7];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_n[snp_indx] = atof(valbuf.c_str());
        else snp_n[snp_indx] = abnormal;
    }
    meta_raw.close();
}

void gcta::read_mtcojofile(string mtcojolist_file, string &target_pheno, vector<string> &covar_pheno, int &ncovar, vector<string> &snplist, vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, vector<vector<double>> &snp_freq,  vector<vector<double>> &snp_b, vector<vector<double>> &snp_se, vector<vector<double>> &snp_pval, vector<vector<double>> &snp_n, vector<double> &popu_prev, vector<double> &smpl_prev) {

    ifstream meta_list(mtcojolist_file.c_str());
    if (!meta_list)
        LOGGER.e(0, "Can not open the file [" + mtcojolist_file + "] to read.");

    string strbuf="", target_pheno_file="", prevbuf1="", prevbuf2="";
    double d_prev1 = 0.0, d_prev2 = 0.0;
    vector<string> vs_buf, covar_pheno_file, snplistbuf;

    // Retrieve the GWAS summary data file
    // The 1st row: the target trait
    int line_number = 1;
    std::getline(meta_list, strbuf);
    std::istringstream linebuf(strbuf);
    std::istream_iterator<string> begin_title(linebuf), end_title;
    vector<string> line_elements(begin_title, end_title);
    if(line_elements.size() != 2 && line_elements.size() != 4) {
        LOGGER.e(0, "Format of file [" + mtcojolist_file + "] is not correct, line " + to_string(line_number) + ".");
    }
    target_pheno=line_elements[0]; target_pheno_file=line_elements[1];
    // prevelance
    d_prev1 = abnormal; d_prev2 = abnormal;
    if(line_elements.size()==4) {
        prevbuf1 = line_elements[2]; prevbuf2 = line_elements[3];
        StrFunc::to_upper(prevbuf1); StrFunc::to_upper(prevbuf2);
        // available, 1 - sample prevelance, 2 - population prevelance
        if(prevbuf1 != "NA"  &&  prevbuf1!= "NAN" && prevbuf1!= ".") {
            d_prev1 = atof(prevbuf1.c_str());
            if(d_prev1 <= 0)
                LOGGER.e(0, "Invalid sample prevelance for [" + target_pheno + "].");
        }
        if(prevbuf2 != "NA"  &&  prevbuf2!= "NAN" && prevbuf2 != ".") {
            d_prev2 = atof(prevbuf2.c_str());
            if(d_prev2 <= 0)
                LOGGER.e(0, "Invalid population prevelance for [" + target_pheno + "].");
        }
    }
    smpl_prev.push_back(d_prev1); popu_prev.push_back(d_prev2);
 
    // The rest rows: the covariate traits
    while(std::getline(meta_list, strbuf)) {
        line_number++;
        //std::istringstream linebuf(strbuf);
        linebuf.clear();
        linebuf.str(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 2 && line_elements.size() != 4) {
            LOGGER.e(0, "Format of  line [" + mtcojolist_file + "] is not correct, line " + to_string(line_number) + ".");
        }
        covar_pheno.push_back(line_elements[0]);
        covar_pheno_file.push_back(line_elements[1]);
        // prevelance
        d_prev1 = abnormal; d_prev2 = abnormal;
        if(line_elements.size()==4) {
            prevbuf1 = line_elements[2]; prevbuf2 = line_elements[3];
            StrFunc::to_upper(prevbuf1); StrFunc::to_upper(prevbuf2);
            // available, 1 - sample prevelance, 2 - population prevelance
            if(prevbuf1 != "NA"  &&  prevbuf1!= "NAN" && prevbuf1!= ".") {
                d_prev1 = atof(prevbuf1.c_str());
                if(d_prev1 <= 0)
                    LOGGER.e(0, "Invalid sample prevelance for [" + line_elements[0] + "].");
            }
            if(prevbuf2 != "NA"  &&  prevbuf2!= "NAN" && prevbuf2 != ".") {
                d_prev2 = atof(prevbuf2.c_str());
                if(d_prev2 <= 0)
                    LOGGER.e(0, "Invalid population prevelance for [" + line_elements[0] + "].");
            }
        }
        smpl_prev.push_back(d_prev1); popu_prev.push_back(d_prev2);
    }
    meta_list.close();
    
    // Read the SNPs
    // Target trait
    snplistbuf=read_snp_metafile(target_pheno_file);
    update_id_map_kp(snplistbuf, _snp_name_map, _include);

    // Covariates
    ncovar = line_number-1;
    int i=0;
    for( i=0; i<ncovar; i++) {
        snplistbuf=read_snp_metafile(covar_pheno_file[i]);
        update_id_map_kp(snplistbuf, _snp_name_map, _include);
    }
    
    // Read the association results
    int nsnp = _snp_name_map.size();
    map<string, int>::iterator iter;
    snplist.clear(); snplist.resize(nsnp);
    for(i=0, iter=_snp_name_map.begin(); iter!=_snp_name_map.end(); i++, iter++) snplist[i] = iter->first;

    // initialization
    snp_a1.clear(); snp_a2.clear(); snp_freq.clear();
    snp_b.clear(); snp_se.clear();
    snp_pval.clear(); snp_n.clear();
    snp_a1.resize(ncovar+1); snp_a2.resize(ncovar+1); snp_freq.resize(ncovar+1);
    snp_b.resize(ncovar+1); snp_se.resize(ncovar+1);
    snp_pval.resize(ncovar+1); snp_n.resize(ncovar+1);
    for( i=0; i<=ncovar; i++) {
        snp_a1[i].resize(nsnp); snp_a2[i].resize(nsnp); snp_freq[i].resize(nsnp, abnormal);
        snp_b[i].resize(nsnp, abnormal); snp_se[i].resize(nsnp, abnormal);
        snp_pval[i].resize(nsnp, abnormal); snp_n[i].resize(nsnp, abnormal);
    }
    
    map<string,int> snplist_map;
    for(i=0; i<nsnp; i++) snplist_map.insert(pair<string,int>(snplist[i], i));
    
    // Target trait
    read_multi_metafile(target_pheno_file, snplist_map,  snp_a1[0], snp_a2[0], snp_freq[0], snp_b[0], snp_se[0], snp_pval[0], snp_n[0]);

    // Covariates
    for(i=0; i<ncovar; i++) {
        read_multi_metafile(covar_pheno_file[i], snplist_map,  snp_a1[i+1], snp_a2[i+1], snp_freq[i+1], snp_b[i+1], snp_se[i+1], snp_pval[i+1], snp_n[i+1]);
    }
}

vector<string> remove_bad_snps(vector<string> snplist, vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, vector<vector<double>> &snp_freq,  vector<vector<double>> &snp_b, vector<vector<double>> snp_se, vector<vector<double>> snp_pval, vector<vector<double>> snp_n, vector<string> snp_ref_a1, vector<string> snp_ref_a2, map<string, int> id_map, string target_pheno, vector<string> covar_pheno, int ncovar, string outfile_name) {
    int i=0, j=0, snpindx=0, nsnp = snplist.size();
    double eps_large = -1e5, eps_small = 1e-6;
    vector<string> badsnps2, allelebuf, bad_ref_a1, bad_ref_a2;
    vector<int> bad_indx;
    map<string, int>::iterator iter;
    bool iterFlag=true;
    
    for( i=0; i<nsnp; i++) {
        iterFlag=true;
        allelebuf.clear(); allelebuf.resize(ncovar*2+4);
        // Remove MAF < 0.01
        iter = id_map.find(snplist[i]);
        snpindx = iter -> second;
//        if(snp_ref_freq[snpindx] <= freq_thresh*2 || snp_ref_freq[snpindx] >= (1-freq_thresh)*2) badsnps1.push_back(snplist[i]);
        // Alleles
        allelebuf[0]=snp_ref_a1[snpindx]; allelebuf[1]=snp_ref_a2[snpindx];
        for(j=0; j<=ncovar; j++) {
            // Alleles
            allelebuf[(j+1)*2]=snp_a1[j][i]; allelebuf[(j+1)*2+1]=snp_a2[j][i];
            // Removing SNPs with missing value
            if(snp_b[j][i] < eps_large | snp_se[j][i] < eps_large | snp_pval[j][i] < eps_large | snp_n[j][i] < eps_large) {
                iterFlag=false; break;
            }
            // Removing SNPs with extremely small SE
            if(snp_se[j][i] < eps_small) {
                iterFlag=false; break;
            }
            // Alignment of b
            if(snp_ref_a1[snpindx] != snp_a1[j][i]) {
                snp_b[j][i] = -1*snp_b[j][i];
                snp_freq[j][i] = 1-snp_freq[j][i];
            }
        }
        
        // Remove SNPs with multiple alleles
        stable_sort(allelebuf.begin(), allelebuf.end());
        allelebuf.erase(unique(allelebuf.begin(), allelebuf.end()), allelebuf.end());
        if(allelebuf.size()!=2) iterFlag=false;
        
        // Collect bad SNPs
        if(iterFlag==false) {
            badsnps2.push_back(snplist[i]);
            bad_indx.push_back(i);
            bad_ref_a1.push_back(snp_ref_a1[snpindx]);
            bad_ref_a2.push_back(snp_ref_a2[snpindx]);
        }
    }
    
    if (!badsnps2.empty()) {
        string badsnpfile = outfile_name + ".badsnps", strbuf="";
        ofstream obadsnp(badsnpfile.c_str());
        obadsnp << "SNP\tRef_A1\tRef_A2";
        obadsnp << "\t" << target_pheno << "_A1\t" << target_pheno << "_A2\t" << target_pheno << "_b\t" << target_pheno << "_se\t" << target_pheno << "_p\t" << target_pheno << "_N";
        for(i=0; i<ncovar; i++) {
            obadsnp << "\t" << covar_pheno[i] << "_A1\t" << covar_pheno[i] << "_A2\t" << covar_pheno[i] << "_b\t" << covar_pheno[i] << "_se\t" << covar_pheno[i] << "_p\t" << covar_pheno[i] << "_N";
        }
        obadsnp << endl;
        int nbadsnps2 = badsnps2.size();
        for (i = 0; i < nbadsnps2; i++) {
            obadsnp << badsnps2[i] <<"\t" << bad_ref_a1[i] << "\t" << bad_ref_a2[i];
            for(j=0; j<=ncovar; j++) {
                obadsnp <<"\t" <<snp_a1[j][bad_indx[i]]<<"\t" <<snp_a2[j][bad_indx[i]]<<"\t" ;
                if(snp_a1[j][bad_indx[i]]==bad_ref_a1[i]) obadsnp<<o_digital_number(snp_b[j][bad_indx[i]]);
                else obadsnp<<o_digital_number(-1*snp_b[j][bad_indx[i]]);
                obadsnp<<"\t" <<o_digital_number(snp_se[j][bad_indx[i]])<<"\t" <<o_digital_number(snp_pval[j][bad_indx[i]])<<"\t" <<o_digital_number(snp_n[j][bad_indx[i]]);
            }
            obadsnp << endl;
        }
        obadsnp.close();
        LOGGER.i(0,  to_string(nbadsnps2) + " SNPs have missing value or mismatch alleles. These SNPs have been saved in [" + badsnpfile + "].");
    }
    
    stable_sort(badsnps2.begin(), badsnps2.end());
    return(badsnps2);
}

vector<string> gcta::clumping_meta(vector<string> snplist, vector<int> include_snp, vector<double> snp_pval, double pval_thresh1, double pval_thresh2, int wind_size, double r2_thresh) {
    
    wind_size = wind_size*1e3;
    
    vector<pair<double, int>> snp_pvalbuf;
    int i=0, j=0, indx = 0, nsnp = include_snp.size();
    double pvalbuf;
    string snpbuf = "", snpbuf_left="", snpbuf_right="";
    
    // Sort the p-value
    for(i=0; i<nsnp; i++) snp_pvalbuf.push_back(make_pair(snp_pval[include_snp[i]], include_snp[i]));
    stable_sort(snp_pvalbuf.begin(), snp_pvalbuf.end());
    
    // Start to clump
    map<string, bool> clumped_snp;
    for(i=0; i<nsnp; i++) {
        pvalbuf = snp_pvalbuf[i].first;
        if( pvalbuf >= pval_thresh1) continue;
        indx = snp_pvalbuf[i].second;
        clumped_snp.insert(pair<string,bool>(snplist[indx],false));
    }

    map<string, int>::iterator iter;
    map<string, bool>::iterator iter_clump;
    int geno_indx=0, geno_indx_j = 0, geno_indx_buf = 0, geno_indx_center = 0, nindi=_keep.size(), m_clumped_snp=clumped_snp.size();
    int left_indx = 0, right_indx = 0;
    double r2_left=0.0, r2_right=0.0;
    vector<string> indices_snp;
    eigenVector x(nindi), x_j(nindi);
    for(i=0; i<m_clumped_snp; i++) {
        indx = snp_pvalbuf[i].second;
        snpbuf = snplist[indx];
       
        if( clumped_snp[snpbuf]) continue;
        iter = _snp_name_map.find(snpbuf);
        if(iter==_snp_name_map.begin() || iter==_snp_name_map.end()) break;
        geno_indx = iter -> second;
        geno_indx_center = std::find(_include.begin(), _include.end(), geno_indx) - _include.begin();
        makex_eigenVector(geno_indx_center, x, false, true);
        
        geno_indx_j = geno_indx_center;
        // Left side
        while(1) {
            r2_left=-1;  geno_indx_j--;

            if(geno_indx_j<0) break;
            snpbuf_left = _snp_name[_include[geno_indx_j]];
            
            if(_chr[geno_indx]==_chr[_include[geno_indx_j]] && abs(_bp[geno_indx] - _bp[_include[geno_indx_j]]) < wind_size) {
                iter_clump = clumped_snp.find(snpbuf_left);
                if(iter_clump==clumped_snp.end()) continue;
                // r2
                makex_eigenVector(geno_indx_j, x_j, false, true);
                r2_left= x.dot(x_j) / sqrt(x.dot(x) * x_j.dot(x_j));
                r2_left = r2_left*r2_left;

                // save the SNP
                if(r2_left >= r2_thresh) iter_clump->second=true;
            } else{
                break;
            }
        }

        // Right side
        geno_indx_j = geno_indx_center;
        while (1) {
            r2_right=-1; geno_indx_j++;
            
            if(geno_indx_j >= nsnp) break;
            snpbuf_right = _snp_name[_include[geno_indx_j]];
           
            if(_chr[geno_indx] == _chr[_include[geno_indx_j]] && abs(_bp[geno_indx] - _bp[_include[geno_indx_j]]) < wind_size ) {
                 iter_clump = clumped_snp.find(snpbuf_right);
                if(iter_clump==clumped_snp.end()) continue;
                // r2
                makex_eigenVector(geno_indx_j, x_j, false, true);
                r2_right= x.dot(x_j) / sqrt(x.dot(x) * x_j.dot(x_j));
                r2_right = r2_right*r2_right;
                // Save the SNP
                if(r2_right >= r2_thresh) iter_clump->second=true;
            } else{
                break;
            }
        }
        indices_snp.push_back(snpbuf);
    }

    return(indices_snp);
}

vector<double> gcta::gsmr_meta(vector<string> snplist, vector<int> include_snp,  map<string,int> snp_sumdata_map, vector<double> bzx, vector<double> bzx_se, vector<double> bzx_pval, vector<double> bzy, vector<double> bzy_se, double pval_thresh1, double pval_thresh2, int wind_size, double r2_thresh, double gwas_thresh, double heidi_thresh, int nsnp_gsmr, int nsnp_heidi, bool flag_heidi) {
    
    int i=0, j=0, nsnp = snplist.size(), nindi=_keep.size();
    vector<string> indices_snp;
    
    if(nsnp < nsnp_gsmr) LOGGER.e(0, "Not enough SNPs to perfrom the GSMR analysis. Only " + to_string(nsnp) + " were retained in the summary data before the clumping analysis.");
    
    // clumping analysis
    indices_snp = clumping_meta(snplist, include_snp, bzx_pval, pval_thresh1, pval_thresh2, wind_size, r2_thresh);
    int n_indices_snp = indices_snp.size();
    
    LOGGER.i(0, to_string(n_indices_snp) + " index SNPs were obtained from the clumping analysis.");
    
    if(n_indices_snp < nsnp_gsmr) LOGGER.e(0, "Not enough SNPs to perfrom the GSMR analysis. At least " + to_string(nsnp_gsmr) + "are required from the GSMR analysis.");
    
    // estimate cov(bxy1, bxy2)
    map<string, int>::iterator iter;
    double x_cov=0.0, x_sd1 = 0.0, x_sd2 = 0.0;
    eigenVector bxy(n_indices_snp);
    for(i=0; i<n_indices_snp; i++) {
        iter = snp_sumdata_map.find(indices_snp[i]);
        bxy(i) = bzy[iter->second]/bzx[iter->second];
     }
    
    MatrixXf x_sub(nindi, n_indices_snp);
    vector<int> snp_sn(n_indices_snp);
    int i_buf = 0;
    for( i = 0; i < n_indices_snp; i++ ) {
        iter = _snp_name_map.find(indices_snp[i]);
        if(iter!=_snp_name_map.end()) {
            i_buf = iter->second;
            snp_sn[i] = find(_include.begin(), _include.end(), i_buf) - _include.begin();
        }
    }
    
    eigenMatrix ld_r_mat(n_indices_snp, n_indices_snp);
    make_XMat_subset(x_sub, snp_sn, true);

    ld_r_mat = MatrixXd::Identity(n_indices_snp, n_indices_snp);
    for(i=0; i<(n_indices_snp-1); i++) {
        x_sd1 = x_sub.col(i).norm();
        for(j=(i+1); j<n_indices_snp; j++) {
            x_cov = x_sub.col(i).dot(x_sub.col(j));
            x_sd2 = x_sub.col(j).norm();
            ld_r_mat(i,j) = ld_r_mat(j,i) = x_cov/(x_sd1*x_sd2);
        }
    }
    
    eigenVector zscore_inv1(n_indices_snp), zscore_inv2(n_indices_snp);
    eigenMatrix cov_bxy(n_indices_snp, n_indices_snp);
    for(i=0; i<n_indices_snp; i++) {
        iter = snp_sumdata_map.find(indices_snp[i]);
        zscore_inv1(i) = bzx_se[iter->second]/bzx[iter->second];
        zscore_inv2(i) = bzy_se[iter->second]/bzx[iter->second];
    }
    for(i=0; i<n_indices_snp; i++) {
        for(j=i; j<n_indices_snp; j++) {
            if(i==j) {
                cov_bxy(i,j) = zscore_inv2(i)*zscore_inv2(j) + zscore_inv1(i)*zscore_inv1(j)*bxy(i)*bxy(j) - bxy(i)*bxy(j)*zscore_inv1(i)*zscore_inv1(j)*zscore_inv1(i)*zscore_inv1(j);
            } else {
                cov_bxy(i,j) = cov_bxy(j,i) = ld_r_mat(i,j)*zscore_inv2(i)*zscore_inv2(j) + ld_r_mat(i,j)*zscore_inv1(i)*zscore_inv1(j)*bxy(i)*bxy(j) - bxy(i)*bxy(j)*zscore_inv1(i)*zscore_inv1(j)*zscore_inv1(i)*zscore_inv1(j);
            }
        }
    }
    
    // heidi-outlier
    string top_indices_snp = "";
    double lower_bounder=0.0, upper_bounder = 0.0, min_bzx_pval=1.0;
    eigenVector bxy_sort(bxy);
    std::sort(bxy_sort.data(), bxy_sort.data()+bxy_sort.size());
    lower_bounder = quantile(bxy_sort, 0.4);
    upper_bounder = quantile(bxy_sort, 0.6);
    vector<int> heidi_ref_snp_id;
    for(i=0; i<n_indices_snp; i++) {
        if(bxy(i) >= lower_bounder && bxy(i) <= upper_bounder) {
            heidi_ref_snp_id.push_back(i);
            iter = snp_sumdata_map.find(indices_snp[i]);
            if(bzx_pval[iter->second] < min_bzx_pval) {
                min_bzx_pval = bzx_pval[iter->second];
                top_indices_snp=indices_snp[i];
            }
        }
    }
    
    // heidi test
    map<string,int> indices_snp_map;
    vector<int> include_gsmr;
    int indxbuf = 0;
    for(i=0; i<n_indices_snp; i++) indices_snp_map.insert(pair<string,int>(indices_snp[i], i));
    iter = indices_snp_map.find(top_indices_snp);
    indxbuf = iter->second;
    double bxy_diff = 0.0, top_bxy = bxy(indxbuf), var_bxy_diff=0.0, pval_heidi=0.0;
    for( i=0; i<n_indices_snp; i++) {
        if(i==indxbuf) {
            include_gsmr.push_back(i);
            continue;
        }
        bxy_diff = bxy(i) - top_bxy;
        var_bxy_diff = cov_bxy(indxbuf, indxbuf) + cov_bxy(i, i) - 2*cov_bxy(indxbuf, i);
        pval_heidi=StatFunc::pchisq(bxy_diff*bxy_diff/var_bxy_diff, 1);
        if(pval_heidi >= heidi_thresh) {
            include_gsmr.push_back(i);
        }
    }
    
    int n_snp_gsmr = include_gsmr.size();
    LOGGER.i(0, to_string(n_snp_gsmr) + " index SNPs were remaining after the HEIDI-outlier test.");
    if(n_snp_gsmr  < nsnp_gsmr) LOGGER.e(0, "Not enough SNPs to perfrom the GSMR analysis. At least " + to_string(nsnp_gsmr) + "are required from the GSMR analysis.");
    
    // update cov(bxy_i, bxy_j) and bxy
    double logdet = 0.0;
    eigenVector bxy_heidi(n_snp_gsmr);
    eigenMatrix cov_bxy_inv(n_snp_gsmr, n_snp_gsmr);
    n_snp_gsmr = include_gsmr.size();
    for(i=0; i<n_snp_gsmr; i++) {
        bxy_heidi(i) = bxy(include_gsmr[i]);
        for(j=i; j<n_snp_gsmr; j++) {
            cov_bxy_inv(i,j) = cov_bxy_inv(j,i) = cov_bxy(include_gsmr[i], include_gsmr[j]);
        }
    }
    if (!comput_inverse_logdet_LU_mkl(cov_bxy_inv, logdet)) LOGGER.e(0, "The variance-covaraince matrix of bxy is not invertible.");

    // gsmr
    if(n_snp_gsmr < nsnp_gsmr) LOGGER.e(0, "Not enough SNPs to perfrom the GSMR analysis. At least " + to_string(nsnp_heidi) + "are required from the GSMR analysis.");
    double bxy_gsmr = 0.0, bxy_gsmr_se = 0.0, bxy_gsmr_pval = 0.0;
    VectorXd vec_1= VectorXd::Ones(n_snp_gsmr);
    eigenMatrix mat_buf;
    mat_buf = vec_1.transpose()*cov_bxy_inv*vec_1;
    bxy_gsmr_se  = 1/mat_buf(0,0);
    mat_buf = vec_1.transpose()*cov_bxy_inv*bxy_heidi;
    bxy_gsmr = bxy_gsmr_se*mat_buf(0,0);
    bxy_gsmr_pval = StatFunc::pchisq( bxy_gsmr*bxy_gsmr/bxy_gsmr_se, 1);
    bxy_gsmr_se = sqrt(bxy_gsmr_se);
    
    vector<double> rst(3);
    rst[0] = bxy_gsmr; rst[1] = bxy_gsmr_se; rst[2] = bxy_gsmr_pval;
    return (rst);
}

double read_ld_marker(string ref_ld_dirt) {
    // Read the number of markers
    int i = 0, chr_num = 22;
    double d_buf = 0.0, ttl_mk_num = 0.0;
    string filestr = "", strbuf = "";
    for(i=0; i<chr_num; i++) {
        filestr = ref_ld_dirt + to_string(i+1)+ ".l2.M_5_50";
        ifstream ref_marker(filestr.c_str());
        if (!ref_marker) LOGGER.e(0, "Can not open the file [" + filestr + "] to read.");
        std::getline(ref_marker, strbuf);
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 1) LOGGER.e(0, "The file format of [" + filestr + "] is not correct.");
        d_buf = atof(line_elements[0].c_str());
        ttl_mk_num += d_buf;
        ref_marker.close();
    }
    return(ttl_mk_num);
}

vector<string> read_ld_score(string ld_dirt, map<string,int> snplist_map, int nsnp, vector<double> &ld_score) {

    // Read the reference / weighted LD score
    int i = 0, chr_num = 22, line_number = 0, n_cm_ld_snp = 0, indxbuf=0;
    double ldscbuf = 0.0;
    string filestr = "", strbuf = "", snpbuf = "";
    vector<string> ld_score_snps;
    vector<bool> flag_ld_score_snp;
    map<string,int>::iterator snp_iter;
    
    ld_score.resize(nsnp); flag_ld_score_snp.resize(nsnp);
    for(i=0; i<nsnp; i++) flag_ld_score_snp[i] = false;
    for(i=0; i<chr_num; i++) {
        line_number = 0;
        filestr = ld_dirt  + to_string(i+1)+ ".l2.ldscore";
        ifstream ldsc_marker(filestr.c_str());
        if (!ldsc_marker) LOGGER.e(0, "Can not open the file [" + filestr + "] to read.");
        while(std::getline(ldsc_marker, strbuf)) {
            line_number ++;
            std::istringstream linebuf(strbuf);
            std::istream_iterator<string> begin(linebuf), end;
            vector<string> line_elements(begin, end);
            if(line_elements.size() != 6) LOGGER.e(0, "The file format of [" + filestr + "] is not correct in line" + to_string(line_number) + ".");
            if(line_number==1) continue;
            snpbuf = line_elements[1];
            ldscbuf  = atof(line_elements[5].c_str());
            // save the data
            snp_iter = snplist_map.find(snpbuf);
            if(snp_iter!=snplist_map.end()) {
                indxbuf = snp_iter->second;
                ld_score[indxbuf] = ldscbuf;
                flag_ld_score_snp[indxbuf]=true;
                ld_score_snps.push_back(snpbuf);
                n_cm_ld_snp++;
            }
        }
        ldsc_marker.close();
    }
    
    // The SNPs in common
    int j=0;
    vector<string> snplist_common(n_cm_ld_snp);
    snp_iter = snplist_map.begin();
    for(i=0, j=0; i<nsnp; i++, snp_iter++) {
        if(flag_ld_score_snp[i]) {
            snplist_common[j] = snp_iter->first;
            j++;
        }
    }

    return(snplist_common);
}

eigenVector update_weights_hsq(double intercept, double h, int ttl_mk_num, int n_ld_snp, eigenVector ref_ld, eigenVector w_ld, eigenVector n) {
    int i=0;
    if(h < 0.0) h = 0.0;
    if(h > 1.0) h = 1.0;
    for(i=0; i<n_ld_snp; i++) {
        if(ref_ld(i) < 1.0) ref_ld(i) = 1.0;
        if(w_ld(i) < 1.0) w_ld(i) = 1.0;
    }
    eigenVector intercept_vec = intercept * eigenVector::Ones(n_ld_snp);
    eigenVector denominator = (intercept_vec + h/(double)ttl_mk_num*n.cwiseProduct(ref_ld));
    denominator = 2*w_ld.cwiseProduct(denominator.cwiseProduct(denominator));
    eigenVector numerator = eigenVector::Ones(n_ld_snp);
    numerator  = numerator.cwiseQuotient(denominator);
    return numerator;
}

eigenVector compute_irls(double &intercept, double &hsq, eigenMatrix x, eigenVector y, eigenVector wt, eigenVector ref_ld, eigenVector w_ld, eigenVector n, int ttl_mk_num, int n_ld_snp, bool intercept_flag, bool x_flag) {
    eigenVector w;
    eigenVector wx(x.col(0));
    // x = x*sqrt(w); y = y*sqrt(w);
    wt = wt.cwiseSqrt();
    wt = wt/wt.sum();
    x.col(0) = x.col(0).cwiseProduct(wt);
    x.col(1) = x.col(1).cwiseProduct(wt);
    y = y.cwiseProduct(wt);
    // b = (x'x)^-1x'y
    if(intercept_flag) {
        eigenMatrix xt_x;
        eigenVector xt_y, b_coeff;
        xt_x = x.transpose()*x;
        xt_y = x.transpose()*y;
        b_coeff = xt_x.ldlt().solve(xt_y);
        // update h2
        hsq=b_coeff(0)*(double)ttl_mk_num/n.mean();
        intercept = b_coeff(1);
    } else {
        double xt_x = 0.0, xt_y = 0.0, b_coeff = 0.0;
        xt_x = x.col(0).transpose()*x.col(0);
        xt_y = x.col(0).transpose()*y;
        b_coeff = xt_y/xt_x;
        // update h2
        hsq = b_coeff*(double)ttl_mk_num/n.mean();
    }
    
    // weights
    if(x_flag) {
        w = update_weights_hsq(intercept, hsq, ttl_mk_num, n_ld_snp, wx, w_ld, n);
    } else {
        w = update_weights_hsq(intercept, hsq, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n);
    }
    
    return w;
}

vector<double> est_hsq_trait_1_step(eigenVector chival, eigenVector n, eigenVector ref_ld, eigenVector w_ld, int n_ld_snp, int ttl_mk_num) {

    // Estimate prior of the weights
    eigenVector denominator = ref_ld.cwiseProduct(n);
    double h_prior = (chival.mean() - 1)*(double)ttl_mk_num/denominator.mean();
    double intercept_prior = 1.0;
    eigenVector wt_ttl = update_weights_hsq(intercept_prior, h_prior, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n);
    eigenMatrix x(n_ld_snp, 2);
    x.col(1).setOnes(); x.col(0) = ref_ld.cwiseProduct(n)/n.mean();
    
    // Iteratively reweighted least squares
    wt_ttl = compute_irls(intercept_prior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, true, false);
    wt_ttl = compute_irls(intercept_prior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, true, false);
    // Estimate intercept and slope
    double intercept_posterior = 0.0, h_posterior = 0.0;
    wt_ttl = compute_irls(intercept_posterior, h_posterior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, true, false);
    
    vector<double> ldsc_est(2);
    ldsc_est[0]  = intercept_posterior; ldsc_est[1] = h_posterior;
    return ldsc_est;
}

vector<double> est_hsq_trait_2_steps(eigenVector chival, eigenVector n, eigenVector ref_ld, eigenVector w_ld, int n_ld_snp, int ttl_mk_num) {
    int i = 0, n_subset_snp = 0;
    double thresh = 30;
    vector<int> subset_indx;
    eigenVector subset_chi, subset_n, subset_ref_ld, subset_w_ld;

    // Estimate prior of the weights
    eigenVector denominator = ref_ld.cwiseProduct(n);
    double h_prior = (chival.mean() - 1)*(double)ttl_mk_num/denominator.mean();
    
    double intercept_prior = 1.0;
    eigenVector wt_ttl = update_weights_hsq(intercept_prior, h_prior, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n);
    eigenVector subset_wt;
    eigenMatrix x(n_ld_snp,2), subset_x;
    x.col(1).setOnes(); x.col(0) = ref_ld.cwiseProduct(n)/n.mean();
    
    // SNPs with chi^2 < 30
    for(i=0; i<n_ld_snp; i++) {
        if(chival(i) < thresh) subset_indx.push_back(i);
    }
    n_subset_snp = subset_indx.size();
    
    subset_chi.resize(n_subset_snp); subset_n.resize(n_subset_snp);
    subset_ref_ld.resize(n_subset_snp);  subset_w_ld.resize(n_subset_snp);
    subset_wt.resize(n_subset_snp); subset_x.resize(n_subset_snp,2);
    for(i=0; i<n_subset_snp; i++) {
        subset_chi(i) = chival(subset_indx[i]);
        subset_n(i) = n(subset_indx[i]);
        subset_ref_ld(i) = ref_ld(subset_indx[i]);
        subset_w_ld(i) = w_ld(subset_indx[i]);
        subset_wt(i) = wt_ttl(subset_indx[i]);
        subset_x.row(i) = x.row(subset_indx[i]);
    }
    
    // Iteratively reweighted least squares
    subset_wt = compute_irls(intercept_prior, h_prior, subset_x, subset_chi, subset_wt, subset_ref_ld, subset_w_ld, subset_n, ttl_mk_num, n_subset_snp, true, true);
    subset_wt = compute_irls(intercept_prior, h_prior, subset_x, subset_chi, subset_wt, subset_ref_ld, subset_w_ld, subset_n, ttl_mk_num, n_subset_snp, true, true);
    
    // Estimate intercept
    double intercept_posterior = 0.0;
    subset_wt = compute_irls(intercept_posterior, h_prior, subset_x, subset_chi, subset_wt, subset_ref_ld, subset_w_ld, subset_n, ttl_mk_num, n_subset_snp, true, true);
    
    // Iteratively reweighted least squares
    chival = chival - eigenVector::Ones(n_ld_snp)*intercept_posterior;
    wt_ttl = compute_irls(intercept_posterior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, false, false);
    wt_ttl = compute_irls(intercept_posterior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, false, false);
    
    // Estimate slope
    double h_posterior = 0.0;
    wt_ttl = compute_irls(intercept_posterior, h_posterior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, false, false);
    
    vector<double> ldsc_est(2);
    ldsc_est[0]  = intercept_posterior; ldsc_est[1] = h_posterior;
    return ldsc_est;
}

eigenVector update_weights_gcov(double intercept1, double h1, double intercept2, double h2, double intercept_gcov, double gcov, int ttl_mk_num, int n_ld_snp, eigenVector ref_ld, eigenVector w_ld, eigenVector n1, eigenVector n2, eigenVector n_gcov ) {
    int i=0;
    if(h1 < 0.0) h1 = 0.0; if(h1 > 1.0) h1 = 1.0;
    if(h2 < 0.0) h2 = 0.0; if(h2 > 1.0) h2 = 1.0;
    if(gcov < -1.0) gcov = -1.0; if(gcov > 1.0) gcov = 1.0;
    for(i=0; i<n_ld_snp; i++) {
        if(ref_ld(i) < 1.0) ref_ld(i) = 1.0;
        if(w_ld(i) < 1.0) w_ld(i) = 1.0;
    }
    
    double d1=0.0, d2=0.0, d3=0.0;
    eigenVector w(n_ld_snp);

    for(i=0; i<n_ld_snp; i++) {
        d1 = n1(i)*h1*ref_ld(i)/ttl_mk_num + intercept1;
        d2 = n2(i)*h2*ref_ld(i)/ttl_mk_num + intercept2;
        d3 = n_gcov(i)*gcov*ref_ld(i)/ttl_mk_num + intercept_gcov;
        w(i) = 1/( w_ld(i) * (d1*d2 + d3*d3) );
    }
    return w;
}

eigenVector compute_irls_gcov(double &intercept_gcov, double &gcov, eigenMatrix x, eigenVector y, eigenVector wt, eigenVector ref_ld, eigenVector w_ld, eigenVector n_gcov, double intercept1,  double hsq1, eigenVector n1, double intercept2,  double hsq2, eigenVector n2, int ttl_mk_num, int n_ld_snp, bool intercept_flag, bool x_flag) {

    eigenVector w;
    // x = x*sqrt(w); y = y*sqrt(w);
    wt = wt.cwiseSqrt();
    wt = wt/wt.sum();
    x.col(0) = x.col(0).cwiseProduct(wt);
    x.col(1) = x.col(1).cwiseProduct(wt);
    y = y.cwiseProduct(wt);
    // b = (x'x)^-1x'y
    if(intercept_flag) {
        eigenMatrix xt_x;
        eigenVector xt_y, b_coeff;
        xt_x = x.transpose()*x;
        xt_y = x.transpose()*y;
        b_coeff = xt_x.ldlt().solve(xt_y);
        // update h2
        gcov=b_coeff(0)*(double)ttl_mk_num/n_gcov.mean();
        intercept_gcov = b_coeff(1);
    } else {
        double xt_x = 0.0, xt_y = 0.0, b_coeff = 0.0;
        xt_x = x.col(0).transpose()*x.col(0);
        xt_y = x.col(0).transpose()*y;
        b_coeff = xt_y/xt_x;
        // update h2
        gcov = b_coeff*(double)ttl_mk_num/n_gcov.mean();
    }
    
    // weights
    if(x_flag) {
        w = update_weights_gcov(intercept1, hsq1, intercept2, hsq2, intercept_gcov, gcov, ttl_mk_num, n_ld_snp, x.col(0), w_ld, n1,  n2, n_gcov);
    } else {
        w = update_weights_gcov(intercept1, hsq1, intercept2, hsq2, intercept_gcov, gcov, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n1,  n2, n_gcov);
    }
    return w;
}

vector<double> est_gcov_trait_1_step(eigenVector zscore, eigenVector n_gcov, eigenVector ref_ld, eigenVector w_ld, double intercept1, double hsq1, eigenVector n1, double intercept2, double hsq2, eigenVector n2, int n_ld_snp, int ttl_mk_num) {
    
    // Estimate prior of the weights
    eigenVector denominator = ref_ld.cwiseProduct(n_gcov);
    double gcov_prior = zscore.mean()*(double)ttl_mk_num/denominator.mean();
    double intercept_gcov_prior = 0.0;
    eigenVector wt_ttl = update_weights_gcov(intercept1, hsq1, intercept2, hsq2, intercept_gcov_prior, gcov_prior, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n1, n2, n_gcov);
    eigenMatrix x(n_ld_snp, 2);
    x.col(1).setOnes(); x.col(0) = ref_ld.cwiseProduct(n_gcov)/n_gcov.mean();
    
    // Estimate intercept and slope by Iteratively reweighted least squares
    double intercept_gcov_posterior = 0.0, gcov_posterior = 0.0;

    wt_ttl = compute_irls_gcov(intercept_gcov_posterior, gcov_posterior, x, zscore, wt_ttl, ref_ld, w_ld, n_gcov, intercept1,  hsq1, n1, intercept2,  hsq2, n2, ttl_mk_num, n_ld_snp, true, false);
    wt_ttl = compute_irls_gcov(intercept_gcov_posterior, gcov_posterior, x, zscore, wt_ttl, ref_ld, w_ld, n_gcov, intercept1,  hsq1, n1, intercept2,  hsq2, n2, ttl_mk_num, n_ld_snp, true, false);
    wt_ttl = compute_irls_gcov(intercept_gcov_posterior, gcov_posterior, x, zscore, wt_ttl, ref_ld, w_ld, n_gcov, intercept1,  hsq1, n1, intercept2,  hsq2, n2, ttl_mk_num, n_ld_snp, true, false);
    
    vector<double> ldsc_est(2);
    ldsc_est[0]  = intercept_gcov_posterior; ldsc_est[1] = gcov_posterior/sqrt(hsq1*hsq2);
    
    return ldsc_est;
}

double transform_hsq_L(double P, double K, double hsq) {
    double t = StatFunc::qnorm(1.0 - K);
    double z = StatFunc::dnorm(t);
    double C = (K * (1 - K) / (z * z))*(K * (1 - K) / (P * (1 - P)));
    return (hsq * C);
}

vector<vector<vector<double>>> mtcojo_ldsc(vector<vector<double>> snp_b, vector<vector<double>> snp_se, vector<vector<double>> snp_n, vector<int> rg_pairs_indx1, vector<int> rg_pairs_indx2, vector<string> snplist, vector<int> include_sumdata, vector<double> popu_prev, vector<double> smpl_prev, string ref_ld_dirt, string w_ld_dirt) {
    int i=0, ntrait = snp_b.size(), nsnp = include_sumdata.size();
    double ttl_mk_num = 0.0;
    vector<double> ref_ld_vec, w_ld_vec;
    vector<string> ref_ld_snps, w_ld_snps;
    map<string,int> snplist_map;
    for(i=0; i<nsnp; i++) {
        snplist_map.insert(pair<string,int>(snplist[include_sumdata[i]], i));
    }
    
    // Read the total number of markers
    ttl_mk_num = read_ld_marker(ref_ld_dirt);
    
    // Read the reference LD scores
    ref_ld_snps = read_ld_score(ref_ld_dirt, snplist_map, nsnp, ref_ld_vec);
    // Read the weighted LD scores
    w_ld_snps = read_ld_score(w_ld_dirt, snplist_map, nsnp, w_ld_vec);
    
    // SNPs in common
    map<string,int> w_ld_snp_map;
    vector<string> cm_ld_snps;
    int n_ref_ld_snps = ref_ld_snps.size(), n_w_ld_snps = w_ld_snps.size();
    for(i=0; i<n_w_ld_snps; i++)
        w_ld_snp_map.insert(pair<string, int>(w_ld_snps[i], i));
    for(i=0; i<n_ref_ld_snps; i++) {
        if(w_ld_snp_map.find(ref_ld_snps[i])!=w_ld_snp_map.end()) {
            cm_ld_snps.push_back(ref_ld_snps[i]);
        }
    }
    
    // Re-order the variables
    int j = 0, n_cm_ld_snps = cm_ld_snps.size(), indxbuf = 0;
    map<string,int>::iterator iter;
    eigenVector ref_ld(n_cm_ld_snps), w_ld(n_cm_ld_snps);
    eigenMatrix bhat_z(n_cm_ld_snps, ntrait), bhat_n(n_cm_ld_snps, ntrait);
    for(i=0; i<n_cm_ld_snps; i++) {
        iter = snplist_map.find(cm_ld_snps[i]);
        if(iter != snplist_map.end()) {
            indxbuf = iter->second;
            ref_ld(i) = ref_ld_vec[indxbuf]; w_ld(i) = w_ld_vec[indxbuf];
            for(j=0; j<ntrait; j++) {
                bhat_z(i,j) = snp_b[j][include_sumdata[indxbuf]]/snp_se[j][include_sumdata[indxbuf]];
                bhat_n(i,j) = snp_n[j][include_sumdata[indxbuf]];
            }
        }
    }

    // Estimate h2 of those traits
    vector<double> rst_ldsc(2);
    vector<vector<vector<double>>>  ldsc_trait(3);
    
    ldsc_trait[0].resize(ntrait); ldsc_trait[1].resize(ntrait); ldsc_trait[2].resize(ntrait);
    for(i=0; i<ntrait; i++) {
        ldsc_trait[0][i].resize(ntrait); ldsc_trait[1][i].resize(ntrait); ldsc_trait[2][i].resize(ntrait);
    }
    for(i=0; i<ntrait; i++) {
        rst_ldsc = est_hsq_trait_2_steps(bhat_z.col(i).cwiseProduct(bhat_z.col(i)), bhat_n.col(i), ref_ld, w_ld, n_cm_ld_snps, ttl_mk_num);
        ldsc_trait[0][i][i] = rst_ldsc[0];

        if(rst_ldsc[1]>0) {
            if(popu_prev[i] > 0 ) {
                ldsc_trait[2][i][i] = rst_ldsc[1];
                ldsc_trait[1][i][i] = transform_hsq_L(smpl_prev[i], popu_prev[i], rst_ldsc[1]);
            } else {
                ldsc_trait[1][i][i] = rst_ldsc[1];
                ldsc_trait[2][i][i] = rst_ldsc[1];
            }
        } else {
            if(i==0) {
                LOGGER.e(0, "The analysis stopped, because of negative heritability estimate of the target trait. Exiting ...");
            } else {
                LOGGER.e(0, "The analysis stopped, because of negative heritability estimate of covariate" + to_string(i) + ". Exiting ...");
            }
        }
    }

    // Estimate cov(g1,g2) of those traits
    eigenVector gcov_n1n2(n_cm_ld_snps);
    for(i=0; i<ntrait; i++) {
        for(j=i+1; j<ntrait; j++) {
            gcov_n1n2 = bhat_n.col(i).cwiseProduct(bhat_n.col(j));
            rst_ldsc = est_gcov_trait_1_step(bhat_z.col(i).cwiseProduct(bhat_z.col(j)), gcov_n1n2.cwiseSqrt(), ref_ld, w_ld, ldsc_trait[0][i][i], ldsc_trait[2][i][i], bhat_n.col(i), ldsc_trait[0][j][j], ldsc_trait[2][j][j], bhat_n.col(j), n_cm_ld_snps, ttl_mk_num);
            ldsc_trait[0][i][j] = ldsc_trait[0][j][i] = rst_ldsc[0];
            ldsc_trait[1][i][j] = ldsc_trait[1][j][i] = rst_ldsc[1];
        }
    }

    return ldsc_trait;
}

eigenMatrix mtcojo_cond_single_covar(eigenVector bzy, eigenVector bzy_se,  eigenMatrix bzx, eigenMatrix bzx_se, double bxy, vector<vector<double>> ldsc_intercept, vector<vector<double>> ldsc_slope, int nsnp) {
    int i=0;
    eigenMatrix mtcojo_est(nsnp, 3);
    
    double var_bzx_buf = 0.0, cov_bzx_bzy = 0.0;
    for(i=0; i<nsnp; i++) {
        mtcojo_est(i,0) = bzy(i) - bzx(i,0)*bxy;
        var_bzx_buf = bxy*bxy*bzx_se(i,0)*bzx_se(i,0);
        cov_bzx_bzy = bxy*ldsc_intercept[0][1]*bzx_se(i,0)*bzy_se(i);
        mtcojo_est(i,1) = sqrt(bzy_se(i)*bzy_se(i) + var_bzx_buf - 2*cov_bzx_bzy);
        mtcojo_est(i,2) = StatFunc::pchisq(mtcojo_est(i,0)*mtcojo_est(i,0)/mtcojo_est(i,1)/mtcojo_est(i,1), 1);
     }
    
    return mtcojo_est;
}

eigenMatrix mtcojo_cond_multiple_covars(eigenVector bzy, eigenVector bzy_se,  eigenMatrix bzx, eigenMatrix bzx_se, vector<double> bxy, vector<vector<double>> ldsc_intercept, vector<vector<double>> ldsc_slope, int nsnp, int ncovar) {
    int i=0, j=0;
    eigenMatrix mtcojo_est(nsnp, 3);
    
    double var_bzx_buf = 0.0, cov_bzx_bzy = 0.0;
    eigenVector bjxy(ncovar), bxy_hat(ncovar);
    MatrixXd d_mat = MatrixXd::Zero(ncovar,ncovar), r_mat = MatrixXd::Identity(ncovar, ncovar);
    
    for(i=0; i<ncovar; i++) {
        d_mat(i,i) = sqrt(ldsc_slope[i+1][i+1]);
        bxy_hat(i) = bxy[i];
        for(j=1; j<ncovar; j++) {
            if(i==j) continue;
            r_mat(i,j) = r_mat(j,i) = ldsc_slope[i+1][j+1];
        }
    }
    
    bjxy = d_mat.ldlt().solve(r_mat.ldlt().solve(d_mat*bxy_hat));
    
    double d_buf = 0.0;
    eigenVector bjxy_buf(ncovar), se_bzx_bzy(ncovar);
    eigenMatrix bzx_se_buf=MatrixXd::Zero(ncovar,ncovar), bzx_intercept=MatrixXd::Identity(ncovar,ncovar);
    for(i=0; i<(ncovar-1); i++) {
        for(j=i+1; j<ncovar; j++) {
            bzx_intercept(i,j) = bzx_intercept(j,i) = ldsc_intercept[i+1][j+1];
        }
    }
    
    for(i=0; i<nsnp; i++) {
        d_buf = 0.0; cov_bzx_bzy = 0.0;
        for(j=0; j<ncovar; j++) {
            d_buf += bzx(i,j) * bjxy(j);
            bzx_se_buf(j,j) = bzx_se(i,j);
            cov_bzx_bzy += bzx_se(i,j)*bzy_se(i)*bjxy(j)*ldsc_intercept[0][j+1];
        }
        var_bzx_buf = (bjxy.transpose()*bzx_se_buf*bzx_intercept*bzx_se_buf*bjxy)(0,0);
        mtcojo_est(i,0) = bzy(i) - d_buf;
        mtcojo_est(i,1) = sqrt(bzy_se(i)*bzy_se(i) + var_bzx_buf - 2*cov_bzx_bzy);
        mtcojo_est(i,2) = StatFunc::pchisq(mtcojo_est(i,0)*mtcojo_est(i,0)/mtcojo_est(i,1)/mtcojo_est(i,1), 1);
    }
    
    return mtcojo_est;
}

void mtcojo_cond_output(string output_file, vector<string> snplist, vector<int> include_snp_meta, map<string, int> snp_name_map, vector<int> snp_chr_g, vector<int> snp_bp_g, vector<string> snp_a1_g, vector<string> snp_a2_g, vector<double> snp_freq, vector<double> snp_freq_g, vector<double> snp_b, vector<double> snp_se, vector<double> snp_pval, eigenMatrix mtcojo_est, vector<double> snp_n, int nsnp) {
    
    int i=0, snpindx=0;
    ofstream ofile(output_file.c_str());
    if (!ofile) LOGGER.e(0, "Can not open the file [" + output_file + "] to write.");
    map<string,int>::iterator iter;
    
    ofile << "Chr\tSNP\tbp\tA1\tA2\tfreq\tb\tse\tp\tN\tfreq_geno\tbC\tbC_se\tbC_pval" <<endl;
    for (i = 0; i < nsnp; i++) {
        iter=snp_name_map.find(snplist[include_snp_meta[i]]);
        snpindx = iter -> second;
        ofile << snp_chr_g[snpindx] <<"\t"<<snplist[include_snp_meta[i]] << "\t" << snp_bp_g[snpindx] << "\t" <<snp_a1_g[snpindx] << "\t" << snp_a2_g[snpindx] << "\t" << o_digital_number(snp_freq[include_snp_meta[i]]) << "\t" << o_digital_number(snp_b[include_snp_meta[i]]) << "\t" << o_digital_number(snp_se[include_snp_meta[i]]) << "\t"  << snp_pval[include_snp_meta[i]] << "\t" << o_digital_number(snp_n[include_snp_meta[i]]) << "\t" << 0.5*snp_freq_g[snpindx] << "\t" << mtcojo_est(i,0) << "\t" << mtcojo_est(i,1) <<"\t"<< mtcojo_est(i,2) << "\t" <<endl;
    }
    ofile.close();
}

void gcta::mtcojo(string mtcojolist_file, string ref_ld_dirt, string w_ld_dirt, double clump_thresh1, double clump_thresh2, int clump_wind_size, double clump_r2_thresh, double gwas_thresh, double heidi_thresh, int nsnp_heidi, int nsnp_gsmr, bool heidi_flag) {
    
    string target_pheno="";
    vector<string> covar_pheno, snplist, badsnps;
    vector<double> popu_prev, smpl_prev;
    vector<vector<string>> snp_a1, snp_a2;
    vector<vector<double>> snp_freq, snp_b, snp_se, snp_pval, snp_n;
    int i=0, j=0, nsnp=0, ncovar=0;

    // Read the summary data
    LOGGER.i(0,"\n");
    LOGGER.i(0, "Reading GWAS summary-level statistics from [" + mtcojolist_file + "] ...");
    read_mtcojofile(mtcojolist_file, target_pheno, covar_pheno, ncovar, snplist, snp_a1, snp_a2, snp_freq,  snp_b, snp_se, snp_pval, snp_n, popu_prev, smpl_prev);
    if(_mu.empty()) calcu_mu();
    
    nsnp = snplist.size();
    map<string, int> snp_sumdata_map;
    vector<int> include_sumdata(nsnp);
    for(i=0; i<nsnp; i++) {
        snp_sumdata_map.insert(pair<string,int>(snplist[i], i));
        include_sumdata[i] = i;
    }
    
    // QC of SNPs
    LOGGER.i(0, to_string(nsnp) + " SNPs were matched with the genotype data ...");
    badsnps = remove_bad_snps(snplist, snp_a1, snp_a2, snp_freq,  snp_b, snp_se, snp_pval, snp_n, _allele1, _allele2, _snp_name_map, target_pheno, covar_pheno, ncovar, _out);
    
    if(badsnps.size()>0) {
        update_id_map_rm(badsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(badsnps, snp_sumdata_map, include_sumdata);
    }
    nsnp = include_sumdata.size();
    if(nsnp<1) LOGGER.e(0, "None SNPs were retained after quality controls ...");
    else LOGGER.i(0, to_string(nsnp) + " SNPs were retained after quality controls ...\n");
    
    // GSMR analysis
    vector<double> gsmr_rst(3), bxy_est(ncovar);
    for(i=1; i<=ncovar; i++) {
        LOGGER.i(0, "GSMR analysis for covariate #" + to_string(i) + " ...");
        gsmr_rst =  gsmr_meta(snplist, include_sumdata, snp_sumdata_map, snp_b[i], snp_se[i], snp_pval[i], snp_b[0], snp_se[0],  clump_thresh1, clump_thresh2, clump_wind_size, clump_r2_thresh, gwas_thresh, heidi_thresh, nsnp_heidi, nsnp_gsmr, heidi_flag);
        bxy_est[i-1] = gsmr_rst[0];
        LOGGER.i(0, "GSMR analysis for covariate #" + to_string(i) + " completed.");
    }

    // LDSC analysis to estimate h2 and rg
    vector<int> rg_pairs_indx1(ncovar*(ncovar+1)), rg_pairs_indx2(ncovar*(ncovar+1));
    vector<vector<vector<double>>> ldsc_est;
    int k = 0;
    for(i=0; i<ncovar; i++) {
        for(j=0; j<=ncovar; j++) {
            rg_pairs_indx1[k] = i+1;
            rg_pairs_indx2[k] = j;
            k++;
        }
    }
    LOGGER.i(0, "\n LD score regression analysis to estimate h2 and rg ...");
    ldsc_est = mtcojo_ldsc(snp_b, snp_se, snp_n, rg_pairs_indx1, rg_pairs_indx2, snplist, include_sumdata, popu_prev, smpl_prev, ref_ld_dirt, w_ld_dirt);
    LOGGER.i(0, "LD score regression analysis completed.\n");
    
    // mtcojo analysis
    eigenVector snp_bzy(nsnp), snp_bzy_se(nsnp);
    eigenMatrix mtcojo_est, snp_bzx(nsnp, ncovar), snp_bzx_se(nsnp, ncovar);
    for(i=0; i<nsnp; i++) {
        snp_bzy(i) = snp_b[0][include_sumdata[i]]; snp_bzy_se(i) = snp_se[0][include_sumdata[i]];
        for(j = 0; j<ncovar; j++) {
            snp_bzx(i,j) = snp_b[j+1][include_sumdata[i]]; snp_bzx_se(i,j) = snp_se[j+1][include_sumdata[i]];
        }
    }
    
    LOGGER.i(0, "mtCOJO analysis to adjust the target trait ...");
    if(ncovar==1) {
        mtcojo_est = mtcojo_cond_single_covar(snp_bzy, snp_bzy_se, snp_bzx, snp_bzx_se, bxy_est[0], ldsc_est[0], ldsc_est[1], nsnp);
    } else {
        mtcojo_est = mtcojo_cond_multiple_covars(snp_bzy, snp_bzy_se, snp_bzx, snp_bzx_se, bxy_est, ldsc_est[0], ldsc_est[1], nsnp, ncovar);
    }
   LOGGER.i(0, "mtCOJO analysis completed ...");
    
    // Output
    string output_filename = _out + ".mtcojo.cma";
    LOGGER.i(0, "Saving the mtCOJO analysis results of " + to_string(nsnp) + " remaining SNPs to [" + output_filename + "] ...");
    mtcojo_cond_output(output_filename, snplist, include_sumdata, _snp_name_map,_chr,  _bp, _allele1, _allele2, snp_freq[0], _mu, snp_b[0], snp_se[0], snp_pval[0], mtcojo_est, snp_n[0], nsnp);
}

