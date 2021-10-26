#include "gcta.h"
#include "Logger.h"
#include <sstream>
#include <iterator>

bool determine_gwas_file(string input_file) {
    int nelements = 0;
    bool file_type = false;
    string strbuf;

    ifstream file_list(input_file);
    if(!file_list) 
        LOGGER.e(0, "cannot open file [" + input_file + "] to read.");

    std::getline(file_list, strbuf);
    std::istringstream linebuf(strbuf);
    vector<string> line_elements((istream_iterator<string>(linebuf)), istream_iterator<string>());
    nelements = line_elements.size();

    if(nelements==1) {
        file_type = true;
    } else if(nelements==8) {
        file_type = false;
    } else {
        LOGGER.e(0, "the GWAS summary data should be in GCTA-COJO format. Please check.");
    }

    return file_type;
}

void update_id_map_kp_by_map(map<string, int> &id_list_map, map<string, int> &id_map, vector<int> &keep) {
    int i = 0;
    map<string, int> id_map_buf(id_map);
    map<string, int>::iterator iter, iter_buf;
    for (iter = id_list_map.begin(); iter != id_list_map.end(); iter++) {
        iter_buf = id_map_buf.find(iter->first);
        if(iter_buf == id_map_buf.end()) continue;
        id_map_buf.erase(iter_buf);
    }
    for (iter = id_map_buf.begin(); iter != id_map_buf.end(); iter++) id_map.erase(iter->first);

    keep.clear();
    for (iter = id_map.begin(); iter != id_map.end(); iter++) keep.push_back(iter->second);
    stable_sort(keep.begin(), keep.end());
}

void read_gsmr_file_list(string gsmr_file_list, vector<string> &pheno_name, vector<string> &pheno_file, vector<double> &popu_prev, vector<double> &smpl_prev) {

    ifstream meta_list(gsmr_file_list.c_str());
    if (!meta_list)
        LOGGER.e(0, "cannot open the file [" + gsmr_file_list + "] to read.");
    
    string strbuf="", prevbuf1="", prevbuf2="";
    // Retrieve the GWAS summary data file
    // The 1st row: the target trait
    int line_number = 1;
    while(std::getline(meta_list, strbuf)) 
    {
        std::istringstream linebuf(strbuf);
        vector<string> line_elements((istream_iterator<string>(linebuf)), istream_iterator<string>());
        
        if(line_elements.size() != 2 && line_elements.size() != 4)
            LOGGER.e(0, "the format of file [" + gsmr_file_list + "] is incorrect, line " + to_string(line_number) + ".");
        
        pheno_name.push_back(line_elements[0]); pheno_file.push_back(line_elements[1]);
        
        // prevelance
        double d_prev1 = nan(""), d_prev2 = nan("");
        if(line_elements.size()==4) {
            prevbuf1 = line_elements[2]; prevbuf2 = line_elements[3];
            StrFunc::to_upper(prevbuf1); StrFunc::to_upper(prevbuf2);
            // available, 1 - sample prevelance, 2 - population prevelance
            if(prevbuf1 != "NA"  &&  prevbuf1!= "NAN" && prevbuf1!= ".") {
                d_prev1 = atof(prevbuf1.c_str());
                if(d_prev1 < 0 || d_prev1 > 1)
                    LOGGER.e(0, "invalid sample prevalence for [" + pheno_name[line_number] + "].");
            }
            if(prevbuf2 != "NA"  &&  prevbuf2!= "NAN" && prevbuf2 != ".") {
                d_prev2 = atof(prevbuf2.c_str());
                if(d_prev2 < 0 || d_prev2 > 1)
                    LOGGER.e(0, "invalid population prevalence for [" + pheno_name[line_number] + "].");
            }
        }
        smpl_prev.push_back(d_prev1); popu_prev.push_back(d_prev2);
        line_number++;
    }
    meta_list.close();
}

void gcta::read_gsmrfile(string expo_file_list, string outcome_file_list, double gwas_thresh, int nsnp_gsmr, int gsmr_so_alg) {
    int i = 0, j = 0;
    double pval_thresh = gsmr_so_alg >= 0 ?  1.0 : gwas_thresh;
    vector<string> expo_gwas_file, outcome_gwas_file, gwas_data_file, pheno_name_buf;
    vector<double> popu_prev_buf, smpl_prev_buf;
    
    // Read the SNPs
    // Exposure
    LOGGER.i(0, "\nReading GWAS summary data for exposure(s) from [" + expo_file_list + "].");
    popu_prev_buf.clear(); smpl_prev_buf.clear(); pheno_name_buf.clear();
    read_gsmr_file_list(expo_file_list, pheno_name_buf, expo_gwas_file, popu_prev_buf, smpl_prev_buf);
    _meta_popu_prev = popu_prev_buf; _meta_smpl_prev = smpl_prev_buf; 
    _gwas_trait_name = pheno_name_buf;
    gwas_data_file = expo_gwas_file; 

    _expo_num = pheno_name_buf.size();

    vector<string> snplist;
    map<string, int> gws_snp_name_map;
    for(i=0; i<_expo_num; i++) {
        if(expo_gwas_file[i].substr(expo_gwas_file[i].length()-3,3)!=".gz")
            snplist=read_snp_metafile_txt(expo_gwas_file[i], gws_snp_name_map, pval_thresh);
        else
            snplist=read_snp_metafile_gz(expo_gwas_file[i], gws_snp_name_map, pval_thresh);
        if(i==0) init_meta_snp_map(snplist, _meta_snp_name_map, _meta_snp_name, _meta_remain_snp);
        else update_meta_snp_map(snplist, _meta_snp_name_map, _meta_snp_name, _meta_remain_snp, true);
    }

    // Outcome
    LOGGER.i(0, "Reading GWAS summary data for outcome(s) from [" + outcome_file_list + "].");
    popu_prev_buf.clear(); smpl_prev_buf.clear(); pheno_name_buf.clear();
    read_gsmr_file_list(outcome_file_list, pheno_name_buf, outcome_gwas_file, popu_prev_buf, smpl_prev_buf);
    _meta_popu_prev.insert(_meta_popu_prev.end(), popu_prev_buf.begin(), popu_prev_buf.end());
    _meta_smpl_prev.insert(_meta_smpl_prev.end(), smpl_prev_buf.begin(), smpl_prev_buf.end());
    _gwas_trait_name.insert(_gwas_trait_name.end(), pheno_name_buf.begin(), pheno_name_buf.end());
    gwas_data_file.insert(gwas_data_file.end(), outcome_gwas_file.begin(), outcome_gwas_file.end());
    
    _outcome_num = pheno_name_buf.size();

    map<string, int> outcome_snp_name_map;
    vector<string> outcome_snp_name;
    vector<int> outcome_remain_snp;  
    for(i=0; i<_outcome_num; i++) {
        if(outcome_gwas_file[i].substr(outcome_gwas_file[i].length()-3,3)!=".gz")
            snplist=read_snp_metafile_txt(outcome_gwas_file[i], gws_snp_name_map, pval_thresh);
        else
            snplist=read_snp_metafile_gz(outcome_gwas_file[i], gws_snp_name_map, pval_thresh);
        update_meta_snp_list(snplist, _meta_snp_name_map);
        if(i==0) init_meta_snp_map(snplist, outcome_snp_name_map, outcome_snp_name, outcome_remain_snp);
        else update_meta_snp_map(snplist, outcome_snp_name_map, outcome_snp_name, outcome_remain_snp, false);
    }  
    // SNPs in common between exposures and outcomes
    update_id_map_kp(outcome_snp_name, _meta_snp_name_map, _meta_remain_snp); 

    // Keep significant SNPs
    update_id_map_kp_by_map(gws_snp_name_map, _meta_snp_name_map, _meta_remain_snp);
    // Initialization of variables
    int nsnp = _meta_snp_name_map.size(), npheno = _expo_num + _outcome_num;
    vector<vector<string>> snp_a1, snp_a2;
    eigenMatrix snp_freq; 

    init_gwas_variable(snp_a1, snp_a2, snp_freq, _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, npheno, nsnp); 

    // reset SNP variables
    update_meta_snp(_meta_snp_name_map, _meta_snp_name, _meta_remain_snp);
    
    LOGGER.i(0, to_string(nsnp) + " genome-wide significant SNPs in common between the exposure(s) and the outcome(s).");
    if(nsnp<nsnp_gsmr) LOGGER.e(0, "not enough SNPs to perform the GSMR analysis.");

    // Reading the summary data
    _meta_vp_trait.resize(npheno);
    _snp_val_flag.clear(); _snp_val_flag.resize(npheno);
    for(i=0; i<npheno; i++) {
        _snp_val_flag[i].resize(nsnp);
        for(j=0; j<nsnp; j++) _snp_val_flag[i][j] = false;
    }

    // Summary data
    for(i=0; i<npheno; i++) {
        eigenVector snp_freq_buf(nsnp), snp_b_buf(nsnp), snp_se_buf(nsnp), snp_pval_buf(nsnp), snp_n_buf(nsnp);
        if(gwas_data_file[i].substr(gwas_data_file[i].length()-3,3)!=".gz")
            _meta_vp_trait[i] = read_single_metafile_txt(gwas_data_file[i], _meta_snp_name_map,  snp_a1[i], snp_a2[i], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i]);
        else _meta_vp_trait[i] = read_single_metafile_gz(gwas_data_file[i], _meta_snp_name_map,  snp_a1[i], snp_a2[i], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i]);
        if(_meta_vp_trait[i] < 0) LOGGER.e(0, "negative phenotypic variance of trait " + _gwas_trait_name[i] + ".");
        snp_freq.col(i) = snp_freq_buf;
        _meta_snp_b.col(i) = snp_b_buf;
        _meta_snp_se.col(i) = snp_se_buf;
        _meta_snp_pval.col(i) = snp_pval_buf;
        _meta_snp_n_o.col(i) = snp_n_buf;
    }

    // QC of SNPs
    LOGGER.i(0, "Filtering out SNPs with multiple alleles or missing value ...");
    vector<string> badsnps;
    badsnps = remove_bad_snps(_meta_snp_name, _meta_remain_snp, _snp_val_flag, snp_a1, snp_a2, snp_freq,  _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, 
                              _snp_name_map, _allele1, _allele2, _outcome_num, _expo_num, _out);
    if(badsnps.size()>0) {
        update_id_map_rm(badsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(badsnps, _meta_snp_name_map, _meta_remain_snp);
    }

    // For output
    _meta_snp_a1 = snp_a1[0]; _meta_snp_a2 = snp_a2[0];
    _meta_snp_freq = snp_freq;

    nsnp = _meta_remain_snp.size();
    if(nsnp<1) LOGGER.e(0, "no SNP is retained for the GSMR analysis.");
    else LOGGER.i(0, to_string(nsnp) + " SNPs are retained after filtering.");

    // Only keep SNPs with p-value < threshold
    vector<string> keptsnps;
    keptsnps = filter_meta_snp_pval(_meta_snp_name, _meta_remain_snp, _meta_snp_pval, 0, npheno, _snp_val_flag, gwas_thresh);
    if(keptsnps.size()>0) {
        update_id_map_kp(keptsnps, _snp_name_map, _include);
    }

    std::stringstream ss;
    ss << std::scientific << std::setprecision(1) << gwas_thresh;
    LOGGER.i(0, to_string(_include.size()) + " genome-wide significant SNPs with p < " + ss.str() + " are in common among the exposure(s), the outcome(s) and the LD reference sample.\n");
}

eigenMatrix gcta::rho_sample_overlap(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, int nexpo, int noutcome, 
                vector<string> snp_name, vector<int> snp_remain, string ref_ld_dirt, string w_ld_dirt, vector<string> trait_name, int gsmr_so_alg) {

    int i = 0, j = 0;
    eigenMatrix ldsc_intercept(nexpo, noutcome);
    for(i=0; i<nexpo; i++) 
        for(j=0; j<noutcome; j++)
            ldsc_intercept(i,j) = 0.0;
    switch(gsmr_so_alg) {
        case 0 : { ldsc_intercept = sample_overlap_ldsc(snp_val_flag, snp_b, snp_se, snp_n, nexpo, noutcome, 
                                  snp_name, snp_remain, ref_ld_dirt, w_ld_dirt, trait_name);
                break; }
        case 1 : { ldsc_intercept = sample_overlap_rb(snp_val_flag, snp_b, snp_se, snp_pval, snp_n, nexpo, noutcome, 
                                  snp_name, snp_remain, trait_name);
                break; }
    }
    return ldsc_intercept;
}

eigenMatrix gcta::sample_overlap_ldsc(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, int nexpo, int noutcome, 
                vector<string> snp_name, vector<int> snp_remain, string ref_ld_dirt, string w_ld_dirt, vector<string> trait_name) {

    int i = 0, j = 0, ttl_mk_num = 0.0, ntrait = nexpo + noutcome, nsnp = snp_remain.size();
    vector<int> nsnp_cm_trait(ntrait);    
    vector<double> ref_ld_vec, w_ld_vec;
    vector<string> cm_ld_snps;
    vector<vector<bool>> snp_flag;
    map<string,int> ldsc_snp_name_map;
    eigenVector ref_ld, w_ld;
    eigenMatrix bhat_z, bhat_n, ldsc_intercept(nexpo, noutcome);

    // Initialize the parameters
    for(i=0; i<nsnp; i++) {
        ldsc_snp_name_map.insert(pair<string,int>(snp_name[snp_remain[i]], i));
    }
    for(i=0; i<ntrait; i++) nsnp_cm_trait[i] = 0;
    // Read the LD scores
    cm_ld_snps = read_snp_ldsc(ldsc_snp_name_map, snp_name, snp_remain, ttl_mk_num, ref_ld_dirt, w_ld_dirt, ref_ld_vec, w_ld_vec);
     // Re-order the SNP effects and LD scores
    reorder_snp_effect(snp_remain, bhat_z, bhat_n, snp_b, snp_se, snp_n, snp_flag, snp_val_flag, nsnp_cm_trait,
                       cm_ld_snps, ldsc_snp_name_map, ref_ld, w_ld, ref_ld_vec, w_ld_vec, ntrait);
    // Estimate sample overlap
    LOGGER.i(0, "LD score regression analysis to estimate sample overlap between each pair of exposure and outcome ...");
    int n_cm_ld_snps = cm_ld_snps.size();
    vector<double> rst_ldsc(2);
    
    // Univariate LDSC analysis
    LOGGER.i(0, "Univariate LD score regression analysis ...");
    eigenMatrix ldsc_var_h2;
    ldsc_var_h2 = ldsc_snp_h2(bhat_z, bhat_n, ref_ld, w_ld, snp_flag, nsnp_cm_trait, n_cm_ld_snps, ttl_mk_num, trait_name, ntrait);
    
    // Bivariate LDSC analysis
    LOGGER.i(0, "Bivariate LD score regression analysis ...");
    int k = 0, nproc = nexpo * noutcome;
    eigenMatrix ldsc_var_rg;
    vector<int> trait_indx1(nproc), trait_indx2(nproc);
    for( i = 0, k = 0; i < nexpo; i++) {
        for( j = 0; j < noutcome; j++, k++) {
            trait_indx1[k] = i; trait_indx2[k] = j+nexpo;
        }
    }
    ldsc_var_rg = ldsc_snp_rg(ldsc_var_h2, bhat_z, bhat_n, ref_ld, w_ld, snp_flag, trait_indx1, trait_indx2, n_cm_ld_snps, ttl_mk_num, trait_name);
    for( i = 0, k = 0; i < nexpo; i++ ) {
        for( j = 0; j < noutcome; j++, k++) {
            ldsc_intercept(i,j) = ldsc_var_rg(i,0);
        }
    }

    // Print the intercept of rg
    stringstream ss;
    LOGGER.i(0, "Intercept:");
    for(i=0; i<nexpo; i++) {
        ss.str("");
        ss << trait_name[i] <<": ";
        for(j=0; j<noutcome; j++)
            ss << ldsc_intercept(i,j) << " ";
        LOGGER.i(0, ss.str());
    }
    LOGGER.i(0, "LD score regression analysis completed.");

    return ldsc_intercept;
}

eigenMatrix gcta::sample_overlap_rb(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, int nexpo, int noutcome, 
                vector<string> snp_name, vector<int> snp_remain, vector<string> trait_name) {

    int i = 0, j = 0, k = 0, nproc = nexpo * noutcome, nsnp = snp_remain.size();
    double pval_thresh = 0.01;
    eigenMatrix ldsc_intercept(nexpo, noutcome);
    vector<int> trait_indx1(nproc), trait_indx2(nproc);

    // Estimate sample overlap
    LOGGER.i(0, "Using correlation of SNP effects to estimate sample overlap between each pair of exposure and outcome ...");
    // Initialize the variable
    for( i = 0, k = 0; i < nexpo; i++) {
        for( j = 0; j < noutcome; j++, k++) {
            trait_indx1[k] = i; trait_indx2[k] = j+nexpo;
        }
    }

    for( i = 0; i < nproc; i++) {
        // SNPs in common
        int n_cm_snps_buf = 0;
        vector<bool> snp_pair_flag(nsnp);
        for(k = 0; k < nsnp; k++) {
            snp_pair_flag[k] = snp_val_flag[trait_indx1[i]][snp_remain[k]] && snp_val_flag[trait_indx2[i]][snp_remain[k]];
            if(!snp_pair_flag[k]) continue;
            if(snp_pval(snp_remain[k], trait_indx1[i]) < pval_thresh || snp_pval(snp_remain[k], trait_indx2[i]) < pval_thresh)
                snp_pair_flag[k] = 0;
            n_cm_snps_buf += snp_pair_flag[k];
        }
        // Estimate cor(b1, b2)
        eigenVector bhat1_zscore(n_cm_snps_buf), bhat2_zscore(n_cm_snps_buf);
        for(k = 0, j = 0; k<nsnp; k++) {
            if(!snp_pair_flag[k]) continue;
            bhat1_zscore(j) = snp_b(snp_remain[k], trait_indx1[i]) / snp_se(snp_remain[k], trait_indx1[i]);
            bhat2_zscore(j) = snp_b(snp_remain[k], trait_indx2[i]) / snp_se(snp_remain[k], trait_indx2[i]);
            j++;
        }

        double average1 = bhat1_zscore.mean(), average2 = bhat2_zscore.mean();
        bhat1_zscore = bhat1_zscore - average1*VectorXd::Ones(n_cm_snps_buf);
        bhat2_zscore = bhat2_zscore - average2*VectorXd::Ones(n_cm_snps_buf);
        double sd1 = bhat1_zscore.norm(), sd2 = bhat2_zscore.norm();
        ldsc_intercept(trait_indx1[i],trait_indx2[i]-nexpo) = bhat1_zscore.dot(bhat2_zscore)/(sd1*sd2);
    }

    // Print the correlation
    stringstream ss;
    LOGGER.i(0, "Correlation:");
    for(i=0; i<nexpo; i++) {
        ss.str("");
        ss << trait_name[i] <<": ";
        for(j=0; j<noutcome; j++)
            ss << ldsc_intercept(i,j) << " ";
        LOGGER.i(0, ss.str());
    }
    LOGGER.i(0, "Correlation computation completed.");

    return ldsc_intercept;
}

void collect_snp_instru(stringstream &ss, map<string,int> &snp_instru_map, int expo_indx, int outcome_indx, vector<string> snp_instru) {
    int i = 0, nsnp_instru = snp_instru.size();

    ss << expo_indx << " " << outcome_indx << endl;
    // With HEIDI-outlier test
    if(nsnp_instru == 0) {
        ss << "NA";
    } else {
        ss << snp_instru[0];
        for( i = 1; i < nsnp_instru; i++ ) 
            ss << " " << snp_instru[i];
    }
    ss << endl;
    // Save the SNP instruments
    for(i=0; i<nsnp_instru; i++) 
        snp_instru_map.insert(pair<string, int>(snp_instru[i], i));
}

void collect_gsmr_trait(stringstream &ss, vector<string> pheno_str, int nexpo, int noutcome) {
    int i = 0;
    
    ss << "#trait_begin" <<endl;
    // Exposure
    ss << pheno_str[0];
    for(i=1; i<nexpo; i++) 
        ss << " " << pheno_str[i];
    ss << endl;
    // Outcome
    ss << pheno_str[nexpo];
    for(i=1; i<noutcome; i++) 
        ss << " " << pheno_str[i+nexpo];
    ss << endl;
    ss << "#trait_end" <<endl;
}

void collect_gsmr_result(stringstream &ss, vector<vector<double>> bxy_est, int gsmr_alg_flag, vector<string> pheno_name, int expo_num_buf, int outcome_num_buf, int gsmr_beta_version, int n_gsmr_rst_item) {
    bool output_forward = false, output_reverse = false;
    int i=0, j=0, k=0, t=0;

    switch(gsmr_alg_flag) {
        case 0 : output_forward=true; break;
        case 1 : output_reverse=true; break;
        case 2 : output_forward=true; output_reverse=true; break;
    }

    if(gsmr_beta_version) {
        ss << "Exposure\tOutcome\tbxy\tse\tp\tnsnp\tmulti_snp_based_heidi_outlier" <<endl;
    } else {
        ss << "Exposure\tOutcome\tbxy\tse\tp\tnsnp" <<endl;
    }
    
    k=0; 
    if(output_forward) {
        for(i=0; i<expo_num_buf; i++) {
            for(j=0; j<outcome_num_buf; j++, k++) {
                ss << pheno_name[i] << "\t" << pheno_name[expo_num_buf+j];
                for(t=0; t<n_gsmr_rst_item; t++) 
                    ss << '\t' << bxy_est[t][k];
                ss << endl;
            }
        }
    }
    if(output_reverse) {
        for(i=0; i<outcome_num_buf; i++) {
            for(j=0; j<expo_num_buf; j++, k++) {
                ss << pheno_name[expo_num_buf+i] << "\t" << pheno_name[j];
                for(t=0; t<n_gsmr_rst_item; t++) 
                    ss << "\t" << bxy_est[t][k];
                ss << endl;
            }
        }
    }
}

void collect_snp_instru_effect(stringstream &ss, vector<vector<bool>> snp_flag, map<string, int> snp_instru_map, map<string, int> meta_snp_map, vector<string> snp_a1, vector<string> snp_a2, eigenVector snp_freq, eigenMatrix snp_b, eigenMatrix snp_se) {
    map<string, int>::iterator iter1, iter2;
    int i = 0, snpindx = 0, npheno = snp_b.cols();
    string snpbuf = "";

    ss << "#effect_begin" << endl;
    for(iter1 = snp_instru_map.begin(); iter1 != snp_instru_map.end(); iter1++) {
        snpbuf = iter1 -> first;
        iter2 = meta_snp_map.find(snpbuf); 
        if(iter2 == meta_snp_map.end()) continue;
        snpindx = iter2 -> second;
        ss << snpbuf << " " << snp_a1[snpindx] << " " << snp_a2[snpindx] << " " << snp_freq(snpindx);
        for( i = 0; i < npheno; i++ ) {
            if(snp_flag[i][snpindx]) ss << " " << snp_b(snpindx, i) << " " << snp_se(snpindx, i);
            else ss << " nan nan";
        }
        ss << endl;
    }
    ss << "#effect_end" << endl;
}

void gcta::gsmr(int gsmr_alg_flag, string ref_ld_dirt, string w_ld_dirt, double freq_thresh, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, bool o_snp_instru_flag, int gsmr_so_alg, int gsmr_beta_version) {
    if(gsmr_beta_version) { _n_gsmr_rst_item = 5; _gsmr_beta_version = 1; }
    else { _n_gsmr_rst_item = 4; _gsmr_beta_version = 0; }
   
    vector<vector<double>> bxy_est;
    map<string, int> snp_instru_map;
    vector<string> afsnps;
    std::stringstream ss, ss_gsmr, ss_pleio;
    ss << "#marker_begin" <<endl;
    int ntrait = _expo_num + _outcome_num;

    // Calculate allele frequency
    if (_mu.empty()) calcu_mu();
    LOGGER.i(0, "Checking allele frequencies among the GWAS summary data and the reference sample...");
    afsnps = remove_freq_diff_snps(_meta_snp_name, _meta_remain_snp, _snp_name_map, _mu, _meta_snp_freq, _snp_val_flag, ntrait, freq_thresh, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(afsnps, _meta_snp_name_map, _meta_remain_snp);
    }

    // Remove monomorphic SNPs
    afsnps = remove_mono_snps(_snp_name_map, _mu, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(afsnps, _meta_snp_name_map, _meta_remain_snp);
    }

    // Only keep the AF for the target trait
    int nsnp_freq = _meta_snp_name_map.size();
    _meta_snp_freq.conservativeResize(nsnp_freq,1);

    // Estimate intercept from LDSC regression
    _r_pheno_sample = rho_sample_overlap(_snp_val_flag, _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, _expo_num, _outcome_num, 
                                        _meta_snp_name, _meta_remain_snp, ref_ld_dirt, w_ld_dirt, _gwas_trait_name, gsmr_so_alg);

   // GSMR analysis
    switch(gsmr_alg_flag) {
        case 0 : { 
            bxy_est = forward_gsmr(ss, snp_instru_map, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, ss_pleio); 
            break;
        }
        case 1 : {
            bxy_est = reverse_gsmr(ss, snp_instru_map, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, ss_pleio);
            break;
        }
        case 2 : {
            vector<vector<double>> bxy_est_buf;
            bxy_est = forward_gsmr(ss, snp_instru_map, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, ss_pleio);
            bxy_est_buf = reverse_gsmr(ss, snp_instru_map, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, ss_pleio);
            int i = 0;
            for(i=0; i<_n_gsmr_rst_item; i++) bxy_est[i].insert(bxy_est[i].end(), bxy_est_buf[i].begin(), bxy_est_buf[i].end());
            break;
        }
    }
    ss << "#marker_end" << endl;
    LOGGER.i(0, "");

    // Save pleiotropic SNPs
    bool pleio_flag = ss_pleio.str().size() > 0 ? true : false;   
    if(pleio_flag) {
        string pleio_snpfile = _out + ".pleio_snps";	
        ofstream o_pleio_snp(pleio_snpfile.c_str());	
        if(!o_pleio_snp) LOGGER.e(0, "cannot open file [" + pleio_snpfile + "] to write pleiotropic SNPs.");
        o_pleio_snp << ss_pleio.str();            
        o_pleio_snp.close();	
        LOGGER.i(0,  "The pleiotropic SNPs filtered by HEIDI-outlier analysis have been saved in [" + pleio_snpfile + "].");
    }

    // Output GSMR result
    collect_gsmr_result(ss_gsmr, bxy_est, gsmr_alg_flag, _gwas_trait_name, _expo_num, _outcome_num, _gsmr_beta_version, _n_gsmr_rst_item);
    // Output the SNP instruments
    if(o_snp_instru_flag) {
        int nsnp_instru = snp_instru_map.size();
        if(nsnp_instru >= nsnp_gsmr) {
            stringstream ss_effect, ss_pheno;
            collect_gsmr_trait(ss_pheno, _gwas_trait_name, _expo_num, _outcome_num);
            collect_snp_instru_effect(ss_effect, _snp_val_flag, snp_instru_map, _meta_snp_name_map, _meta_snp_a1, _meta_snp_a2, _meta_snp_freq, _meta_snp_b, _meta_snp_se);
            string output_filename = _out + ".eff_plot.gz";
            LOGGER.i(0, "Saving the SNP instruments for the GSMR plots to [" + output_filename + "] ...");
            gzofstream zofile(output_filename.c_str());
            if (!zofile) LOGGER.e(0, "cannot open the file [" + output_filename + "] to write.");
            zofile << ss_pheno.str() << "#gsmr_begin" << endl << ss_gsmr.str() 
                   << "#gsmr_end" << endl << ss_effect.str() << ss.str();
            zofile.close();
        } else {
            LOGGER.w(0, "Not enough SNP instruments to be saved in the compressed text file.");
        }
    }

    // Output GSMR results
    string output_filename = _out + ".gsmr";
    LOGGER.i(0, "Saving the GSMR analyses results of " + to_string(_expo_num) + " exposure(s) and "
                 + to_string(_outcome_num) + " outcome(s) to [" + output_filename + "] ...");
    ofstream ofile(output_filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + output_filename + "] to write.");             
    ofile << ss_gsmr.str();
    ofile.close();
    LOGGER.i(0, "\nGSMR analyses completed.");
}

vector<vector<double>> gcta::forward_gsmr(stringstream &ss, map<string,int> &snp_instru_map, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, stringstream &ss_pleio) {
    int i=0, j=0, k=0, t=0, m=_expo_num*_outcome_num, nsnp = _meta_remain_snp.size();
    vector<bool> snp_pair_flag(nsnp);
    vector<vector<double>> bxy_est;
     
    bxy_est.resize(_n_gsmr_rst_item);
    for(i=0; i<_n_gsmr_rst_item; i++) bxy_est[i].resize(m);
    
    // GSMR analysis
    vector<double> gsmr_rst(_n_gsmr_rst_item);
    for(i=0, t=0; i<_expo_num; i++) {
        for(j=0; j<_outcome_num; j++, t++) {
            string err_msg = "";
            string pleio_snps = "";
            vector<string> snp_instru;
            for(k=0; k<nsnp; k++) snp_pair_flag[k] = _snp_val_flag[i][_meta_remain_snp[k]] && _snp_val_flag[j+_expo_num][_meta_remain_snp[k]];
            LOGGER.i(0, "\nForward GSMR analysis for exposure #" + to_string(i+1) + " and outcome #" + to_string(j+1) + " ...");         
            gsmr_rst =  gsmr_meta(snp_instru, _meta_snp_b.col(i), _meta_snp_se.col(i), _meta_snp_pval.col(i), 
                                  _meta_snp_b.col(j+_expo_num), _meta_snp_se.col(j+_expo_num), _meta_snp_pval.col(j+_expo_num), _r_pheno_sample(i,j), snp_pair_flag, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, pleio_snps, err_msg);
            if(std::isnan(gsmr_rst[3]))
                LOGGER.w(0, err_msg);
            else
                LOGGER.i(0, "Forward GSMR analysis for exposure #" + to_string(i+1) + " and outcome #" + to_string(j+1) + " completed.");
            for(k=0; k<_n_gsmr_rst_item; k++) bxy_est[k][t] = gsmr_rst[k];
            // Saving pleiotropic SNPs
            if(pleio_snps.size() > 0) {
                ss_pleio << _gwas_trait_name[i] << " " << _gwas_trait_name[j+_expo_num] << " " << pleio_snps << endl;
            }
            // Saving the SNP instruments
            collect_snp_instru(ss, snp_instru_map, i+1, j+_expo_num+1, snp_instru); 
        }
    }


    return bxy_est;
}

vector<vector<double>> gcta::reverse_gsmr(stringstream &ss, map<string,int> &snp_instru_map, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, stringstream &ss_pleio) {
     int i=0, j=0, k=0, t=0, m=_expo_num*_outcome_num, nsnp = _meta_remain_snp.size();
     vector<bool> snp_pair_flag(nsnp);
     vector<vector<double>> bxy_est;
     
     bxy_est.resize(_n_gsmr_rst_item);
     for(i=0; i<_n_gsmr_rst_item; i++) bxy_est[i].resize(m);
    
    // GSMR analysis
    vector<double> gsmr_rst(_n_gsmr_rst_item);
    for(i=0, t=0; i<_outcome_num; i++) {
        for(j=0; j<_expo_num; j++, t++) {
            string err_msg = "";
            string pleio_snps = "";
            vector<string> snp_instru;
            for(k=0; k<nsnp; k++) snp_pair_flag[k] = _snp_val_flag[i+_expo_num][_meta_remain_snp[k]] && _snp_val_flag[j][_meta_remain_snp[k]];
            LOGGER.i(0, "\nReverse GSMR analysis for exposure #" + to_string(j+1) + " and outcome #" + to_string(i+1) + " ...");
            gsmr_rst =  gsmr_meta(snp_instru, _meta_snp_b.col(i+_expo_num), _meta_snp_se.col(i+_expo_num), _meta_snp_pval.col(i+_expo_num), 
                                  _meta_snp_b.col(j), _meta_snp_se.col(j), _meta_snp_pval.col(j), _r_pheno_sample(j,i), snp_pair_flag, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, pleio_snps, err_msg);
            if(std::isnan(gsmr_rst[3])) 
                LOGGER.w(0, err_msg);
            else
                LOGGER.i(0, "Reverse GSMR analysis for exposure #" + to_string(j+1) + " and outcome #" + to_string(i+1) + " completed.");
            for(k=0; k<_n_gsmr_rst_item; k++) bxy_est[k][t] = gsmr_rst[k]; 
            // Saving pleiotropic SNPs
            if(pleio_snps.size() > 0) {
                ss_pleio << _gwas_trait_name[i+_expo_num] << " " << _gwas_trait_name[j] << " " << pleio_snps << endl;
            }            
            // Saving the SNP instruments
            collect_snp_instru(ss, snp_instru_map, i+_expo_num+1, j+1, snp_instru);
        }
    }

    return bxy_est;
}
