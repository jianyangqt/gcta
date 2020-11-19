#include "gcta.h"
#include "Logger.h"
#include <iterator>

void gcta::read_pc_adjust_file(string pcadjust_list_file, string pc_file) {

    int ncovar=0, i=0, j=0;
    vector<string> pheno_file, snplist;
    string strbuf = "", eigenvalue_file = pc_file + ".eigenval";
    std::map<string,int>::iterator iter;

    // Read the summary data
    LOGGER.i(0, "\nReading eigenvalues from [" + eigenvalue_file + "]...");
    // Read eigenvalue
    int line_number = 0;
    std::istringstream linebuf;
    ifstream in_lambda(eigenvalue_file.c_str());
    if (!in_lambda) LOGGER.e(0, "Cannot open the file [" + eigenvalue_file + "] to read.");
    while(std::getline(in_lambda, strbuf)) {
        line_number++;
        linebuf.clear();
        linebuf.str(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() > 1) {
            LOGGER.w(0, "Only the first element would be accepted. File [" + eigenvalue_file + "], line " + to_string(line_number) + ".");
        }
        _eigen_value.push_back(atof(line_elements[0].c_str()));
    }
    in_lambda.close();
    
    // Read the file list for PC adjust analysis
    line_number = 0;
    ifstream meta_list(pcadjust_list_file.c_str());
    if (!meta_list) LOGGER.e(0, "Cannot open the file [" + pcadjust_list_file + "] to read.");
    while(std::getline(meta_list, strbuf)) {
        line_number++;
        linebuf.clear();
        linebuf.str(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() > 1) {
            LOGGER.w(0, "Only the first element would be accepted. File [" + pcadjust_list_file + "], line " + to_string(line_number) + ".");
        }
        pheno_file.push_back(line_elements[0]);
    }
    meta_list.close();

    // Read the SNPs to know the total number of markers
    map<string,int> gws_snp_name_map;
    ncovar = pheno_file.size() - 1;
    if(ncovar == 0)
       LOGGER.e(0, "At least 1 PC loading is required.");
    if(ncovar > _eigen_value.size())
        LOGGER.e(0, "There are " + to_string(ncovar) + " covariates summary data. " + to_string(_eigen_value.size()) + " eigenvalues are provided.");

    // Covariates
    for( i=1; i<=ncovar; i++) {
        if(pheno_file[i].substr(pheno_file[i].length()-3,3)!=".gz")
            snplist=read_snp_metafile_txt(pheno_file[i], gws_snp_name_map, -9);
        else 
            snplist=read_snp_metafile_gz(pheno_file[i], gws_snp_name_map, -9);
        if( i == 1 ) init_meta_snp_map(snplist, _meta_snp_name_map, _meta_snp_name, _meta_remain_snp);
        else update_id_map_kp(snplist, _meta_snp_name_map, _meta_remain_snp);
        update_id_map_kp(snplist, _snp_name_map, _include);
    }
    _ttl_mk_num = _meta_snp_name_map.size();

    // Target trait
    if(pheno_file[0].substr(pheno_file[0].length()-3,3)!=".gz")
        snplist=read_snp_metafile_txt(pheno_file[0], gws_snp_name_map, -9);
    else 
        snplist=read_snp_metafile_gz(pheno_file[0], gws_snp_name_map, -9);
    update_id_map_kp(snplist, _meta_snp_name_map, _meta_remain_snp);
    update_id_map_kp(snplist, _snp_name_map, _include);

    // Initialization of variables
    int nsnp = _meta_snp_name_map.size();
    eigenMatrix snp_freq;
    vector<vector<string>> snp_a1, snp_a2;
    init_gwas_variable(snp_a1, snp_a2, snp_freq, _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, ncovar+1, nsnp); 

    // reset SNP variables
    update_meta_snp(_meta_snp_name_map, _meta_snp_name, _meta_remain_snp);
    LOGGER.i(0, to_string(nsnp) + " SNPs in common between the summary data and the PC loading(s).");

    _snp_val_flag.clear(); _snp_val_flag.resize(ncovar+1);
    for(i=0; i<ncovar+1; i++) {
        _snp_val_flag[i].resize(nsnp);
        for(j=0; j<nsnp; j++) _snp_val_flag[i][j] = false;
    }
    
    // Read the meta analysis
    _meta_vp_trait.resize(ncovar+1);
    // Target trait and Covariates
    for(i=0; i<=ncovar; i++) {
        eigenVector snp_freq_buf(nsnp), snp_b_buf(nsnp), snp_se_buf(nsnp), snp_pval_buf(nsnp), snp_n_buf(nsnp);
        if(pheno_file[i].substr(pheno_file[i].length()-3,3)!=".gz")
            _meta_vp_trait(i) = read_single_metafile_txt(pheno_file[i], _meta_snp_name_map, snp_a1[i], snp_a2[i], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i]);
        else _meta_vp_trait(i) = read_single_metafile_gz(pheno_file[i], _meta_snp_name_map, snp_a1[i], snp_a2[i], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i]);
        if(_meta_vp_trait(i) < 0) {
            if(i==0) LOGGER.e(0, "Negative phenotypic variance of the target trait.");
            else LOGGER.e(0, "Negative phenotypic variance of the covariate #" + to_string(i+1) + ".");
        }
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
                             _snp_name_map, _allele1, _allele2, 1, ncovar, _out);                             
    if(badsnps.size()>0) {
        update_id_map_rm(badsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(badsnps, _meta_snp_name_map, _meta_remain_snp);
    }

    // For output
    _meta_snp_a1 = snp_a1[0]; _meta_snp_a2 = snp_a2[0];
    _meta_snp_freq = snp_freq;  
   
    if(nsnp<1) LOGGER.e(0, "None SNPs are retained after filtering.");
    else LOGGER.i(0, to_string(nsnp) + " SNPs are retained after filtering.");
    LOGGER.i(0, to_string(_include.size()) + " SNPs are in common between the summary data and the LD reference sample.");
}

vector<string> update_snp_freq(vector<string> meta_snp_name, vector<int> meta_snp_remain, map<string,int> snp_name_map, vector<double> ref_freq, eigenMatrix meta_freq, vector<vector<bool>> snp_flag, int npheno, string outfile_name) {
    int i = 0, j = 0, nsnp = meta_snp_remain.size();
    string snpbuf="";
    vector<string> afsnps;
    
    for( i=0; i<nsnp; i++ ) {
        snpbuf = meta_snp_name[meta_snp_remain[i]];
        // AF from the reference sample
        int refsnp_index = 0;
        double a1_freq = 0.0;
        bool freq_flag = false;
        map<string,int>::iterator iter_ref;
        iter_ref = snp_name_map.find(snpbuf);
        if( iter_ref != snp_name_map.end()) {
            refsnp_index = iter_ref -> second;
            a1_freq = ref_freq[refsnp_index]/2;
            freq_flag = true;
        }
            
        // AF from the GWAS summary
        for( j=0; j<npheno; j++) {
            if(!snp_flag[j][meta_snp_remain[i]]) continue;
            if(!std::isnan(meta_freq(meta_snp_remain[i],j))) {
                if(!freq_flag) a1_freq = meta_freq(meta_snp_remain[i],j);
                freq_flag = true;
            }
            meta_freq(meta_snp_remain[i],j) = a1_freq;
        }

        // Update AF
        for( j=0; (j<npheno) & freq_flag; j++) {
            if(!snp_flag[j][meta_snp_remain[i]]) continue;
            if(std::isnan(meta_freq(meta_snp_remain[i],j))) 
                meta_freq(meta_snp_remain[i],j) = a1_freq;
        }

        if(!freq_flag) {
            afsnps.push_back(snpbuf);
        }
    }  

    if (!afsnps.empty()) {
        string afsnpfile = outfile_name + ".miss_freq.badsnps", strbuf="";
        ofstream oafsnp(afsnpfile.c_str());
        if(!oafsnp) LOGGER.e(0, "Cannot open file [" + afsnpfile + "] to write bad SNPs.");
        int nafsnps = afsnps.size();
        for (i = 0; i < nafsnps; i++) oafsnp << afsnps[i] << endl;
        oafsnp.close();
        LOGGER.i(0,  to_string(nafsnps) + " SNP(s) do not have allele frequency. These SNPs have been saved in [" + afsnpfile + "].");
        if(nafsnps > nsnp*0.05) 
            LOGGER.e(0, "There are too many SNPs where allele frequencies are not available. Please check your summary datasets.");
    }
    return(afsnps);
}

double init_pc_meta(vector<int> snp_remain, eigenVector snp_freq, eigenVector snp_b, eigenVector snp_se, eigenVector &snp_n_o, double vp_trait) {
    
    int i = 0, m1 = snp_remain.size(), m2 = snp_n_o.size();

    snp_n_o.resize(m2); snp_n_o.setZero(m2);
    // Update sample size per SNP
    for( i = 0; i < m1; i ++) {
        double vf = 2*snp_freq(snp_remain[i])*(1-snp_freq(snp_remain[i]));
        snp_n_o(snp_remain[i]) = (vp_trait - vf*snp_b(snp_remain[i])*snp_b(snp_remain[i]))/(vf*snp_se(snp_remain[i])*snp_se(snp_remain[i])) + 1.0;
    }

    // Estimate total sample size
    eigenVector snp_n_sort(snp_n_o);
    std::stable_sort(snp_n_sort.data(), snp_n_sort.data()+m2);
    double n_o = CommFunc::quantile(snp_n_sort, 0.50);

    return n_o;
}

double est_bxy_pc(vector<int> snp_remain, eigenVector bzy, eigenVector bzy_freq, eigenVector &bzx, eigenVector bzx_freq, eigenVector bzx_n, double n_o, int ttl_mk_num, double eigen_value)
{
    int i = 0, m = snp_remain.size();
    double bxy = 0.0;
    eigenVector bzy_v;

    for(i=0; i<m; i++) {
        // Standardise SNP effects
        double vf1 = 2*bzx_freq(snp_remain[i])*(1-bzx_freq(snp_remain[i]));
        double vf2 = 2*bzy_freq(snp_remain[i])*(1-bzy_freq(snp_remain[i]));
        double b1 = sqrt(vf1)*bzx(snp_remain[i]), b2 = sqrt(vf2)*bzy(snp_remain[i]);
        bxy += b1*b2*bzx_n(snp_remain[i]);
        bzx(snp_remain[i]) = b1;
    }
    bxy *= n_o/((double)ttl_mk_num*eigen_value);

    return bxy;
}

void restrict_snp_effect(map<string,int> meta_snp_name_map, vector<string> snp_name, vector<int> snp_remain, 
                         eigenMatrix &snp_freq, eigenMatrix &snp_b, eigenMatrix &snp_se, eigenMatrix &snp_n_o, 
                         eigenVector &bzy, eigenVector &bzy_se, eigenVector &bzy_n, int ncovar) {
    int m = snp_remain.size();
    eigenMatrix snp_freq_tmp(m, 1), snp_b_tmp(m, ncovar), snp_se_tmp(m, ncovar), snp_n_tmp(m, ncovar);

    bzy.resize(m); bzy_se.resize(m); bzy_n.resize(m);
    #pragma omp parallel for
    for(int i = 0; i < m; i++) {
        map<string, int>::iterator iter = meta_snp_name_map.find(snp_name[snp_remain[i]]);
        snp_freq_tmp(i, 0) = snp_freq(iter->second, 1);
        bzy(i) = snp_b(iter->second, 0); bzy_se(i) = snp_se(iter->second, 0);
        bzy_n(i) = snp_n_o(iter->second, 0);
        for(int j = 0; j < ncovar; j++) {
            snp_b_tmp(i, j) = snp_b(iter->second, j+1);
            snp_se_tmp(i, j) = snp_se(iter->second, j+1);
            snp_n_tmp(i, j) = snp_n_o(iter->second, j+1);
        }
    }

    snp_freq.resize(m, 1); snp_b.resize(m, ncovar); snp_se.resize(m, ncovar); snp_n_o.resize(m,ncovar);
    snp_freq = snp_freq_tmp; snp_b = snp_b_tmp; snp_se = snp_se_tmp; snp_n_o = snp_n_tmp;
}

void init_ld_snp_index(int snp_start, int &snp_end, vector<int> snp_remain, vector<int> snp_chr, vector<int> snp_bp, int wind_size) {
    int snp_index = snp_start, snp_test_start = snp_remain[snp_start], snp_test_end = 0, m = snp_remain.size();

    while(1) {
        snp_index++;
        if(snp_index>=m) break;
        snp_test_end = snp_remain[snp_index];
        if(snp_chr[snp_test_end] == snp_chr[snp_test_start] && abs(snp_bp[snp_test_end] - snp_bp[snp_test_start]) < wind_size) {
            continue;
        }
        else break;
    }
    snp_end = --snp_index;
}

double calcu_ztz_product_b(map<string,int> meta_snp_name_map, vector<string> snp_name, vector<int> snp_remain, 
                           eigenVector snp_delta, eigenVector ztz_vec, double ztz_target, 
                           int snp_target, int snp_wind_start, int snp_wind_end, vector<int> snp_chr, vector<int> snp_bp, 
                           vector<double> lambda, int wind_size, double nsnps, int ncovar) {
    int i = 0, j = 0, snp_test_start = -1, snp_test_end = -1;

    // Start SNP
    for(i=snp_wind_start; i<=snp_wind_end; i++) {
        if(abs(snp_bp[snp_remain[snp_target]] - snp_bp[snp_remain[i]]) >= wind_size) continue;
        snp_test_start = i; break;
    }
    // End SNP
    for(i=snp_wind_end; i>=snp_wind_start; i--) {
        if(abs(snp_bp[snp_remain[snp_target]] - snp_bp[snp_remain[i]]) >= wind_size) continue;
        snp_test_end = i; break;
    }
    
    int nsnp_test = snp_test_end - snp_test_start + 1;
    eigenVector ztz_tmp = ztz_vec.segment(snp_test_start, nsnp_test), delta = snp_delta.segment(snp_test_start, nsnp_test);
    double d_buf = ztz_tmp.dot(delta);
    d_buf /= (ztz_target*(double)nsnps);

    return d_buf;
}

void gcta::adjust_snp_effect_for_pc(eigenVector &bzy_adj, eigenVector &bzx_hat, eigenVector bzy, eigenVector bxy_hat, int wind_size) {
    int snp_index = 0, m = _include.size(), n = _keep.size(), ncovar = _meta_snp_b.cols();
//cout<<"ncovar "<<ncovar<<endl;

    // Estimate var(z)
    eigenVector msx(m);
    for(int i=0; i<m; i++) {
        eigenVector z_buf(n);
        makex_eigenVector(_include[i], z_buf, false, true);
        msx(i) = z_buf.squaredNorm()/(double)(n-1);
    }

    // Open the window
    int snp_start = -1, snp_end = -1, snp_buf1_start = -1, snp_buf1_end = -1, snp_buf2_start = -1, snp_buf2_end = -1; 
    snp_start = 0;
    init_ld_snp_index(snp_start, snp_end, _include, _chr, _bp, wind_size);
    snp_buf2_start = snp_end + 1;
    init_ld_snp_index(snp_buf2_start, snp_buf2_end, _include, _chr, _bp, wind_size);
    // adjust SNP effects
    eigenVector bxy_adj(ncovar);
    for(int i=0; i<ncovar; i++) bxy_adj(i) = bxy_hat(i)/_eigen_value[i];
    eigenVector snp_delta = (_meta_snp_b.array()*_meta_snp_n_o.array()).matrix()*bxy_adj;

eigenVector bxy_adj2(ncovar);
for(int i=0; i<ncovar; i++) bxy_adj2(i) = 1/_eigen_value[i];
eigenVector snp_delta2 = (_meta_snp_b.array()*_meta_snp_n_o.array()).matrix()*bxy_adj2;

    // determine start and end of each chromosome
    vector<int> chr_start_index, chr_end_index;
    chr_start_index.push_back(_chr[0]);
    for(int i=1, j=0; i<m; i++) {
        if(_chr[i]==chr_start_index[j]) continue;
        chr_start_index.push_back(i);
        chr_end_index.push_back(i-1);
    }
    chr_end_index.push_back(m-1);
    int nchr = chr_start_index.size();

//for(int i=0; i<nchr; i++) {
//    cout<<"i "<< i+1 << " " << chr_start_index[i] << " "<< chr_end_index[i] << endl;
//}    
    eigenMatrix x_sub, x_sub_buf1, x_sub_buf2;
    eigenMatrix ztz_buf1, ztz_in, ztz_buf2; 
    while(1) {
        int m_buf1 = snp_buf1_end - snp_buf1_start + 1, m_in = snp_end - snp_start + 1, m_buf2 = snp_buf2_end - snp_buf2_start + 1;

        eigenVector z_buf(n);
        
        if(m_buf1==1) {
            // the first iteration
            // re-code for SNPs in the window
            x_sub.resize(n, m_in);
            for(int i=snp_start; i<=snp_end; i++) {
                makex_eigenVector(_include[i], z_buf, false, true);
                z_buf = z_buf/sqrt(msx(i));
                x_sub.col(i-snp_start) = z_buf;
            }
        } else {
            x_sub = x_sub_buf2;
        }      
        // buf 1
        if(m_buf1>1) {
            ztz_buf1 = ztz_buf2.transpose();
        }
        // within window
        ztz_in = x_sub.transpose()*x_sub/(double)(n-1);
        // buf 2
        if(m_buf2>1) {          
            // re-code for SNPs in the left buf window
            x_sub_buf2.resize(n, m_buf2);
            for(int i = snp_buf2_start; i <= snp_buf2_end; i++) {
                makex_eigenVector(_include[i], z_buf, false, true);
                z_buf = z_buf/sqrt(msx(i));
                x_sub_buf2.col(i-snp_buf2_start) = z_buf;
            }         
            ztz_buf2 = x_sub.transpose()*x_sub_buf2/(double)(n-1);                 
        }
//LOGGER<<"snp "<<snp_buf1_start << " "<<m_buf1 << " "<< snp_start<<" "<< m_in << " "<< snp_buf2_start<<" "<<m_buf2<<endl;
        // adjust SNP effect
        eigenVector dbuf1(m_in), dbuf2(m_in), dbuf3(m_in);
        dbuf1.setZero(m_in); dbuf2.setZero(m_in); dbuf3.setZero(m_in);
        if(m_buf1>1) dbuf1 = (ztz_buf1*snp_delta.segment(snp_buf1_start, m_buf1)).array()/msx.segment(snp_start, m_in).array().sqrt();             
        dbuf2 = (ztz_in*snp_delta.segment(snp_start, m_in)).array()/msx.segment(snp_start, m_in).array().sqrt();
        if(m_buf2>1) dbuf3 = (ztz_buf2*snp_delta.segment(snp_buf2_start, m_buf2)).array()/msx.segment(snp_start, m_in).array().sqrt();
        eigenVector bzy_delta(m_in);
        int nsnp_region = m_in;
        if(m_buf1>1) nsnp_region += m_buf1;
        if(m_buf2>1) nsnp_region += m_buf2;      
        bzy_delta = bzy.segment(snp_start, m_in) - (dbuf1 + dbuf2 + dbuf3)/(double)nsnp_region;               
        bzy_adj.segment(snp_start, m_in) = bzy_delta;
            
// Test bzx
//cout<<"nsnp_region "<<nsnp_region<<endl;
//LOGGER<<"one window complete"<<endl;        
eigenVector t1(m_in), t2(m_in), t3(m_in);
t1.setZero(m_in); t2.setZero(m_in); t3.setZero(m_in);
if(m_buf1>1) t1 = (ztz_buf1*snp_delta2.segment(snp_buf1_start, m_buf1)).array()/msx.segment(snp_start, m_in).array().sqrt();             
t2 = (ztz_in*snp_delta2.segment(snp_start, m_in)).array()/msx.segment(snp_start, m_in).array().sqrt();
if(m_buf2>1) t3 = (ztz_buf2*snp_delta2.segment(snp_buf2_start, m_buf2)).array()/msx.segment(snp_start, m_in).array().sqrt();
eigenVector bzx_delta(m_in);    
bzx_delta = (t1 + t2 + t3)/(double)nsnp_region;             
bzx_hat.segment(snp_start, m_in) = bzx_delta;

        // move the window
        if(snp_end == m-1) break;
        bool end_chr_flag1 = false, end_chr_flag2 = false;
        for(int i=0; i<nchr; i++) {
            if(snp_end == chr_end_index[i]) 
                end_chr_flag1=true;
            if(snp_buf2_end == chr_end_index[i])
                end_chr_flag2=true;
        }
        if(end_chr_flag1) {
            snp_buf1_start = snp_buf1_end = -1;
            snp_start = snp_end + 1;
            init_ld_snp_index(snp_start, snp_end, _include, _chr, _bp, wind_size);
            snp_buf2_start = snp_end + 1;
            init_ld_snp_index(snp_buf2_start, snp_buf2_end, _include, _chr, _bp, wind_size);
        } else if(end_chr_flag2) {
            snp_buf1_start = snp_start; snp_buf1_end = snp_end;
            snp_start = snp_buf2_start; snp_end = snp_buf2_end;
            snp_buf2_start = snp_buf2_end = -1;           
        } else {
            snp_buf1_start = snp_start; snp_buf1_end = snp_end;
            snp_start = snp_buf2_start; snp_end = snp_buf2_end;
            snp_buf2_start = snp_end + 1;
            init_ld_snp_index(snp_buf2_start, snp_buf2_end, _include, _chr, _bp, wind_size);
        }
    }
}

void output_snp_effect_for_pc(string output_file, vector<string> meta_snp_name, vector<int> meta_snp_remain, 
                              map<string,int> snp_name_map, vector<string> snp_name, vector<string> snp_a1, vector<string> snp_a2,
                              eigenVector snp_freq, eigenVector snp_b, eigenVector snp_se, eigenVector snp_pval, eigenVector snp_n, 
                              eigenVector snp_b_adj, eigenVector bzx_hat) {
    int i=0, meta_nsnp = meta_snp_remain.size();
    
    output_file = output_file + ".pcadj.cma";
    ofstream ofile(output_file.c_str());
    if (!ofile) LOGGER.e(0, "Cannot open the file [" + output_file + "] to write.");
    
    ofile << "SNP\tA1\tA2\tfreq\tb\tse\tp\tN\tbC\tbzx" <<endl;
    for (i = 0; i < meta_nsnp; i++) {
        string snpbuf = meta_snp_name[meta_snp_remain[i]];
        map<string,int>::iterator iter = snp_name_map.find(snpbuf);
        if(iter==snp_name_map.end()) continue;
        int i_buf = iter->second;
        ofile << snp_name[i_buf] << "\t" <<snp_a1[meta_snp_remain[i]] << "\t" << snp_a2[meta_snp_remain[i]] << "\t" << snp_freq(i_buf)
              << "\t" << snp_b(i_buf) << "\t" << snp_se(i_buf) << "\t"  << snp_pval(meta_snp_remain[i]) << "\t" << snp_n(i_buf)
              << "\t" << snp_b_adj(i_buf) << "\t" << bzx_hat(i_buf)<<endl;
    }
    ofile.close();
}

void gcta::pc_adjust(string pcadjust_list_file, string pc_file, double freq_thresh, int wind_size) {

    // Read the summary data
    read_pc_adjust_file(pcadjust_list_file, pc_file);

    int i = 0, j = 0, npheno = _meta_snp_b.cols();

    // Calculate allele frequency
    if (_mu.empty()) calcu_mu();

    // Check allele frequency
    vector<string> afsnps;
    LOGGER.i(0, "Checking allele frequencies among the GWAS summary data and the reference sample...");
    afsnps = remove_freq_diff_snps(_meta_snp_name, _meta_remain_snp, _snp_name_map, _mu, _meta_snp_freq, _snp_val_flag, npheno, freq_thresh, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(afsnps, _meta_snp_name_map, _meta_remain_snp);
    }
    // Check missing allele frequency
    LOGGER.i(0, "Update allele frequencies for the GWAS summary data ...");
    afsnps = update_snp_freq(_meta_snp_name, _meta_remain_snp, _snp_name_map, _mu, _meta_snp_freq, _snp_val_flag, npheno, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(afsnps, _meta_snp_name_map, _meta_remain_snp);
    }

    // Re-estimate sample sizes
    eigenVector n_o(npheno);
    for(i=0; i<npheno; i++) {
       eigenVector snp_n = _meta_snp_n_o.col(i);
       n_o(i) = init_pc_meta(_meta_remain_snp, _meta_snp_freq.col(i), _meta_snp_b.col(i), _meta_snp_se.col(i), snp_n, _meta_vp_trait(i));
       _meta_snp_n_o.col(i) = snp_n;
    }
    
    // bxy (PC -> phenotype)
    eigenVector bxy_hat(npheno-1);
    for(i=1; i<npheno; i++) {
        eigenVector bzx_tmp = _meta_snp_b.col(i);
        bxy_hat(i-1) = est_bxy_pc(_meta_remain_snp, _meta_snp_b.col(0), _meta_snp_freq.col(0), bzx_tmp, _meta_snp_freq.col(i), _meta_snp_n_o.col(i), n_o(i), _ttl_mk_num, _eigen_value[i-1]);
        _meta_snp_b.col(i) = bzx_tmp;
        LOGGER << "PC" << i << ", bxy = " << bxy_hat(i-1) << endl;
    }

    // Restrict effects to SNPs in the reference sample
    eigenVector bzy, bzy_se, bzy_pval, bzy_n;
    restrict_snp_effect(_meta_snp_name_map, _snp_name, _include, _meta_snp_freq, _meta_snp_b, _meta_snp_se, _meta_snp_n_o, 
                        bzy, bzy_se, bzy_n, npheno-1);

    // Estimate delta
    int m = _include.size();
    eigenVector bzy_adj(m);
    eigenVector bzx_hat(m);
    bzy_adj.setZero(m);
    wind_size = wind_size*1000;
    
    adjust_snp_effect_for_pc(bzy_adj, bzx_hat, bzy, bxy_hat, wind_size);
    
    output_snp_effect_for_pc(_out, _meta_snp_name, _meta_remain_snp, _snp_name_map, _snp_name, _meta_snp_a1, _meta_snp_a2, 
                             _meta_snp_freq.col(0), bzy, bzy_se, _meta_snp_pval.col(0), bzy_n, bzy_adj, bzx_hat);
}
