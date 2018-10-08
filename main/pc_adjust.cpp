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
cout<<"_ttl_mk_num " << _ttl_mk_num<<endl;
cout<<"_meta_remain_snp.size()1 "<<_meta_remain_snp.size()<<endl;
    // Target trait
    if(pheno_file[0].substr(pheno_file[0].length()-3,3)!=".gz")
        snplist=read_snp_metafile_txt(pheno_file[0], gws_snp_name_map, -9);
    else 
        snplist=read_snp_metafile_gz(pheno_file[0], gws_snp_name_map, -9);
    update_id_map_kp(snplist, _snp_name_map, _include);
cout<<"_meta_remain_snp.size()2 "<<_meta_remain_snp.size()<<endl;
    // Initialization of variables
    int nsnp = _snp_name_map.size();
    eigenMatrix snp_freq;
    vector<vector<string>> snp_a1, snp_a2;
cout<<"_meta_remain_snp.size()3 "<<_meta_remain_snp.size()<<endl;
    init_gwas_variable(snp_a1, snp_a2, snp_freq, _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, ncovar+1, nsnp); 
cout<<"_meta_remain_snp.size()4 "<<_meta_remain_snp.size()<<endl;
    _snp_val_flag.clear(); _snp_val_flag.resize(ncovar+1);
    for(i=0; i<ncovar+1; i++) {
        _snp_val_flag[i].resize(nsnp);
        for(j=0; j<nsnp; j++) _snp_val_flag[i][j] = false;
    }
cout<<"_meta_remain_snp.size()5 "<<_meta_remain_snp.size()<<endl;
    LOGGER.i(0, to_string(_include.size()) + " SNPs in common between the summary data and the LD reference sample.");
    // Read the meta analysis
    eigenVector snp_freq_buf(nsnp), snp_b_buf(nsnp), snp_se_buf(nsnp), snp_pval_buf(nsnp), snp_n_buf(nsnp);
cout<<"_meta_remain_snp.size()6 "<<_meta_remain_snp.size()<<endl;
    _meta_vp_trait.resize(ncovar+1);
cout<<"_meta_remain_snp.size()7 "<<_meta_remain_snp.size()<<endl;
    // Target trait and Covariates
    for(i=0; i<=ncovar; i++) {
        if(pheno_file[i].substr(pheno_file[i].length()-3,3)!=".gz")
            _meta_vp_trait(i) = read_single_metafile_txt(pheno_file[i], _snp_name_map, snp_a1[i], snp_a2[i], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i]);
        else _meta_vp_trait(i) = read_single_metafile_gz(pheno_file[i], _snp_name_map, snp_a1[i], snp_a2[i], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i]);
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
cout<<"_meta_remain_snp.size()8 "<<_meta_remain_snp.size()<<endl;     
    // QC of SNPs
    LOGGER.i(0, "Filtering out SNPs with multiple alleles or missing value ...");
    
    vector<string> badsnps;
    badsnps = remove_bad_snps(_snp_name, _include, _snp_val_flag, snp_a1, snp_a2, snp_freq,  _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, 
                             _snp_name_map, _allele1, _allele2, 1, ncovar, _out);
    if(badsnps.size()>0) {
        update_id_map_rm(badsnps, _snp_name_map, _include);
    }
cout<<"_meta_remain_snp.size()9 "<<_meta_remain_snp.size()<<endl;
    // For output
    _meta_snp_a1 = snp_a1[0]; _meta_snp_a2 = snp_a2[0];
    _meta_snp_freq = snp_freq;  
cout<<"_meta_remain_snp.size()10 "<<_meta_remain_snp.size()<<endl;    
    nsnp = _include.size();    
    if(nsnp<1) LOGGER.e(0, "None SNPs are retained after filtering.");
    else LOGGER.i(0, to_string(nsnp) + " SNPs are retained after filtering." + to_string(_include.size()) + " SNPs are in common with LD reference sample.");
}

double init_pc_meta(eigenVector msx, eigenVector &snp_b, eigenVector snp_se, eigenVector &snp_n_o, double vp_trait) {
    
    int i = 0, m = snp_b.size();

    // Update sample size per SNP
    snp_n_o = (vp_trait*eigenVector::Ones(m) - (msx.array()*snp_b.array()*snp_b.array()).matrix()).array() / (msx.array()*snp_se.array()*snp_se.array()) + 1;
    
    // Estimate total sample size
    eigenVector snp_n_sort(snp_n_o);
    std::stable_sort(snp_n_sort.data(), snp_n_sort.data()+m);
    double n_o = CommFunc::quantile(snp_n_sort, 0.50);
    
    // Standardise SNP effects
    snp_b = msx.array().sqrt()*snp_b.array();

    return n_o;
}

double est_bxy_pc(eigenVector bzy, eigenVector bzx, eigenVector bzx_n, double n_o, int ttl_mk_num, double eigen_value)
{
    double bxy = 0.0;
    bxy = n_o*bzy.dot((bzx.array()*bzx_n.array()).matrix())/((double)ttl_mk_num*eigen_value);

    return bxy;
}

void init_ld_snp_index(int &snp_start, int &snp_end, vector<int> snp_slct_in, vector<int> snp_remain, vector<int> snp_chr, vector<int> snp_bp, int wind_size) {
    int snp_index = snp_start;

    while(1) {
        snp_index++;
        snp_end = snp_remain[snp_index];
        if(snp_chr[snp_end] == snp_chr[snp_start] && abs(snp_bp[snp_end] - snp_bp[snp_start]) < wind_size) 
            snp_slct_in.push_back(snp_end);
        else break;
    }
}

double calcu_ztz_product_b(eigenVector bzx, eigenVector bzx_n, eigenVector ztz_vec, double ztz_target, int snp_target, int snp_wind_start, int snp_wind_end, vector<int> snp_chr, vector<int> snp_bp, double lambda, int wind_size, double nsnps) {
    int i = 0, snp_test_start = -1, snp_test_end = -1;
    double d_buf = 0.0;

    for(i=snp_wind_start; i<=snp_wind_end; i++) {
        if(abs(snp_bp[snp_target] - snp_bp[i]) >= wind_size) continue;
        if( snp_test_start < 0 ) snp_test_start = i;
        if( i > snp_test_end ) snp_test_end = i;
    }

    for(i=snp_test_start; i<=snp_test_end; i++) {
        d_buf += ztz_vec(i-snp_wind_start)*bzx(i)*bzx_n(i);
    }
    d_buf /= ztz_target*lambda*(double)nsnps;

    return d_buf;
}

void gcta::adjust_snp_effect_for_pc(eigenVector bzy_delta, eigenVector msx, int wind_size) {    
    int snp_index = 0, m = _include.size(), n = _keep.size(), ncovar = _meta_snp_b.size() - 1;

    int snp_start = 0, snp_end = 0, snp_buf1_start = 0, snp_buf1_end = 0, snp_buf2_start = 0, snp_buf2_end = 0; 
    vector<int> snp_slct_in, snp_in_buf1, snp_in_buf2;
    snp_start = _include[0];
    init_ld_snp_index(snp_start, snp_end, snp_slct_in, _include, _chr, _bp, wind_size);
    snp_buf2_start = snp_end + 1;
    init_ld_snp_index(snp_buf2_start, snp_buf2_end, snp_in_buf2, _include, _chr, _bp, wind_size);

    // adjust SNP effects
    while(1) {
        int m_buf1 = snp_buf1_end - snp_buf1_start + 1, m_in = snp_end - snp_start + 1, m_buf2 = snp_buf2_end - snp_buf2_start + 1;
        eigenVector z_buf(n);
        eigenMatrix ztz_buf1, ztz_in, ztz_buf2;
        if(m_in == 0) break;
        // re-code
        int i = 0, j = 0;
        eigenVector z_i;
        eigenMatrix x_sub;
        for(i=snp_start; i<=snp_end; i++) {
            makex_eigenVector(_include[i], z_i, false, true);
            x_sub.col(i-snp_start) = z_i;
        }
        // buf 1
        if(m_buf1>1) {
            ztz_buf1.resize(m_buf1, m_in);
            for( i = snp_buf1_start; i <= snp_buf1_end; i++ ) {
                makex_eigenVector(_include[i], z_buf, false, true);
                for( j = snp_start; j <= snp_end; j++ ) {
                    ztz_buf1(i-snp_buf1_start,j-snp_start) = z_buf.dot(x_sub.col(j-snp_start))*min(_meta_snp_n_o(_include[i]), _meta_snp_n_o(_include[j])) * sqrt(msx(_include[i]) * msx(_include[j]) / (msx(_include[i]) * msx(_include[j])));
                }
            }
        }
        // within window
        ztz_in.resize(m_in, m_in);
        for( i = snp_start; i <= snp_end; i++ ) {
            for( j = snp_start; j <= snp_end; j++ ) {
                ztz_in(i-snp_start,j-snp_start) = x_sub.col(i-snp_start).dot(x_sub.col(j-snp_start))*min(_meta_snp_n_o(_include[i]), _meta_snp_n_o(_include[j])) * sqrt(msx(_include[i]) * msx(_include[j]) / (msx(_include[i]) * msx(_include[j])));
            }
        }
        // buf 2
        if(m_buf2>1) {
            ztz_buf2.resize(m_in, m_buf2);
            for( i = snp_start; i <= snp_end; i++ ) {
                makex_eigenVector(_include[i], z_i, false, true);
                for( j = snp_buf2_start; j <= snp_buf2_end; j++ ) {
                    ztz_buf2(i-snp_start,j-snp_buf2_start) = x_sub.col(i-snp_start).dot(z_buf)*min(_meta_snp_n_o(_include[i]), _meta_snp_n_o(_include[j])) * sqrt(msx(_include[i]) * msx(_include[j]) / (msx(_include[i]) * msx(_include[j])));
                }
            }            
        }
        // ztz at the target SNP
        double ztz_target = _meta_snp_n_o(_include[snp_start])*msx(_include[snp_start]);
        // adjust SNP effect
        for( i = snp_start; i <= snp_end; i++ ) {
            for( j = 0; j < ncovar; j++) {
                // check the distance
                double dbuf1 = 0.0, dbuf2 = 0.0, dbuf3 = 0.0;
                dbuf1 = calcu_ztz_product_b(_meta_snp_b.col(j+1), _meta_snp_n_o.col(j+1), ztz_buf1.row(i), ztz_target, i, snp_buf1_start, snp_buf1_end, _chr, _bp, _eigen_value[j], wind_size, _ttl_mk_num);
                dbuf2 = calcu_ztz_product_b(_meta_snp_b.col(j+1), _meta_snp_n_o.col(j+1), ztz_in.col(i), ztz_target, i, snp_buf1_start, snp_buf1_end, _chr, _bp, _eigen_value[j], wind_size, _ttl_mk_num);
                dbuf3 = calcu_ztz_product_b(_meta_snp_b.col(j+1), _meta_snp_n_o.col(j+1), ztz_buf2.col(i), ztz_target, i, snp_buf1_start, snp_buf1_end, _chr, _bp, _eigen_value[j], wind_size, _ttl_mk_num);
                bzy_delta(i) += dbuf1 + dbuf2 + dbuf3;
            }
        } 
    }
}

void gcta::pc_adjust(string pcadjust_list_file, string pc_file, double freq_thresh, int wind_size) {

    // Read the summary data
    read_pc_adjust_file(pcadjust_list_file, pc_file);

    int i = 0, j = 0, m = _include.size(), npheno = _meta_snp_b.size();

    // Calculate allele frequency
    if (_mu.empty()) calcu_mu();
cout<<"_snp_num "<<_snp_num<<endl;    
    eigenVector msx(m);
    for(i=0; i<m; i++)
        msx(i) = _mu[_include[i]]*(1.0 - _mu[_include[i]]/2);
ofstream odata("snp_flag.txt");
for(i=0; i<_snp_val_flag[0].size(); i++)
    odata <<_snp_val_flag[0][i]<<" " <<_snp_val_flag[1][i]<<endl;
odata.close();
cout<<"output complete"<<endl;
cin.get();
    // Check allele frequency
    vector<string> afsnps;
    LOGGER.i(0, "Checking allele frequencies among the GWAS summary data and the reference sample...");
    afsnps = remove_freq_diff_snps(_snp_name, _include, _snp_name_map, _mu, _meta_snp_freq, _snp_val_flag, npheno, freq_thresh, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
    }
    
    eigenVector n_o(npheno);

    // Update summary statistics
    for(i=0; i<npheno; i++) {
       eigenVector snp_b = _meta_snp_b.col(i), snp_se = _meta_snp_se.col(i), snp_n = _meta_snp_se.col(i);
       n_o(i) = init_pc_meta(msx, snp_b, snp_se, snp_n, _meta_vp_trait(i));
    }
    
    // bxy (PC -> phenotype)
    eigenVector bxy_hat(npheno-1);
    LOGGER << "bxy" << endl;
    for(i=1; i<npheno; i++) {
        bxy_hat(i-1) = est_bxy_pc(_meta_snp_b.col(0), _meta_snp_b.col(i), _meta_snp_n_o.col(i), n_o(i), _ttl_mk_num, _eigen_value[i]);
        LOGGER << i << " " << bxy_hat(i-1) << endl;
    }

    // Estimate delta
    eigenVector bzy_delta(m);
    bzy_delta.setZero(m);
    adjust_snp_effect_for_pc(bzy_delta, msx, wind_size);
}