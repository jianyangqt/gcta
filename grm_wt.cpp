//
//  grm_wt.cpp
//  gcta
//
//  Created by Jian Yang on 14/03/13.
//
//

#include "gcta.h"

// LD smoothing approach
void gcta::calcu_lds(string i_ld_file, eigenVector &wt, int ldwt_wind, int ldwt_seg, double rsq_cutoff)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size(), size = 0;
    double rsq_thre = 0.01;
    if(rsq_cutoff > 0.0) rsq_thre = rsq_cutoff;
    
    // calculate mean LD
    cout << "Calculating mean LD rsq between SNPs within a region of at least " << ldwt_wind / 1000 << "Kb long (LD rsq threshold = " << rsq_thre << ") ... " << endl;
    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, ldwt_wind);
    eigenVector mean_rsq = eigenVector::Zero(m), snp_num = eigenVector::Zero(m), max_rsq = eigenVector::Zero(m);
    eigenVector ssx_sqrt_i;
    calcu_ssx_sqrt_i(ssx_sqrt_i);
    calcu_ld_blk(ssx_sqrt_i, brk_pnt1, brk_pnt3, mean_rsq, snp_num, max_rsq, false, rsq_thre, false);
    if (brk_pnt2.size() > 1) calcu_ld_blk(ssx_sqrt_i, brk_pnt2, brk_pnt3, mean_rsq, snp_num, max_rsq, true, rsq_thre, false);

    wt = eigenVector::Ones(m);
    get_lds_brkpnt(brk_pnt1, brk_pnt2, ldwt_seg);
    int mean_size = 0, count = 0;
    for (i = 0; i < brk_pnt1.size() - 1; i++) {
        size = brk_pnt1[i + 1] - brk_pnt1[i] + 1;
        if (size > 2){
            mean_size += size; 
            count++;
        }
    }
    mean_size /= count;
    get_lds_brkpnt(brk_pnt1, brk_pnt2, 0, mean_size);

    cout << "Segment-based LD smoothing (segment length =  " << ldwt_seg / 1000 << "Kb) ... " << endl;
    eigenVector ld_sum = (mean_rsq.array() * snp_num.array()).array() + 1.0;
    double wt_buf = 0.0;
    for (i = 0; i < brk_pnt1.size() - 1; i++) {
        size = brk_pnt1[i + 1] - brk_pnt1[i] + 1;
        if (size < 3) continue;
        wt_buf = ld_sum.segment(brk_pnt1[i],size).sum();
        if(wt_buf < 1e-6) wt_buf = (double)size;
        wt_buf = (double)size / wt_buf; 
        for (j = 0, k = brk_pnt1[i]; j < size; j++, k++) wt[k] = wt_buf;
    }
    
    for (i = 0; i < brk_pnt2.size() - 1; i++) {
        size = brk_pnt2[i + 1] - brk_pnt2[i] + 1;
        if (size < 3) continue;
        wt_buf = ld_sum.segment(brk_pnt2[i],size).sum();
        if(wt_buf < 1e-6) wt_buf = (double)size;
        wt_buf = (double)size / wt_buf; 
        for (j = 0, k = brk_pnt2[i]; j < size; j++, k++) wt[k] = 0.5*(wt[k] + wt_buf);
    }

    // debug
    double wt_m = wt.mean();
    cout<<"wt mean = " << wt_m << endl;
    cout<<"wt variance = " << (wt - eigenVector::Constant(m, wt_m)).squaredNorm() / (m - 1.0) <<endl;
    cout<<"wt range = "<<wt.minCoeff()<<" ~ "<<wt.maxCoeff()<<endl;
    cout << "wt: " << wt.segment(0,10).transpose() << endl;

    if(_maf.size() < 1) calcu_maf();

    // debug
    eigenVector wt_o(wt);

    // adjust wt for maf
    adj_wt_4_maf(wt);

    // debug
    string wt_file = _out + ".ldwt";
    ofstream owt(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << _bp[_include[i]] << " " << wt_o[i] << " " << wt[i] << " " << _maf[i] << " " << mean_rsq[i] << " " << snp_num[i] << endl;
    owt << endl;
    cout<<"Adjusted LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();
}

void gcta::get_lds_brkpnt(vector<int> &brk_pnt1, vector<int> &brk_pnt2, int ldwt_seg, int wind_snp_num)
{
    unsigned long i = 0, j = 0, k = 0, m = _include.size();

    brk_pnt1.clear();
    brk_pnt1.push_back(0);
    bool chr_start = true;
    int prev_length = 0;
    for (i = 1, j = 0; i < m; i++) {
        if (i == (m - 1)){
            if(!chr_start && prev_length < 0.5 * wind_snp_num) brk_pnt1[j - 1] = brk_pnt1[j] = m - 1;
            else brk_pnt1.push_back(m - 1);
        }
        else if (_chr[_include[i]] != _chr[_include[brk_pnt1[j]]]) {
            if(!chr_start && prev_length < 0.5 * wind_snp_num){                
                brk_pnt1[j - 1] = i - 1;
                brk_pnt1[j] = i;
            }
            else{
                brk_pnt1.push_back(i - 1);
                j++;
                brk_pnt1.push_back(i);
                j++;
            }
            chr_start = true;
        }
        else if ((_bp[_include[i]] - _bp[_include[brk_pnt1[j]]] > ldwt_seg) && (i - brk_pnt1[j] > wind_snp_num)) {
            prev_length  = i - brk_pnt1[j];
            chr_start = false;
            brk_pnt1.push_back(i - 1);
            j++;
            brk_pnt1.push_back(i);
            j++;
        }
    }
    stable_sort(brk_pnt1.begin(), brk_pnt1.end());
    brk_pnt1.erase(unique(brk_pnt1.begin(), brk_pnt1.end()), brk_pnt1.end());

    brk_pnt2.clear();
    for (i = 1; i < brk_pnt1.size() && brk_pnt1.size() > 2; i++) {
        if ((_chr[_include[brk_pnt1[i - 1]]] == _chr[_include[brk_pnt1[i]]]) && (brk_pnt1[i] - brk_pnt1[i - 1] > 1)) {
            int i_buf = (brk_pnt1[i - 1] + brk_pnt1[i]) / 2;
            brk_pnt2.push_back(i_buf);
            brk_pnt2.push_back(i_buf + 1);
        }
    }
}

// Weighting GRM using the Speed et al. 2013 AJHG method
void gcta::calcu_ldak(eigenVector &wt, int ldwt_seg, double rsq_cutoff)
{

    unsigned long i = 0, j = 0, k = 0, l = 0, n = _keep.size(), m = _include.size();
    double rsq_thre = 0.01; //3.8416 / (double)n;
    if(rsq_thre < rsq_cutoff) rsq_thre = rsq_cutoff;

    eigenVector ssx_sqrt_i;
    calcu_ssx_sqrt_i(ssx_sqrt_i);
    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, ldwt_seg);


    wt = eigenVector::Zero(m);
    eigenVector mrsq;
    cout << "Calculating the LD based SNP weights (block size of " << ldwt_seg / 1000 << "Kb with an overlap of "<<ldwt_seg/2000<<"Kb between windows, and a least 3000 SNPs within a window) ..." << endl;
    calcu_ldak_blk(wt, mrsq, ssx_sqrt_i, brk_pnt1, false, rsq_cutoff);
    if (brk_pnt2.size() > 1) calcu_ldak_blk(wt, mrsq, ssx_sqrt_i, brk_pnt2, true, rsq_cutoff);
    
     // debug
    double wt_m = wt.mean();
    cout<<"wt mean = " << wt_m << endl;
    cout<<"wt variance = " << (wt - eigenVector::Constant(m, wt_m)).squaredNorm() / (m - 1.0) <<endl;
    cout<<"wt range = "<<wt.minCoeff()<<" ~ "<<wt.maxCoeff()<<endl;
    cout << "wt: " << wt.segment(0,10).transpose() << endl;

    if(_maf.size() < 1) calcu_maf();

    // debug
    string wt_file = _out + ".ldwt";
    ofstream owt(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << " " << _maf[i] << " " << mrsq[i] << endl;
    owt << endl;
    cout<<"LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();

    // adjust wt for maf
    /*adj_wt_4_maf(wt);

    // debug
    wt_file = _out + ".adj.ldwt";
    owt.open(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << " " << _maf[i] << endl;
    owt << endl;
    cout<<"Adjusted LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();*/
}

void gcta::calcu_ldak_blk(eigenVector &wt, eigenVector &sum_rsq, eigenVector &ssx_sqrt_i, vector<int> &brk_pnt, bool second, double rsq_cutoff)
{
    unsigned long i = 0, n = _keep.size(), m = _include.size(), size = 0;
    double rsq_thre = 0.01; //3.8416 / (double)n;
    if(rsq_thre < rsq_cutoff) rsq_thre = rsq_cutoff;

    for (i = 0; i < brk_pnt.size() - 1; i++) {
        int j = 0, k = 0;
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        // debug
        cout << "size = " << size << endl;

        MatrixXf rsq_sub(size, size+1);
        rsq_sub.block(0,0,size,size) = _geno.block(0,brk_pnt[i],n,size).transpose()*_geno.block(0,brk_pnt[i],n,size);
        eigenVector ssx_sqrt_i_sub_sqrt = ssx_sqrt_i.segment(brk_pnt[i],size); 
        #pragma omp parallel for private(k)
        for (j = 0; j < size; j++) {
            rsq_sub(j,j) = 1.01;
            for (k = j + 1; k < size; k++) {
                rsq_sub(j,k) *= (ssx_sqrt_i_sub_sqrt[j] * ssx_sqrt_i_sub_sqrt[k]);
                rsq_sub(j,k) = rsq_sub(j,k)*rsq_sub(j,k);
                if (rsq_sub(j,k) <= rsq_cutoff) rsq_sub(j,k) = 0.0;
                rsq_sub(k,j) = rsq_sub(j,k);
            }
        }

        
        VectorXf c = VectorXf::Ones(size);

        FullPivLU<MatrixXf> lu(rsq_sub.block(0,0,size,size));
        VectorXf wt_sub = lu.solve(c);

        //VectorXf wt_sub = rsq_sub.lu().solve(c);
        if(!(wt_sub.minCoeff() > -1e10 && wt_sub.minCoeff() < 1e10)) throw("Error: LU decomposition error!");
        for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
            if (second) wt[k] = 0.5*(wt[k] + wt_sub(j));
            else wt[k] = wt_sub(j);
        }
    }
}

void gcta::calcu_ldwt(string i_ld_file, eigenVector &wt, int wind_size, double rsq_cutoff)
{
    int i = 0, m = _include.size();
    double rsq_thre = 0.01;
    if(rsq_cutoff > 0.0) rsq_thre = rsq_cutoff;
    
    // calculate mean LD
    cout << "Calculating mean LD rsq between SNPs within a region of at least " << wind_size / 1000 << "Kb long (LD rsq threshold = " << rsq_thre << ") ... " << endl;
    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_size);
    eigenVector mean_rsq = eigenVector::Zero(m), snp_num = eigenVector::Zero(m), max_rsq = eigenVector::Zero(m);
    eigenVector ssx_sqrt_i;
    calcu_ssx_sqrt_i(ssx_sqrt_i);
    calcu_ld_blk(ssx_sqrt_i, brk_pnt1, brk_pnt3, mean_rsq, snp_num, max_rsq, false, rsq_thre, false);
    if (brk_pnt2.size() > 1) calcu_ld_blk(ssx_sqrt_i, brk_pnt2, brk_pnt3, mean_rsq, snp_num, max_rsq, true, rsq_thre, false);

    wt = 1.0 / ((mean_rsq.array() * snp_num.array()).array() + 1.0);
    double min_coeff = wt.minCoeff();
    for(i = 0; i < m; i++){
        if(wt[i] > 0.99) wt[i] = min_coeff;
    }

    // debug
    double wt_m = wt.mean();
    cout<<"wt mean = " << wt_m << endl;
    cout<<"wt variance = " << (wt - eigenVector::Constant(m, wt_m)).squaredNorm() / (m - 1.0) <<endl;
    cout<<"wt range = "<<wt.minCoeff()<<" ~ "<<wt.maxCoeff()<<endl;
    cout << "wt: " << wt.segment(0,10).transpose() << endl;
    
    // debug
    eigenVector wt_o(wt);

    adj_wt_4_maf(wt);

    // debug
    string wt_file = _out + ".ldwt";
    ofstream owt(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << _bp[_include[i]] << " " << wt_o[i] << " " << wt[i] << " " << _maf[i] << " " << mean_rsq[i] << " " << snp_num[i] << endl;
    owt << endl;
    cout<<"Adjusted LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();
}

void gcta::read_mrsq_mb(string i_ld_file, vector<float> &seq, vector<double> &mrsq_mb, eigenVector &wt, eigenVector &snp_m)
{
    ifstream ild(i_ld_file.c_str());
    if (!ild) throw ("Error: can not open the file [" + i_ld_file + "] to read.");

    int i = 0, j = 0, k = 0;
    string snp_name_buf, str_buf;
    double fbuf = 0.0;
    cout << "Reading LD mean rsq for SNPs from [" + i_ld_file + "] ..." << endl;
    vector<string> snp_name;
    vector<float> snp_num, mrsq, freq;
    int m = 0;
    getline(ild, str_buf); // get the header
    while (ild) {
        ild >> snp_name_buf;
        if (ild.eof()) break;
        snp_name.push_back(snp_name_buf);
        ild >> str_buf;
        fbuf = atof(str_buf.c_str());
        if (fbuf < 0.0 || fbuf > 1.0) throw ("Error: invalid value of \"allele frequency\" for the SNP " + snp_name_buf + ".");
        freq.push_back(fbuf);
        ild >> str_buf;
        fbuf = atof(str_buf.c_str());
        if (fbuf > 1.0 || fbuf < 0.0) {
            throw ("Warning: invalid value of \"mean LD rsq\" for the SNP " + snp_name_buf + ", which will be set to zero.");
            fbuf = 0.0;
        }
        mrsq.push_back(fbuf);
        ild >> str_buf;
        fbuf = atof(str_buf.c_str());
        if (fbuf < 1.0) throw ("Error: invalid value of \"number of SNPs\" for the SNP " + snp_name_buf + ".");
        snp_num.push_back(fbuf);
        m++;
        getline(ild, str_buf);
    }
    ild.close();
    cout << "LD mean rsq for " << m << " SNPs read from [" + i_ld_file + "]." << endl;

    // calculate weights
    #pragma omp parallel for
    for (i = 0; i < m; i++) mrsq[i] = 1.0 / (mrsq[i] * snp_num[i] + 1.0);

    // extract the mean LD rsq for SNPs included in this run of analysis
    map<string, int> snp_map;
    for (i = 0; i < _include.size(); i++) snp_map.insert(pair<string, int>(_snp_name[_include[i]], i));
    map<string, int>::iterator iter, end = snp_map.end();
    int icount = 0;
    #pragma omp parallel for
    for (i = 0; i < m; i++) {
        iter = snp_map.find(snp_name[i]);
        if (iter != end) {
            wt[iter->second] = mrsq[i];
            snp_m[iter->second] = snp_num[i];
            icount++;
        }
    }
    cout << _include.size() << " SNPs are included in the analysis and mean LD rsq for " << icount << " SNPs are updated from the file [" << i_ld_file << "]." << endl;
}

void gcta::adj_wt_4_maf(eigenVector &wt)
{
    if(_maf.size() < 1) calcu_maf();

    vector<float> seq;
    vector< vector<int> > maf_bin_pos;
    assign_snp_2_mb(seq, maf_bin_pos, 100);
    int seq_size = seq.size();

    int i = 0, j = 0;
    eigenVector wt_mb_mean = eigenVector::Zero(seq_size - 1), wt_mb_sd = eigenVector::Zero(seq_size - 1);
    double mean = 0.0, sd = 0.0;
    int max_snp_num = 0;
    for (i = 0; i < seq_size - 1; i++) {
        int size = maf_bin_pos[i].size();
        if (size > 1) {

     /*       VectorXf y(size);
            for(j = 0; j < size; j++) y(j) = wt[maf_bin_pos[i][j]];
            wt_mb_mean[i] = y.mean();
            wt_mb_sd[i] = sqrt((y - VectorXf::Constant(size, wt_mb_mean[i])).squaredNorm() / (size - 1.0));
            eigen_func::inverse_norm_rank_transform(y);
            y = y.array()*wt_mb_sd[i];
            y = y.array() + abs(y.minCoeff());
            for(j = 0; j < size; j++) wt[maf_bin_pos[i][j]] = y(j);
*/
            
            for (j = 0; j < size; j++)  wt_mb_mean[i] += wt[maf_bin_pos[i][j]];
            wt_mb_mean[i] /= (double) size;
            for (j = 0; j < size; j++)  wt_mb_sd[i] += (wt[maf_bin_pos[i][j]] - wt_mb_mean[i]) * (wt[maf_bin_pos[i][j]] - wt_mb_mean[i]);
            wt_mb_sd[i] /= (double) (size - 1);
            wt_mb_sd[i] = sqrt(wt_mb_sd[i]);
            
            cout << size << " SNPs with " << setprecision(4) << seq[i] << " <= MAF < " << seq[i + 1] << ", mean wt = " << wt_mb_mean[i] << ", sd wt = " << wt_mb_sd[i] << endl;
            if(size > max_snp_num){
                max_snp_num = size;
                mean = wt_mb_mean[i];
                sd = wt_mb_sd[i];
            }
        }
    }

    // debug
    cout << "Max mean = " << mean << "; max sd = " << sd <<endl;

    //int m = wt.size();
    //mean = wt.mean();
    //sd = sqrt((wt - eigenVector::Constant(m, mean)).squaredNorm() / (m - 1.0));
    for (i = 0; i < seq_size - 1; i++) {
        if (maf_bin_pos[i].size() > 30){
            for (j = 0; j < maf_bin_pos[i].size(); j++) wt[maf_bin_pos[i][j]] = sd * (wt[maf_bin_pos[i][j]] - wt_mb_mean[i]) / wt_mb_sd[i] + mean;
        }
    }
}

void gcta::cal_sum_rsq_mb(eigenVector &sum_rsq_mb)
{
    vector<float> seq;
    vector< vector<int> > maf_bin_pos;
    assign_snp_2_mb(seq, maf_bin_pos, 100);
    int seq_size = seq.size();

    int i = 0, j = 0;
    for (i = 0; i < seq_size - 1; i++) {
        if (maf_bin_pos[i].size() > 0) {
            double sum_rsq_mb_buf = 0.0;
            for (j = 0; j < maf_bin_pos[i].size(); j++)  sum_rsq_mb_buf += sum_rsq_mb[maf_bin_pos[i][j]];
            sum_rsq_mb_buf /= (double) maf_bin_pos[i].size();
            for (j = 0; j < maf_bin_pos[i].size(); j++)  sum_rsq_mb[maf_bin_pos[i][j]] = sum_rsq_mb_buf;
        }
    }
}

void gcta::assign_snp_2_mb(vector<float> &seq, vector< vector<int> > &maf_bin_pos, int mb_num)
{
    int i = 0, j = 0, m = _include.size();

     // create frequency bin
    seq.clear();
    double x = 0.5;
    seq.push_back(0.501);
    for (i = 0; i < mb_num; i++) {
        x = exp(log(x) - 0.1053605);
        seq.push_back(x);
    }
    int seq_size = seq.size();

    // assign SNPs to MAF bins
    if(_maf.size() < 1) calcu_maf();
    maf_bin_pos.clear();
    maf_bin_pos.resize(seq_size - 1);
    for (j = 1; j < seq_size; j++) {
        for (i = 0; i < m; i++) {
            if (_maf[i] < seq[j - 1] && _maf[i] >= seq[j]) maf_bin_pos[j - 1].push_back(i);
        }
    }
}

/////////////
// not working

void gcta::make_grm_pca(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, double ldwt_seg, bool mlmassoc)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();

    if(grm_d_flag) throw("Error: --domiance not supported in this analysis.");
    if (grm_xchr_flag) throw("Error: --make-grm-xchr not supported in this analysis.");
    if (grm_mtd > 0) throw("Error: --make-grm-alg not supported in this analysis.");
    if (mlmassoc) throw("Error: --mlma not supported in this analysis.");
    check_autosome();

    if (!mlmassoc) cout << "\nCalculating the PCA-based" << ((grm_d_flag) ? " dominance" : "") << " genetic relationship matrix (GRM)" << (grm_xchr_flag ? " for the X chromosome" : "") << (_dosage_flag ? " using imputed dosage data" : "") << " ... (Note: default speed-optimized mode, may use huge RAM)" << endl;
    else cout << "\nCalculating the PCA-based genetic relationship matrix (GRM) ... " << endl;
    cout << "(block size of " << ldwt_seg / 1000 << " SNPs)" << endl;

    vector<float> seq;
    vector< vector<int> > maf_bin_pos; 
    assign_snp_2_mb(seq, maf_bin_pos, 100);
    double trace = 0.0;
    for (i = 0; i < seq.size() - 1; i++) {
        if (maf_bin_pos[i].size() > 0){
            
            // debug
            cout<<seq[i + 1] << " < MAF <= "<<seq[i]<<endl;

            make_grm_pca_blk(maf_bin_pos[i], ldwt_seg / 1000, trace);
        }
    }
    _grm = _grm.array() / trace;
    _grm_N = MatrixXf::Constant(n,n,trace);

    // debug
    double diag_m = 0.0, diag_v = 0.0, off_m = 0.0, off_v = 0.0;
    calcu_grm_var(diag_m, diag_v, off_m, off_v);
    cout<<"\nMean of diagonals = "<<diag_m<<endl;
    cout<<"Variance of diagonals = "<<diag_v<<endl;
    cout<<"Mean of off-diagonals = " << off_m <<endl;
    cout<<"Variance of off-diagonals = " << off_v <<endl;

    // Output A_N and A
    string out_buf = _out;
    output_grm(output_bin);
    _out = out_buf;
}

void gcta::make_grm_pca_blk(vector<int> & maf_bin_pos_i, int ldwt_seg, double &trace)
{
    int i = 0, j = 0, k = 0, n = _keep.size();

    int length = maf_bin_pos_i.size() / (1 + (maf_bin_pos_i.size() / ldwt_seg));
    vector<int> brk_pnt;
    brk_pnt.push_back(0);
    for(i = length; i < maf_bin_pos_i.size(); i+=length){
        brk_pnt.push_back(i);
        brk_pnt.push_back(i+1);
    }
    brk_pnt.push_back(maf_bin_pos_i.size()-1);

    for (i = 0; i < brk_pnt.size() - 1; i++) {
        int size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        vector<int> snp_indx(size);
        for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++) snp_indx[k] = maf_bin_pos_i[j];
        MatrixXf X;
        make_XMat_subset(X, snp_indx, true);

        VectorXf ssx_sqrt_i(size);
        for (j = 0; j < size; j++){
            ssx_sqrt_i[j] = X.col(j).squaredNorm();
            if (ssx_sqrt_i[j] < 1.0e-50) ssx_sqrt_i[j] = 0.0;
            else ssx_sqrt_i[j] = 1.0 / sqrt(ssx_sqrt_i[j]);
        }

        MatrixXf rsq = X.transpose() * X;

        #pragma omp parallel for private(k)
        for (j = 0; j < size; j++) {
            rsq(j,j) = 1.0;
            for (k = j + 1; k < size; k++){
                rsq(j,k) *= (ssx_sqrt_i[j] * ssx_sqrt_i[k]);
                rsq(k,j) = rsq(j,k);
            }
        }

        SelfAdjointEigenSolver<MatrixXf> eigensolver(rsq);
        VectorXf eval = eigensolver.eigenvalues().array() / eigensolver.eigenvalues().sum();

        double VE = 0.0;
        int eff_size = 0;
        for(j = size - 1; j >= 0; j--) {
            VE += eval[j];
            eff_size++;
            if(VE>0.99) break;
        }

        // debug
        cout<< "size = " << size <<", " << "eff_size = "<<eff_size<<endl;

        X = X * eigensolver.eigenvectors().block(0, size - eff_size, size, eff_size);

        col_std(X);
        MatrixXf grm_buf = X * X.transpose();
        grm_buf = grm_buf.array() * (double)size / (double) eff_size;

        trace += grm_buf.diagonal().mean();

        // debug
        cout<<"trace = "<<grm_buf.diagonal().mean()<<endl;

        #ifdef SINGLE_PRECISION
        _grm = _grm + grm_buf;
        #else
        _grm = _grm + grm_buf.cast<double>();
        #endif
    }
}

void gcta::col_std(MatrixXf &X)
{
    int i = 0, col_num = X.cols(), row_num = X.rows();
    double sd;
    VectorXf m(col_num);
    for(i = 0; i < col_num; i++) m(i) = X.col(i).mean();
    for(i = 0; i < col_num; i++) X.col(i) = X.col(i).array() - m(i);
    for(i = 0; i < col_num; i++){
        sd = sqrt(X.col(i).squaredNorm() / (row_num - 1.0));
        if(CommFunc::FloatEqual(sd, 0.0)) X.col(i) = VectorXf::Zero(row_num);
        else X.col(i) = X.col(i).array() / sd;
    }
}



