//
//  grm_wt.cpp
//  gcta
//
//  Created by Jian Yang on 14/03/13.
//
//

#include "gcta.h"

// LD smoothing approach
void gcta::calcu_lds(eigenVector &wt, int wind_size)
{
    unsigned long i = 0, n = _keep.size(), m = _include.size();

    cout << "Calculating SNP weights using LD smoothing approach (block size of " << wind_size / 1000 << "Kb with an overlap of "<<wind_size/2000<<"Kb between blocks) ..." << endl;
    eigenVector ssx_sqrt_i;
    calcu_ssx_sqrt_i(ssx_sqrt_i);

    vector<int> brk_pnt1, brk_pnt2;
    get_lds_brkpnt(brk_pnt1, brk_pnt2, wind_size);
    int mean_size = 0, count = 0;
    for (i = 0; i < brk_pnt1.size() - 1; i++) {
        int size = brk_pnt1[i + 1] - brk_pnt1[i] + 1;
        if (size > 2){
            mean_size += size; 
            count++;
        }
    }
    mean_size /= count;
    get_lds_brkpnt(brk_pnt1, brk_pnt2, 0, mean_size);

    eigenVector m_maf = eigenVector::Zero(m);
    calcu_maf();
    calcu_lds_blk(wt, m_maf, ssx_sqrt_i, brk_pnt1, false);
    if (brk_pnt2.size() > 1) calcu_lds_blk(wt, m_maf, ssx_sqrt_i, brk_pnt2, true);

     // debug
    double wt_m = wt.mean();
    cout<<"wt mean = " << wt_m << endl;
    cout<<"wt variance = " << (wt - eigenVector::Constant(m, wt_m)).squaredNorm() / (m - 1.0) <<endl;
    cout<<"wt range = "<<wt.minCoeff()<<" ~ "<<wt.maxCoeff()<<endl;
    cout << "wt: " << wt.segment(0,10).transpose() << endl;

    // debug
    string wt_file = _out + ".ldwt";
    ofstream owt(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << " " << m_maf[i] << endl;
    owt << endl;
    cout<<"LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();

   // adjust wt for log(maf)
    eigenMatrix X(m, 3);
    X.col(0) = eigenVector::Ones(m);
    X.col(1) = m_maf;
    eigenVector e = wt - X * ((X.transpose() * X).inverse() * (X.transpose() * wt));
    wt = e.array() + wt.mean();

    // debug
    cout << "wt adj: " << wt.segment(0,10).transpose() << endl;

    // debug
    wt_file = _out + ".adj.ldwt";
    owt.open(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << " " << m_maf[i] << endl;
    owt << endl;
    cout<<"Adjusted LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();

}

void gcta::get_lds_brkpnt(vector<int> &brk_pnt1, vector<int> &brk_pnt2, int wind_size, int wind_snp_num)
{
    unsigned long i = 0, j = 0, k = 0, m = _include.size();

    brk_pnt1.clear();
    brk_pnt1.push_back(0);
    bool chr_start = true;
    for (i = 1, j = 0; i < m; i++) {
        if (i == (m - 1)) brk_pnt1[j - 1] = brk_pnt1[j] = m - 1;
        else if (_chr[_include[i]] != _chr[_include[brk_pnt1[j]]]) {
            if(chr_start){
                brk_pnt1.push_back(i - 1);
                j++;
                brk_pnt1.push_back(i);
                j++;
            }
            else{                
                brk_pnt1[j - 1] = i - 1;
                brk_pnt1[j] = i;
            }
            chr_start = true;
        }
        else if ((_bp[_include[i]] - _bp[_include[brk_pnt1[j]]] > wind_size) && (i - brk_pnt1[j] > wind_snp_num)) {
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

void gcta::calcu_lds_blk(eigenVector &wt, eigenVector &m_maf, eigenVector &ssx_sqrt_i, vector<int> &brk_pnt, bool second)
{
    unsigned long i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size(), size = 0;

    for (i = 0; i < brk_pnt.size() - 1; i++) {
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        MatrixXf rsq_sub(size, size);
        rsq_sub = _geno.block(0,brk_pnt[i],n,size).transpose()*_geno.block(0,brk_pnt[i],n,size);
        eigenVector ssx_sqrt_i_sub_sqrt = ssx_sqrt_i.segment(brk_pnt[i],size);
        #pragma omp parallel for private(k)
        for (j = 0; j < size; j++) {
            rsq_sub(j,j) = 1.0;
            for (k = j + 1; k < size; k++){
                rsq_sub(j,k) *= (ssx_sqrt_i_sub_sqrt[j] * ssx_sqrt_i_sub_sqrt[k]);
                rsq_sub(k,j) = rsq_sub(j,k);
            }
        }

        SelfAdjointEigenSolver<MatrixXf> eigensolver(rsq_sub);       
        VectorXf eval = eigensolver.eigenvalues();
        double sum = eval.sum();
        double eff_m = (sum*sum) / eval.squaredNorm();
        double wt_buf = eff_m / (double)size;
        // debug
        /*MatrixXf grm = _geno.block(0,brk_pnt[i],n,size) * _geno.block(0,brk_pnt[i],n,size).transpose();
        grm = grm.array() / double (size);
        double off_num = 0.5*n*(n - 1.0), off_m = 0.0, off_v = 0.0;
        for (j = 1; j < n; j++) off_m += grm.row(j).segment(0, j).sum();
        off_m /= off_num;
        for (j = 1; j < n; j++) off_v += (grm.row(j).segment(0, j) -  VectorXf::Constant(j, off_m).transpose()).squaredNorm();
        off_v /= (off_num - 1.0);
        double eff_m = 1.0 / off_v;
        double wt_buf = eff_m / (double)size;*/

        // debug
        cout<<"size = "<<size<<"; eff_m = "<<eff_m <<endl;

        double m_maf_buf = 0.0;
        for (k = brk_pnt[i]; k <= brk_pnt[i+1]; k++) m_maf_buf += _maf[k];
        m_maf_buf /= (double) size;

        for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
            if (second){
                wt[k] = 0.5*(wt[k] + wt_buf);
                m_maf[k] = 0.5*(m_maf[k] + m_maf_buf);
            }
            else{
                wt[k] = wt_buf;
                m_maf[k] = m_maf_buf;
            }
        }
    }
}

// Weighting GRM using the Speed et al. 2013 AJHG method
void gcta::calcu_ldak(eigenVector &wt, int wind_size, double rsq_cutoff)
{
    unsigned long i = 0, n = _keep.size(), m = _include.size();

    cout << "Calculating the LD based SNP weights (block size of " << wind_size / 1000 << "Kb with an overlap of "<<wind_size/2000<<"Kb between windows, and a least 3000 SNPs within a window) ..." << endl;
    eigenVector ssx_sqrt_i;
    calcu_ssx_sqrt_i(ssx_sqrt_i);

    vector<int> brk_pnt1, brk_pnt2;
    get_lds_brkpnt(brk_pnt1, brk_pnt2, wind_size, 3000);
    calcu_ldak_blk(wt, ssx_sqrt_i, brk_pnt1, false, rsq_cutoff);
    if (brk_pnt2.size() > 1) calcu_ldak_blk(wt, ssx_sqrt_i, brk_pnt2, true, rsq_cutoff);
    if(_maf.size() < 1) calcu_maf();
    
    // adjust wt for maf
    eigenMatrix X(m, 3);
    X.col(0) = eigenVector::Ones(m);
    for(i = 0; i < m; i++){
        X(i,1) = log(_maf[i]);
        X(i,2) = X(i,1) * X(i,1);
    }
    eigenVector e = wt - X * ((X.transpose() * X).inverse() * (X.transpose() * wt));
    wt = e.array() + wt.mean();

    // debug
    double wt_m = wt.mean();
    cout<<"wt mean = " << wt_m << endl;
    cout<<"wt variance = " << (wt - eigenVector::Constant(m, wt_m)).squaredNorm() / (m - 1.0) <<endl;
    cout<<"wt range = "<<wt.minCoeff()<<" ~ "<<wt.maxCoeff()<<endl;
    cout << "wt: " << wt.segment(0,10).transpose() << endl;

    // debug
    string wt_file = _out + ".ldwt";
    ofstream owt(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << " " << _maf[i] << endl;
    owt << endl;
    cout<<"LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    owt.close();

    // debug
    /*wt_file = _out + ".adj.ldwt";
    owt.open(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << " " << _maf[i] << endl;
    owt << endl;
    cout<<"Adjusted LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
    */
}

void gcta::calcu_ldak_blk(eigenVector &wt, eigenVector &ssx_sqrt_i, vector<int> &brk_pnt, bool second, double rsq_cutoff)
{
    unsigned long i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size(), size = 0;
    double rsq_adj = 1.0 / (double)n;

    for (i = 0; i < brk_pnt.size() - 1; i++) {
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        MatrixXf rsq_sub(size, size);
        rsq_sub = _geno.block(0,brk_pnt[i],n,size).transpose()*_geno.block(0,brk_pnt[i],n,size);
        eigenVector ssx_sqrt_i_sub_sqrt = ssx_sqrt_i.segment(brk_pnt[i],size);        
        #pragma omp parallel for private(k)
        for (j = 0; j < size; j++) {
            rsq_sub(j,j) = 1.001;
            for (k = j + 1; k < size; k++) {
                rsq_sub(j,k) *= (ssx_sqrt_i_sub_sqrt[j] * ssx_sqrt_i_sub_sqrt[k]);
                rsq_sub(j,k) = rsq_sub(j,k)*rsq_sub(j,k) - rsq_adj;
                if (rsq_sub(j,k) <= rsq_cutoff) rsq_sub(j,k) = 0.0;
                rsq_sub(k,j) = rsq_sub(j,k);
            }
        }
        VectorXf c = VectorXf::Ones(size);
        VectorXf wt_sub = rsq_sub.lu().solve(c);
        if(!(wt_sub.minCoeff() > -1e10 && wt_sub.minCoeff() < 1e10)) throw("Error: LU decomposition error!");
        for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
            if (second) wt[k] = 0.5*(wt[k] + wt_sub(j));
            else wt[k] = wt_sub(j);
        }
    }
}

void gcta::make_grm_pca(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, double wind_size, bool mlmassoc)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();

    if(grm_d_flag) throw("Error: --domiance not supported in this analysis.");
    if (grm_xchr_flag) throw("Error: --make-grm-xchr not supported in this analysis.");
    if (grm_mtd > 0) throw("Error: --make-grm-alg not supported in this analysis.");
    if (mlmassoc) throw("Error: --mlma not supported in this analysis.");
    check_autosome();

    if (!mlmassoc) cout << "\nCalculating the PCA-based" << ((grm_d_flag) ? " dominance" : "") << " genetic relationship matrix (GRM)" << (grm_xchr_flag ? " for the X chromosome" : "") << (_dosage_flag ? " using imputed dosage data" : "") << " ... (Note: default speed-optimized mode, may use huge RAM)" << endl;
    else cout << "\nCalculating the PCA-based genetic relationship matrix (GRM) ... " << endl;
    cout << "(block size of " << wind_size / 1000 << " SNPs)" << endl;

    vector<float> seq;
    vector< vector<int> > maf_bin_pos; 
    assign_snp_2_mb(seq, maf_bin_pos, 100);
    double trace = 0.0;
    for (i = 0; i < seq.size() - 1; i++) {
        if (maf_bin_pos[i].size() > 0){
            
            // debug
            cout<<seq[i] << " < MAF <= "<<seq[i + 1]<<endl;

            make_grm_pca_blk(maf_bin_pos[i], wind_size / 1000, trace);
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

void gcta::make_grm_pca_blk(vector<int> & maf_bin_pos_i, int wind_size, double &trace)
{
    int i = 0, j = 0, k = 0, n = _keep.size();

    int length = maf_bin_pos_i.size() / (1 + (maf_bin_pos_i.size() / wind_size));
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

void gcta::assign_snp_2_mb(vector<float> &seq, vector< vector<int> > &maf_bin_pos, int mb_num)
{
    int i = 0, j = 0, m = _include.size();

     // create frequency bin
    seq.clear();
    double x = 0.5;
    seq.push_back(0.501);
    for (i = 0; i < mb_num; i++) {
        x = exp(log(x) - 0.1);
        seq.insert(seq.begin(), x);
    }
    int seq_size = seq.size();

    // assign SNPs to MAF bins
    if(_maf.size() < 1) calcu_maf();
    maf_bin_pos.clear();
    maf_bin_pos.resize(seq_size - 1);
    for (j = 1; j < seq_size; j++) {
        for (i = 0; i < m; i++) {
            if (_maf[i] >= seq[j - 1] && _maf[i] < seq[j]) maf_bin_pos[j - 1].push_back(i);
        }
    }
}

void gcta::calcu_ldwt(string i_ld_file, eigenVector &wt, int wind_size, double rsq_cutoff)
{
    unsigned long i = 0, j = 0, k = 0, l = 0, n = _keep.size(), m = _include.size();
    double rsq_thre = 3.92 / (double)n;
    if(rsq_thre < rsq_cutoff) rsq_thre = rsq_cutoff;

    // assign SNPs to MAF bins
    vector<float> seq;
    vector< vector<int> > maf_bin_pos;
    assign_snp_2_mb(seq, maf_bin_pos, 100);
    int seq_size = seq.size();

    vector<double> mrsq_mb;
    wt = eigenVector::Zero(m);
    eigenVector snp_num = eigenVector::Zero(m);
    eigenVector max_rsq = eigenVector::Zero(m);
    // Read mean LD from file and calculate the mean LD in each MAF bin
    if (!i_ld_file.empty()) read_mrsq_mb(i_ld_file, seq, mrsq_mb, wt, snp_num);
    else {
        // calcualte LD between SNPs
        cout << "Calculating mean LD rsq (block size = " << wind_size / 1000 << "Kb; LD rsq threshold = " << rsq_thre << ") ... " << endl;
        eigenVector ssx_sqrt_i;
        calcu_ssx_sqrt_i(ssx_sqrt_i);
        vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
        get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_size);
        calcu_ld_blk(ssx_sqrt_i, brk_pnt1, brk_pnt3, wt, snp_num, max_rsq, false, rsq_thre);
        if (brk_pnt2.size() > 1) calcu_ld_blk(ssx_sqrt_i, brk_pnt2, brk_pnt3, wt, max_rsq, snp_num, true, rsq_thre);

        // adjust LD due to chance correlation
        double n_r = 1.0 / n;
        #pragma omp parallel for
        for (i = 0; i < m; i++){
            wt[i] = wt[i] * snp_num[i] + 1.0;
            snp_num[i]++;
            wt[i] /= snp_num[i]; 
        }
        
        // MAF bin weighting
        cout << "Calculating mean LD rsq in MAF bins ... " << endl;
        mrsq_mb.resize(seq_size - 1);
        for (i = 0; i < seq_size - 1; i++) {
            if (maf_bin_pos[i].size() > 0) {
                long double sum_snp_num = 0.0;
                for (j = 0; j < maf_bin_pos[i].size(); j++) {
                    k = maf_bin_pos[i][j];
                    mrsq_mb[i] += wt[k] * snp_num[k];
                    sum_snp_num += snp_num[k];
                }
                mrsq_mb[i] /= sum_snp_num;
                cout << maf_bin_pos[i].size() << " SNPs with " << setprecision(4) << seq[i] << " <= MAF < " << seq[i + 1] << ", mean LD rsq = " << mrsq_mb[i] << endl;
            }
        }
    }

    for (i = 0; i < seq_size - 1; i++) {
        if (maf_bin_pos[i].size() > 0 && mrsq_mb[i] > 0.0) {
            #pragma omp parallel for private(k)
            for (j = 0; j < maf_bin_pos[i].size(); j++) {
                k = maf_bin_pos[i][j];
                wt[k] = mrsq_mb[i] / wt[k];
            }
        }
    }

    // debug
    double wt_m = wt.mean();
    cout<<"wt mean = " << wt_m << endl;
    cout<<"wt variance = " << (wt - eigenVector::Constant(m, wt_m)).squaredNorm() / (m - 1.0) <<endl;
    cout<<"wt range = "<<wt.minCoeff()<<" ~ "<<wt.maxCoeff()<<endl;
    cout << "wt: " << wt.segment(0,10).transpose() << endl;

    // debug
    string wt_file = _out + ".ldwt";
    ofstream owt(wt_file.data());
    if(!owt) throw("Error: can not open [" + wt_file + "] to read.");
    for (i = 0; i < m; i++)  owt << _snp_name[_include[i]] << " " << wt[i] << endl;
    owt << endl;
    cout<<"LD weights for all SNPs have bene saved in [" + wt_file + "]."<<endl;
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

    // adjusting LD due to chance correlation
    #pragma omp parallel for
    for (i = 0; i < m; i++) {
        mrsq[i] = mrsq[i] * snp_num[i] + 1.0;
        snp_num[i]++;
        mrsq[i] /= snp_num[i];
    }

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

    // assign SNPs to MAF bins
    #pragma omp parallel for
    for (i = 0; i < m; i++) {
        if (freq[i] > 0.5) freq[i] = 1.0 - freq[i];
    }
    int seq_size = seq.size();
    vector< vector<int> > maf_bin_pos(seq_size - 1);
    for (j = 1; j < seq_size; j++) {
        for (i = 0; i < m; i++) {
            if (freq[i] >= seq[j - 1] && freq[i] < seq[j]) maf_bin_pos[j - 1].push_back(i);
        }
    }

    // calculate mean LD rsq in MAF bins
    cout << "Calculating mean LD rsq in MAF bins ... " << endl;
    mrsq_mb.clear();
    mrsq_mb.resize(seq_size - 1);
    for (i = 0; i < seq_size - 1; i++) {
        if (maf_bin_pos[i].size() > 0) {
            long double sum_snp_num = 0.0;
            for (j = 0; j < maf_bin_pos[i].size(); j++) {
                k = maf_bin_pos[i][j];
                mrsq_mb[i] += mrsq[k] * snp_num[k];
                sum_snp_num += snp_num[k];
            }
            mrsq_mb[i] /= sum_snp_num;
            cout << maf_bin_pos[i].size() << " SNPs with " << setprecision(4) << seq[i] << " <= MAF < " << seq[i + 1] << ", mean LD rsq = " << mrsq_mb[i] << endl;
        }
    }
}

void gcta::adj_wt_4_maf(eigenVector &wt)
{
    int i = 0, m = _include.size();
    eigenVector log_maf(m);
    for(i = 0; i < m; i++){
        log_maf[i] = _mu[_include[i]] * 0.5;
        if(log_maf[i] > 0.5) log_maf[i] = 1.0 - log_maf[i];
        log_maf[i] = log(log_maf[i]);
    }

    eigenVector y = wt;
    y = y.array() - y.mean();
    log_maf = log_maf.array() - log_maf.mean();
    double beta = y.dot(log_maf) / log_maf.dot(log_maf);
    wt = wt.array() - log_maf.array()*beta;

    // debug
    cout << "beta = " << beta;
}


