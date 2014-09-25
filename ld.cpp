/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for estimating the LD structure
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::read_LD_target_SNPs(string snplistfile)
{
    // Read snplist file
    _ld_target_snp.clear();
    read_snplist(snplistfile, _ld_target_snp);

    int i = 0, prev_target_snp_num = _ld_target_snp.size(), prev_include_size = _include.size();
    map<string, int> snp_map_buf;
    map<string, int>::iterator iter;
    for (i = 0; i < _snp_num; i++) snp_map_buf.insert(pair<string, int>(_snp_name[i], i));
    for (i = _ld_target_snp.size() - 1; i >= 0; i--) {
        iter = snp_map_buf.find(_ld_target_snp[i]);
        if (iter != snp_map_buf.end()) _snp_name_map.insert(*iter);
        else _ld_target_snp.erase(_ld_target_snp.begin() + i);
    }
    _include.clear();
    for (iter = _snp_name_map.begin(); iter != _snp_name_map.end(); iter++) _include.push_back(iter->second);
    stable_sort(_include.begin(), _include.end());
    if (_ld_target_snp.size() == 0) throw ("Error: no target SNPs are retained to estimate the LD structure.");
    else cout << prev_target_snp_num << " target SNPs read from [" + snplistfile + "], " << _ld_target_snp.size() << " of which exist in the data." << endl;
}

void gcta::LD_Blocks(int stp, double wind_size, double alpha, bool IncldQ, bool save_ram)
{
    int i = 0, j = 0;
    if (_mu.empty()) calcu_mu();

    // Read snplist file
    vector<int> smpl_buf, smpl;
    vector<string> uni_snp;
    for (i = 0; i < _include.size(); i++) uni_snp.push_back(_snp_name[_include[i]]);
    StrFunc::match(_ld_target_snp, uni_snp, smpl_buf);
    for (i = 0; i < smpl_buf.size(); i++) {
        if (smpl_buf[i]>-1) smpl.push_back(smpl_buf[i]);
    }
    int SNP_SmplNum = smpl.size(); // smpl is the position of _include

    // Calculate LD structure
    cout << "Estimating LD structure..." << endl;
    vector<int> K(SNP_SmplNum);
    vector<double> r2(SNP_SmplNum), md_r2(SNP_SmplNum), max_r2(SNP_SmplNum), dL(SNP_SmplNum), dR(SNP_SmplNum);
    vector<string> max_r2_snp(SNP_SmplNum);
    vector< vector<double> > r(SNP_SmplNum);
    vector<string> L_SNP(SNP_SmplNum), R_SNP(SNP_SmplNum);
    vector< vector<string> > snp_ls(SNP_SmplNum);
    EstLD(smpl, wind_size, snp_ls, r, r2, md_r2, max_r2, max_r2_snp, dL, dR, K, L_SNP, R_SNP, alpha, IncldQ);

    // Save result
    string SavFileName = _out + ".rsq.ld";
    ofstream SavFile(SavFileName.c_str());
    if (!SavFile) throw ("Error: can not open the file [" + SavFileName + "] to save result!");
    SavFile << "target_SNP\tfreq\tL_region\tR_region\tL_snp\tR_snp\tnSNPs\tmean_rsq\tmedian_rsq\tmax_rsq\tmax_rsq_snp" << endl;
    for (i = 0; i < SNP_SmplNum; i++) SavFile << _snp_name[_include[smpl[i]]] << "\t" << 0.5 * _mu[_include[smpl[i]]] << "\t" << dL[i] << "\t" << dR[i] << "\t" << L_SNP[i] << "\t" << R_SNP[i] << "\t" << K[i] << "\t" << r2[i] << "\t" << md_r2[i] << "\t" << max_r2[i] << "\t" << max_r2_snp[i] << endl;
    SavFile.close();
    SavFileName = _out + ".r.ld";
    SavFile.open(SavFileName.c_str());
    if (!SavFile) throw ("Error: can not open the file [" + SavFileName + "] to save result.");
    for (i = 0; i < SNP_SmplNum; i++) {
        for (j = 0; j < r[i].size(); j++) SavFile << r[i][j] << " ";
        SavFile << endl;
    }
    SavFile.close();
    SavFileName = _out + ".snp.ld";
    SavFile.open(SavFileName.c_str());
    if (!SavFile) throw ("Can not open the file [" + SavFileName + "] to save result.");
    for (i = 0; i < SNP_SmplNum; i++) {
        for (j = 0; j < snp_ls[i].size(); j++) SavFile << snp_ls[i][j] << " ";
        SavFile << endl;
    }
    SavFile.close();
    cout << "Results have been saved in [" + _out + ".rsq.ld]" + ", [" + _out + ".r.ld]" + " and [" + _out + ".snp.ld].\n" << endl;
}

void gcta::EstLD(vector<int> &smpl, double wind_size, vector< vector<string> > &snp, vector< vector<double> > &r, vector<double> &r2, vector<double> &md_r2, vector<double> &max_r2, vector<string> &max_r2_snp, vector<double> &dL, vector<double> &dR, vector<int> &K, vector<string> &L_SNP, vector<string> &R_SNP, double alpha, bool IncldQ)
{
    int i = 0, j = 0, L = 0, R = 0, maxL = 0, maxR = 0, i_buf = 0;
    map<int, int> smpl_snp_map;
    for (i = 0; i < smpl.size(); i++) smpl_snp_map.insert(pair<int, int>(_include[smpl[i]], i));

    cout << "Parameters used to search SNPs in LD with the given SNPs: window size=" << (int) (wind_size * 0.001) << "Kb, significant level=" << alpha << endl;
    vector<double> rst, y, x;
    for (i = 0; i < smpl.size(); i++) {
        vector<int> buf;
        vector<double> r_buf, rsq_buf;
        maxL = maxR = L = R = smpl[i];
        makex(L, y);
        if (IncldQ) {
            buf.push_back(L);
            rsq_buf.push_back(1.0);
            r_buf.push_back(1.0);
        }
        while (1) {
            if (R == _include.size() - 1) break;
            if (_chr[_include[R]] != _chr[_include[R + 1]]) break;
            if (_bp[_include[R + 1]] - _bp[_include[smpl[i]]] > wind_size) break;
            R++;
            if (smpl_snp_map.find(_include[R]) != smpl_snp_map.end()) continue;
            makex(R, x);
            reg(y, x, rst);
            if (rst[2] < alpha) {
                maxR = R;
                buf.push_back(R);
                rsq_buf.push_back(rst[3]);
                r_buf.push_back(rst[4]);
            }
        }
        while (1) {
            if (L == 0) break;
            if (_chr[_include[L]] != _chr[_include[L - 1]]) break;
            if (_bp[_include[smpl[i]]] - _bp[_include[L - 1]] > wind_size) break;
            L--;
            if (smpl_snp_map.find(_include[L]) != smpl_snp_map.end()) continue;
            makex(L, x);
            reg(y, x, rst);
            if (rst[2] < alpha) {
                maxL = L;
                buf.insert(buf.begin(), L);
                rsq_buf.insert(rsq_buf.begin(), rst[3]);
                r_buf.insert(r_buf.begin(), rst[4]);
            }
        }
        if (buf.size() == 0) {
            K[i] = 0;
            dL[i] = 0;
            dR[i] = 0;
            L_SNP[i] = "NA";
            R_SNP[i] = "NA";
            r[i].push_back(0.0);
            snp[i].push_back("NA");
            r2[i] = 0.0;
            md_r2[i] = 0.0;
            max_r2[i] = 0.0;
            max_r2_snp[i] = "NA";
        } else {
            K[i] = buf.size();
            dL[i] = _bp[_include[smpl[i]]] - _bp[_include[maxL]];
            dR[i] = _bp[_include[maxR]] - _bp[_include[smpl[i]]];
            L_SNP[i] = _snp_name[_include[maxL]];
            R_SNP[i] = _snp_name[_include[maxR]];
            for (j = 0; j < K[i]; j++) {
                r[i].push_back(r_buf[j]);
                snp[i].push_back(_snp_name[_include[buf[j]]]);
            }
            r2[i] = CommFunc::mean(rsq_buf);
            md_r2[i] = CommFunc::median(rsq_buf);
            i_buf = max_element(rsq_buf.begin(), rsq_buf.end()) - rsq_buf.begin();
            max_r2[i] = rsq_buf[i_buf];
            max_r2_snp[i] = snp[i][i_buf];
        }
        cout << i + 1 << " of " << smpl.size() << " target SNPs.\r";
    }
}

eigenMatrix gcta::reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table) {
    int N = x.size();
    if (N != y.size() || N < 1) throw ("Error: The lengths of x and y do not match.");

    int i = 0;
    double d_buf = 0.0, y_mu = 0.0, x_mu = 0.0, x_var = 0.0, y_var = 0.0, cov = 0.0;
    for (i = 0; i < N; i++) {
        x_mu += x[i];
        y_mu += y[i];
    }
    x_mu /= (double) N;
    y_mu /= (double) N;
    for (i = 0; i < N; i++) {
        d_buf = (x[i] - x_mu);
        x_var += d_buf*d_buf;
        d_buf = (y[i] - y_mu);
        y_var += d_buf*d_buf;
    }
    x_var /= (double) (N - 1.0);
    y_var /= (double) (N - 1.0);
    for (i = 0; i < N; i++) cov += (x[i] - x_mu)*(y[i] - y_mu);
    cov /= (double) (N - 1);
    double a = 0.0, b = 0.0, sse = 0.0, a_se = 0.0, b_se = 0.0, p = 0.0, rsq = 0.0, r = 0.0;
    if (x_var > 0.0) b = cov / x_var;
    a = y_mu - b*x_mu;
    for (i = 0; i < N; i++) {
        d_buf = y[i] - a - b * x[i];
        sse += d_buf*d_buf;
    }
    if (x_var > 0.0) {
        a_se = sqrt((sse / (N - 2.0))*(1.0 / N + x_mu * x_mu / (x_var * (N - 1.0))));
        b_se = sqrt(sse / x_var / (N - 1.0) / (N - 2.0));
    }
    if (x_var > 0.0 && y_var > 0.0) {
        r = cov / sqrt(y_var * x_var);
        rsq = r*r;
    }
    double t = 0.0;
    if (b_se > 0.0) t = fabs(b / b_se);
    p = StatFunc::t_prob(N - 2.0, t, true);
    rst.clear();
    rst.push_back(b);
    rst.push_back(b_se);
    rst.push_back(p);
    rst.push_back(rsq);
    rst.push_back(r);

    eigenMatrix reg_sum(3, 3);
    if (table) {
        reg_sum(2, 0) = rsq;
        reg_sum(1, 0) = b;
        reg_sum(1, 1) = b_se;
        reg_sum(1, 2) = p;
        if (a_se > 0.0) t = fabs(a / a_se);
        p = StatFunc::t_prob(N - 2.0, t, true);
        reg_sum(0, 0) = a;
        reg_sum(0, 1) = a_se;
        reg_sum(0, 2) = p;
        return (reg_sum);
    }
    return (reg_sum);
}

void gcta::rm_cor_snp(int m, int start, float *rsq, double rsq_cutoff, vector<int> &rm_snp_ID1)
{
    int i = 0, j = 0, i_buf = 0;

    rm_snp_ID1.clear();
    vector<int> rm_snp_ID2;
    for (i = 0; i < m; i++) {
        for (j = 0; j < i; j++) {
            if (rsq[i * m + j] > rsq_cutoff) {
                rm_snp_ID1.push_back(i + start);
                rm_snp_ID2.push_back(j + start);
            }
        }
    }
    vector<int> rm_uni_ID(rm_snp_ID1);
    rm_uni_ID.insert(rm_uni_ID.end(), rm_snp_ID2.begin(), rm_snp_ID2.end());
    stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    map<int, int> rm_uni_ID_count;
    for (i = 0; i < rm_uni_ID.size(); i++) {
        i_buf = count(rm_snp_ID1.begin(), rm_snp_ID1.end(), rm_uni_ID[i]) + count(rm_snp_ID2.begin(), rm_snp_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
    }
    map<int, int>::iterator iter1, iter2;
    #pragma omp parallel for
    for (i = 0; i < rm_snp_ID1.size(); i++) {
        iter1 = rm_uni_ID_count.find(rm_snp_ID1[i]);
        iter2 = rm_uni_ID_count.find(rm_snp_ID2[i]);
        if (iter1->second < iter2->second) {
            i_buf = rm_snp_ID1[i];
            rm_snp_ID1[i] = rm_snp_ID2[i];
            rm_snp_ID2[i] = i_buf;
        }
    }
    stable_sort(rm_snp_ID1.begin(), rm_snp_ID1.end());
    rm_snp_ID1.erase(unique(rm_snp_ID1.begin(), rm_snp_ID1.end()), rm_snp_ID1.end());
}

void gcta::calcu_mean_rsq(int wind_size, double rsq_cutoff, bool dominance_flag)
{
    check_autosome();

    int i = 0, m = _include.size();

    cout << "Calculating mean LD rsq between SNPs (block size of " << wind_size / 1000 << "Kb with an overlap of "<<wind_size/2000<<"Kb between blocks); LD rsq threshold = " << rsq_cutoff << ") ... " << endl;
    if(dominance_flag) cout<<"(SNP genotypes are coded for dominance effects)"<<endl;
    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_size);

    eigenVector mean_rsq = eigenVector::Zero(m), snp_num = eigenVector::Zero(m), max_rsq = eigenVector::Zero(m);
    calcu_ld_blk(brk_pnt1, brk_pnt3, mean_rsq, snp_num, max_rsq, false, rsq_cutoff, dominance_flag);
    if (brk_pnt2.size() > 1) calcu_ld_blk(brk_pnt2, brk_pnt3, mean_rsq, snp_num, max_rsq, true, rsq_cutoff, dominance_flag);

    string mrsq_file = "";
    if(dominance_flag) mrsq_file = _out + ".d.mrsq.ld";
    else mrsq_file = _out + ".mrsq.ld";
    ofstream o_mrsq(mrsq_file.data());
    o_mrsq<<"SNP chr bp freq mean_rsq snp_num max_rsq"<<endl;
    for (i = 0; i < m; i++) o_mrsq << _snp_name[_include[i]] << " " << _chr[_include[i]] << " " << _bp[_include[i]] << " " << 0.5 * _mu[_include[i]] << " " << mean_rsq[i] << " " << snp_num[i] << " " << max_rsq[i] << endl;
    o_mrsq << endl;
    cout << "Mean and maximum LD rsq for " << m << " SNPs have been saved in the file [" + mrsq_file + "]." << endl;
}

void gcta::calcu_ssx_sqrt_i(eigenVector &ssx_sqrt_i)
{
    int i = 0, m = _include.size();
    ssx_sqrt_i.resize(m);
    for (i = 0; i < m; i++){
        ssx_sqrt_i[i] = _geno.col(i).squaredNorm();
        if (ssx_sqrt_i[i] < 1.0e-50) ssx_sqrt_i[i] = 0.0;
        else ssx_sqrt_i[i] = 1.0 / sqrt(ssx_sqrt_i[i]);
    }
}

void gcta::get_ld_blk_pnt(vector<int> &brk_pnt1, vector<int> &brk_pnt2, vector<int> &brk_pnt3, int wind_bp, int wind_snp)
{
    unsigned long i = 0, j = 0, k = 0, m = _include.size();

    brk_pnt1.clear();
    brk_pnt1.push_back(0);
    bool chr_start = true;
    for (i = 1, j = 0; i < m; i++) {
        if (i == (m - 1)){
            if(chr_start 
                || ((_bp[_include[i]] - _bp[_include[brk_pnt1[j]]] > 0.5 * wind_bp)
                && (i - brk_pnt1[j] > 0.5 * wind_snp))) brk_pnt1.push_back(m - 1);
            else brk_pnt1[j - 1] = brk_pnt1[j] = m - 1;
        }
        else if (_chr[_include[i]] != _chr[_include[brk_pnt1[j]]] || _bp[_include[i]] - _bp[_include[i-1]] > 1e6) {
            if(chr_start 
                || ((_bp[_include[i-1]] - _bp[_include[brk_pnt1[j]]] > 0.5 * wind_bp)
                && (i - 1 - brk_pnt1[j] > 0.5 * wind_snp))){                
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
        else if ((_bp[_include[i]] - _bp[_include[brk_pnt1[j]]] > wind_bp) && (i - brk_pnt1[j] >= wind_snp)) {
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
    brk_pnt3.clear();
    for (i = 1; i < brk_pnt1.size() && brk_pnt1.size() > 2; i++) {
        if ((_chr[_include[brk_pnt1[i - 1]]] == _chr[_include[brk_pnt1[i]]]) && (brk_pnt1[i] - brk_pnt1[i - 1] > 1)) {
            int i_buf = (brk_pnt1[i - 1] + brk_pnt1[i]) / 2;
            brk_pnt2.push_back(i_buf);
            brk_pnt2.push_back(i_buf + 1);
            brk_pnt3.push_back(brk_pnt1[i]);
            brk_pnt3.push_back(brk_pnt1[i]);
        }
    }
}

void gcta::calcu_ld_blk(vector<int> &brk_pnt, vector<int> &brk_pnt3, eigenVector &mean_rsq, eigenVector &snp_num, eigenVector &max_rsq, bool second, double rsq_cutoff, bool dominance_flag)
{
    int i = 0, j = 0, k = 0, s1 = 0, s2 = 0, n = _keep.size(), m = _include.size(), size = 0, size_limit = 10000;

    for (i = 0; i < brk_pnt.size() - 1; i++) {
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        // debug
        cout << "size = " << size <<endl;

        if (second) {
            s1 = brk_pnt3[i] - brk_pnt[i];
            s2 = s1 + 1;
        }
        else {
            s1 = 0;
            s2 = size - 1;
        }

        eigenVector rsq_size(size), mean_rsq_sub(size), max_rsq_sub = eigenVector::Constant(size, -1.0);

        // make genotype matrix
        vector<int> snp_indx(size);
        for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++) snp_indx[k] = j;
        MatrixXf X_sub;
        if(dominance_flag) make_XMat_d_subset(X_sub, snp_indx, true);
        else make_XMat_subset(X_sub, snp_indx, true);
        eigenVector ssx_sqrt_i_sub(size);
        for (j = 0; j < size; j++){
            ssx_sqrt_i_sub[j] = X_sub.col(j).squaredNorm();
            if (ssx_sqrt_i_sub[j] < 1.0e-50) ssx_sqrt_i_sub[j] = 0.0;
            else ssx_sqrt_i_sub[j] = 1.0 / sqrt(ssx_sqrt_i_sub[j]);
        }

        if (size > size_limit) calcu_ld_blk_split(size, size_limit, X_sub, ssx_sqrt_i_sub, rsq_cutoff, rsq_size, mean_rsq_sub, max_rsq_sub, s1, s2, second);
        else {
            MatrixXf rsq_sub = X_sub.transpose() * X_sub;
            #pragma omp parallel for private(k)
            for (j = 0; j < size; j++) {
                rsq_size[j] = 0.0;
                mean_rsq_sub[j] = 0.0;
                for (k = 0; k < size; k++) {
                    if (second) {
                        if (j <= s1 && k <= s1) continue;
                        if (j >= s2 && k >= s2) continue;
                    }
                    if (k == j) continue;
                    rsq_sub(j,k) *= (ssx_sqrt_i_sub[j] * ssx_sqrt_i_sub[k]);
                    rsq_sub(j,k) = rsq_sub(j,k) * rsq_sub(j,k);
                    if (rsq_sub(j,k) >= rsq_cutoff) {
                        mean_rsq_sub[j] += rsq_sub(j,k);
                        rsq_size[j] += 1.0;
                    }
                    if (rsq_sub(j,k) > max_rsq_sub[j]) max_rsq_sub[j] = rsq_sub(j,k);
                }
                if (rsq_size[j] > 0.0) mean_rsq_sub[j] /= rsq_size[j];
            }
        }

        for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
            if (second) {
                if (rsq_size[j] > 0.0) {
                    mean_rsq[k] = (mean_rsq[k] * snp_num[k] + mean_rsq_sub[j] * rsq_size[j]) / (snp_num[k] + rsq_size[j]);
                    snp_num[k] = (snp_num[k] + rsq_size[j]);
                    if(max_rsq[k] < max_rsq_sub[j]) max_rsq[k] = max_rsq_sub[j];
                }
            }
            else {
                mean_rsq[k] = mean_rsq_sub[j];
                snp_num[k] = rsq_size[j];
                max_rsq[k] = max_rsq_sub[j];
            }
        }
    }
}

void gcta::calcu_ld_blk_split(int size, int size_limit, MatrixXf &X_sub, eigenVector &ssx_sqrt_i_sub, double rsq_cutoff, eigenVector &rsq_size, eigenVector &mean_rsq_sub, eigenVector &max_rsq_sub, int s1, int s2, bool second)
{
    int i = 0, j = 0, k = 0, m = 0, n = _keep.size();
    vector<int> brk_pnt_sub;
    brk_pnt_sub.push_back(0);
    for (i = size_limit; i < size - size_limit; i += size_limit) {
        brk_pnt_sub.push_back(i - 1);
        brk_pnt_sub.push_back(i);
        j = i;
    }
    j = (size - j) / 2 + j;
    brk_pnt_sub.push_back(j - 1);
    brk_pnt_sub.push_back(j);
    brk_pnt_sub.push_back(size - 1);

    for (i = 0; i < brk_pnt_sub.size() - 1; i++) {
        int size_sub = brk_pnt_sub[i + 1] - brk_pnt_sub[i] + 1;
        if (size_sub < 3) continue;

        eigenVector ssx_sqrt_i_sub_sub = ssx_sqrt_i_sub.segment(brk_pnt_sub[i], size_sub);
        MatrixXf rsq_sub_sub = X_sub.block(0,brk_pnt_sub[i],n,size_sub).transpose() * X_sub;
        eigenVector rsq_size_sub(size_sub), mean_rsq_sub_sub(size_sub), max_rsq_sub_sub = eigenVector::Constant(size_sub, -1.0);

        #pragma omp parallel for private(k)
        for (j = 0; j < size_sub; j++) {
            unsigned long s = j + brk_pnt_sub[i];
            rsq_size_sub[j] = 0.0;
            mean_rsq_sub_sub[j] = 0.0;
            for (k = 0; k < size; k++) {
                if (second) {
                    if (s <= s1 && k <= s1) continue;
                    if (s >= s2 && k >= s2) continue;
                }
                if (k == s) continue;
                rsq_sub_sub(j,k) *= (ssx_sqrt_i_sub_sub[j] * ssx_sqrt_i_sub[k]);
                rsq_sub_sub(j,k) = rsq_sub_sub(j,k) * rsq_sub_sub(j,k);
                if (rsq_sub_sub(j,k) >= rsq_cutoff) {
                    mean_rsq_sub_sub[j] += rsq_sub_sub(j,k);
                    rsq_size_sub[j] += 1.0;
                }
                if(rsq_sub_sub(j,k) > max_rsq_sub_sub[j]) max_rsq_sub_sub[j] = rsq_sub_sub(j,k);
            }

            if (rsq_size_sub[j] > 0.0) mean_rsq_sub_sub[j] /= rsq_size_sub[j];
        }

        for (j = 0, k = brk_pnt_sub[i]; j < size_sub; j++, k++) {
            mean_rsq_sub[k] = mean_rsq_sub_sub[j];
            rsq_size[k] = rsq_size_sub[j];
            max_rsq_sub[k] = max_rsq_sub_sub[j];
        }
    }
}

void gcta::ld_seg(string i_ld_file, int seg_size, int wind_size, double rsq_cutoff, bool dominance_flag)
{
    int i = 0, j = 0, k = 0, m = 0;
    vector<float> mrsq, snp_num, max_rsq;
    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;

    if(!i_ld_file.empty()){
        ifstream ild(i_ld_file.c_str());
        if (!ild) throw ("Error: can not open the file [" + i_ld_file + "] to read.");

        vector<string> snp_name;
        vector<int> chr, bp;
        vector<double> freq;
        string snp_name_buf, str_buf;
        int i_buf = 0;
        double fbuf = 0.0;
        cout << "Reading LD mean rsq for SNPs from [" + i_ld_file + "] ..." << endl;
        getline(ild, str_buf); // get the header
        while (ild) {
            ild >> snp_name_buf;
            if (ild.eof()) break;
            snp_name.push_back(snp_name_buf);
            if(!(ild >> i_buf)) throw("Error: in the file [" + i_ld_file + "].");
            if (i_buf < 1 || i_buf > 50) throw ("Error: invalid value of \"chr\" for the SNP " + snp_name_buf + ".");
            chr.push_back(i_buf);
            if(!(ild >> i_buf)) throw("Error: in the file [" + i_ld_file + "].");
            if (i_buf < 0) throw ("Error: invalid value of \"bp\" for the SNP " + snp_name_buf + ".");
            bp.push_back(i_buf);
            if(!(ild >> fbuf)) throw("Error: in the file [" + i_ld_file + "].");
            if (fbuf < 0.0 || fbuf > 1.0) throw ("Error: invalid value of \"allele frequency\" for the SNP " + snp_name_buf + ".");
            freq.push_back(fbuf*2.0);
            if(!(ild >> fbuf)) throw("Error: in the file [" + i_ld_file + "].");
            if (fbuf < -1.0 || fbuf > 1.0) throw ("Error: invalid value of \"mean LD rsq\" for the SNP " + snp_name_buf + ".");
            mrsq.push_back(fbuf);
            if(!(ild >> fbuf)) throw("Error: in the file [" + i_ld_file + "].");
            if (fbuf < 0) throw ("Error: invalid value of \"number of SNPs in LD\" for the SNP " + snp_name_buf + ".");
            snp_num.push_back(fbuf);
            if(!(ild >> fbuf)) throw("Error: in the file [" + i_ld_file + "].");
            if (fbuf < 0 || fbuf > 2.0) throw ("Error: invalid value of \"max LD rsq\" for the SNP " + snp_name_buf + ".");
            max_rsq.push_back(fbuf);
            m++;
            getline(ild, str_buf);
        }
        ild.close();
        cout << "LD mean rsq for " << m << " SNPs read from [" + i_ld_file + "]." << endl;

        _snp_num = m;
        _snp_name = snp_name;
        _chr = chr;
        _bp = bp;
        _mu = freq;
        init_include();
    }
    else {
        if(_snp_name.size() < 1) throw("Error: need to input the LD file or PLINK files for LD calculation.");
        check_autosome();
        cout << "Calculating mean LD rsq between SNPs (block size of " << wind_size / 1000 << "Kb with an overlap of "<<wind_size / 2000<<"Kb between blocks); LD rsq threshold = " << rsq_cutoff << ") ... " << endl;
        if(dominance_flag) cout<<"(SNP genotypes are coded for dominance effects)"<<endl;
        get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_size);
        eigenVector mrsq_buf = eigenVector::Zero(m), snp_num_buf = eigenVector::Zero(m), max_rsq_buf = eigenVector::Zero(m);
        calcu_ld_blk(brk_pnt1, brk_pnt3, mrsq_buf, snp_num_buf, max_rsq_buf, false, rsq_cutoff, dominance_flag);
        if (brk_pnt2.size() > 1) calcu_ld_blk(brk_pnt2, brk_pnt3, mrsq_buf, snp_num_buf, max_rsq_buf, true, rsq_cutoff, dominance_flag);
        mrsq.resize(m);
        snp_num.resize(m);
        max_rsq.resize(m);
        for(i = 0; i < m; i++){
            mrsq[i] = mrsq_buf(i);
            snp_num[i] = snp_num_buf(i);
            max_rsq[i] = max_rsq_buf(i);
        }
    }

    get_lds_brkpnt(brk_pnt1, brk_pnt2, seg_size);
    int size = 0, mean_size = 0, count = 0;
    for (i = 0; i < brk_pnt1.size() - 1; i++) {
        size = brk_pnt1[i + 1] - brk_pnt1[i] + 1;
        if (size > 2){
            mean_size += size; 
            count++;
        }
    }
    mean_size /= count;
    get_lds_brkpnt(brk_pnt1, brk_pnt2, 0, mean_size);

    vector<double> lds(m);
    for (i = 0; i < brk_pnt1.size() - 1; i += 2){
        size = brk_pnt1[i + 1] - brk_pnt1[i] + 1;
        double ld_score = 0.0;
        for(j = brk_pnt1[i]; j <= brk_pnt1[i + 1]; j++) ld_score += mrsq[j] * snp_num[j] + 1.0; 
        ld_score /= (double)size;
        for(j = brk_pnt1[i]; j <= brk_pnt1[i + 1]; j++) lds[j] = ld_score; 
    }

    string lds_file;
    if(dominance_flag) lds_file = _out + ".d.mrsq.ld";
    else lds_file = _out + ".mrsq.ld";
    ofstream o_lds(lds_file.data());
    o_lds<<"SNP chr bp freq mean_rsq snp_num max_rsq mean_lds"<<endl;
    for (i = 0; i < m; i++) o_lds << _snp_name[_include[i]] << " " << _chr[_include[i]] << " " << _bp[_include[i]] << " " << 0.5 * _mu[_include[i]] << " " << mrsq[i] << " " << snp_num[i] << " " << max_rsq[i] << " " << lds[i] << endl;
    o_lds << endl;
    cout << "Mean and maximum LD rsq for " << m << " SNPs have been saved in the file [" + lds_file + "]." << endl;
}

// calculate maximum LD rsq between SNPs
void gcta::calcu_max_ld_rsq(int wind_size, double rsq_cutoff, bool dominance_flag)
{
    int i = 0, m = _include.size(), max_size = 0.9 * _keep.size();

    cout << "Calculating maximum LD rsq between SNPs (block size of " << wind_size / 1000 << "Kb with an overlap of "<<wind_size/2000<<"Kb between blocks; LD rsq threshold = " << rsq_cutoff << ") ... " << endl;
    cout << "(Maximum number of SNPs allowed in a block = " << max_size << " due to computatinal limitation)" << endl;
    if(dominance_flag) cout<<"(SNP genotypes are coded for dominance effects)"<<endl;

    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt_max_limit(brk_pnt1, brk_pnt2, brk_pnt3, wind_size, max_size);

    eigenVector multi_rsq = eigenVector::Constant(m, -1.0), max_rsq = eigenVector::Constant(m, -999), multi_rsq_adj = eigenVector::Constant(m, -999);
    vector<int> max_pos(m);
    calcu_max_ld_rsq_blk(multi_rsq, multi_rsq_adj, max_rsq, max_pos, brk_pnt1, rsq_cutoff, dominance_flag);
    if (brk_pnt2.size() > 1) calcu_max_ld_rsq_blk(multi_rsq, multi_rsq_adj, max_rsq, max_pos, brk_pnt2, rsq_cutoff, dominance_flag);

    string max_rsq_file = "";
    if(dominance_flag) max_rsq_file = _out + ".d.max_rsq.ld";
    else max_rsq_file = _out + ".max_rsq.ld";
    ofstream o_max_rsq(max_rsq_file.data());
    o_max_rsq<<"SNP freq max_rsq max_snp multi_rsq multi_rsq_adj"<<endl;
    for (i = 0; i < m; i++){
        if(max_rsq[i] > 0.0 || multi_rsq[i] > -998) o_max_rsq << _snp_name[_include[i]] << " " << 0.5 * _mu[_include[i]] << " " << max_rsq[i] << " " << _snp_name[_include[max_pos[i]]] << " " << multi_rsq[i] << " " <<multi_rsq_adj[i] << endl;
        else o_max_rsq << _snp_name[_include[i]] << " " << 0.5 * _mu[_include[i]] << " NA NA NA NA" << endl;
    }
    o_max_rsq << endl;
    cout << "Maximum LD rsq for " << m << " SNPs have been saved in the file [" + max_rsq_file + "]." << endl;    
}

void gcta::get_ld_blk_pnt_max_limit(vector<int> &brk_pnt1, vector<int> &brk_pnt2, vector<int> &brk_pnt3, int wind_bp, int wind_snp)
{
    unsigned long i = 0, j = 0, k = 0, m = _include.size();

    brk_pnt1.clear();
    brk_pnt1.push_back(0);
    for (i = 1, j = 0; i < m; i++) {
        if (i == (m - 1)) brk_pnt1.push_back(m - 1);
        else if (_chr[_include[i]] != _chr[_include[brk_pnt1[j]]] || _bp[_include[i]] - _bp[_include[i-1]] > 1e6) {
            brk_pnt1.push_back(i - 1);
            j++;
            brk_pnt1.push_back(i);
            j++;
        }
        else if ((_bp[_include[i]] - _bp[_include[brk_pnt1[j]]] > wind_bp) || (i - brk_pnt1[j] >= wind_snp)) {
            brk_pnt1.push_back(i - 1);
            j++;
            brk_pnt1.push_back(i);
            j++;
        }
    }
    stable_sort(brk_pnt1.begin(), brk_pnt1.end());
    brk_pnt1.erase(unique(brk_pnt1.begin(), brk_pnt1.end()), brk_pnt1.end());

    brk_pnt2.clear();
    brk_pnt3.clear();
    for (i = 1; i < brk_pnt1.size() && brk_pnt1.size() > 2; i++) {
        if ((_chr[_include[brk_pnt1[i - 1]]] == _chr[_include[brk_pnt1[i]]]) && (brk_pnt1[i] - brk_pnt1[i - 1] > 1)) {
            int i_buf = (brk_pnt1[i - 1] + brk_pnt1[i]) / 2;
            brk_pnt2.push_back(i_buf);
            brk_pnt2.push_back(i_buf + 1);
            brk_pnt3.push_back(brk_pnt1[i]);
            brk_pnt3.push_back(brk_pnt1[i]);
        }
    }
}

void gcta::calcu_max_ld_rsq_blk(eigenVector &multi_rsq, eigenVector &multi_rsq_adj, eigenVector &max_rsq, vector<int> &max_pos, vector<int> &brk_pnt, double rsq_cutoff, bool dominance_flag)
{
	int i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size(), size = 0;

    for (i = 0; i < brk_pnt.size() - 1; i++)
    {
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        // make genotype matrix
        vector<int> snp_indx(size);
        for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++) snp_indx[k] = j;
        MatrixXf X_sub;
        if(dominance_flag) make_XMat_d_subset(X_sub, snp_indx, true);
        else make_XMat_subset(X_sub, snp_indx, true);
        eigenVector ssx_sqrt_i_sub(size);
        for (j = 0; j < size; j++){
            ssx_sqrt_i_sub[j] = X_sub.col(j).squaredNorm();
            if (ssx_sqrt_i_sub[j] < 1.0e-50) ssx_sqrt_i_sub[j] = 0.0;
            else ssx_sqrt_i_sub[j] = 1.0 / sqrt(ssx_sqrt_i_sub[j]);
        }

        MatrixXf rsq_sub = X_sub.transpose() * X_sub;
        #pragma omp parallel for private(k)
        for (j = 0; j < size; j++) {
            rsq_sub(j,j) = 1.0;
            for (k = j + 1; k < size; k++){
                rsq_sub(j,k) *= (ssx_sqrt_i_sub[j] * ssx_sqrt_i_sub[k]);
                rsq_sub(k,j) = rsq_sub(j,k);
            }
        }


        SelfAdjointEigenSolver<MatrixXf> pca(rsq_sub.array());

                // debug
       // ofstream tmp("tmp_R.txt");
        //tmp<< rsq_sub.col(1998) << endl;


        VectorXf d_i = pca.eigenvalues();
        double eff_m = 0;
        for(j = 0; j < size; j ++){
        	if(d_i(j) < 1e-5) d_i(j) = 0.0;
        	else{
        		d_i(j) = 1.0 / d_i(j);
        		eff_m++;
        	}
        }

        // debug
        //ofstream tmp2("tmp_eval.txt");
        //tmp2<<pca.eigenvalues()<< endl;



        // debug
        cout << "size = " << size << "; eff_m = " << eff_m << endl;
        MatrixXf R_i = pca.eigenvectors() * d_i.asDiagonal() * pca.eigenvectors().transpose();
        R_i = R_i.array();


        //JacobiSVD<MatrixXf> svd;
        /*svd.compute(rsq_sub, ComputeThinV);
        VectorXf d_i = svd.singularValues();
        double eff_m = 0;
        for(j = 0; j < size; j ++){
        	if(d_i(j) < 1e-6) d_i(j) = 0.0;
        	else{
        		d_i(j) = 1.0 / d_i(j);
        		eff_m++;
        	}
        }

        // debug
        cout << "size = " << size << "; eff_m = " << eff_m << endl;

        MatrixXf R_i = svd.matrixV() * d_i.asDiagonal() * svd.matrixV().transpose();*/
        VectorXf Q_diag(size);
        for(j = 0; j < size; j ++) Q_diag(j) = R_i.col(j).dot(rsq_sub.row(j).transpose());
        VectorXf multi_rsq_buf(size);
    	for(j = 0; j < size; j ++){
    		if(fabs(Q_diag[j] - 1.0) < 0.01) multi_rsq_buf[j] = 1.0 - 1.0 / R_i(j,j);
    		else multi_rsq_buf[j] = 1.0;
            if(multi_rsq_buf[j] > 1.0) multi_rsq_buf[j] = 1.0;
    	}
        VectorXf multi_rsq_buf_adj = multi_rsq_buf.array() - (1.0 - multi_rsq_buf.array()) * (eff_m / ((double)n - eff_m - 1.0));

        rsq_sub.diagonal() = VectorXf::Zero(size);
        rsq_sub = rsq_sub.array() * rsq_sub.array();
        VectorXf max_rsq_buf(size);
        vector<int> max_pos_buf(size);
        VectorXf::Index max_pos_pnt;
        for(j = 0; j < size; j++){
            rsq_sub.col(j).maxCoeff(&max_pos_pnt);
            max_pos_buf[j] = (int)max_pos_pnt;
            max_rsq_buf[j] = rsq_sub.col(j)[max_pos_buf[j]];
            if(multi_rsq_buf[j] < max_rsq_buf[j]) multi_rsq_buf[j] = max_rsq_buf[j];
        }

/*      for(j = 0; j < size; j++){
            if(multi_rsq_buf[j] > 1.0) multi_rsq_buf[j] = 1.0;
            if(max_rsq_buf[j] > 1.0) max_rsq_buf[j] = 1.0;
        }
        */
        for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
            if(multi_rsq_adj[k] <= multi_rsq_buf_adj[j]){
            	multi_rsq[k] = multi_rsq_buf[j];
            	multi_rsq_adj[k] = multi_rsq_buf_adj[j];
            }
            if(max_rsq[k] < max_rsq_buf[j]){
                max_rsq[k] = max_rsq_buf[j];
                max_pos[k] = max_pos_buf[j] + brk_pnt[i];
            }
        }
    }
}

bool gcta::bending_eigenval_Xf(VectorXf &eval)
{
    int j = 0;
    double eval_m = eval.mean();
    if (eval.minCoeff() > 0.0) return false;
    double S = 0.0, P = 0.0;
    for (j = 0; j < eval.size(); j++) {
        if (eval[j] >= 0) continue;
        S += eval[j];
        P = -eval[j];
    }
    double W = S * S * 100.0 + 1;
    for (j = 0; j < eval.size(); j++) {
        if (eval[j] >= 0) continue;
        eval[j] = P * (S - eval[j])*(S - eval[j]) / W;
    }
    eval *= eval_m / eval.mean();
    return true;
}

