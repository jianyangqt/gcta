/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for estimating the genetic relationship matrix
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::enable_grm_bin_flag() {
    _grm_bin_flag = true;
}

void gcta::check_autosome() {
    for (int i = 0; i < _include.size(); i++) {
        if (_chr[_include[i]] > _autosome_num) throw ("Error: this option is for the autosomal SNPs only. Please check the option --autosome.");
    }
}

void gcta::check_chrX() {
    for (int i = 0; i < _include.size(); i++) {
        if (_chr[_include[i]] != (_autosome_num + 1)) throw ("Error: this option is for SNPs on the X chromosome only.");
    }
}

void gcta::check_sex() {
    for (int i = 0; i < _keep.size(); i++) {
        if (_sex[_keep[i]] != 1 && _sex[_keep[i]] != 2) throw ("Error: Sex information of the individual \"" + _fid[_keep[i]] + " " + _pid[_keep[i]] + "\" is missing.\nUse --update-sex option to update the sex information of the individuals.");
    }
}

void gcta::make_grm(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, double make_grm_scl, bool mlmassoc, bool diag_f3_flag, string subpopu_file)
{
    bool have_mis = false;

    if (!grm_d_flag && grm_xchr_flag) check_chrX();
    else check_autosome();

    unsigned long i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();
    if(grm_d_flag) have_mis = make_XMat_d(_geno);
    else  have_mis = make_XMat(_geno);

    eigenVector sd_SNP;
    if (grm_d_flag) std_XMat_d(_geno, sd_SNP, false, make_grm_scl);
    else{
        if(subpopu_file.empty()) std_XMat(_geno, sd_SNP, grm_xchr_flag, false, make_grm_scl);
        else std_XMat_subpopu(subpopu_file, _geno, sd_SNP, grm_xchr_flag, false, make_grm_scl);
    }

    if (!mlmassoc) cout << "\nCalculating the" << ((grm_d_flag) ? " dominance" : "") << " genetic relationship matrix (GRM)" << (grm_xchr_flag ? " for the X chromosome" : "") << (_dosage_flag ? " using imputed dosage data" : "") << " ... (Note: default speed-optimized mode, may use huge RAM)" << endl;
    else cout << "\nCalculating the genetic relationship matrix (GRM) ... " << endl;

    // count the number of missing genotypes
    vector< vector<int> > miss_pos;
    vector< vector<bool> > X_bool;
    if(have_mis){
        miss_pos.resize(n);
        X_bool.resize(n);
        #pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            X_bool[i].resize(m);
            for (j = 0; j < m; j++) {
                if (_geno(i,j) < 1e5) X_bool[i][j] = true;
                else {
                    _geno(i,j) = 0.0;
                    miss_pos[i].push_back(j);
                    X_bool[i][j] = false;
                }
            }
        }
    }

    // Calculate A_N matrix
    if(have_mis){
        _grm_N.resize(n, n);
        #pragma omp parallel for private(j, k)
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                int miss_j = 0;
                for (k = 0; k < miss_pos[j].size(); k++) miss_j += (int) X_bool[i][miss_pos[j][k]];
                _grm_N(i,j) = m - miss_pos[i].size() - miss_j;
                _grm_N(j,i) = _grm_N(i,j);
            }
        }
    }
    else _grm_N = MatrixXf::Constant(n,n,m);

    // Calcuate WW'
    #ifdef SINGLE_PRECISION
    _grm = _geno * _geno.transpose();
    #else
    _grm = (_geno * _geno.transpose()).cast<double>();
    #endif

    // Calculate A matrix
    for (j = 0; j < m; j++) {
        sd_SNP[j] = powf(sd_SNP[j]*sd_SNP[j], make_grm_scl-1.0);
    }
    _grm_N = _grm_N.array() * sd_SNP.mean();
    
    #pragma omp parallel for private(j)
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            if(_grm_N(i,j) > 0) _grm(i,j) /= _grm_N(i,j);
            else _grm(i,j) = 0.0;
            _grm(j,i) = _grm(i,j);
        }
    }

    // GRM summary
    double diag_m = 0.0, diag_v = 0.0, off_m = 0.0, off_v = 0.0;
    calcu_grm_var(diag_m, diag_v, off_m, off_v);
    cout<<"\nSummary of the GRM:" << endl;
    cout<<"Mean of diagonals = "<<diag_m<<endl;
    cout<<"Variance of diagonals = "<<diag_v<<endl;
    cout<<"Mean of off-diagonals = " << off_m <<endl;
    cout<<"Variance of off-diagonals = " << off_v <<endl;

    // re-calcuate the diagonals (Fhat3+1)
    if (diag_f3_flag) {
        #pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            _grm(i,i) = 0.0;
            double non_missing = 0.0;
            for (j = 0; j < m; j++) {
                if (_geno(i,j) < 1e5){
                    _grm(i,i) += _geno(i,j)*(_geno(i,j)+(_mu[_include[j]] - 1.0) * sd_SNP[j]);
                    non_missing += 1.0;
                } 
            }
            _grm(i,i) /= non_missing; 
        }
    }

    if (inbred) {
        #pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) _grm(i,j) *= 0.5;
        }
    }

    if (mlmassoc && grm_mtd == 0) {
        for (j = 0; j < m; j++) {
            if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
            else sd_SNP[j] = 1.0 / sd_SNP[j];
        }
        #pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                if (_geno(i,j) < 1e5) _geno(i,j) *= sd_SNP[j];
                else _geno(i,j) = 0.0;
            }
        }
    } 
    else {
        // Output A_N and A
        string out_buf = _out;
        if (grm_d_flag) _out += ".d";
        output_grm(output_bin);
        _out = out_buf;
    }
}

void gcta::calcu_grm_var(double &diag_m, double &diag_v, double &off_m, double &off_v)
{
    int i = 0, n = _keep.size();
    diag_m = _grm.diagonal().mean();
    diag_v = (_grm.diagonal() - eigenVector::Constant(n, diag_m)).squaredNorm() / ((double)n - 1.0);
    double off_num = 0.5*n*(n - 1.0);
    off_m = 0.0;
    for (i = 1; i < n; i++) off_m += _grm.row(i).segment(0, i).sum();
    off_m /= off_num;
    off_v = 0.0;
    for (i = 1; i < n; i++) off_v += (_grm.row(i).segment(0, i) -  eigenVector::Constant(i, off_m).transpose()).squaredNorm();
    off_v /= (off_num - 1.0);    
}

void gcta::output_grm(bool output_grm_bin)
{
    int i = 0, j = 0;
    string grm_file;
    if (output_grm_bin) {
        // Save matrix A in binary file
        grm_file = _out + ".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out | ios::binary);
        if (!A_Bin) throw ("Error: can not open the file [" + grm_file + "] to write.");
        float f_buf = 0.0;
        int size = sizeof (float);
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j <= i; j++) {
                f_buf = (float) (_grm(i, j));
                A_Bin.write((char*) &f_buf, size);
            }
        }
        A_Bin.close();
        cout << "GRM of " << _keep.size() << " individuals has been saved in the file [" + grm_file + "] (in binary format)." << endl;

        string grm_N_file = _out + ".grm.N.bin";
        fstream N_Bin(grm_N_file.c_str(), ios::out | ios::binary);
        if (!N_Bin) throw ("Error: can not open the file [" + grm_N_file + "] to write.");
        f_buf = 0.0;
        size = sizeof (int);
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j <= i; j++) {
                f_buf = (float) (_grm_N(i, j));
                N_Bin.write((char*) &f_buf, size);
            }
        }
        N_Bin.close();
        cout << "Number of SNPs to calcuate the genetic relationship between each pair of individuals has been saved in the file [" + grm_N_file + "] (in binary format)." << endl;
    } 
    else {
        // Save A matrix in txt format
        grm_file = _out + ".grm.gz";
        gzofstream zoutf;
        zoutf.open(grm_file.c_str());
        if (!zoutf.is_open()) throw ("Error: can not open the file [" + grm_file + "] to write.");
        cout << "Saving the genetic relationship matrix to the file [" + grm_file + "] (in compressed text format)." << endl;
        zoutf.setf(ios::scientific);
        zoutf.precision(6);
        for (i = 0; i < _keep.size(); i++) {
            if (_grm_N.rows() > 0){
                zoutf.setf(ios::scientific);
                zoutf.precision(6);
                for (j = 0; j <= i; j++) zoutf << i + 1 << '\t' << j + 1 << '\t' << _grm_N(i, j) << '\t' << _grm(i, j) << endl;
            }
            else{ 
                for (j = 0; j <= i; j++) zoutf << i + 1 << '\t' << j + 1 << "\t0\t" << _grm(i, j) << endl;
            }
        }
        zoutf.close();
        cout << "The genetic relationship matrix has been saved in the file [" + grm_file + "] (in compressed text format)." << endl;
    }

    string famfile = _out + ".grm.id";
    ofstream Fam(famfile.c_str());
    if (!Fam) throw ("Error: can not open the file [" + famfile + "] to write.");
    for (i = 0; i < _keep.size(); i++) Fam << _fid[_keep[i]] + "\t" + _pid[_keep[i]] << endl;
    Fam.close();
    cout << "IDs for the GRM file [" + grm_file + "] have been saved in the file [" + famfile + "]." << endl;
}

int gcta::read_grm_id(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only)
{
    // read GRM IDs
    string grm_id_file = grm_file + ".grm.id";
    if (out_id_log) cout << "Reading IDs of the GRM from [" + grm_id_file + "]." << endl;
    ifstream i_grm_id(grm_id_file.c_str());
    if (!i_grm_id) throw ("Error: can not open the file [" + grm_id_file + "] to read.");
    string str_buf, id_buf;
    vector<string> fid, pid;
    grm_id.clear();
    while (i_grm_id) {
        i_grm_id >> str_buf;
        if (i_grm_id.eof()) break;
        fid.push_back(str_buf);
        id_buf = str_buf + ":";
        i_grm_id >> str_buf;
        pid.push_back(str_buf);
        id_buf += str_buf;
        grm_id.push_back(id_buf);
        getline(i_grm_id, str_buf);
    }
    i_grm_id.close();
    int n = grm_id.size();
    if (out_id_log) cout << n << " IDs read from [" + grm_id_file + "]." << endl;

    if (_id_map.empty()) {
        _fid = fid;
        _pid = pid;
        _indi_num = _fid.size();
        _sex.resize(_fid.size());
        init_keep();
    }

    return (n);
}

void gcta::read_grm(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only, bool dont_read_N)
{
    if (_grm_bin_flag) read_grm_bin(grm_file, grm_id, out_id_log, read_id_only, dont_read_N);
    else read_grm_gz(grm_file, grm_id, out_id_log, read_id_only);
}

void gcta::read_grm_gz(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only) {
    int n = read_grm_id(grm_file, grm_id, out_id_log, read_id_only);

    if (read_id_only) return;

    string grm_gzfile = grm_file + ".grm.gz", str_buf;
    const int MAX_LINE_LENGTH = 1000;
    char buf[MAX_LINE_LENGTH];
    gzifstream zinf;
    zinf.open(grm_gzfile.c_str());
    if (!zinf.is_open()) throw ("Error: can not open the file [" + grm_gzfile + "] to read.");

    int indx1 = 0, indx2 = 0, nline = 0;
    double grm_buf = 0.0, grm_N_buf;
    string errmsg = "Error: failed to read [" + grm_gzfile + "]. The format of the GRM file has been changed?\nError occurs in line:\n";
    cout << "Reading the GRM from [" + grm_gzfile + "]." << endl;
    _grm.resize(n, n);
    _grm_N.resize(n, n);
    while (1) {
        zinf.getline(buf, MAX_LINE_LENGTH, '\n');
        if (zinf.fail() || !zinf.good()) break;
        stringstream ss(buf);
        if (!(ss >> indx1)) throw (errmsg + buf);
        if (!(ss >> indx2)) throw (errmsg + buf);
        if (!(ss >> grm_N_buf)) throw (errmsg + buf);
        if (!(ss >> grm_buf)) throw (errmsg + buf);
        if (indx1 < indx2 || indx1 > n || indx2 > n) throw (errmsg + buf);
        if (grm_N_buf == 0) cout << "Warning: " << buf << endl;
        _grm_N(indx1 - 1, indx2 - 1) = _grm_N(indx2 - 1, indx1 - 1) = grm_N_buf;
        _grm(indx1 - 1, indx2 - 1) = _grm(indx2 - 1, indx1 - 1) = grm_buf;
        nline++;
        if (ss >> str_buf) throw (errmsg + buf);
    }
    zinf.close();
    if (!_within_family && nline != (int) (n * (n + 1)*0.5)){
        stringstream errmsg_tmp;
        errmsg_tmp << "Error: there are " << nline << " lines in the [" << grm_gzfile << "] file. The expected number of lines is " << (int) (n * (n + 1)*0.5) << "." << endl;
        throw(errmsg_tmp.str());
        //throw ("Error: incorrect number of lines in the grm file. *.grm.gz file and *.grm.id file are mismatched?");
    }
    cout << "GRM for " << n << " individuals are included from [" + grm_gzfile + "]." << endl;
}

void gcta::read_grm_bin(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only, bool dont_read_N)
{
    int i = 0, j = 0, n = read_grm_id(grm_file, grm_id, out_id_log, read_id_only);

    if (read_id_only) return;

    string grm_binfile = grm_file + ".grm.bin";
    ifstream A_bin(grm_binfile.c_str(), ios::in | ios::binary);
    if (!A_bin.is_open()) throw ("Error: can not open the file [" + grm_binfile + "] to read.");
    _grm.resize(n, n);
    cout << "Reading the GRM from [" + grm_binfile + "]." << endl;
    int size = sizeof (float);
    float f_buf = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            if (!(A_bin.read((char*) &f_buf, size))) throw ("Error: the size of the [" + grm_binfile + "] file is incomplete?");
            _grm(j, i) = _grm(i, j) = f_buf;
        }
    }
    A_bin.close();

    if(!dont_read_N){
        string grm_Nfile = grm_file + ".grm.N.bin";
        ifstream N_bin(grm_Nfile.c_str(), ios::in | ios::binary);
        if (!N_bin.is_open()) throw ("Error: can not open the file [" + grm_Nfile + "] to read.");
        _grm_N.resize(n, n);
        cout << "Reading the number of SNPs for the GRM from [" + grm_Nfile + "]." << endl;
        size = sizeof (float);
        f_buf = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                if (!(N_bin.read((char*) &f_buf, size))) throw ("Error: the size of the [" + grm_Nfile + "] file is incomplete?");
                _grm_N(j, i) = _grm_N(i, j) = f_buf;
            }
        }
        N_bin.close();
    }

    cout << "GRM for " << n << " individuals are included from [" + grm_binfile + "]." << endl;
}

void gcta::rm_cor_indi(double grm_cutoff) {
    cout << "Pruning the GRM with a cutoff of " << grm_cutoff << " ..." << endl;

    int i = 0, j = 0, i_buf = 0;

    // identify the positions where you see a value > than the threshold
    vector<int> rm_grm_ID1, rm_grm_ID2;
    for (i = 0; i < _keep.size(); i++) {
        for (j = 0; j < i; j++) {
            if (_grm(_keep[i], _keep[j]) > grm_cutoff) {
                rm_grm_ID1.push_back(_keep[i]);
                rm_grm_ID2.push_back(_keep[j]);
            }
        }
    }

    // count the number of appearance of each "position" in the vector, which involves a few steps
    vector<int> rm_uni_ID(rm_grm_ID1);
    rm_uni_ID.insert(rm_uni_ID.end(), rm_grm_ID2.begin(), rm_grm_ID2.end());
    stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    map<int, int> rm_uni_ID_count;
    for (i = 0; i < rm_uni_ID.size(); i++) {
        i_buf = count(rm_grm_ID1.begin(), rm_grm_ID1.end(), rm_uni_ID[i]) + count(rm_grm_ID2.begin(), rm_grm_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
    }

    // swapping
    map<int, int>::iterator iter1, iter2;
    for (i = 0; i < rm_grm_ID1.size(); i++) {
        iter1 = rm_uni_ID_count.find(rm_grm_ID1[i]);
        iter2 = rm_uni_ID_count.find(rm_grm_ID2[i]);
        if (iter1->second < iter2->second) {
            i_buf = rm_grm_ID1[i];
            rm_grm_ID1[i] = rm_grm_ID2[i];
            rm_grm_ID2[i] = i_buf;
        }
    }
    
    stable_sort(rm_grm_ID1.begin(), rm_grm_ID1.end());
    rm_grm_ID1.erase(unique(rm_grm_ID1.begin(), rm_grm_ID1.end()), rm_grm_ID1.end());
    vector<string> removed_ID;
    for (i = 0; i < rm_grm_ID1.size(); i++) removed_ID.push_back(_fid[rm_grm_ID1[i]] + ":" + _pid[rm_grm_ID1[i]]);

    // update _keep and _id_map
    update_id_map_rm(removed_ID, _id_map, _keep);

    cout << "After pruning the GRM, there are " << _keep.size() << " individuals (" << removed_ID.size() << " individuals removed)." << endl;
}

void gcta::adj_grm(double adj_grm_fac) {
    cout << "Adjusting the GRM for sampling errors ..." << endl;
    int i = 0, j = 0, n = _keep.size();
    double off_mean = 0.0, diag_mean = 0.0, off_var = 0.0, diag_var = 0.0, d_buf = 0.0;
    for (i = 0; i < n; i++) {
        diag_mean += _grm(_keep[i], _keep[i]);
        for (j = 0; j < i; j++) off_mean += _grm(_keep[i], _keep[j]);
    }
    diag_mean /= n;
    off_mean /= 0.5 * n * (n - 1.0);
    for (i = 0; i < n; i++) {
        d_buf = _grm(_keep[i], _keep[i]) - diag_mean;
        diag_var += d_buf*d_buf;
        for (j = 0; j < i; j++) {
            d_buf = _grm(_keep[i], _keep[j]) - off_mean;
            off_var += d_buf*d_buf;
        }
    }
    diag_var /= n - 1.0;
    off_var /= 0.5 * n * (n - 1.0) - 1.0;
    for (i = 0; i < _keep.size(); i++) {
        d_buf = 1.0 - (adj_grm_fac + 1.0 / _grm_N(_keep[i], _keep[i])) / diag_var;
        if (_grm(_keep[i], _keep[i]) > 0) _grm(_keep[i], _keep[i]) = 1.0 + d_buf * (_grm(_keep[i], _keep[i]) - 1.0);
        for (j = 0; j < i; j++) {
            if (_grm_N(_keep[i], _keep[j]) > 0) _grm(_keep[i], _keep[j]) *= 1.0 - (adj_grm_fac + 1.0 / _grm_N(_keep[i], _keep[j])) / off_var;
        }
    }
}

void gcta::dc(int dosage_compen) {
    cout << "Parameterizing the GRM under the assumption of ";
    if (dosage_compen == 1) cout << "full dosage compensation ..." << endl;
    else if (dosage_compen == 0) cout << "no dosage compensation ..." << endl;

    int i = 0, j = 0, i_buf = 0;
    double c1 = 1.0, c2 = 1.0;
    if (dosage_compen == 1) {
        c1 = 2.0;
        c2 = sqrt(2.0);
    }// full dosage compensation
    else if (dosage_compen == 0) {
        c1 = 0.5;
        c2 = sqrt(0.5);
    } // on dosage compensation
    for (i = 0; i < _keep.size(); i++) {
        for (j = 0; j <= i; j++) {
            i_buf = _sex[_keep[i]] * _sex[_keep[j]];
            if (i_buf == 1) _grm(i, j) *= c1;
            else if (i_buf == 2) _grm(i, j) *= c2;
        }
    }
}

void gcta::manipulate_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool dont_read_N)
{
    int i = 0, j = 0;

    vector<string> grm_id;
    if (merge_grm_flag) merge_grm(grm_file);
    else read_grm(grm_file, grm_id, true, false, dont_read_N);

    if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if (!remove_indi_file.empty()) remove_indi(remove_indi_file);
    if (grm_cutoff>-1.0) rm_cor_indi(grm_cutoff);
    if (!sex_file.empty()) update_sex(sex_file);
    if (adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
    if (dosage_compen>-1) dc(dosage_compen);
    if (grm_cutoff>-1.0 || !keep_indi_file.empty() || !remove_indi_file.empty()) {
        eigenMatrix grm_buf(_grm);
        _grm.resize(_keep.size(), _keep.size());
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j <= i; j++) _grm(i, j) = grm_buf(_keep[i], _keep[j]);
        }
        grm_buf.resize(0,0);
        MatrixXf grm_N_buf = _grm_N;
        _grm_N.resize(_keep.size(), _keep.size());
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j <= i; j++) _grm_N(i, j) = grm_N_buf(_keep[i], _keep[j]);
        }
    }
}

void gcta::save_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool output_grm_bin) {
    if (dosage_compen>-1) check_sex();
    manipulate_grm(grm_file, keep_indi_file, remove_indi_file, sex_file, grm_cutoff, adj_grm_fac, dosage_compen, merge_grm_flag);
    output_grm(output_grm_bin);
}

void gcta::merge_grm(string merge_grm_file) {
    vector<string> grm_files, grm_id;
    read_grm_filenames(merge_grm_file, grm_files);

    int f = 0, i = 0, j = 0;
    for (f = 0; f < grm_files.size(); f++) {
        read_grm(grm_files[f], grm_id, false, true);
        update_id_map_kp(grm_id, _id_map, _keep);
    }
    vector<string> uni_id;
    for (i = 0; i < _keep.size(); i++) uni_id.push_back(_fid[_keep[i]] + ":" + _pid[_keep[i]]);
    _n = uni_id.size();
    if (_n == 0) throw ("Error: no individual is in common in the GRM files.");
    else cout << _n << " individuals in common in the GRM files." << endl;

    vector<int> kp;
    eigenMatrix grm = eigenMatrix::Zero(_n, _n);
    eigenMatrix grm_N = eigenMatrix::Zero(_n, _n);
    for (f = 0; f < grm_files.size(); f++) {
        cout << "Reading the GRM from the " << f + 1 << "th file ..." << endl;
        read_grm(grm_files[f], grm_id);
        StrFunc::match(uni_id, grm_id, kp);
        for (i = 0; i < _n; i++) {
            for (j = 0; j <= i; j++) {
                if (kp[i] >= kp[j]) {
                    grm(i, j) += _grm(kp[i], kp[j]) * _grm_N(kp[i], kp[j]);
                    grm_N(i, j) += _grm_N(kp[i], kp[j]);
                } else {
                    grm(i, j) += _grm(kp[j], kp[i]) * _grm_N(kp[j], kp[i]);
                    grm_N(i, j) += _grm_N(kp[j], kp[i]);
                }
            }
        }
    }
    for (i = 0; i < _n; i++) {
        for (j = 0; j <= i; j++) {
            if (grm_N(i, j) == 0) _grm(i, j) = 0;
            else _grm(i, j) = grm(i, j) / grm_N(i, j);
            _grm_N(i, j) = grm_N(i, j);
        }
    }
    grm.resize(0, 0);
    grm_N.resize(0, 0);
    cout << "\n" << grm_files.size() << " GRMs have been merged together." << endl;
}

void gcta::read_grm_filenames(string merge_grm_file, vector<string> &grm_files, bool out_log) {
    ifstream merge_grm(merge_grm_file.c_str());
    if (!merge_grm) throw ("Error: can not open the file [" + merge_grm_file + "] to read.");
    string str_buf;
    grm_files.clear();
    vector<string> vs_buf;
    while (getline(merge_grm, str_buf)) {
        if (!str_buf.empty()) {
            if (StrFunc::split_string(str_buf, vs_buf) == 1) grm_files.push_back(vs_buf[0]);
        }
    }
    if (out_log) cout << "There are " << grm_files.size() << " GRM file names specified in [" + merge_grm_file + "]." << endl;
    if (grm_files.size() > 1000) throw ("Error: too many GRM file names specified in [" + merge_grm_file + "]. Maximum is 1000.");
    if (grm_files.size() < 1) throw ("Error: no GRM file name is found in [" + merge_grm_file + "].");
}

void gcta::grm_bK(string grm_file, string keep_indi_file, string remove_indi_file, double threshold, bool grm_out_bin_flag)
{
    int i = 0, j = 0;
    vector<string> grm_id;
    read_grm(grm_file, grm_id);
    if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if (!remove_indi_file.empty()) remove_indi(remove_indi_file);
    if (!keep_indi_file.empty() || !remove_indi_file.empty()) {
        eigenMatrix grm_buf(_grm);
        _grm.resize(_keep.size(), _keep.size());
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j <= i; j++) _grm(i, j) = grm_buf(_keep[i], _keep[j]);
        }
        grm_buf.resize(0,0);
        MatrixXf grm_N_buf = _grm_N;
        _grm_N.resize(_keep.size(), _keep.size());
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j <= i; j++) _grm_N(i, j) = grm_N_buf(_keep[i], _keep[j]);
        }
    }

    cout << "\nThe off-diagonals that are < " << threshold << " are set to zero.\n" << endl;
    for (i = 0; i < _keep.size(); i++) {
        for (j = 0; j < i; j++){
            if(_grm(i, j) < threshold) _grm(i, j) = 0.0;
        }
    }

    output_grm(grm_out_bin_flag);
}

void gcta::pca(string grm_file, string keep_indi_file, string remove_indi_file, double grm_cutoff, bool merge_grm_flag, int out_pc_num)
{
    manipulate_grm(grm_file, keep_indi_file, remove_indi_file, "", grm_cutoff, -2.0, -2, merge_grm_flag, true);
    _grm_N.resize(0, 0);
    int i = 0, j = 0, n = _keep.size();
    cout << "\nPerforming principal component analysis ..." << endl;

    SelfAdjointEigenSolver<MatrixXd> eigensolver(_grm.cast<double>());
    MatrixXd evec = (eigensolver.eigenvectors());
    VectorXd eval = eigensolver.eigenvalues();

    string eval_file = _out + ".eigenval";
    ofstream o_eval(eval_file.c_str());
    if (!o_eval) throw ("Error: can not open the file [" + eval_file + "] to read.");
    for (i = n - 1; i >= 0; i--) o_eval << eval(i) << endl;
    o_eval.close();
    cout << "Eigenvalues of " << n << " individuals have been saved in [" + eval_file + "]." << endl;
    string evec_file = _out + ".eigenvec";
    ofstream o_evec(evec_file.c_str());
    if (!o_evec) throw ("Error: can not open the file [" + evec_file + "] to read.");
    if (out_pc_num > n) out_pc_num = n;
    for (i = 0; i < n; i++) {
        o_evec << _fid[_keep[i]] << " " << _pid[_keep[i]];
        for (j = n - 1; j >= (n - out_pc_num); j--) o_evec << " " << evec(i, j);
        o_evec << endl;
    }
    o_evec.close();
    cout << "The first " << out_pc_num << " eigenvectors of " << n << " individuals have been saved in [" + evec_file + "]." << endl;
}

void gcta::snp_pc_loading(string pc_file, int grm_N)
{
    // read eigenvectors and eigenvalues
    string eigenval_file = pc_file + ".eigenval";
    ifstream in_eigenval(eigenval_file.c_str());
    if (!in_eigenval) throw ("Error: can not open the file [" + eigenval_file + "] to read.");
    string eigenvec_file = pc_file + ".eigenvec";
    ifstream in_eigenvec(eigenvec_file.c_str());
    if (!in_eigenvec) throw ("Error: can not open the file [" + eigenvec_file + "] to read.");
  
    cout << "Reading eigenvectors from [" + eigenvec_file + "]." << endl;
    vector<string> eigenvec_ID;
    vector< vector<string> > eigenvec_str;
    int eigenvec_num = read_fac(in_eigenvec, eigenvec_ID, eigenvec_str);
    cout << eigenvec_num << " eigenvectors of " << eigenvec_ID.size() << " individuals are included from [" + eigenvec_file + "]." << endl;
    update_id_map_kp(eigenvec_ID, _id_map, _keep);

    cout << "\nReading eigenvalues from [" + eigenval_file + "]." << endl;
    vector<double> eigenval_buf;
    double d_buf = 0.0;
    int eigenval_num = 0;
    while(in_eigenval && eigenval_num < eigenvec_num){
        in_eigenval >> d_buf;
        if(d_buf > 1e10 || d_buf < 1e-10) throw("Error: invalid eigenvalue in the file [" + eigenval_file + "].");
        eigenval_buf.push_back(d_buf);
        eigenval_num++;
    }
    if(eigenvec_num != eigenval_num) throw("Error: inconsistent numbers of eigenvalues and eigenvectors in the files [" + eigenval_file + "] and [" + eigenvec_file + "]");
    cout << eigenval_num << " eigenvalues read from [" + eigenval_file + "]" << endl;  

    int i = 0, j = 0;
    vector<string> uni_id;
    map<string, int> uni_id_map;
    map<string, int>::iterator iter;
    for(i=0; i<_keep.size(); i++){
        uni_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
        uni_id_map.insert(pair<string,int>(_fid[_keep[i]]+":"+_pid[_keep[i]], i));
    }
    _n=_keep.size();
    if(_n < 1) throw("Error: no individual is in common between the input files.");
    cout << _n << " individuals in common between the input files are included in the analysis."<<endl;
    
    eigenMatrix eigenvec(eigenvec_num, _n);
    for(i = 0; i < eigenvec_ID.size(); i++){
        iter = uni_id_map.find(eigenvec_ID[i]);
        if(iter == uni_id_map.end()) continue;
        for(j = 0; j < eigenvec_num; j++) eigenvec(j, iter->second) = atof(eigenvec_str[i][j].c_str());
    }

    eigenVector inv_eigenval(eigenval_num);
    for(i = 0; i < eigenval_num; i++)  inv_eigenval(i) = 1.0 / (eigenval_buf[i] * grm_N);

    // calculating SNP loading
    if (_mu.empty()) calcu_mu();
    cout << "\nCalculating SNP loading ..." << endl;
    int m = _include.size();
    eigenMatrix snp_loading(m, eigenvec_num);
    eigenVector x(_n);
    for(j = 0; j < m ; j++) {
        makex_eigenVector(j, x, false, true);
        x = x.array() / sqrt(_mu[_include[j]]*(1.0 - 0.5*_mu[_include[j]]));
        snp_loading.row(j) = (eigenvec * x).array() * inv_eigenval.array();
    }

    string filename = _out + ".pcl";
    cout << "\nSaving the PC loading of " << m << " SNPs to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if(!ofile) throw("Can not open the file [" + filename + "] to write.");
    ofile << "SNP\trefA";
    for(i = 0; i < eigenval_num; i++) ofile << "\tpc" << i+1 << "_loading";
    ofile << endl;
    for(i = 0; i < m; i++){
        ofile << _snp_name[_include[i]] << "\t" << _ref_A[_include[i]];
        for(j = 0; j < eigenvec_num; j++) ofile << "\t" << snp_loading(i, j);
        ofile << "\n";
    }
    ofile.close();
}



