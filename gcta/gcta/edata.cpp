/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions to build a ERM
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::read_efile(string efile)
{
    ifstream einf;
    einf.open(efile.c_str());
    if (!einf.is_open()) throw ("Error: can not open the file [" + efile + "] to read.");
    cout << "Reading gene expression / methylation data from [" + efile + "] ..." << endl;
    
    string str_buf="";
    vector<string> vs_buf;
    getline(einf, str_buf); // reading the header
    int col_num = StrFunc::split_string(str_buf, vs_buf, " \t\n");
    if(col_num < 3) throw ("Error: there needs be at least 3 columns in the file [" + efile + "].");
    _probe_num = col_num - 2;
    _probe_name.resize(_probe_num);
    int i=0;
    for(i=0; i<_probe_num; i++) _probe_name[i]=vs_buf[i+2];
    _indi_num = 0;
    while(getline(einf,str_buf)) _indi_num++;
    einf.close();

    einf.open(efile.c_str());
    getline(einf, str_buf);
    i=0;
    int j=0;
    stringstream errmsg;
    _fid.resize(_indi_num);
    _pid.resize(_indi_num);
    _probe_data.resize(_indi_num, _probe_num);
    string id_buf="";
    while (getline(einf, str_buf)) {
        stringstream ss(str_buf);
        if (!(ss >> id_buf)){ errmsg<<"Error: in line "<<i+2<<"."; throw(errmsg.str()); }
        _fid[i]=id_buf;
        if (!(ss >> id_buf)){ errmsg<<"Error: in line "<<i+2<<"."; throw(errmsg.str()); }
        _pid[i]=id_buf;
        for(j=0; j<_probe_num; j++){
            if (!(ss >> id_buf)){ errmsg<<"Error: in line "<<i+2<<"."; throw(errmsg.str()); }
            if(id_buf=="-9") _probe_data(i,j)=1e10;
            else _probe_data(i,j)=atof(id_buf.c_str());
        }
        i++;
    }
    einf.close();
    cout<<"Expression data for "<<_probe_num<<" probes of "<<_indi_num<<" individuals have been included from the file [" + efile + "]."<<endl; 

    // Initialize _keep and _e_include
    init_keep();
    init_e_include();
}

void gcta::init_e_include() {
    _e_include.clear();
    _e_include.resize(_probe_num);
    _probe_name_map.clear();
    int i = 0, size = 0;
    for (i = 0; i < _probe_num; i++) {
        _e_include[i] = i;
        _probe_name_map.insert(pair<string, int>(_probe_name[i], i));
        if (size == _probe_name_map.size()) throw ("Error: Duplicate probe name found: \"" + _probe_name[i] + "\".");
        size = _probe_name_map.size();
    }
}

void gcta::std_probe(vector< vector<bool> > &X_bool, bool divid_by_std)
{
    eigenMatrix X(_probe_data);
    _probe_data.resize(_keep.size(), _e_include.size());

    int i=0, j=0;
    #pragma omp parallel for private(j)
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<_e_include.size(); j++) _probe_data(i,j)=X(_keep[i], _e_include[j]);
    }
    eigenVector mu(_e_include.size()), nonmiss(_e_include.size());

    X_bool.resize(_keep.size());
    for(i=0; i<_keep.size(); i++) X_bool[i].resize(_e_include.size());
    #pragma omp parallel for private(i)
    for(j=0; j<_e_include.size(); j++){
        mu(j)=0.0;
        nonmiss(j)=0.0;
        for(i=0; i<_keep.size(); i++){
            if(_probe_data(i,j)<1e9){
                mu(j)+=_probe_data(i,j);
                nonmiss(j)+=1.0;
                X_bool[i][j] = true;
            }
            else X_bool[i][j] = false;
        }
        mu(j)/=nonmiss(j);
    }

    #pragma omp parallel for private(j)
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<_e_include.size(); j++){
            if(_probe_data(i,j)<1e9) _probe_data(i,j) -= mu(j);
            else _probe_data(i,j) = 0.0;
        }
    }

    if(divid_by_std){
        eigenVector sd(_e_include.size());
        #pragma omp parallel for
        for(j=0; j<_e_include.size(); j++){
            sd(j) = sqrt((_probe_data.col(j).dot(_probe_data.col(j))) / (nonmiss(j) - 1.0));
        }

        #pragma omp parallel for private(j)
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<_e_include.size(); j++) _probe_data(i,j) /= sd(j);
        }
    }
}

void gcta::std_probe_ind(vector< vector<bool> > &X_bool, bool divid_by_std)
{
     int i = 0, j = 0, n = _keep.size(), m = _e_include.size();
    eigenMatrix X(_probe_data);
    _probe_data.resize(n, m);

    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        for(j=0; j<m; j++) _probe_data(i,j)=X(_keep[i], _e_include[j]);
    }
    X.resize(0,0);

    eigenVector mu(n), nonmiss(n);

    X_bool.resize(n);
    for(i=0; i<n; i++) X_bool[i].resize(m);
    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        mu(i)=0.0;
        nonmiss(i)=0.0;
        for(j=0; j<m; j++){
            if(_probe_data(i,j)<1e9){
                mu(i)+=_probe_data(i,j);
                nonmiss(i)+=1.0;
                X_bool[i][j] = true;
            }
            else X_bool[i][j] = false;
        }
        mu(i)/=nonmiss(i);
    }

    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            if(_probe_data(i,j)<1e9) _probe_data(i,j) -= mu(i);
            else _probe_data(i,j) = 0.0;
        }
    }

    if(divid_by_std){
        eigenVector sd(n);
        #pragma omp parallel for
        for(i=0; i<n; i++){
            sd(i) = sqrt((_probe_data.row(i).dot(_probe_data.row(i))) / (nonmiss(i) - 1.0));
        }

        #pragma omp parallel for private(j)
        for(i=0; i<n; i++){
            for(j=0; j<m; j++) _probe_data(i,j) /= sd(i);
        }
    }
}

void gcta::make_erm(int erm_mtd, bool output_bin)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _e_include.size();
    
    cout << "Recoding gene expression / methylation data ..." << endl;
    bool divid_by_std = false;
    vector< vector<bool> > X_bool;
    if(erm_mtd < 2){
        if(erm_mtd == 0) divid_by_std = true;
        else if(erm_mtd == 1) divid_by_std = false;
        std_probe(X_bool, divid_by_std);
    }
    else std_probe_ind(X_bool, false);

    cout << "\nCalculating expression relationship matrix (ERM) ... " << endl;

    // count the number of missing genotypes
    vector< vector<int> > miss_pos(n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (X_bool[i][j] == false) miss_pos[i].push_back(j);
        }
    }

    // Calculate A_N matrix
    _grm_N.resize(n, n);
    if(erm_mtd == 0){
        #pragma omp parallel for private(j, k)
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                int miss_j = 0;
                for (k = 0; k < miss_pos[j].size(); k++) miss_j += (int) X_bool[i][miss_pos[j][k]];
                _grm_N(i,j) = m - miss_pos[i].size() - miss_j;
            }
        }
    }
    else if (erm_mtd > 0){
        if (erm_mtd == 1){
            eigenVector nonmiss(m);
            #pragma omp parallel for private(i)
            for (j = 0; j < m; j++) {
                nonmiss[j] = 0.0;
                for (i = 0; i < n; i++) {
                    if (X_bool[i][j] == true) nonmiss[j] += 1.0;
                }
            }
            eigenVector var(m);
            #pragma omp parallel for
            for(j=0; j<m; j++){
                if((nonmiss(j) - 1.0) > 0.0) var(j) = (_probe_data.col(j).dot(_probe_data.col(j))) / (nonmiss(j) - 1.0);
                else var(j) = 0.0;

            }
            double sum_var = var.sum();
            #pragma omp parallel for private(j, k)
            for (i = 0; i < n; i++) {
                double i_miss_sum_var = 0.0;
                for (k = 0; k < miss_pos[i].size(); k++) i_miss_sum_var += var[miss_pos[i][k]];
                for (j = 0; j <= i; j++) {
                    double j_miss_sum_var = 0.0;
                    for (k = 0; k < miss_pos[j].size(); k++){
                        if (X_bool[i][miss_pos[j][k]] == true) j_miss_sum_var += var[miss_pos[j][k]];
                    }
                    _grm_N(i,j) = sum_var - i_miss_sum_var - j_miss_sum_var;
                }
            }
        }
        else{
            eigenVector ssx(n);
            #pragma omp parallel for
            for(i=0; i<n; i++) ssx(i) = _probe_data.row(i).dot(_probe_data.row(i));
            #pragma omp parallel for private(j, k)
            for (i = 0; i < n; i++) {
                for (j = 0; j <= i; j++) {
                    double ssx_i = ssx(i);
                    double ssx_j = ssx(j);
                    for (k = 0; k < miss_pos[j].size(); k++){
                        int l = miss_pos[j][k];
                        if(X_bool[i][l] == true) ssx_i -= _probe_data(i,l) * _probe_data(i,l);
                    }
                    for (k = 0; k < miss_pos[i].size(); k++){
                        int l = miss_pos[i][k];
                        if(X_bool[j][l] == true) ssx_j -= _probe_data(j,l) * _probe_data(j,l);
                    }
                    _grm_N(i,j) = sqrt(ssx_i*ssx_j);            
                }
            }
        }
    }

    // Calculate A matrix
    _grm = _probe_data * _probe_data.transpose();
    #pragma omp parallel for private(j)
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            if (_grm_N(i,j) > 0.0) _grm(i,j) /= _grm_N(i,j);
            else _grm(i,j) = 0.0;
        }
    }
    _grm = _grm.array() / _grm.diagonal().mean();

    // Output A_N and A
    string out_buf = _out;
    _out += ".E";
    output_grm(output_bin);
    _out = out_buf;
}
