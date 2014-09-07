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
#include "CommFunc.h"

// local function without declaration
eigenMatrix init_X(eigenMatrix &X, vector<int> keep, vector<int> include) {
    eigenMatrix xbuf ( X );
    int i=0, j=0;
    int nrow = keep.size(), ncol = include.size();
    
    X.resize(nrow, ncol);
    for( i=0; i<nrow; i++ ) {
        for( j=0; j<ncol; j++ ) {
            X(i,j)=xbuf(keep[i], include[j]);
        }
    }
    return X;
}
eigenVector est_geno_mean(eigenMatrix geno_data, int *miss_indx, vector<missValue> miss_vale_pos, int mean_center) {
    int i=0, ncol=geno_data.cols();
    eigenVector geno_mu(ncol), genobuf(ncol);
    
    if(mean_center) geno_mu.setZero();
    else {
        for( i=0; i<ncol; i++) {
            genobuf = geno_data.col(i);
            geno_mu(i) = CommFunc::eigenMean(genobuf, miss_indx[i], miss_vale_pos[miss_indx[i]]);
        }
    }
    return geno_mu;
}

eigenVector est_geno_var(eigenMatrix geno_data, int *miss_indx, vector<missValue> miss_vale_pos, int var_center) {
    int i=0, indxbuf =0, nmiss=0;
    int nrow = geno_data.rows(), ncol=geno_data.cols();
    double dbuf=0.0;
    eigenVector genobuf, geno_var;
    
    // var_center: 0 - method 2, 1 - method 1, 2 - method 3,
    if(!var_center) {
        // method 2
        geno_var.resize(ncol);
        for(i=0; i<ncol; i++) {
            genobuf = geno_data.col(i);
            indxbuf = miss_indx[i];
            dbuf = CommFunc::eigenSd(genobuf, indxbuf, miss_vale_pos[indxbuf]);
            dbuf = pow(dbuf, 2);
            geno_var(i) = dbuf;
        }
    } else if( var_center == 1 ) {
        // method 1
        geno_var.resize(ncol);
        for( i=0; i<ncol; i++) geno_var(i) = 1;
    } else {
        // method 3
        geno_var.resize(nrow);
        for( i=0; i<nrow; i++) {
            genobuf = geno_data.row(i);
            indxbuf = miss_indx[i];
            dbuf = CommFunc::eigenSd(genobuf, indxbuf, miss_vale_pos[indxbuf]);
            nmiss = indxbuf ? miss_vale_pos[indxbuf].num : 0;
            dbuf = pow(dbuf, 2);
            geno_var(i) = dbuf*(ncol-nmiss-1);
        }
    }
    return geno_var;
}

// main function declared in gcta.h
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
            StrFunc::to_upper(id_buf);
            if(id_buf=="-9" || id_buf=="NA") _probe_data(i,j)=1e10;
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

void gcta::read_tefile(string efile)
{
    ifstream einf;
    einf.open(efile.c_str());
    if (!einf.is_open()) throw ("Error: can not open the file [" + efile + "] to read.");
    cout << "Reading gene expression / methylation data from [" + efile + "] ..." << endl;
    
    string str_buf="";
    vector<string> vs_buf;
    getline(einf, str_buf); // reading the first line, family ID
    int col_num = StrFunc::split_string(str_buf, vs_buf, " \t\n");
    if(col_num < 2) throw ("Error: there needs be at least 2 columns in the file [" + efile + "].");
    _indi_num = col_num - 1;
    _fid.resize(_indi_num);
    _pid.resize(_indi_num);
    int i = 0;
    for(i=0; i<_indi_num; i++) _fid[i] = vs_buf[i+1];
    
    getline(einf, str_buf); // reading the second line, individual ID
    col_num = StrFunc::split_string(str_buf, vs_buf, " \t\n");
    if(col_num!=_indi_num+1) throw ("Error: number of family ID and individual ID is different.");
    for(i =0; i< _indi_num; i++) _pid[i] = vs_buf[i+1];
    
    _probe_num = 0; // read the number of probes
    while(getline(einf,str_buf)) _probe_num++;
    einf.close();
    
    einf.open(efile.c_str());  // re-open the file to read the gene expression and methylation data
    for( i = 0; i < 2; i++ ) getline(einf, str_buf);  // skip family ID and individual ID
    i=0;
    int j=0;
    stringstream errmsg;
    _probe_name.resize(_probe_num);
    _probe_data.resize(_indi_num, _probe_num);
    string id_buf="";
    while (getline(einf, str_buf)) {
        stringstream ss(str_buf);
        if (!(ss >> id_buf)){ errmsg<<"Error: in line "<<i+2<<"."; throw(errmsg.str()); }
        _probe_name[i]=id_buf;
        for(j=0; j<_indi_num; j++){
            if (!(ss >> id_buf)){ errmsg<<"Error: in line "<<i+2<<"."; throw(errmsg.str()); }
            StrFunc::to_upper(id_buf);
            if(id_buf=="-9" || id_buf=="NA") _probe_data(j,i)=1e10;
            else _probe_data(j,i)=atof(id_buf.c_str());
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

void gcta::make_erm(int erm_mtd, int erm_indi, bool output_bin)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _e_include.size();
    bool est_var = true;
    
    cout << "Recoding gene expression / methylation data ..." << endl;
    MatrixXi X_bool;
    
    // default: 0, method 2, estimate mean and var for each probe / marker
    // (xij - meani) * (xik - meani) / sum(var)
    int mean_center = 0, var_center = 0;
    // method 1: mean = 0, var = 1, don't estimate mean and var for each probe / makrer
    // wij * wik / m*1
    if(!erm_mtd) {
        mean_center = 1; var_center = 1;
    } else if( erm_mtd ==2 ) {
        //method 3:  ( xij - meanj ) * ( xik - meank)
        mean_center = 1; var_center = 2;
    }
    // initialize the X
    _probe_data = init_X(_probe_data, _keep, _e_include);
    // get the missing value
    X_bool = CommFunc::find_missing_X(_probe_data);
    // initialize the positions of miss value
    _miss_row_indx = CommFunc::find_missing_pos(X_bool, _miss_row_indx, _miss_row, true);
    _miss_col_indx = CommFunc::find_missing_pos(X_bool, _miss_col_indx, _miss_col);
    
    // pre-normalise the matrix
    // method 1: standise the x to w
    if(!erm_mtd) _probe_data = CommFunc::scale(_probe_data, _miss_col_indx, _miss_col);
    // method 3: xij - meanj
    else if(erm_mtd==2) _probe_data = CommFunc::scale(_probe_data, _miss_row_indx, _miss_row, true,  true);
    // --make-erm-indi: re-scale the individual-wised value
    if(erm_indi) {
        _probe_data = CommFunc::scale(_probe_data, _miss_row_indx, _miss_row, true);
        mean_center = var_center = 0;
    }
    
    eigenVector geno_mu, geno_var;
    
    // method 1 and method 3: set the mean to zero
    geno_mu = est_geno_mean( _probe_data, _miss_col_indx, _miss_col, mean_center );
    // method 1: var = m
    // method 2: var = sum(var(probe))
    // method 3: var = ssx(x.k), k - one individual
    if(erm_mtd != 2) {
        geno_var = est_geno_var(_probe_data, _miss_col_indx, _miss_col, var_center);
    } else {
        geno_var = est_geno_var(_probe_data, _miss_row_indx, _miss_row, var_center);
    }

    // estimate the _grm_N
    _grm_N.resize(n, n);
    double sumvar = 0.0, sumbuf=0.0, sumbuf2=0.0;
    missValue missbuf;
    // method 1 and 2, _grm_N = sum(var) - missingness
    if(  erm_mtd != 2 ) {
        sumvar = geno_var.sum();
        for( i=0; i<n; i++ ) {
            for( j=0; j<=i; j++) {
                sumbuf = sumvar;
                for(k=0; k<m; k++) {
                    // substract var, if one of xij or xik is equal to zero
                    if( (X_bool(i,k) & X_bool(j,k)) == 0 ) sumbuf -= geno_var(k);
                }
                _grm_N(i,j) = sumbuf;
            }
        }
    }  else {
        // method 3
        for( i=0; i<n; i++ ) {
            for( j=0; j<=i; j++) {
                _grm_N(i,j) = sqrt(geno_var(i) * geno_var(j));
            }
        }
    }
    
    // geno_data - mu for each column
    int indxbuf = 0, nmiss = 0, *posbuf = 0;
    eigenVector vecbuf;
    for( i = 0; i< m; i++) {
        vecbuf = _probe_data.col(i);
        for( j = 0; j < n; j++) vecbuf(j) -= geno_mu(i);
        // set the missing value to zero
        indxbuf = _miss_col_indx[i];
        if(indxbuf) {
            nmiss = _miss_col[indxbuf].num;
            posbuf = CommFunc::split_miss_pos(_miss_col[indxbuf], posbuf);
            vecbuf = CommFunc::assignVale(vecbuf, posbuf, nmiss);
        }
        _probe_data.col(i) = vecbuf;
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
    
    free(_miss_row_indx);
    free(_miss_col_indx);
}

// data management for gene expression and methylation
void gcta::extract_probe(string probelistfile)
{
    vector<string> probelist;
    read_snplist(probelistfile, probelist, "probes");
    update_id_map_kp(probelist, _probe_name_map, _e_include);
    cout << _e_include.size() << " probes are extracted from [" + probelistfile + "]." << endl;
}

void gcta::extract_single_probe(string probename)
{
    vector<string> probelist;
    probelist.push_back(probename);
    update_id_map_kp(probelist, _probe_name_map, _e_include);
    if (_e_include.empty()) throw ("Error: can not find the probe [" + probename + "] in the data.");
    else cout << "Only the probe [" + probename + "] is included in the analysis." << endl;
}

void gcta::exclude_probe(string probelistfile)
{
    vector<string> probelist;
    read_snplist(probelistfile, probelist);
    int prev_size = _e_include.size();
    update_id_map_rm(probelist, _probe_name_map, _e_include);
    cout << prev_size - _e_include.size() << " probes are excluded from [" + probelistfile + "] and there are " << _e_include.size() << " probes remaining." << endl;
}

void gcta::exclude_single_probe(string probename)
{
    vector<string> probelist;
    probelist.push_back(probename);
    int include_size = _e_include.size();
    update_id_map_rm(probelist, _probe_name_map, _e_include);
    if (_e_include.size() == include_size) throw ("Error: can not find the probe [" + probename + "] in the data.");
    else cout << "The probe [" + probename + "] has been excluded from the analysis." << endl;
}