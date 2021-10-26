/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for BLUP analysis using sumamry data from EWAS
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */


#include "gcta.h"

void gcta::calcu_eR()
{
    // _probe_data rotated
    eigenMatrix X(_probe_data);
    _probe_data.resize(_e_include.size(), _keep.size());

    int i=0, j=0;
    #pragma omp parallel for private(j)
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<_e_include.size(); j++) _probe_data(j,i)=X(_keep[i], _e_include[j]);
    }
    eigenVector m(_e_include.size()), nonmiss(_e_include.size());

    #pragma omp parallel for private(i)
    for(j=0; j<_e_include.size(); j++){
        m(j)=0.0;
        nonmiss(j)=0.0;
        for(i=0; i<_keep.size(); i++){
            if(_probe_data(j,i)<1e9){
                m(j)+=_probe_data(j,i);
                nonmiss(j)+=1.0;
            }
        }
        m(j)/=nonmiss(j);
    }

    #pragma omp parallel for private(j)
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<_e_include.size(); j++){
            if(_probe_data(j,i)<1e9) _probe_data(j,i)-=m(j);
            else _probe_data(j,i)=0.0;
        }
    }
    eigenVector d(_e_include.size());

    #pragma omp parallel for
    for(j=0; j<_e_include.size(); j++){
        d(j)=_probe_data.row(j).dot(_probe_data.row(j));
    }

    _ecojo_wholeR=_probe_data*_probe_data.transpose();

    #pragma omp parallel for private(j)
    for(i=0; i<_e_include.size(); i++){
        _ecojo_wholeR(i,i)=1.0;
        for(j=i+1; j<_e_include.size(); j++) _ecojo_wholeR(i,j)=_ecojo_wholeR(j,i)=_ecojo_wholeR(i,j)/sqrt(d(i)*d(j));
    }
    _probe_data.resize(0,0);
}

void gcta::read_eR(string eR_file)
{
    ifstream eR_inf(eR_file.c_str());
    if (!eR_inf.is_open()) LOGGER.e(0, "cannot open the file [" + eR_file + "] to read.");
    LOGGER << "Reading correlation matrix of gene expression from [" + eR_file + "] ..." << endl;
    
    string str_buf="";
    getline(eR_inf, str_buf); // reading the probe names
    _probe_num = StrFunc::split_string(str_buf, _probe_name, " \t\n");
    LOGGER<<_probe_num<<" probes found in the file. \nReading correlation matrix ..."<<endl; 
    int i = 0, j = 0;
    _ecojo_wholeR.resize(_probe_num, _probe_num);
    for(i = 0; i < _probe_num; i++) {
        for(j = 0; j < _probe_num; j++) {
            if(!(eR_inf >> _ecojo_wholeR(i,j))) LOGGER.e(0, "incorrect format of [" + eR_file + "].");
        }
    }
    eR_inf.close();
    LOGGER<<"Correlation matrix for "<<_probe_num<<" probes have been included from the file [" + eR_file + "]."<<endl; 

    init_e_include();
}

void gcta::read_e_metafile(string e_metafile)
{
    LOGGER << "\nReading expression-trait association summary-level statistics from [" + e_metafile + "] ..." << endl;
    ifstream e_meta(e_metafile.c_str());
    if (!e_meta) LOGGER.e(0, "cannot open the file [" + e_metafile + "] to read.");

    string str_buf="";
    double d_buf=0.0;
    vector<string> vs_buf, probe_buf;
    vector<double> z_buf, n_buf;
    
    getline(e_meta, str_buf); // the header line
    if (StrFunc::split_string(str_buf, vs_buf) < 3) LOGGER.e(0, "there needs to be at least 3 columns in the file [" + e_metafile + "].");
    stringstream errmsg;
    int line=1;
    while(getline(e_meta, str_buf)){
        stringstream iss(str_buf);
        if(!(iss >> str_buf)){ errmsg<<"in line "<<line<<"."; LOGGER.e(0, errmsg.str()); }
        if (_probe_name_map.find(str_buf) == _probe_name_map.end()) continue;
        probe_buf.push_back(str_buf);
        if(!(iss >> d_buf)){ errmsg<<"in line "<<line<<"."; LOGGER.e(0, errmsg.str()); }
        z_buf.push_back(d_buf);
        if(!(iss >> d_buf)){ errmsg<<"in line "<<line<<"."; LOGGER.e(0, errmsg.str()); }
        n_buf.push_back(d_buf);
        line++;
    }
    e_meta.close();
    if(probe_buf.size()<1) LOGGER.e(0, "no probe remains in the analysis.");
    LOGGER << "GWAS summary statistics of " << probe_buf.size() << " probs read from [" + e_metafile + "]." << endl;

    LOGGER << "Matching the summary data to the genotype data ..." << endl;
    update_id_map_kp(probe_buf, _probe_name_map, _e_include);
    vector<int> indx(_e_include.size());
    map<string, int> id_map;
    int i=0;
    for (i = 0; i < probe_buf.size(); i++) id_map.insert(pair<string, int>(probe_buf[i], i));
    map<string, int>::iterator iter;
    for (i = 0; i < _e_include.size(); i++) {
        iter = id_map.find(_probe_name[_e_include[i]]);
        indx[i] = iter->second;
    }
    _ecojo_z.resize(_e_include.size());
    _ecojo_b.resize(_e_include.size());
    _ecojo_se.resize(_e_include.size());
    _ecojo_n.resize(_e_include.size());
    _ecojo_pval.resize(_e_include.size());

    //#pragma omp parallel for private(d_buf)
    for (i = 0; i < _e_include.size(); i++) {
        _ecojo_z[i]=z_buf[indx[i]];
        _ecojo_pval[i] = StatFunc::pchisq(_ecojo_z[i]*_ecojo_z[i], 1);
        _ecojo_n[i]=n_buf[indx[i]];
        d_buf=sqrt(_ecojo_z[i]*_ecojo_z[i] + _ecojo_n[i]);
        if(d_buf<1e-30){
            _ecojo_b[i]=0.0;
            _ecojo_se[i]=-1.0;    
        }
        else{
            _ecojo_b[i]=_ecojo_z[i]/d_buf; 
            _ecojo_se[i]=1.0/d_buf;
        }
    }
 }

void gcta::run_ecojo_slct(string e_metafile, double p_cutoff, double collinear)
{
    bool joint_only=false, backward=false;
    _ecojo_p_cutoff = p_cutoff;
    _ecojo_collinear = collinear;
    read_e_metafile(e_metafile);
    calcu_eR();

    int i = 0, j = 0;
    vector<int> slct, remain;
    eigenVector bC, bC_se, pC;
    LOGGER << endl;
    if (!joint_only && !backward) {
        LOGGER << "Performing stepwise model selection on " << _e_include.size() << " probes to select association signals ... (p-value cutoff = " << _ecojo_p_cutoff << "; ";
        LOGGER << "collinearity cutoff = " << _ecojo_collinear << ")"<< endl;
        ecojo_slct(slct, remain, bC, bC_se, pC);
        if (slct.empty()) {
            LOGGER << "No probe has been selected." << endl;
            return;
        }
    }
    else {
        for (i = 0; i < _e_include.size(); i++) slct.push_back(i);
        if (backward) {
            LOGGER << "Performing backward selection on " << _e_include.size() << " probes at p-value cutoff = " << _ecojo_p_cutoff << " ..." << endl;
            ecojo_slct_stay(slct, bC, bC_se, pC);
        }
    }

    // joint analysis
    eigenVector bJ, bJ_se, pJ;
    LOGGER << "Performing joint analysis on all the " << slct.size();
    if (joint_only) LOGGER << " probes ..." << endl;
    else LOGGER << " selected signals ..." << endl;
    if (slct.size() >= _keep.size()) LOGGER.e(0, "too many probes. The number of probes in a joint analysis should not be larger than the sample size.");
    ecojo_joint(slct, bJ, bJ_se, pJ);
    ecojo_slct_output(joint_only, slct, bJ, bJ_se, pJ);
}

void gcta::ecojo_slct_output(bool joint_only, vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ)
{
    string filename = _out + ".slct.ecojo";
    if (joint_only) LOGGER << "Saving the joint analysis result of " << slct.size() << " probes to [" + filename + "] ..." << endl;
    else LOGGER << "Saving the " << slct.size() << " independent signals to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Probe\tb\tse\tz\tn\tbJ\tbJ_se\tzJ\tpJ"<< endl;
    int i = 0, j = 0;
    for (i = 0; i < slct.size(); i++) {
        j = slct[i];
        ofile << _probe_name[_e_include[j]] << "\t" << _ecojo_b[j] << "\t" <<_ecojo_se[j] << "\t" << _ecojo_z[j] << "\t" << _ecojo_n[j] << "\t"<< bJ[i] << "\t" << bJ_se[i] << "\t" << bJ[i]/bJ_se[i] << "\t" << pJ[i] << "\t" << endl;
    }
    ofile.close();
}

void gcta::ecojo_slct(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC)
{
    int i = 0, i_buf = 0;
    vector<double> p_buf;
    eigenVector2Vector(_ecojo_pval, p_buf);
    int m = min_element(p_buf.begin(), p_buf.end()) - p_buf.begin();
    if (p_buf[m] >= _ecojo_p_cutoff) return;
    slct.push_back(m);
    for (i = 0; i < _e_include.size(); i++) {
        if (i != m) remain.push_back(i);
    }
    int prev_num = 0;
    ecojo_init_R(slct);
    ecojo_init_RC(slct, remain);
    if (_ecojo_p_cutoff > 1e-3) LOGGER << "Performing forward model selection because the significance level is too low..." << endl;

    while (!remain.empty()) {
        if (ecojo_slct_entry(slct, remain, bC, bC_se, pC)) {
            if (_ecojo_p_cutoff <= 1e-3) ecojo_slct_stay(slct, bC, bC_se, pC);
            ecojo_init_RC(slct, remain);
        }
        else break;        
        if (slct.size() % 5 == 0 && slct.size() > prev_num) LOGGER << slct.size() << " associated probes have been selected." << endl;
        if (slct.size() > prev_num) prev_num = slct.size();
    }
    if (_ecojo_p_cutoff > 1e-3) {
        LOGGER << "Performing backward elimination..." << endl;
        ecojo_slct_stay(slct, bC, bC_se, pC);
    }
    LOGGER << "Finally, " << slct.size() << " associated probes are selected." << endl;
}

bool gcta::ecojo_slct_entry(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC)
{
    int i = 0, m = 0;
    ecojo_cond(slct, remain, bC, bC_se, pC);
    vector<double> pC_buf;
    eigenVector2Vector(pC, pC_buf);

    while (true) {
        m = min_element(pC_buf.begin(), pC_buf.end()) - pC_buf.begin();
        if (pC_buf[m] >= _ecojo_p_cutoff){

            ecojo_init_RC(slct, remain);
            ecojo_cond(slct, remain, bC, bC_se, pC);


            // debug
/*            LOGGER<<"here"<<endl;
            ofstream tmp("cond.txt");
            for(int j=0; j<remain.size(); j++){
                tmp<<_probe_name[_e_include[remain[j]]]<<" "<<bC[j]<<" "<<bC_se[j]<<" "<<pC[j]<<endl;
            }
            tmp.close();
            ofstream oR("R.txt");
            for(int j=0; j<_ecojo_R.rows(); j++){
                for(int k=0; k<_ecojo_R.cols(); k++) oR<<_ecojo_R(j,k)<<" ";
                oR<<endl;
            }
            oR.close();
            ofstream oRC("RC.txt");
            for(int j=0; j<_ecojo_RC.rows(); j++){
                for(int k=0; k<_ecojo_RC.cols(); k++) oRC<<_ecojo_RC(j,k)<<" ";
                oRC<<endl;
            }
            oRC.close();*/

            return false;
        }
        if (ecojo_insert_R(slct, remain[m])){
            slct.push_back(remain[m]);
            stable_sort(slct.begin(), slct.end());
            remain.erase(remain.begin() + m);           
            return (true);
        }
        pC_buf.erase(pC_buf.begin() + m);
        remain.erase(remain.begin() + m);

        // debug
        //LOGGER<<"remain = "<<remain.size()<<endl;
    }
}

void gcta::ecojo_slct_stay(vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ)
{
    vector<double> pJ_buf;
    while(!slct.empty()){
        ecojo_joint(slct, bJ, bJ_se, pJ);
        eigenVector2Vector(pJ, pJ_buf);
        int m = max_element(pJ_buf.begin(), pJ_buf.end()) - pJ_buf.begin();
        if(pJ[m] > _ecojo_p_cutoff){
            slct.erase(slct.begin() + m);
            ecojo_erase_R(slct);
        }
        else break;
    }
}

void gcta::ecojo_joint(const vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ)
{
    int i = 0, size = slct.size();
    eigenVector b(size), n(size);
    for(i = 0; i < size; i++){
        b[i] = _ecojo_b[slct[i]];
        n[i] = _ecojo_n[slct[i]];
    }
    bJ = _ecojo_R*b;
    bJ_se = eigenVector::Constant(size, -1);
    pJ = eigenVector::Constant(size, 2);
    double chisq=0.0;
    for(i = 0; i < size; i++){
        bJ_se[i] = _ecojo_R(i,i) / n[i];
        if (bJ_se[i] > 1.0e-30) {
            bJ_se[i] = sqrt(bJ_se[i]);
             chisq = bJ[i] / bJ_se[i];
            pJ[i] = StatFunc::pchisq(chisq*chisq, 1);
        }
    }
}

void gcta::ecojo_cond(const vector<int> &slct, const vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC)
{
    int i=0, j=0;
    eigenVector b1(slct.size()), b2(remain.size());
    for(i=0; i<slct.size(); i++) b1[i]=_ecojo_b[slct[i]];
    for(i=0; i<remain.size(); i++) b2[i]=_ecojo_b[remain[i]];

    bC = eigenVector::Constant(remain.size(), 0);
    bC_se = eigenVector::Constant(remain.size(), -1);
    pC = eigenVector::Constant(remain.size(), 2);
    double chisq = 0.0;
    eigenVector RCRi;
    //#pragma omp parallel for private(RCRi)
    for (i = 0; i < remain.size(); i++) {
        RCRi=_ecojo_R*_ecojo_RC.col(i);
        bC[i]=b2[i]-RCRi.dot(b1);
        bC_se[i] = (1 - _ecojo_RC.col(i).dot(RCRi)) / _ecojo_n[remain[i]];
        if (bC_se[i] > 1e-30) {
            bC_se[i] = sqrt(bC_se[i]);
            chisq = bC[i] / bC_se[i];
            pC[i] = StatFunc::pchisq(chisq*chisq, 1);
        }
    }
}

bool gcta::ecojo_init_R(const vector<int> &slct)
{
    int i=0, j=0, size=slct.size();
    _ecojo_R.resize(size, size);

    #pragma omp parallel for private(j)
    for(i=0; i<size; i++){
        for(j=0; j<size; j++) _ecojo_R(i,j)=_ecojo_wholeR(slct[i],slct[j]);
    }
    ecojo_inv_R();
    if ((1 - eigenVector::Constant(size, 1).array() / _ecojo_R.diagonal().array()).maxCoeff() > _ecojo_collinear) return false;

    return true;
}

void gcta::ecojo_init_RC(const vector<int> &slct, const vector<int> &remain) {
    int i = 0, j = 0;
    _ecojo_RC.resize(slct.size(), remain.size());

    #pragma omp parallel for private(j)
    for (i = 0; i < slct.size(); i++){
        for (j = 0; j < remain.size(); j++) _ecojo_RC(i,j)=_ecojo_wholeR(slct[i], remain[j]);
    }
}

bool gcta::ecojo_insert_R(const vector<int> &slct, int insert_indx)
{
    eigenMatrix R_buf(_ecojo_R);
    vector<int> ix(slct);
    ix.push_back(insert_indx);
    stable_sort(ix.begin(), ix.end());
    int i = 0, j = 0;
    _ecojo_R.resize(ix.size(), ix.size());

    #pragma omp parallel for private(j)
    for(i=0; i<ix.size(); i++){
        for(j=0; j<ix.size(); j++) _ecojo_R(i,j)=_ecojo_wholeR(ix[i],ix[j]);
    }

    ecojo_inv_R();
    if((1 - eigenVector::Constant(ix.size(), 1).array() / _ecojo_R.diagonal().array()).maxCoeff() > _ecojo_collinear){
        _ecojo_R=R_buf;
        return false;
    }

    return true;
}

void gcta::ecojo_erase_R(const vector<int> &slct)
{
    int i = 0, j = 0;
    _ecojo_R.resize(slct.size(), slct.size());

    #pragma omp parallel for private(j)
    for(i=0; i<slct.size(); i++){
        for(j=0; j<slct.size(); j++) _ecojo_R(i,j)=_ecojo_wholeR(slct[i],slct[j]);
    }
    ecojo_inv_R();
}

void gcta::run_ecojo_blup_efile(string e_metafile, double lambda)
{
    read_e_metafile(e_metafile);
    LOGGER << "Recoding gene expression data ..." << endl;
    calcu_eR();
    ecojo_blup(lambda);
}

void gcta::run_ecojo_blup_eR(string e_metafile, double lambda)
{
    read_e_metafile(e_metafile);
    ecojo_blup(lambda);
}

void gcta::ecojo_blup(double lambda)
{
    LOGGER << "\nPerforming joint analysis on all the " << _e_include.size() << " probes ..." << endl;
    int i = 0, j=0;
    double d_n = _ecojo_n.mean();
    double diag_val=1.0+_e_include.size()*(1.0/lambda-1.0)/d_n;
    for(i=0; i<_e_include.size(); i++) _ecojo_wholeR(i,i)=diag_val;
    double logdet=0.0;
    if (!comput_inverse_logdet_LDLT_mkl(_ecojo_wholeR, logdet)) {
        LOGGER<<"Note: no solution to LDLT decomposition. Switching to LU decomposition."<<endl;
        _ecojo_wholeR = _ecojo_wholeR.lu().solve(eigenMatrix::Identity(_e_include.size(), _e_include.size()));
        //if (!comput_inverse_logdet_LU_mkl(_ecojo_wholeR, logdet)) LOGGER.e(0, "\n  the correlation matrix is not invertible.");
    }
    eigenVector bJ=_ecojo_wholeR*_ecojo_b;

    string filename = _out + ".blup.ecojo";
    LOGGER << "Saving the BLUP analysis result of " << _e_include.size() << " probes to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Probe\tb\tse\tz\tn\tb_blup"<< endl;
    for (i = 0; i < _e_include.size(); i++) {
        ofile << _probe_name[_e_include[i]] << "\t" << _ecojo_b[i] << "\t" <<_ecojo_se[i] << "\t" << _ecojo_z[i] << "\t" << _ecojo_n[i] << "\t"<< bJ[i] << endl;
    }
    ofile.close();
}

void gcta::ecojo_inv_R() {
    int i = 0, j = 0, k = 0;
    string errmsg = "\n  the correlation matrix is not invertible.";

    double logdet=0.0;
    if (!comput_inverse_logdet_LDLT_mkl(_ecojo_R, logdet)) {
        LOGGER<<"Note: no solution to LDLT decomposition. Switching to LU decomposition."<<endl;
        _ecojo_R = _ecojo_R.lu().solve(eigenMatrix::Identity(_ecojo_R.cols(), _ecojo_R.cols()));
    }
}

