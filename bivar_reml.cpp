/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for bivariate REML analysis
 *
 * 2012 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::fit_bivar_reml(string grm_file, string phen_file, string qcovar_file, string covar_file, string keep_indi_file, string remove_indi_file, string sex_file, int mphen, int mphen2, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool m_grm_flag, bool pred_rand_eff, bool est_fix_eff, int reml_mtd, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, vector<int> drop, bool no_lrt, double prevalence, double prevalence2, bool no_constrain, bool ignore_Ce, vector<double> &fixed_rg_val, bool bivar_no_constrain) {
    _bivar_reml = true;
    _bivar_no_constrain = bivar_no_constrain;
    no_lrt = true;
    _fixed_rg_val = fixed_rg_val;
    _reml_mtd = reml_mtd;
    _reml_max_iter = MaxIter;
    int i = 0, j = 0, k = 0;
    bool grm_flag = (!grm_file.empty());
    bool qcovar_flag = (!qcovar_file.empty());
    bool covar_flag = (!covar_file.empty());
    if (m_grm_flag) grm_flag = false;

    // Read data
    stringstream errmsg;
    int qcovar_num = 0, covar_num = 0;
    vector<string> phen_ID, qcovar_ID, covar_ID, grm_id, grm_files;
    vector< vector<string> > phen_buf, qcovar, covar; // save individuals by column

    if (grm_flag) {
        read_grm(grm_file, grm_id, true, false, !(adj_grm_fac > -1.0));
        update_id_map_kp(grm_id, _id_map, _keep);
        grm_files.push_back(grm_file);
    }
    else if (m_grm_flag) {
        read_grm_filenames(grm_file, grm_files, false);
        for (i = 0; i < grm_files.size(); i++) {
            read_grm(grm_files[i], grm_id, false, true, !(adj_grm_fac > -1.0));
            update_id_map_kp(grm_id, _id_map, _keep);
        }
    }
    read_phen(phen_file, phen_ID, phen_buf, mphen, mphen2);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if (qcovar_flag) {
        qcovar_num = read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if (covar_flag) {
        covar_num = read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if (!remove_indi_file.empty()) remove_indi(remove_indi_file);
    if (grm_flag) {
        if (grm_cutoff>-1.0) rm_cor_indi(grm_cutoff);
        if (!sex_file.empty()) update_sex(sex_file);
        if (adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
        if (dosage_compen>-1) dc(dosage_compen);
        _grm_N.resize(0, 0);
    }

    vector<string> uni_id;
    map<string, int> uni_id_map;
    map<string, int>::iterator iter;
    for (i = 0; i < _keep.size(); i++) {
        uni_id.push_back(_fid[_keep[i]] + ":" + _pid[_keep[i]]);
        uni_id_map.insert(pair<string, int>(_fid[_keep[i]] + ":" + _pid[_keep[i]], i));
    }
    _n = _keep.size();
    if (_n < 1) throw ("Error: no individuals are in common in the input files.");
    cout << _n << " individuals are in common in these files." << endl;

    // construct model terms
    int n1 = 0, n2 = 0;
    vector<string> ystr1(_n), ystr2(_n);
    mphen--;
    mphen2--;
    for (i = 0; i < phen_ID.size(); i++) {
        iter = uni_id_map.find(phen_ID[i]);
        if (iter == uni_id_map.end()) continue;
        if (phen_buf[i][mphen] != "NA" && phen_buf[i][mphen] != "-9") n1++;
        if (phen_buf[i][mphen2] != "NA" && phen_buf[i][mphen2] != "-9") n2++;
        ystr1[iter->second] = phen_buf[i][mphen];
        ystr2[iter->second] = phen_buf[i][mphen2];
    }

    _n = n1 + n2;
    _y = eigenVector::Zero(_n);
    vector<int> nms1, nms2;
    int two_tr_comm = 0;
    for (i = 0, j = 0, k = 0; i < _keep.size(); i++) {
        bool tr1_miss = true;
        if (ystr1[i] != "NA" && ystr1[i] != "-9") {
            (_y.segment(0, n1))(j) = atof(ystr1[i].c_str());
            nms1.push_back(i);
            j++;
            tr1_miss = false;
        }
        bool tr2_miss = true;
        if (ystr2[i] != "NA" && ystr2[i] != "-9") {
            (_y.segment(n1, n2))(k) = atof(ystr2[i].c_str());
            nms2.push_back(i);
            k++;
            tr2_miss = false;
        }
        if (!tr1_miss && !tr2_miss) two_tr_comm++;
    }
    eigenVector y1_tmp = (_y.segment(0, n1)).array() - (_y.segment(0, n1)).mean();
    _y_Ssq = y1_tmp.squaredNorm() / (n1 - 1.0);
    if (!(fabs(_y_Ssq) < 1e30)) throw ("Error: the phenotypic variance for trait 1 is infinite. Please check the missing data in your phenotype file. Missing values should be represented by \"NA\" or \"-9\".");
    eigenVector y2_tmp = (_y.segment(n1, n2)).array() - (_y.segment(n1, n2)).mean();
    _y2_Ssq = y2_tmp.squaredNorm() / (n2 - 1.0);
    if (!(fabs(_y2_Ssq) < 1e30)) throw ("Error: the phenotypic variance for trait 2 is infinite. Please check the missing data in your phenotype file. Missing values should be represented by \"NA\" or \"-9\".");
    cout << nms1.size() << " non-missing phenotypes for trait #1 and " << nms2.size() << " for trait #2" << endl;
    if (!ignore_Ce) {
        if (two_tr_comm == 0) {
            ignore_Ce = true;
            cout << "Note: the residual covariance component is ignored because no individuals were measured for both traits." << endl;
        } else if ((double) two_tr_comm / (double) _keep.size() < 0.1) {
            ignore_Ce = true;
            cout << "Note: the residual covariance component is ignored because < 10% of individuals were measured for both traits." << endl;
        }
    }

    _ncase = 0.0;
    _ncase2 = 0.0;
    eigenVector y1 = _y.segment(0, n1), y2 = _y.segment(n1, n2);
    _flag_CC = check_case_control(_ncase, y1);
    if (_flag_CC) cout << "for trait #1" << endl;
    else prevalence = -1.0;
    _flag_CC2 = check_case_control(_ncase2, y2);
    if (_flag_CC2) cout << "for trait #2" << endl;
    else prevalence2 = -1.0;
    //if(flag_CC2!=_flag_CC) throw("Error: for a bivariate analysis, the two traits should be both quantitative or both binary.");
    if ((_flag_CC && prevalence<-1) || (_flag_CC2 && prevalence2<-1)) cout << "Note: we can specify the disease prevalence by the option --reml-bivar-prevalence so that GCTA can transform the variance explained to the underlying liability scale." << endl;

    int pos = 0;
    _r_indx.clear();
    _bivar_pos.resize(3);
    if (grm_flag) {
        for (i = 0; i < 3 + 3 - ignore_Ce; i++) _r_indx.push_back(i);
        _Asp.resize(_r_indx.size());
        for (i = 0; i < _r_indx.size(); i++) (_Asp[i]).resize(_n, _n);
        if (!no_lrt) drop_comp(drop);
        _bivar_pos[0].push_back(pos);
        for (j = 0; j < n1; j++) {
            (_Asp[pos]).startVec(j);
            for (i = 0; i < n1; i++) (_Asp[pos]).insertBack(i, j) = _grm(_keep[nms1[i]], _keep[nms1[j]]);
        }
        pos++;

        _bivar_pos[1].push_back(pos);
        for (j = 0; j < n2; j++) {
            (_Asp[pos]).startVec(j + n1);
            for (i = 0; i < n2; i++) (_Asp[pos]).insertBack(i + n1, j + n1) = _grm(_keep[nms2[i]], _keep[nms2[j]]);
        }
        pos++;

        _bivar_pos[2].push_back(pos);
        for (j = 0; j < n1; j++) {
            (_Asp[pos]).startVec(j);
            for (i = 0; i < n2; i++) (_Asp[pos]).insertBack(i + n1, j) = _grm(_keep[nms2[i]], _keep[nms1[j]]);
        }
        for (j = 0; j < n2; j++) {
            (_Asp[pos]).startVec(j + n1);
            for (i = 0; i < n1; i++) (_Asp[pos]).insertBack(i, j + n1) = _grm(_keep[nms1[i]], _keep[nms2[j]]);
        }
        pos++;

        for (j = 0; j < pos; j++) (_Asp[j]).finalize();
        _grm.resize(0, 0);
    } 
    else if (m_grm_flag) {
        if (!sex_file.empty()) update_sex(sex_file);
        for (i = 0; i < 3 * grm_files.size() + 3 - ignore_Ce; i++) _r_indx.push_back(i);
        _Asp.resize(_r_indx.size());
        for (i = 0; i < _r_indx.size(); i++) (_Asp[i]).resize(_n, _n);
        if (!no_lrt) drop_comp(drop);
        string prev_file = grm_files[0];
        vector<string> prev_grm_id(grm_id);
        cout << "There are " << grm_files.size() << " GRM file names specified in the file [" + grm_file + "]." << endl;
        vector<int> kp;
        for (k = 0; k < grm_files.size(); k++) {
            cout << "Reading the GRM from the " << k + 1 << "th file ..." << endl;
            read_grm(grm_files[k], grm_id, true, false, !(adj_grm_fac > -1.0));
            if (adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
            if (dosage_compen>-1) dc(dosage_compen);
            StrFunc::match(uni_id, grm_id, kp);
            int pos0 = pos;

            _bivar_pos[0].push_back(pos);
            for (j = 0; j < n1; j++) {
                (_Asp[pos]).startVec(j);
                for (i = 0; i < n1; i++) (_Asp[pos]).insertBack(i, j) = _grm(kp[nms1[i]], _keep[nms1[j]]);
            }
            pos++;

            _bivar_pos[1].push_back(pos);
            for (j = 0; j < n2; j++) {
                (_Asp[pos]).startVec(j + n1);
                for (i = 0; i < n2; i++) (_Asp[pos]).insertBack(i + n1, j + n1) = _grm(kp[nms2[i]], _keep[nms2[j]]);
            }
            pos++;

            _bivar_pos[2].push_back(pos);
            for (j = 0; j < n1; j++) {
                (_Asp[pos]).startVec(j);
                for (i = 0; i < n2; i++) (_Asp[pos]).insertBack(i + n1, j) = _grm(kp[nms2[i]], _keep[nms1[j]]);
            }
            for (j = 0; j < n2; j++) {
                (_Asp[pos]).startVec(j + n1);
                for (i = 0; i < n1; i++) (_Asp[pos]).insertBack(i, j + n1) = _grm(kp[nms1[i]], _keep[nms2[j]]);
            }
            pos++;

            for (j = pos0; j < pos; j++) (_Asp[j]).finalize();
            prev_file = grm_files[k];
            prev_grm_id = grm_id;
        }
        _grm_N.resize(0, 0);
        _grm.resize(0, 0);
    }

    _bivar_pos[0].push_back(pos);
    for (i = 0; i < n1; i++) {
        (_Asp[pos]).startVec(i);
        (_Asp[pos]).insertBack(i, i) = 1.0;
    }
    (_Asp[pos]).finalize();
    pos++;

    _bivar_pos[1].push_back(pos);
    for (i = 0; i < n2; i++) {
        (_Asp[pos]).startVec(i + n1);
        (_Asp[pos]).insertBack(i + n1, i + n1) = 1.0;
    }
    (_Asp[pos]).finalize();
    pos++;

    if (!ignore_Ce) {
        _bivar_pos[2].push_back(pos);
        for (j = 0; j < n1; j++) {
            (_Asp[pos]).startVec(j);
            for (i = 0; i < n2; i++) {
                if (nms2[i] == nms1[j]) (_Asp[pos]).insertBack(i + n1, j) = 1;
            }
        }
        for (j = 0; j < n2; j++) {
            (_Asp[pos]).startVec(j + n1);
            for (i = 0; i < n1; i++) {
                if (nms1[i] == nms2[j]) (_Asp[pos]).insertBack(i, j + n1) = 1;
            }
        }
        pos++;
    }

    // construct X matrix
    vector<eigenMatrix> E_float;
    eigenMatrix qE_float;
    construct_X(_keep.size(), uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);
    eigenMatrix X(_X);
    _X = eigenMatrix::Zero(_n, _X_c * 2);
    for (i = 0; i < n1; i++) (_X.block(0, 0, n1, _X_c)).row(i) = X.row(nms1[i]);
    for (i = 0; i < n2; i++) (_X.block(n1, _X_c, n2, _X_c)).row(i) = X.row(nms2[i]);
    _X_c *= 2;

    // names of variance component
    for (i = 0; i < grm_files.size(); i++) {
        stringstream strstrm;
        if (grm_files.size() == 1) strstrm << "";
        else strstrm << i + 1;
        _var_name.push_back("V(G" + strstrm.str() + ")_tr1");
        _var_name.push_back("V(G" + strstrm.str() + ")_tr2");
        _var_name.push_back("C(G" + strstrm.str() + ")_tr12");
        _hsq_name.push_back("V(G" + strstrm.str() + ")/Vp_tr1");
        _hsq_name.push_back("V(G" + strstrm.str() + ")/Vp_tr2");
    }
    _var_name.push_back("V(e)_tr1");
    _var_name.push_back("V(e)_tr2");
    if (!ignore_Ce) _var_name.push_back("C(e)_tr12");
    _ignore_Ce = ignore_Ce;

    // run REML algorithm
    reml(pred_rand_eff, est_fix_eff, reml_priors, reml_priors_var, prevalence, prevalence2, no_constrain, no_lrt, false);
}

bool gcta::calcu_Vi_bivar(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter) {
    int i = 0, n = Vi.cols();
    double d_buf = 0.0;
    logdet = 0.0;
    string errmsg = "\nError: the V (variance-covariance) matrix is not invertible.";

    Vi = eigenMatrix::Zero(_n, _n);
    for (i = 0; i < _r_indx.size(); i++) Vi += (_Asp[_r_indx[i]]) * prev_varcmp[i];

    if (_V_inv_mtd == 0) {
        if (!comput_inverse_logdet_LDLT_mkl(Vi, logdet)) {
            //cout<<"Note: the variance-covaraince matrix V is non-positive definite. Switching to Cholesky to LU decomposition approach."<<endl;
            _V_inv_mtd = 1;
        }
    }
    if (_V_inv_mtd == 1) {
        if (!comput_inverse_logdet_LU_mkl(Vi, logdet)) throw ("Error: the variance-covaraince matrix V is not invertible.");
    }

    /*    if(!comput_inverse_logdet_LDLT_mkl(Vi, logdet)){
            if(_reml_have_bend_A) throw(errmsg);
            cout<<"Warning: the variance-covaraince matrix V is negative-definite."<<endl;
            bend_A();
            _reml_have_bend_A=true;
            iter=-1;
            cout<<"Restarting iterations ..."<<endl;
            return false;
        }
     */

}

void gcta::constrain_rg(eigenVector &varcmp) {
    static int count = 0;
    int v_pos = 0, c_pos = 0, d = _bivar_pos[0].size() - _bivar_pos[2].size();
    if (_ignore_Ce) d--;
    eigenMatrix G(2, 2);

    for (c_pos = _bivar_pos[2].size() - 1; c_pos >= 0; c_pos--) {
        v_pos = c_pos + d;
        G(0, 0) = varcmp[_bivar_pos[0][v_pos]];
        G(1, 1) = varcmp[_bivar_pos[1][v_pos]];
        G(0, 1) = G(1, 0) = varcmp[_bivar_pos[2][c_pos]];

        SelfAdjointEigenSolver<eigenMatrix> eigensolver(G);
        eigenVector eval = eigensolver.eigenvalues();
        if (eval.minCoeff() <= 0.0) {
            if (count == 0) {
                cout << "Note: to constrain the correlation being from -1 to 1, a genetic (or residual) variance-covariance matrix is bended to be positive definite. In this case, the SE is not reliable." << endl;
                count++;
            }
            bending_eigenval(eval);
            G = eigensolver.eigenvectors() * eigenDiagMat(eval) * eigensolver.eigenvectors().transpose();
            varcmp[_bivar_pos[0][v_pos]] = G(0, 0);
            varcmp[_bivar_pos[1][v_pos]] = G(1, 1);
            varcmp[_bivar_pos[2][c_pos]] = G(0, 1);
        }
    }
}

void gcta::calcu_rg(eigenVector &varcmp, eigenMatrix &Hi, eigenVector &rg, eigenVector &rg_var, vector<string> &rg_name) {
    int i = 0, j = 0;
    double V1 = 0, V2 = 0, C = 0, VarV1 = 0, VarV2 = 0, VarC = 0.0, CovV1V2 = 0.0, CovV1C = 0.0, CovV2C = 0.0;

    rg = eigenVector::Zero(_bivar_pos[0].size() - 1);
    rg_var = rg;
    for (i = 0; i < _bivar_pos[0].size() - 1; i++) {
        V1 = varcmp[_bivar_pos[0][i]];
        V2 = varcmp[_bivar_pos[1][i]];
        C = varcmp[_bivar_pos[2][i]];
        VarV1 = Hi(_bivar_pos[0][i], _bivar_pos[0][i]);
        VarV2 = Hi(_bivar_pos[1][i], _bivar_pos[1][i]);
        VarC = Hi(_bivar_pos[2][i], _bivar_pos[2][i]);
        CovV1V2 = Hi(_bivar_pos[0][i], _bivar_pos[1][i]);
        CovV1C = Hi(_bivar_pos[0][i], _bivar_pos[2][i]);
        CovV2C = Hi(_bivar_pos[1][i], _bivar_pos[2][i]);
        if (V1 * V2 > 0) {
            rg[i] = sqrt(V1 * V2);
            if (rg[i] > 0) rg[i] = C / rg[i];
            rg_var[i] = rg[i] * rg[i]*(VarV1 / (4 * V1 * V1) + VarV2 / (4 * V2 * V2) + VarC / (C * C) + CovV1V2 / (2 * V1 * V2) - CovV1C / (V1 * C) - CovV2C / (V2 * C));
        }

        if (_bivar_pos[0].size() == 2) rg_name.push_back("rG");
        else {
            stringstream strstrm;
            strstrm << "rG" << i + 1;
            rg_name.push_back(strstrm.str());
        }
    }
}

double gcta::lgL_fix_rg(eigenVector &prev_varcmp, bool no_constrain) {
    int i = 0, j = 0;
    if (_fixed_rg_val.size() > _bivar_pos[0].size() - 1) {
        vector<double> rg_val_buf(_fixed_rg_val);
        _fixed_rg_val.clear();
        for (i = 0; i < _bivar_pos[0].size() - 1; i++) _fixed_rg_val.push_back(rg_val_buf[i]);
    }

    cout << "\nCalculating the logLikelihood for the model with the genetic correlation" << (_fixed_rg_val.size() > 1 ? "s" : "") << " being fixed at ";
    for (int i = 0; i < _fixed_rg_val.size() - 1; i++) cout << _fixed_rg_val[i] << "\t";
    cout << _fixed_rg_val[_fixed_rg_val.size() - 1] << endl;

    vector<int> r_indx_buf(_r_indx);
    _bivar_pos_prev = _bivar_pos;
    for (i = _fixed_rg_val.size() - 1; i >= 0; i--) {
        for (j = i + 1; j < _bivar_pos_prev[0].size(); j++) _bivar_pos[0][j]--;
        for (j = i + 1; j < _bivar_pos_prev[1].size(); j++) _bivar_pos[1][j]--;
        for (j = i + 1; j < _bivar_pos_prev[2].size(); j++) _bivar_pos[2][j]--;
        _r_indx.erase(_r_indx.begin() + _bivar_pos[2][i]);
        _bivar_pos[2].erase(_bivar_pos[2].begin() + i);
    }

    _Asp_prev = _Asp;
    eigenMatrix Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(_r_indx.size(), _r_indx.size());
    eigenVector Py(_n);
    eigenVector varcmp(_r_indx.size());
    for (i = 0; i < _r_indx.size(); i++) {
        varcmp[i] = fabs(prev_varcmp[_r_indx[i]]);
        if (varcmp[i] < 1.0e-30) varcmp[i] = 0.1;
    }
    double lgL = reml_iteration(Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, false, no_constrain, true);
    _r_indx = r_indx_buf;
    _bivar_pos = _bivar_pos_prev;

    return lgL;
}

void gcta::update_A(eigenVector &prev_varcmp) {
    int i = 0;
    double g1 = 0.0, g2 = 0.0;
    for (i = 0; i < _fixed_rg_val.size(); i++) {
        g1 = prev_varcmp[_bivar_pos[0][i]];
        g2 = prev_varcmp[_bivar_pos[1][i]];
        _Asp[_bivar_pos_prev[0][i]] = _Asp_prev[_bivar_pos_prev[0][i]] + _Asp_prev[_bivar_pos_prev[2][i]]*(0.5 * _fixed_rg_val[i] * sqrt(g2 / g1));
        _Asp[_bivar_pos_prev[1][i]] = _Asp_prev[_bivar_pos_prev[1][i]] + _Asp_prev[_bivar_pos_prev[2][i]]*(0.5 * _fixed_rg_val[i] * sqrt(g1 / g2));
    }

}

/*
void gcta::init_rg(eigenVector &varcmp)
{
    int i=0;
    double d_buf=0.0, ratio=0.0, delta_0=0.0, delta_1=0.0;
    eigenVector varcmp_buf=varcmp;
    
    for(i=0; i<_fixed_rg_val.size(); i++){
        d_buf=varcmp[_bivar_pos[0][i]]*varcmp[_bivar_pos[1][i]];
        if(d_buf>0.0){
            ratio=sqrt(abs(_fixed_rg_val[i]/((varcmp[_bivar_pos[2][i]]/sqrt(d_buf)))));
            if(d_buf*_fixed_rg_val[i]>0.0) varcmp[_bivar_pos[2][i]]*=ratio;
            else varcmp[_bivar_pos[2][i]]*=(-1.0*ratio);
            varcmp[_bivar_pos[0][i]]/=ratio;
            varcmp[_bivar_pos[1][i]]/=ratio;
            delta_0=varcmp_buf[_bivar_pos[0][i]]-varcmp[_bivar_pos[0][i]];
        }
        else throw("Error: unable to calcuate the genetic correlation because one of the genetic variance components is negative.");
    }
}

void gcta::fix_rg(eigenVector &varcmp, int pos)
{
    int i=0;
    double d_buf=0.0;
    string errmsg="Error: unable to calcuate the genetic correlation because one of the genetic variance components is zero.";
    
    for(i=0; i<_fixed_rg_val.size(); i++){
        if(pos==2){
            d_buf=varcmp[_bivar_pos[0][i]]*varcmp[_bivar_pos[1][i]];
            if(CommFunc::FloatNotEqual(d_buf, 0.0)) varcmp[_bivar_pos[2][i]]=_fixed_rg_val[i]*sqrt(d_buf);
            else throw(errmsg);
        }
        else if(pos==1){
            if(CommFunc::FloatNotEqual(varcmp[_bivar_pos[1][i]], 0.0)) varcmp[_bivar_pos[0][i]]=varcmp[_bivar_pos[2][i]]*varcmp[_bivar_pos[2][i]]/(varcmp[_bivar_pos[1][i]]*_fixed_rg_val[i]*_fixed_rg_val[i]);
            else throw(errmsg);
        }
        else if(pos==0){
            if(CommFunc::FloatNotEqual(varcmp[_bivar_pos[0][i]], 0.0)) varcmp[_bivar_pos[1][i]]=varcmp[_bivar_pos[2][i]]*varcmp[_bivar_pos[2][i]]/(varcmp[_bivar_pos[0][i]]*_fixed_rg_val[i]*_fixed_rg_val[i]);
            else throw(errmsg);
        }
    }
}
 */