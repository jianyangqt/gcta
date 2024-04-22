/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Interface to all the GCTA functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file COPYING for more
 * details
 */

#ifndef _GCTA_H
#define _GCTA_H

#ifndef EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif

#include "cpu.h"
#include <cstdio>
#include "CommFunc.h"
#include "StrFunc.h"
#include "StatFunc.h"
#include "eigen_func.h"
#include <fstream>
#include <iomanip>
#include <bitset>
#include <map>
#include <Eigen/StdVector>
#include "zfstream.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/IterativeSolvers>
#include <omp.h>
#include "Logger.h"
#include "Matrix.hpp"

#ifdef SINGLE_PRECISION
typedef Eigen::SparseMatrix<float, Eigen::ColMajor, long long> eigenSparseMat;
#else
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long long> eigenSparseMat;
#endif
//To avoid potential alignment problem. 
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(eigenSparseMat);
using namespace Eigen;
using namespace std;

#ifdef SINGLE_PRECISION
typedef DiagonalMatrix<float, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;
typedef DynamicSparseMatrix<float> eigenDynSparseMat;
#else
typedef DiagonalMatrix<double, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;
typedef DynamicSparseMatrix<double> eigenDynSparseMat;
#endif

class gcta {
public:
    gcta(int autosome_num, double rm_ld_cutoff, string out);
    gcta();
    virtual ~gcta();

    void read_famfile(string famfile);
    void read_bimfile(string bimfile);
    void read_bedfile(string bedfile);
    vector<string> read_bfile_list(string bfile_list);
    void read_multi_famfiles(vector<string> multi_bfiles);
    void read_multi_bimfiles(vector<string> multi_bfiles);
    void read_multi_bedfiles(vector<string> multi_bfiles);
    void read_imp_info_mach_gz(string zinfofile);
    void read_imp_info_mach(string infofile);
    void read_imp_dose_mach_gz(string zdosefile, string kp_indi_file, string rm_indi_file, string blup_indi_file);
    void read_imp_dose_mach(string dosefile, string kp_indi_file, string rm_indi_file, string blup_indi_file);
    void read_imp_info_beagle(string zinfofile);
    void read_imp_dose_beagle(string zdosefile, string kp_indi_file, string rm_indi_file, string blup_indi_file);
    void update_ref_A(string ref_A_file);
    void update_impRsq(string zinfofile);
    void update_freq(string freq);
    void save_freq(bool ssq_flag);
    void extract_snp(string snplistfile);
    void extract_single_snp(string snpname);
    void extract_region_snp(string snpname, int wind_size);
    void extract_region_bp(int chr, int bp, int wind_size);
    void exclude_snp(string snplistfile);
    void exclude_single_snp(string snpname);
    void exclude_region_snp(string snpname, int wind_size);
    void exclude_region_bp(int chr, int bp, int wind_size);
    void extract_chr(int chr_start, int chr_end);
    void filter_snp_maf(double maf);
    void filter_snp_max_maf(double max_maf);
    void filter_impRsq(double rsq_cutoff);
    void keep_indi(string indi_list_file);
    void remove_indi(string indi_list_file);
    void update_sex(string sex_file);
    void read_indi_blup(string blup_indi_file);
    void save_XMat(bool miss_with_mu, bool std);

    void make_grm(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, bool mlmassoc, bool diag_f3_flag = false, string subpopu_file = "");
    //void make_grm_pca(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, double wind_size, bool mlmassoc);
    void save_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool output_grm_bin);
    void align_grm(string m_grm_file);
    void pca(string grm_file, string keep_indi_file, string remove_indi_file, double grm_cutoff, bool merge_grm_flag, int out_pc_num);
    void snp_pc_loading(string pc_file);
    void project_loading(string pc_load, int N); 

    // bigK + smallK method
    void grm_bK(string grm_file, string keep_indi_file, string remove_indi_file, double threshold, bool grm_out_bin_flag);

    void enable_grm_bin_flag();
    void fit_reml(string grm_file, string phen_file, string qcovar_file, string covar_file, string qGE_file, string GE_file, string keep_indi_file, string remove_indi_file, string sex_file, int mphen, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool m_grm_flag, bool pred_rand_eff, bool est_fix_eff, int reml_mtd, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, vector<int> drop, bool no_lrt, double prevalence, bool no_constrain, bool mlmassoc = false, bool within_family = false, bool reml_bending = false, bool reml_diag_one = false, string weight_file = "");
    //void HE_reg(string grm_file, string phen_file, string keep_indi_file, string remove_indi_file, int mphen); // old HE regression method
    void HE_reg(string grm_file, bool m_grm_flag, string phen_file, string keep_indi_file, string remove_indi_file, int mphen); // allow multiple regression
    void HE_reg_bivar(string grm_file, bool m_grm_flag, string phen_file, string keep_indi_file, string remove_indi_file, int mphen, int mphen2); // estimate genetic covariance between two traits
    void blup_snp_geno();
    void blup_snp_dosage();
    void set_reml_force_inv();
    void set_reml_force_converge();
    void set_cv_blup(bool cv_blup);
    void set_reml_no_converge();
    void set_reml_fixed_var();
    void set_reml_mtd(int reml_mtd);
    void set_reml_allow_constrain_run();
    void set_reml_diag_mul(double value);
    void set_reml_diagV_adj(int method);
    void set_reml_inv_method(int method);


    // bivariate REML analysis
    void fit_bivar_reml(string grm_file, string phen_file, string qcovar_file, string covar_file, string keep_indi_file, string remove_indi_file, string sex_file, int mphen, int mphen2, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool m_grm_flag, bool pred_rand_eff, bool est_fix_eff, int reml_mtd, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, vector<int> drop, bool no_lrt, double prevalence, double prevalence2, bool no_constrain, bool ignore_Ce, vector<double> &fixed_rg_val, bool bivar_no_constrain);

    // LD
    void read_LD_target_SNPs(string snplistfile);
    void LD_Blocks(int stp, double wind_size, double alpha, bool IncldQ = true, bool save_ram = false);
    void calcu_mean_rsq(int wind_size, double rsq_cutoff, bool dominance_flag);
    void calcu_mean_rsq_multiSet(string snpset_filenames_file, int wind_size, double rsq_cutoff, bool dominance_flag);
    void calcu_max_ld_rsq(int wind_size, double rsq_cutoff, bool dominance_flag);
    void ld_seg(string i_ld_file, int seg_size, int wind_size, double rsq_cutoff, bool dominance_flag);
    void set_ldscore_adj_flag(bool ldscore_adj);

    void genet_dst(string bfile, string hapmap_genet_map);

    void GWAS_simu(string bfile, int simu_num, string qtl_file, int case_num, int control_num, double hsq, double K, int seed, bool output_causal, bool simu_emb_flag, int eff_mod=0);
    //void simu_geno_unlinked(int N, int M, double maf);

    void set_diff_freq(double freq_diff);
    void run_massoc_slct(string metafile, int wind_size, double p_cutoff, double collinear, int64_t top_SNPs, bool joint_only, bool GC, double GC_val, bool actual_geno, int mld_slct_alg);
    void run_massoc_cond(string metafile, string snplistfile, int wind_size, double collinear, bool GC, double GC_val, bool actual_geno);
    void run_massoc_sblup(string metafile, int wind_size, double lambda);
    void set_massoc_pC_thresh(double thresh);

    void save_plink();
    void dose2bed();

    void read_IRG_fnames(string snp_info_file, string fname_file, double GC_cutoff);

    // population genetics
    void Fst(string filename);
    void paa(string aa_file);
    void ibc(bool ibc_all);

    // mkl
    void make_grm_mkl(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, bool mlmassoc, bool diag_f3_flag = false);
    void calcu_mean_rsq_mkl(int wind_size, double rsq_cutoff);
    void LD_pruning_mkl(double rsq_cutoff, int wind_size);
    //void make_grm_wt_mkl(string i_ld_file, int wind_m, double wt_ld_cut, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd=0, int wt_mtd=0, bool mlmassoc=false, bool impData_flag=false, int ttl_snp_num=-1);

    // mlma
    void mlma(string grm_file, bool m_grm_flag, string subtract_grm_file, string phen_file, string qcovar_file, string covar_file, int mphen, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, bool no_constrain, bool within_family, bool inbred, bool no_adj_covar);
    void mlma_loco(string phen_file, string qcovar_file, string covar_file, int mphen, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, bool no_constrain, bool inbred, bool no_adj_covar);

    // gene based association test
    void sbat_gene(string sAssoc_file, string gAnno_file, int sbat_wind, double sbat_ld_cutoff, bool sbat_write_snpset, bool GC, double GC_val);
    void sbat(string sAssoc_file, string snpset_file, double sbat_ld_cutoff, bool sbat_write_snpset,bool GC, double GC_val);
    void sbat_seg(string sAssoc_file, int seg_size, double sbat_ld_cutoff, bool sbat_write_snpset,bool GC, double GC_val);

    // gene based association
    //////////////////////////////
    void mbat_gene(string mbat_sAssoc_file, string mbat_gAnno_file, int mbat_wind, double mbat_svd_gamma, double sbat_ld_cutoff, bool mbat_write_snpset, bool GC, double GC_val,bool mbat_print_all_p);
    void svdDecomposition( MatrixXf &X,double &prop, int &eigenvalueNum, VectorXd &eigenvalueUsed,MatrixXd &U_prop);
    void mbat_ACATO(double &mbat_svd_pvalue,double &fastbat_pvalue, double &mbat_pvalue);
    void mbat_calcu_lambda(vector<int> &snp_indx, MatrixXf &rval, VectorXd &eigenval, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx);
    void mbat(string mbat_sAssoc_file, string snpset_file, double mbat_svd_gamma, double sbat_ld_cutoff, bool mbat_write_snpset, bool GC, double GC_val,bool mbat_print_all_p); 

    ////////////////////////////////
    // GSMR
    void read_gsmrfile(string expo_file_list, string outcome_file_list, double gwas_thresh, int nsnp_gsmr, int gsmr_so_alg);
    void gsmr(int gsmr_alg_flag, string ref_ld_dirt, string w_ld_dirt, double freq_thresh, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, bool o_snp_instru_flag, int gsmr_so_alg, int gsmr_beta_version);
    vector<vector<double>> forward_gsmr(stringstream &ss, map<string,int> &snp_instru_map, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, stringstream &ss_pleio);
    vector<vector<double>> reverse_gsmr(stringstream &ss, map<string,int> &snp_instru_map, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, stringstream &ss_pleio);
    eigenMatrix rho_sample_overlap(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, int nexpo, int noutcome, vector<string> snp_name, vector<int> snp_remain, string ref_ld_dirt, string w_ld_dirt, vector<string> trait_name, int gsmr_so_alg);

    eigenMatrix sample_overlap_ldsc(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, int nexpo, int noutcome, vector<string> snp_name, vector<int> snp_remain, string ref_ld_dirt, string w_ld_dirt, vector<string> trait_name);
    eigenMatrix sample_overlap_rb(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, int nexpo, int noutcome, vector<string> snp_name, vector<int> snp_remain, vector<string> trait_name);

    // mtCOJO
    void mtcojo(string mtcojo_bxy_file, string ref_ld_dirt, string w_ld_dirt, double freq_thresh, double gwas_thresh, int clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, int gsmr_beta_version);
    bool mtcojo_ldsc(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, int ntrait, vector<string> snp_name, vector<int> snp_remain, string ref_ld_dirt, string w_ld_dirt, vector<string> trait_name, eigenMatrix &ldsc_intercept, eigenMatrix &ldsc_slope);
    int read_mtcojofile(string mtcojolist_file, double gwas_thresh, int nsnp_gsmr);
    double read_single_metafile_txt(string metafile, map<string, int> id_map, vector<string> &snp_a1, vector<string> &snp_a2, 
                                eigenVector &snp_freq, eigenVector &snp_b, eigenVector &snp_se, eigenVector &snp_pval, eigenVector &snp_n, vector<bool> &snpflag);
    double read_single_metafile_gz(string metafile, map<string, int> id_map, vector<string> &snp_a1, vector<string> &snp_a2, 
                                eigenVector &snp_freq, eigenVector &snp_b, eigenVector &snp_se, eigenVector &snp_pval, eigenVector &snp_n, vector<bool> &snpflag);
    vector<string> read_snp_metafile_txt(string metafile, map<string,int> &gws_snp_name_map, double thresh);
    vector<string> read_snp_metafile_gz(string metafile, map<string,int> &gws_snp_name_map, double thresh);
    
    
    // Adjusted for PC
    void pc_adjust(string pcadjust_list_file, string eigenvalue_file, double freq_thresh, int wind_size);
    void read_pc_adjust_file(string pcadjust_list_file, string pc_file);

    /////////////////////////
    // gene expresion data
    void read_efile(string efile);

    // ecojo
    void read_eR(string eR_file);
    void run_ecojo_slct(string e_metafile, double p_cutoff, double collinear);
    void run_ecojo_blup_efile(string e_metafile, double lambda);
    void run_ecojo_blup_eR(string e_metafile, double lambda);

    // ERM
    void make_erm(int erm_mtd, bool output_bin); 

    // vif threshold
    void set_vif_threshold(double value);



private:
    void init_keep();
    void init_include();
    void get_rsnp(vector<int> &rsnp);
    void get_rindi(vector<int> &rindi);

    void save_famfile();
    void save_bimfile();
    void save_bedfile();

    void update_bim(vector<int> &rsnp);
    void update_fam(vector<int> &rindi);

    void update_include(vector<int> chr_buf, vector<string> snpid_buf, vector<double> gd_buf, vector<int> bp_buf, vector<string> a1_buf, vector<string> a2_buf, int file_indx);
    void update_keep(vector<string> fid_buf, vector<string> pid_buf, vector<string> fa_id_buf, vector<string> mo_id_buf, vector<int> sex_buf, vector<double> pheno_buf, string famfile);
    
    void update_id_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep);
    void update_id_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep);
    void read_snplist(string snplistfile, vector<string> &snplist, string msg = "SNPs");
    void read_indi_list(string indi_list_file, vector<string> &indi_list);

    bool make_XMat(MatrixXf &X);
    bool make_XMat_d(MatrixXf &X);
    //void make_XMat_SNPs(vector< vector<float> > &X, bool miss_with_mu);
    void std_XMat(MatrixXf &X, eigenVector &sd_SNP, bool grm_xchr_flag, bool miss_with_mu, bool divid_by_std);
    void std_XMat_subpopu(string subpopu_file, MatrixXf &X, eigenVector &sd_SNP, bool grm_xchr_flag, bool miss_with_mu, bool divid_by_std);
    void std_XMat_d(MatrixXf &X, eigenVector &sd_SNP, bool miss_with_mu, bool divid_by_std);
    //void std_XMat(vector< vector<float> > &X, vector<double> &sd_SNP, bool grm_xchr_flag, bool divid_by_std = true);
    void makex_eigenVector(int j, eigenVector &x, bool resize = true, bool minus_2p = false);
    void makex_eigenVector_std(int j, eigenVector &x, bool resize = true, double snp_std = 1.0);
    //void make_XMat_eigenMatrix(MatrixXf &X);
    bool make_XMat_subset(MatrixXf &X, vector<int> &snp_indx, bool divid_by_std);
    bool make_XMat_d_subset(MatrixXf &X, vector<int> &snp_indx, bool divid_by_std);


    void calcu_mu(bool ssq_flag = false);
    void calcu_maf();
    void mu_func(int j, vector<double> &fac);
    void check_autosome();
    void check_chrX();
    void check_sex();

    // grm
    void calcu_grm_var(double &diag_m, double &diag_v, double &off_m, double &off_v);
    int read_grm_id(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only);
    void read_grm(string grm_file, vector<string> &grm_id, bool out_id_log = true, bool read_id_only = false, bool dont_read_N = false);
    void read_grm_gz(string grm_file, vector<string> &grm_id, bool out_id_log = true, bool read_id_only = false);
    void read_grm_bin(string grm_file, vector<string> &grm_id, bool out_id_log = true, bool read_id_only = false, bool dont_read_N = false);
    void read_grm_filenames(string merge_grm_file, vector<string> &grm_files, bool out_log = true);
    void merge_grm(string merge_grm_file);
    void rm_cor_indi(double grm_cutoff);
    void adj_grm(double adj_grm_fac);
    void dc(int dosage_compen);
    void manipulate_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool dont_read_N = false);
    void output_grm_vec(vector< vector<float> > &A, vector< vector<int> > &A_N, bool output_grm_bin);
    void output_grm(bool output_grm_bin);

    // reml
    void read_phen(string phen_file, vector<string> &phen_ID, vector< vector<string> > &phen_buf, int mphen, int mphen2 = 0);
    int read_fac(ifstream &ifstrm, vector<string> &ID, vector< vector<string> > &fac);
    int read_covar(string covar_file, vector<string> &covar_ID, vector< vector<string> > &covar, bool qcovar_flag);
    int read_GE(string GE_file, vector<string> &GE_ID, vector< vector<string> > &GE, bool qGE_flag = false);
    bool check_case_control(double &ncase, eigenVector &y);
    double transform_hsq_L(double P, double K, double hsq);
    int constrain_varcmp(eigenVector &varcmp);
    void drop_comp(vector<int> &drop);
    void construct_X(int n, map<string, int> &uni_id_map, bool qcovar_flag, int qcovar_num, vector<string> &qcovar_ID, vector< vector<string> > &qcovar, bool covar_flag, int covar_num, vector<string> &covar_ID, vector< vector<string> > &covar, vector<eigenMatrix> &E_float, eigenMatrix &qE_float);
    void coeff_mat(const vector<string> &vec, eigenMatrix &coeff_mat, string errmsg1, string errmsg2);
    void reml(bool pred_rand_eff, bool est_fix_eff, vector<double> &reml_priors, vector<double> &reml_priors_var, double prevalence, double prevalence2, bool no_constrain, bool no_lrt, bool mlmassoc = false);
    double reml_iteration(eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp, bool prior_var_flag, bool no_constrain, bool reml_bivar_fix_rg = false);
    void init_varcomp(vector<double> &reml_priors_var, vector<double> &reml_priors, eigenVector &varcmp);
    bool calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter);
    bool inverse_H(eigenMatrix &H);
    bool comput_inverse_logdet_LDLT(eigenMatrix &Vi, double &logdet);
    void bend_V(eigenMatrix &Vi);
    void bend_A();
    bool bending_eigenval(eigenVector &eval);
    bool comput_inverse_logdet_PLU(eigenMatrix &Vi, double &logdet);
    bool comput_inverse_logdet_LU(eigenMatrix &Vi, double &logdet);
    double calcu_P(eigenMatrix &Vi, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &P);
    void calcu_Hi(eigenMatrix &P, eigenMatrix &Hi);
    void reml_equation(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp);
    double lgL_reduce_mdl(bool no_constrain);
    void em_reml(eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp);
    void ai_reml(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, double dlogL);
    void calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA);
    void calcu_Vp(double &Vp, double &Vp2, double &VarVp, double &VarVp2, eigenVector &varcmp, eigenMatrix &Hi);
    void calcu_hsq(int i, double Vp, double Vp2, double VarVp, double VarVp2, double &hsq, double &var_hsq, eigenVector &varcmp, eigenMatrix &Hi);
    void calcu_sum_hsq(double Vp, double VarVp, double &sum_hsq, double &var_sum_hsq, eigenVector &varcmp, eigenMatrix &Hi);
    void output_blup_snp(eigenMatrix &b_SNP);

    // within-family reml analysis
    void detect_family();
    bool calcu_Vi_within_family(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter);

    // bivariate REML analysis
    void calcu_rg(eigenVector &varcmp, eigenMatrix &Hi, eigenVector &rg, eigenVector &rg_var, vector<string> &rg_name);
    void update_A(eigenVector &prev_varcmp);
    void constrain_rg(eigenVector &varcmp);
    double lgL_fix_rg(eigenVector &prev_varcmp, bool no_constrain);
    bool calcu_Vi_bivar(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter);

    // GWAS simulation
    void kosambi();
    int read_QTL_file(string qtl_file, vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, vector<int> &have_eff);
    void output_simu_par(vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, double Vp);
    void save_phenfile(vector< vector<double> > &y);
    // not usefully any more
    void GenerCases(string bfile, string qtl_file, int case_num, int control_num, double hsq, double K, bool curr_popu = false, double gnrt = 100);

    // LD
    void EstLD(vector<int> &smpl, double wind_size, vector< vector<string> > &snp, vector< vector<double> > &r, vector<double> &r2, vector<double> &md_r2, vector<double> &max_r2, vector<string> &max_r2_snp, vector<double> &dL, vector<double> &dR, vector<int> &K, vector<string> &L_SNP, vector<string> &R_SNP, double alpha, bool IncldQ);
    eigenMatrix reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table = false);
    void rm_cor_snp(int m, int start, float *rsq, double rsq_cutoff, vector<int> &rm_snp_ID1);
    void get_ld_blk_pnt(vector<int> &brk_pnt1, vector<int> &brk_pnt2, vector<int> &brk_pnt3, int wind_bp, int wind_snp = 0);
    void get_ld_blk_pnt_max_limit(vector<int> &brk_pnt1, vector<int> &brk_pnt2, vector<int> &brk_pnt3, int wind_bp, int wind_snp);
    void get_lds_brkpnt(vector<int> &brk_pnt1, vector<int> &brk_pnt2, int ldwt_seg, int wind_snp_num=0);
    void calcu_ld_blk(vector<int> &brk_pnt, vector<int> &brk_pnt3, eigenVector &mean_rsq, eigenVector &snp_num, eigenVector &max_rsq, bool second, double rsq_cutoff, bool dominance_flag = false);
    void calcu_ld_blk_split(int size, int size_limit, MatrixXf &X_sub, eigenVector &ssx_sqrt_i_sub, double rsq_cutoff, eigenVector &rsq_size, eigenVector &mean_rsq_sub, eigenVector &max_rsq_sub, int s1, int s2, bool second);
    void calcu_ssx_sqrt_i(eigenVector &ssx_sqrt_i);
    void calcu_max_ld_rsq_blk(eigenVector &multi_rsq, eigenVector &multi_rsq_adj, eigenVector &max_rsq, vector<int> &max_pos, vector<int> &brk_pnt, double rsq_cutoff, bool dominance_flag);
    bool bending_eigenval_Xf(VectorXf &eval);
    void calcu_ld_blk_multiSet(vector<int> &brk_pnt, vector<int> &brk_pnt3, vector< vector<bool> > &set_flag, vector<eigenVector> &mean_rsq, vector<eigenVector> &snp_num, vector<eigenVector> &max_rsq, bool second, double rsq_cutoff, bool dominance_flag);
    void calcu_ld_blk_split_multiSet(int size, int size_limit, MatrixXf &X_sub, MatrixXf &X_sub2, vector<int> &used_in_this_set, eigenVector &ssx_sqrt_i_sub, double rsq_cutoff, eigenVector &rsq_size, eigenVector &mean_rsq_sub, eigenVector &max_rsq_sub, int s1, int s2, bool second);

    // Joint analysis of GWAS MA results
    void read_metafile(string metafile, bool GC, double GC_val);
    void init_massoc(string metafile, bool GC, double GC_val);
    void read_fixed_snp(string snplistfile, string msg, vector<int> &pgiven, vector<int> &remain);
    void eigenVector2Vector(eigenVector &x, vector<double> &y);
    //double crossprod(int indx1, int indx2);
    void get_x_vec(float *x, int indx);
    void get_x_mat(float *x, vector<int> &indx);
    void vec_t_mat(float *vec, int nrow, float *mat, int ncol, eigenVector &out); // 1 x n vector multiplied by m x n matrix
    void stepwise_slct(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, int mld_slct_alg, uint64_t top_SNPs);
    bool slct_entry(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC);
    void slct_stay(vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ);
    double massoc_calcu_Ve(const vector<int> &slct, eigenVector &bJ, eigenVector &b);
    void massoc_cond(const vector<int> &slct, const vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC);
    void massoc_joint(const vector<int> &indx, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ);
    bool init_B(const vector<int> &indx);
    void init_Z(const vector<int> &indx);
    bool insert_B_and_Z(const vector<int> &indx, int insert_indx);
    void erase_B_and_Z(const vector<int> &indx, int erase_indx);
    void LD_rval(const vector<int> &indx, eigenMatrix &rval);
    bool massoc_sblup(double lambda, eigenVector &bJ);
    void massoc_slct_output(bool joint_only, vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, eigenMatrix &rval);
    void massoc_cond_output(vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC);

    // read raw genotype data (Illumina)
    char flip_allele(char a);
    void read_one_IRG(ofstream &oped, int ind, string IRG_fname, double GC_cutoff);

    // mkl 
    void make_XMat_mkl(float* X, bool grm_d_flag);
    void std_XMat_mkl(float* X, vector<double> &sd_SNP, bool grm_xchr_flag, bool miss_with_mu = false, bool divid_by_std = true);
    void std_XMat_d_mkl(float *X, vector<double> &sd_SNP, bool miss_with_mu = false, bool divid_by_std = true);
    void output_grm_mkl(float* A, bool output_grm_bin);
    bool comput_inverse_logdet_LDLT_mkl(eigenMatrix &Vi, double &logdet);
    bool comput_inverse_logdet_LU_mkl(eigenMatrix &Vi, double &logdet);
    bool comput_inverse_logdet_LU_mkl_array(int n, float *Vi, double &logdet);
    void LD_pruning_blk_mkl(float *X, vector<int> &brk_pnt, double rsq_cutoff, vector<int> &rm_snp_ID1);
    void calcu_ssx_sqrt_i_mkl(float *X_std, vector<double> &ssx);
    void calcu_ld_blk_mkl(float *X, vector<double> &ssx, vector<int> &brk_pnt, vector<int> &brk_pnt3, eigenVector &mean_rsq, eigenVector &snp_num, eigenVector &max_rsq, bool second, double rsq_cutoff);
    void calcu_ld_blk_split_mkl(int size, int size_limit, float *X_sub, vector<double> &ssx_sub, double rsq_cutoff, vector<double> &rsq_size, vector<double> &mean_rsq_sub, vector<double> &max_rsq_sub, int s1, int s2, bool second);
    

    // mlma
    void mlma_calcu_stat(float *y, float *geno_mkl, unsigned long n, unsigned long m, eigenVector &beta, eigenVector &se, eigenVector &pval);
    void mlma_calcu_stat_covar(float *y, float *geno_mkl, unsigned long n, unsigned long m, eigenVector &beta, eigenVector &se, eigenVector &pval);
    void grm_minus_grm(float *grm, float *sub_grm);

    // population
    void read_subpopu(string filename, vector<string> &subpopu, vector<string> &subpopu_name);

    /*
    // weighting GRM: ldwt_wind = window size for mean LD calculation; ld_seg = block size;
    void calcu_ld_blk_ldwt(eigenVector &ssx_sqrt_i, vector<int> &brk_pnt, vector<int> &brk_pnt3, eigenVector &mean_rsq, eigenVector &snp_num, eigenVector &max_rsq, bool second, double rsq_cutoff, bool adj);
    void calcu_ld_blk_split_ldwt(int size, int size_limit, int s_pnt, eigenVector &ssx_sqrt_i_sub, double rsq_cutoff, eigenVector &rsq_size, eigenVector &mean_rsq_sub, eigenVector &max_rsq_sub, int s1, int s2, bool second, bool adj);
    void calcu_lds(string i_ld_file, eigenVector &wt, int ldwt_wind, int ldwt_seg, double rsq_cutoff);
    void calcu_ldak(eigenVector &wt, int ldwt_seg, double rsq_cutoff);
    void calcu_ldak_blk(eigenVector &wt, eigenVector &sum_rsq, eigenVector &ssx_sqrt_i, vector<int> &brk_pnt, bool second, double rsq_cutoff);
    void calcu_ldwt(string i_ld_file, eigenVector &wt, int wind_size, double rsq_cutoff);
    void read_mrsq_mb(string i_ld_file, vector<float> &seq, vector<double> &mrsq_mb, eigenVector &wt, eigenVector &snp_m);
    void adj_wt_4_maf(eigenVector &wt);
    void cal_sum_rsq_mb(eigenVector &sum_rsq_mb);
    void col_std(MatrixXf &X);
    void assign_snp_2_mb(vector<float> &seq, vector< vector<int> > &maf_bin_pos, int mb_num);
    void make_grm_pca_blk(vector<int> & maf_bin_pos_i, int ldwt_seg, double &trace);
    void glpk_simplex_solver(MatrixXf &rsq, eigenVector &mrsq, eigenVector &wt, int maxiter);
    void calcu_grm_wt_mkl(string i_ld_file, float *X, vector<double> &sd_SNP, eigenVector &wt, int wind_size, double rsq_cutoff, int wt_mtd, int ttl_snp_num);
    */

    // gene based association test
    void sbat_read_snpAssoc(string snpAssoc_file, vector<string> &snp_name, vector<int> &snp_chr, vector<int> &snp_bp, vector<double> &snp_pval);
    void sbat_read_geneAnno(string gAnno_file, vector<string> &gene_name, vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2);
    void sbat_read_snpset(string snpset_file, vector<string> &set_name, vector< vector<string> > &snpset);
    void sbat_calcu_lambda(vector<int> &snp_indx, VectorXd &eigenval, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx);
    void get_sbat_seg_blk(int seg_size, vector< vector<int> > &snp_set_indx, vector<int> &set_chr, vector<int> &set_start_bp, vector<int> &set_end_bp);
    void rm_cor_sbat(MatrixXf &R, double R_cutoff, int m, vector<int> &rm_ID1);

	// gene based association test 
	void gbat_read_snpAssoc(string snpAssoc_file, vector<string>& snp_name, vector<int>& snp_chr, vector<int>& snp_bp, vector<double>& snp_pval);
	void gbat_read_geneAnno(string gAnno_file, vector<string>& gene_name, vector<int>& gene_chr, vector<int>& gene_bp1, vector<int>& gene_bp2);
	void gbat_calcu_ld(MatrixXf & X, eigenVector & sumsq_x, int snp1_indx, int snp2_indx, MatrixXf & C);
	void gbat(string sAssoc_file, string gAnno_file, int wind, int simu_num);
	double gbat_simu_p(int & seed, int size, eigenMatrix & L, int simu_num, double chisq_o);


    //////////////////////
    // gene expresion data
    void init_e_include();
    void std_probe(vector< vector<bool> > &X_bool, bool divid_by_std);
    void std_probe_ind(vector< vector<bool> > &X_bool, bool divid_by_std);
    
    // ecojo
    void calcu_eR();
    void read_e_metafile(string e_metafile);
    void ecojo_slct_output(bool joint_only, vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ);
    void ecojo_slct(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC);
    bool ecojo_slct_entry(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC);
    void ecojo_slct_stay(vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ);
    void ecojo_joint(const vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ);
    void ecojo_cond(const vector<int> &slct, const vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC);
    bool ecojo_init_R(const vector<int> &slct);
    void ecojo_init_RC(const vector<int> &slct, const vector<int> &remain);
    bool ecojo_insert_R(const vector<int> &slct, int insert_indx);
    void ecojo_erase_R(const vector<int> &slct);
    void ecojo_inv_R();
    void ecojo_blup(double lambda);   

    // mtCOJO and GSMR
    void init_meta_snp_map(vector<string> snplist, map<string, int> &snp_name_map, vector<string> &snp_name, vector<int> &remain_snp);
    void init_gwas_variable(vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, eigenMatrix &snp_freq, eigenMatrix &snp_b, eigenMatrix &snp_se, eigenMatrix &snp_pval, eigenMatrix &n, int npheno, int nsnp);
    void update_meta_snp_list(vector<string> &snplist, map<string, int> snp_id_map);
    void update_meta_snp_map(vector<string> snplist, map<string, int> &snp_id_map, vector<string> &snp_id, vector<int> &snp_indx, bool indx_flag);
    void update_meta_snp(map<string,int> &snp_name_map, vector<string> &snp_name, vector<int> &snp_remain);
    vector<string> remove_bad_snps(vector<string> snp_name, vector<int> snp_remain, vector<vector<bool>> snp_flag, vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, eigenMatrix &snp_freq,  eigenMatrix &snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, map<string,int> plink_snp_name_map, vector<string> snp_ref_a1, vector<string> snp_ref_a2, int ntarget, int ncovar, string outfile_name);
    vector<string> remove_freq_diff_snps(vector<string> meta_snp_name, vector<int> meta_snp_remain, map<string,int> snp_name_map, vector<double> ref_freq, eigenMatrix meta_freq, vector<vector<bool>> snp_flag, int ntrait, double freq_thresh, string outfile_name);
    vector<string> remove_mono_snps(map<string,int> snp_name_map, vector<double> ref_snpfreq, string outfile_name);
    vector<string> filter_meta_snp_pval(vector<string> snp_name, vector<int> remain_snp_indx,  eigenMatrix snp_pval, int start_indx, int end_indx, vector<vector<bool>> snp_flag, double pval_thresh);
    vector<double> gsmr_meta(vector<string> &snp_instru, eigenVector bzx, eigenVector bzx_se, eigenVector bzx_pval, eigenVector bzy, eigenVector bzy_se, eigenVector bzy_pval, double rho_pheno, vector<bool> snp_flag, double gwas_thresh, int wind_size, double r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, string &pleio_snps, string &err_msg);
    vector<string> clumping_meta(eigenVector snp_chival, vector<bool> snp_flag, double pval_thresh, int wind_size, double r2_thresh);
    void update_mtcojo_snp_rm(vector<string> adjsnps, map<string,int> &snp_id_map, vector<int> &remain_snp_indx);
    vector<string> read_snp_ldsc(map<string,int> ldsc_snp_name_map, vector<string> snp_name, vector<int> snp_remain, int &ttl_mk_num, string ref_ld_dirt, string w_ld_dirt, vector<double> &ref_ld_vec, vector<double> &w_ld_vec);
    void reorder_snp_effect(vector<int> snp_remain, eigenMatrix &bhat_z, eigenMatrix &bhat_n, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, vector<vector<bool>> &snp_flag, vector<vector<bool>> snp_val_flag, vector<int> &nsnp_cm_trait, vector<string> cm_ld_snps, map<string,int> ldsc_snp_name_map, eigenVector &ref_ld, eigenVector &w_ld, vector<double> ref_ld_vec, vector<double> w_ld_vec, int ntrait);
    eigenMatrix ldsc_snp_h2(eigenMatrix bhat_z, eigenMatrix bhat_n, eigenVector ref_ld, eigenVector w_ld, vector<vector<bool>> snp_flag, vector<int> nsnp_cm_trait, int n_cm_ld_snps, int ttl_mk_num, vector<string> trait_name, int ntrait);
    eigenMatrix ldsc_snp_rg(eigenMatrix ldsc_var_h2, eigenMatrix bhat_z, eigenMatrix bhat_n, eigenVector ref_ld, eigenVector w_ld, vector<vector<bool>> snp_flag, vector<int> trait_indx1, vector<int> trait_indx2, int n_cm_ld_snps, int ttl_mk_num, vector<string> trait_name);

    // Ajust summarydata for PC
    void adjust_snp_effect_for_pc(eigenVector &bzy_adj, eigenVector &bzx_hat, eigenVector bzy, eigenVector bxy_hat, int wind_size);

    // inline functions
    template<typename ElemType>
    void makex(int j, vector<ElemType> &x, bool minus_2p = false) {
        int i = 0;
        x.resize(_keep.size());
        for (i = 0; i < _keep.size(); i++) {
            if (!_snp_1[_include[j]][_keep[i]] || _snp_2[_include[j]][_keep[i]]) {
                if (_allele1[_include[j]] == _ref_A[_include[j]]) x[i] = (_snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]]);
                else x[i] = 2.0 - (_snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]]);
            } else x[i] = _mu[_include[j]];
            if (minus_2p) x[i] -= _mu[_include[j]];
        }
    }

private:
    // read in plink files
    // bim file
    int _autosome_num;
    vector<int> _chr;
    vector<string> _snp_name;
    map<string, int> _snp_name_map;
    map<string, string> _snp_name_per_chr;
    vector<double> _genet_dst;
    vector<int> _bp;
    vector<string> _allele1;
    vector<string> _allele2;
    vector<string> _ref_A; // reference allele
    vector<string> _other_A; // the other allele
    int _snp_num;
    vector<double> _rc_rate;
    vector<int> _include; // initialized in the read_bimfile()
    eigenVector _maf;

    // fam file
    vector<string> _fid;
    vector<string> _pid;
    map<string, int> _id_map;
    vector<string> _fa_id;
    vector<string> _mo_id;
    vector<int> _sex;
    vector<double> _pheno;
    int _indi_num;
    vector<int> _keep; // initialized in the read_famfile()
    eigenMatrix _varcmp_Py; // BLUP solution to the total genetic effects of individuals

    // bed file
    vector< vector<bool> > _snp_1;
    vector< vector<bool> > _snp_2;

    // imputed data
    bool _dosage_flag;
    vector< vector<float> > _geno_dose;
    vector<double> _impRsq;

    // genotypes
    MatrixXf _geno;

    // QC
    double _rm_ld_cutoff;

    // grm 
    MatrixXf _grm_N;
    eigenMatrix _grm;
    float * _grm_mkl;
    float * _geno_mkl;
    bool _grm_bin_flag;

    // reml
    int _n;
    int _X_c;
    vector<int> _r_indx;
    vector<int> _r_indx_drop;
    int _reml_max_iter;
    int _reml_mtd;
    int _reml_inv_mtd;
    double _reml_diag_mul;
    int _reml_diagV_adj;
    bool _cv_blup;
    eigenMatrix _X;
    vector<eigenMatrix> _A;
    eigenVector _y;
    eigenMatrix _Vi;
    eigenMatrix _P;
    eigenVector _b;
    vector<string> _var_name;
    vector<double> _varcmp;
    vector<string> _hsq_name;
    double _y_Ssq;
    double _ncase;
    bool _flag_CC;
    bool _flag_CC2;
    bool _reml_diag_one;
    bool _reml_have_bend_A;
    int _V_inv_mtd;
    bool _reml_force_inv;
    bool _reml_AI_not_invertible;
    bool _reml_force_converge;
    bool _reml_no_converge;
    bool _reml_fixed_var;
    bool _reml_allow_constrain_run = false;

    // within-family reml analysis
    bool _within_family;
    vector<int> _fam_brk_pnt;

    // bivariate reml
    bool _bivar_reml;
    bool _ignore_Ce;
    bool _bivar_no_constrain;
    double _y2_Ssq;
    double _ncase2;
    vector< vector<int> > _bivar_pos;
    vector< vector<int> > _bivar_pos_prev;
    vector< eigenSparseMat > _Asp;
    vector< eigenSparseMat > _Asp_prev;
    vector<eigenMatrix> _A_prev;
    vector<double> _fixed_rg_val;

    vector<double> _mu;
    string _out;
    bool _save_ram;

    // LD
    vector<string> _ld_target_snp;
    bool _ldscore_adj;

    // joint analysis of META
    bool _jma_actual_geno;
    int _jma_wind_size;
    double _jma_p_cutoff;
    double _jma_collinear;
    double _jma_Vp;
    double _jma_Ve;
    double _GC_val;
    int _jma_snpnum_backward;
    int _jma_snpnum_collienar;
    eigenVector _freq;
    eigenVector _beta;
    eigenVector _beta_se;
    eigenVector _chisq;
    eigenVector _pval;
    eigenVector _N_o;
    eigenVector _Nd;
    eigenVector _MSX;
    eigenVector _MSX_B;
    eigenSparseMat _B_N;
    eigenSparseMat _B;
    eigenSparseMat _B_N_i;
    eigenSparseMat _B_i;
    eigenVector _D_N;
    eigenSparseMat _Z_N;
    eigenSparseMat _Z;
    double g_massoc_out_thresh = -1.0;
    double _diff_freq = 0.2;
    
    // GSMR analysis
    int _expo_num;
    int _outcome_num;
    int _n_gsmr_rst_item = 5;
    int _gsmr_beta_version = 0;
    eigenMatrix _r_pheno_sample;
    vector<string> _gwas_trait_name;
    vector<vector<bool>> _snp_val_flag;
    
    // mtCOJO analysis
    string _target_pheno_name="";
    vector<string> _meta_snp_name;
    vector<int> _meta_remain_snp;
    map<string,int> _meta_snp_name_map;
    vector<string> _covar_pheno_name;
    vector<string> _meta_snp_a1;
    vector<string> _meta_snp_a2;
    eigenMatrix _meta_snp_freq;
    eigenMatrix _meta_snp_b;
    eigenMatrix _meta_snp_se;
    eigenMatrix _meta_snp_pval;
    eigenMatrix _meta_snp_n_o;
    eigenVector _meta_vp_trait;
    vector<double> _meta_popu_prev;
    vector<double> _meta_smpl_prev;
    
    // gene-trait association
    map<string, int> _probe_name_map;
    vector<string> _probe_name;
    int _probe_num;
    eigenMatrix _probe_data;
    vector<int> _e_include;
    eigenVector _ecojo_z;
    eigenVector _ecojo_b;
    eigenVector _ecojo_se;
    eigenVector _ecojo_n;
    eigenVector _ecojo_pval;
    eigenMatrix _ecojo_R;
    eigenMatrix _ecojo_RC;
    eigenMatrix _ecojo_wholeR;
    double _ecojo_p_cutoff;
    double _ecojo_collinear;

    // PC adjustment
    double _ttl_mk_num;
    vector<double> _eigen_value;

    //vif threshold
    double _vif_threshold;

};

class locus_bp {
public:
    string locus_name;
    int chr;
    int bp;

    locus_bp(string locus_name_buf, int chr_buf, int bp_buf) {
        locus_name = locus_name_buf;
        chr = chr_buf;
        bp = bp_buf;
    }

    bool operator()(const locus_bp & other) {
        return (chr == other.chr && bp <= other.bp);
    }
};

#endif
