/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * GCTA options
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include <stdlib.h>
#include "gcta.h"

void option(int option_num, char* option_str[]);

int main(int argc, char* argv[])
{
    cout << "*******************************************************************" << endl;
    cout << "* Genome-wide Complex Trait Analysis (GCTA)" << endl;
    cout << "* version 1.24.5" << endl;
    cout << "* (C) 2010-2013 Jian Yang, Hong Lee, Michael Goddard and Peter Visscher" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;

    long int time_used = 0, start = time(NULL);
    time_t curr = time(0);
    cout << "Analysis started: " << ctime(&curr) << endl;
    cout << "Options:" << endl;
    try {
        option(argc, argv);
    } catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    } catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }
    curr = time(0);
    cout << "\nAnalysis finished: " << ctime(&curr);
    time_used = time(NULL) - start;
    cout << "Computational time: " << time_used / 3600 << ":" << (time_used % 3600) / 60 << ":" << time_used % 60 << endl;

    return 0;
}

void option(int option_num, char* option_str[])
{
    int i = 0, j = 0;

    // OpenMP
    bool thread_flag = false;
    int thread_num = 1;

    // raw genotype data
    string RG_fname_file = "", RG_summary_file = "";
    double GC_cutoff = 0.7;

    // data management
    string bfile = "", bfile2 = "", update_sex_file = "", update_freq_file = "", update_refA_file = "", kp_indi_file = "", rm_indi_file = "", extract_snp_file = "", exclude_snp_file = "", extract_snp_name = "", exclude_snp_name = "", out = "gcta";
    bool SNP_major = false, bfile_flag = false, make_bed_flag = false, dose_mach_flag = false, dose_mach_gz_flag = false, dose_beagle_flag = false, bfile2_flag = false, out_freq_flag = false, out_ssq_flag = false;
    bool ref_A = false, recode = false, recode_nomiss = false, save_ram = false, autosome_flag = false;
    int autosome_num = 22, extract_chr_start = 0, extract_chr_end = 0;
    string dose_file = "", dose_info_file = "", update_impRsq_file = "";
    double maf = 0.0, max_maf = 0.0, dose_Rsq_cutoff = 0.0;

    // GRM
    bool ibc = false, ibc_all = false, grm_flag = false, grm_bin_flag = true, m_grm_flag = false, m_grm_bin_flag = true, make_grm_flag = false, make_grm_inbred_flag = false, dominance_flag = false, make_grm_xchar_flag = false, grm_out_bin_flag = true, make_grm_ldwt_flag = false, make_grm_wt_impRsq_flag = false, make_grm_f3_flag = false;
    bool grm_pca_flag = false, pca_flag = false;
    double grm_adj_fac = -2.0, grm_cutoff = -2.0, rm_high_ld_cutoff = -1.0, ldwt_wind = 2e5;
    int dosage_compen = -2, out_pc_num = 20, make_grm_mtd = 0, make_grm_ldwt_mtd = -1, ttl_snp_num = -1;
    string grm_file = "", paa_file = "";

    // LD
    string LD_file = "", i_ld_file = "";
    bool LD = false, LD_search = false, LD_i = false, ld_mean_rsq_flag = false, ld_max_rsq_flag = false;
    int LD_step = 10;
    double LD_wind = 1e7, LD_sig = 0.05, LD_prune_rsq = -1.0, LD_rsq_cutoff = 0.0;

    // initialize paramters for simulation based on real genotype data
    bool simu_qt_flag = false, simu_cc = false, simu_emb_flag = false, simu_output_causal = false;
    int simu_rep = 1, simu_case_num = 0, simu_control_num = 0;
    double simu_h2 = 0.1, simu_K = 0.1, simu_gener = 100, simu_seed = -CommFunc::rand_seed();
    string simu_causal = "";

    // simulate unlinked SNPs
    bool simu_unlinked_flag = false;
    int simu_unlinked_n = 1, simu_unlinked_m = 1;
    double simu_unlinked_maf = 0.0;

    // estimate genetic distance based on hapmap_data
    bool hapmap_genet_dst = false;
    string hapmap_genet_dst_file = "";

    // REML analysis
    bool prevalence_flag = false;
    int mphen = 1, mphen2 = 2, reml_mtd = 0, MaxIter = 100;
    double prevalence = -2.0, prevalence2 = -2.0;
    bool reml_flag = false, pred_rand_eff = false, est_fix_eff = false, blup_snp_flag = false, no_constrain = false, reml_lrt_flag = false, no_lrt = false, bivar_reml_flag = false, ignore_Ce = false, within_family = false, reml_bending = false, HE_reg_flag = false, reml_diag_one = false, bivar_no_constrain = false;
    string phen_file = "", qcovar_file = "", covar_file = "", qgxe_file = "", gxe_file = "", blup_indi_file = "";
    vector<double> reml_priors, reml_priors_var, fixed_rg_val;
    vector<int> reml_drop;
    reml_drop.push_back(1);

    // Joint analysis of GWAS MA
    string massoc_file = "", massoc_init_snplist = "", massoc_cond_snplist = "";
    int massoc_wind = 1e7, massoc_top_SNPs = -1, massoc_mld_slct_alg = 0;
    double massoc_p = 5e-8, massoc_collinear = 0.9, massoc_sblup_fac = -1, massoc_gc_val = -1;
    bool massoc_slct_flag = false, massoc_joint_flag = false, massoc_sblup_flag = false, massoc_gc_flag = false, massoc_actual_geno_flag = false;

    // mixed linear model association 
    bool mlma_flag = false, mlma_loco_flag = false, mlma_no_adj_covar = false;

    // Fst
    string subpopu_file = "";

    // gene-based association test
    string sbat_sAssoc_file = "", sbat_gAnno_file = "", sbat_snpset_file = "";
    int sbat_wind = 50000;

    // gene expression data
    string efile="", eR_file = "", ecojo_ma_file="";
    int make_erm_mtd = 1;
    double ecojo_p = 5e-6, ecojo_collinear = 0.9, ecojo_lambda = -1;
    bool efile_flag=false, eR_file_flag = false, ecojo_slct_flag = false, ecojo_blup_flag = false, make_erm_flag = false;

    int argc = option_num;
    vector<char *> argv(option_num + 2);
    for (i = 0; i < option_num; i++) argv[i] = option_str[i];
    argv[option_num] = "gcta";
    argv[option_num + 1] = "gcta";
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--thread-num") == 0) {
            thread_num = atoi(argv[++i]);
            cout << "--thread-num " << thread_num << endl;
            if (thread_num < 1 || thread_num > 1000) throw ("\nError: --thread-num should be from 1 to 1000.\n");
        }// raw genotype data
        else if (strcmp(argv[i], "--raw-files") == 0) {
            RG_fname_file = argv[++i];
            cout << "--raw-files " << argv[i] << endl;
        } else if (strcmp(argv[i], "--raw-summary") == 0) {
            RG_summary_file = argv[++i];
            cout << "--raw-summary " << argv[i] << endl;
        } else if (strcmp(argv[i], "--gencall") == 0) {
            GC_cutoff = atof(argv[++i]);
            cout << "--gencall " << GC_cutoff << endl;
            if (GC_cutoff < 0.0 || GC_cutoff > 1.0) throw ("\nError: --gencall should be within the range from 0 to 1.\n");
        }            // data management
        else if (strcmp(argv[i], "--bfile") == 0) {
            bfile_flag = true;
            bfile = argv[++i];
            cout << "--bfile " << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-bed") == 0) {
            make_bed_flag = true;
            cout << "--make-bed " << endl;
        } else if (strcmp(argv[i], "--bfile2") == 0) {
            bfile2_flag = true;
            bfile2 = argv[++i];
            cout << "--bfile2 " << argv[i] << endl;
        } else if (strcmp(argv[i], "--dosage-mach") == 0) {
            dose_mach_flag = true;
            dose_beagle_flag = false;
            dose_file = argv[++i];
            dose_info_file = argv[++i];
            cout << "--dosage-mach " << dose_file << " " << dose_info_file << endl;
        } else if (strcmp(argv[i], "--dosage-mach-gz") == 0) {
            dose_mach_gz_flag = true;
            dose_beagle_flag = false;
            dose_file = argv[++i];
            dose_info_file = argv[++i];
            cout << "--dosage-mach-gz " << dose_file << " " << dose_info_file << endl;
        } else if (strcmp(argv[i], "--dosage-beagle") == 0) {
            dose_beagle_flag = true;
            dose_mach_flag = false;
            dose_mach_gz_flag = false;
            dose_file = argv[++i];
            dose_info_file = argv[++i];
            cout << "--dosage-beagle " << dose_file << " " << dose_info_file << endl;
        } else if (strcmp(argv[i], "--imput-rsq") == 0) {
            dose_Rsq_cutoff = atof(argv[++i]);
            cout << "--imput-rsq " << dose_Rsq_cutoff << endl;
            if (dose_Rsq_cutoff < 0.0 || dose_Rsq_cutoff > 1.0) throw ("\nError: --imput-rsq should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--update-imput-rsq") == 0) {
            update_impRsq_file = argv[++i];
            cout << "--update-imput-rsq " << update_impRsq_file << endl;
            CommFunc::FileExist(update_impRsq_file);
        } else if (strcmp(argv[i], "--update-freq") == 0) {
            update_freq_file = argv[++i];
            cout << "--update-freq " << update_freq_file << endl;
            CommFunc::FileExist(update_freq_file);
        } else if (strcmp(argv[i], "--update-ref-allele") == 0) {
            update_refA_file = argv[++i];
            cout << "--update-ref-allele " << update_refA_file << endl;
            CommFunc::FileExist(update_refA_file);
        } else if (strcmp(argv[i], "--keep") == 0) {
            kp_indi_file = argv[++i];
            cout << "--keep " << kp_indi_file << endl;
            CommFunc::FileExist(kp_indi_file);
        } else if (strcmp(argv[i], "--remove") == 0) {
            rm_indi_file = argv[++i];
            cout << "--remove " << rm_indi_file << endl;
            CommFunc::FileExist(rm_indi_file);
        } else if (strcmp(argv[i], "--update-sex") == 0) {
            update_sex_file = argv[++i];
            cout << "--update-sex " << update_sex_file << endl;
            CommFunc::FileExist(update_sex_file);
        } else if (strcmp(argv[i], "--chr") == 0) {
            extract_chr_start = extract_chr_end = atoi(argv[++i]);
            cout << "--chr " << extract_chr_start << endl;
            if (extract_chr_start < 1 || extract_chr_start > 100) throw ("\nError: --chr should be within the range from 1 to 100.\n");
        } else if (strcmp(argv[i], "--autosome-num") == 0) {
            autosome_num = atoi(argv[++i]);
            cout << "--autosome-num " << autosome_num << endl;
            if (autosome_num < 1 || autosome_num > 100) throw ("\nError: invalid number specified after the option --autosome-num.\n");
        } else if (strcmp(argv[i], "--autosome") == 0) {
            autosome_flag = true;
            cout << "--autosome" << endl;
        } else if (strcmp(argv[i], "--extract") == 0) {
            extract_snp_file = argv[++i];
            cout << "--extract " << extract_snp_file << endl;
            CommFunc::FileExist(extract_snp_file);
        } else if (strcmp(argv[i], "--exclude") == 0) {
            exclude_snp_file = argv[++i];
            cout << "--exclude " << exclude_snp_file << endl;
            CommFunc::FileExist(exclude_snp_file);
        } else if (strcmp(argv[i], "--extract-snp") == 0) {
            extract_snp_name = argv[++i];
            cout << "--extract-snp " << extract_snp_name << endl;
        } else if (strcmp(argv[i], "--exclude-snp") == 0) {
            exclude_snp_name = argv[++i];
            cout << "--exclude-snp " << exclude_snp_name << endl;
        } else if (strcmp(argv[i], "--maf") == 0) {
            maf = atof(argv[++i]);
            cout << "--maf " << maf << endl;
            if (maf < 0 || maf > 0.5) throw ("\nError: --maf should be within the range from 0 to 0.5.\n");
        } else if (strcmp(argv[i], "--max-maf") == 0) {
            max_maf = atof(argv[++i]);
            cout << "--max-maf " << max_maf << endl;
            if (max_maf <= 0) throw ("\nError: --max-maf should be > 0.\n");
        } else if (strcmp(argv[i], "--out") == 0) {
            out = argv[++i];
            cout << "--out " << out << endl;
        } else if (strcmp(argv[i], "--freq") == 0) {
            out_freq_flag = true;
            thread_flag = true;
            cout << "--freq" << endl;
        } else if (strcmp(argv[i], "--ssq") == 0) {
            out_ssq_flag = true;
            cout << "--ssq" << endl;
        } else if (strcmp(argv[i], "--recode") == 0) {
            recode = true;
            thread_flag = true;
            cout << "--recode" << endl;
        } else if (strcmp(argv[i], "--recode-nomiss") == 0) {
            recode_nomiss = true;
            thread_flag = true;
            cout << "--recode-nomiss" << endl;
        } else if (strcmp(argv[i], "--save-ram") == 0) {
            save_ram = true;
            cout << "--save-ram" << endl;
        }// GRM
        else if (strcmp(argv[i], "--paa") == 0) {
            paa_file = argv[++i];
            cout << "--paa " << paa_file << endl;
            CommFunc::FileExist(paa_file);
        } else if (strcmp(argv[i], "--ibc") == 0) {
            ibc = true;
            cout << "--ibc" << endl;
        } else if (strcmp(argv[i], "--ibc-all") == 0) {
            ibc = ibc_all = true;
            cout << "--ibc-all" << endl;
        } else if (strcmp(argv[i], "--mgrm") == 0 || strcmp(argv[i], "--mgrm-bin") == 0) {
            m_grm_flag = true;
            grm_file = argv[++i];
            cout << argv[i - 1] << " " << grm_file << endl;
        } else if (strcmp(argv[i], "--mgrm-gz") == 0) {
            m_grm_flag = true;
            m_grm_bin_flag = false;
            grm_bin_flag = false;
            grm_file = argv[++i];
            cout << "--mgrm-gz " << grm_file << endl;
        } else if (strcmp(argv[i], "--grm") == 0 || strcmp(argv[i], "--grm-bin") == 0) {
            grm_flag = true;
            grm_file = argv[++i];
            cout << argv[i - 1] << " " << grm_file << endl;
        } else if (strcmp(argv[i], "--grm-gz") == 0) {
            grm_flag = true;
            m_grm_bin_flag = false;
            grm_bin_flag = false;
            grm_file = argv[++i];
            cout << "--grm-gz " << grm_file << endl;
        } else if (strcmp(argv[i], "--rm-high-ld") == 0) {
            rm_high_ld_cutoff = atof(argv[++i]);
            cout << "--rm-high-ld " << rm_high_ld_cutoff << endl;
            if (rm_high_ld_cutoff <= 0 || rm_high_ld_cutoff >= 1) throw ("\nError: the value to be specified after --rm-high-ld should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--make-grm") == 0 || strcmp(argv[i], "--make-grm-bin") == 0) {
            make_grm_flag = true;
            thread_flag = true;
            cout << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-gz") == 0) {
            make_grm_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            cout << "--make-grm-gz" << endl;
        } else if (strcmp(argv[i], "--make-grm-alg") == 0) {
            make_grm_flag = true;
            make_grm_mtd = atoi(argv[++i]);
            thread_flag = true;
            cout << "--make-grm-alg " << make_grm_mtd << endl;
            if (make_grm_mtd < 0 || make_grm_mtd > 1) throw ("\nError: --make-grm-alg should be 0 or 1.\n");
        } else if (strcmp(argv[i], "--make-grm-f3") == 0) {
            make_grm_f3_flag = true;
            grm_out_bin_flag = true;
            thread_flag = true;
            cout << "--make-grm-f3" << endl;
        } else if (strcmp(argv[i], "--make-grm-ld") == 0) {
            make_grm_flag = true;
            make_grm_ldwt_flag = true;
            make_grm_ldwt_mtd = 0;
            thread_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                ldwt_wind = 200;
                i--;
            } else ldwt_wind = atoi(argv[i]);
            cout << "--make-grm-ld " << ldwt_wind << endl;
            if (ldwt_wind < 0 || ldwt_wind > 20000) throw ("\nError: block size for --make-grm-ld should be between 0Kb to 20Mb.\n");
            ldwt_wind *= 1000;
        } else if (strcmp(argv[i], "--make-grm-ld-alg") == 0) {
            make_grm_flag = true;
            make_grm_ldwt_flag = true;
            make_grm_ldwt_mtd = atoi(argv[++i]);
            cout << "--make-grm-ld-alg " << make_grm_ldwt_mtd << endl;
            if (make_grm_ldwt_mtd < 0 || make_grm_ldwt_mtd > 2) throw ("\nError: --make-grm-ld-alg should be 0 or 2.\n");
        } else if (strcmp(argv[i], "--make-grm-d") == 0 || strcmp(argv[i], "--make-grm-d-bin") == 0) {
            make_grm_flag = true;
            dominance_flag = true;
            thread_flag = true;
            cout << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-d-gz") == 0) {
            make_grm_flag = true;
            dominance_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            cout << "--make-grm-d-gz" << endl;
        } else if (strcmp(argv[i], "--dominance") == 0) {
            dominance_flag = true;
            thread_flag = true;
            cout <<"--dominance"<< endl;
        } else if (strcmp(argv[i], "--make-grm-xchr") == 0 || strcmp(argv[i], "--make-grm-xchr-bin") == 0) {
            make_grm_flag = true;
            make_grm_xchar_flag = true;
            thread_flag = true;
            cout << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-xchr-gz") == 0) {
            make_grm_flag = true;
            make_grm_xchar_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            cout << "--make-grm-xchr-gz" << endl;
        } else if (strcmp(argv[i], "--make-grm-wt-imp") == 0) {
            make_grm_flag = true;
            make_grm_wt_impRsq_flag = true;
            update_impRsq_file = argv[++i];
            thread_flag = true;
            cout << "--make-grm-wt-imp " << update_impRsq_file << endl;
            CommFunc::FileExist(update_impRsq_file);
        } else if (strcmp(argv[i], "--make-grm-inbred") == 0 || strcmp(argv[i], "--make-grm-inbred-bin") == 0) {
            make_grm_flag = true;
            make_grm_inbred_flag = true;
            thread_flag = true;
            cout << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-inbred-gz") == 0) {
            make_grm_flag = true;
            grm_out_bin_flag = false;
            make_grm_inbred_flag = true;
            thread_flag = true;
            cout << "--make-grm-inbred-gz" << endl;
        } else if (strcmp(argv[i], "--grm-adj") == 0) {
            grm_adj_fac = atof(argv[++i]);
            cout << "--grm-adj " << grm_adj_fac << endl;
            if (grm_adj_fac < 0 || grm_adj_fac > 1) throw ("\nError: the value to be specified after --grm-adj should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--dc") == 0) {
            dosage_compen = atoi(argv[++i]);
            cout << "--dc " << dosage_compen << endl;
            if (dosage_compen != 0 && dosage_compen != 1) throw ("\nError: the value to be specified after --dc should be 0 or 1.\n");
        } else if (strcmp(argv[i], "--grm-cutoff") == 0) {
            grm_cutoff = atof(argv[++i]);
            if (grm_cutoff >= -1 && grm_cutoff <= 2) cout << "--grm-cutoff " << grm_cutoff << endl;
            else grm_cutoff = -2;
        } else if (strcmp(argv[i], "--pca") == 0) {
            pca_flag = true;
            thread_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                out_pc_num = 20;
                i--;
            } else out_pc_num = atoi(argv[i]);
            cout << "--pca " << out_pc_num << endl;
            if (out_pc_num < 1) throw ("\nError: the value to be specified after --pca should be positive.\n");
        }// estimation of LD structure
        else if (strcmp(argv[i], "--ld") == 0) {
            LD = true;
            LD_file = argv[++i];
            cout << "--ld " << LD_file << endl;
            CommFunc::FileExist(LD_file);
        } else if (strcmp(argv[i], "--ld-step") == 0) {
            LD_search = true;
            LD_step = atoi(argv[++i]);
            cout << "--ld-step " << LD_step << endl;
            if (LD_step < 1 || LD_step > 20) throw ("\nError: --ld-step should be within the range from 1 to 20.\n");
        } else if (strcmp(argv[i], "--ld-wind") == 0 || strcmp(argv[i], "--ld-pruning-wind") == 0 || strcmp(argv[i], "--make-grm-wt-wind") == 0) {
            LD_wind = atof(argv[++i]);
            cout << argv[i - 1] << " " << LD_wind << endl;
            LD_wind *= 1000;
            if (LD_wind < 1e3 || LD_wind > 2e7) {
                stringstream err_msg;
                err_msg << "\nError: " << argv[i - 1] << " should be 1Kb or 20Mb.\n";
                throw (err_msg.str());
            }
        } else if (strcmp(argv[i], "--ld-sig") == 0) {
            LD_sig = atof(argv[++i]);
            cout << "--ld-sig " << LD_sig << endl;
            if (LD_sig <= 0) throw ("\nError: --ld-sig should be > 0.\n");
        } else if (strcmp(argv[i], "--ld-i") == 0) {
            LD_i = true;
            cout << "--ld-i" << endl;
        } else if (strcmp(argv[i], "--ld-pruning") == 0) {
            thread_flag = true;
            LD_prune_rsq = atof(argv[++i]);
            cout << "--ld-pruning " << LD_prune_rsq << endl;
            if (LD_prune_rsq < 0.0001 || LD_prune_rsq > 0.9999) throw ("\nError: --ld-pruning should be within the range from 0.0001 to 0.9999.\n");
        } else if (strcmp(argv[i], "--ld-mean-rsq") == 0) {
            ld_mean_rsq_flag = true;
            thread_flag = true;
            cout << "--ld-mean-rsq" << endl;
        } else if (strcmp(argv[i], "--ld-rsq-cutoff") == 0 || strcmp(argv[i], "--make-grm-wt-rsq-cutoff") == 0) {
            LD_rsq_cutoff = atof(argv[++i]);
            cout << "--ld-rsq-cutoff " << LD_rsq_cutoff << endl;
            if (LD_rsq_cutoff < 0.0 || LD_rsq_cutoff > 1.0) {
                stringstream err_msg;
                err_msg << "\nError: " << argv[i - 1] << " should be within the range from 0 to 1.\n";
                throw (err_msg.str());
            }
        } else if (strcmp(argv[i], "--ld-max-rsq") == 0) {
            ld_max_rsq_flag = true;
            thread_flag = true;
            cout << "--ld-max-rsq" << endl;
        }
        // simulation based on real genotype data
        else if (strcmp(argv[i], "--simu-qt") == 0) {
            simu_qt_flag = true;
            cout << "--simu-qt" << endl;
        } else if (strcmp(argv[i], "--simu-cc") == 0) {
            simu_cc = true;
            simu_case_num = atoi(argv[++i]);
            simu_control_num = atoi(argv[++i]);
            cout << "--simu-cc " << simu_case_num << " " << simu_control_num << endl;
            if (simu_case_num < 10) throw ("Error: --simu-cc, Invalid number of cases. Minimun number 10.");
            if (simu_control_num < 10) throw ("Error: --simu-cc, Invalid number of controls. Minimum number 10.");
        } else if (strcmp(argv[i], "--simu-rep") == 0) {
            simu_rep = atoi(argv[++i]);
            cout << "--simu-rep " << simu_rep << endl;
            if (simu_rep < 1 || simu_rep > 10000) throw ("Error: --simu-rep should be within the range from 1 to 10000.");
        } else if (strcmp(argv[i], "--simu-hsq") == 0) {
            simu_h2 = atof(argv[++i]);
            cout << "--simu-hsq " << simu_h2 << endl;
            if (simu_h2 > 1.0 || simu_h2 < 0.0) throw ("Error: --simu-h2 should be within the range from 0 to 1.");
        } else if (strcmp(argv[i], "--simu-k") == 0) {
            simu_K = atof(argv[++i]);
            cout << "--simu-k " << simu_K << endl;
            if (simu_K > 0.5 || simu_K < 0.0001) throw ("Error: --simu-K should be within the range from 0.0001 to 0.5.");
        } else if (strcmp(argv[i], "--simu-causal-loci") == 0) {
            simu_causal = argv[++i];
            cout << "--simu-causal-loci " << simu_causal << endl;
            CommFunc::FileExist(simu_causal);
        } else if (strcmp(argv[i], "--simu-embayesb") == 0) { // internal
            simu_emb_flag = true;
            cout << "--simu-embayesb" << endl;
        } else if (strcmp(argv[i], "--simu-ouput-causal") == 0) { // internal
            simu_output_causal = true;
            cout << "--simu-output-causal" << endl;
        } else if (strcmp(argv[i], "--simu-seed") == 0) {
            simu_seed = atof(argv[++i]);
            cout << "--simu-seed " << simu_seed << endl;
            if (simu_seed <= 100) throw ("Error: --simu-seed should be >100.");
        }
        else if (strcmp(argv[i], "--hapmap-genet-dst") == 0) { // calculate genetic dst based on HapMap data
            hapmap_genet_dst = true;
            hapmap_genet_dst_file = argv[++i];
            cout << "--hapmap-genet-dst " << hapmap_genet_dst_file << endl;
        }// estimate variance explained by all SNPs
        else if (strcmp(argv[i], "--HEreg") == 0) {
            HE_reg_flag = true;
            thread_flag = true;
            cout << "--HEreg" << endl;
        } else if (strcmp(argv[i], "--reml") == 0) {
            reml_flag = true;
            thread_flag = true;
            cout << "--reml" << endl;
            if (m_grm_flag) no_lrt = true;
        } else if (strcmp(argv[i], "--prevalence") == 0) {
            prevalence_flag = true;
            prevalence = atof(argv[++i]);
            cout << "--prevalence " << prevalence << endl;
            if (prevalence <= 0 || prevalence >= 1) throw ("\nError: --prevalence should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--reml-pred-rand") == 0) {
            pred_rand_eff = true;
            cout << "--reml-pred-rand" << endl;
        } else if (strcmp(argv[i], "--reml-est-fix") == 0) {
            est_fix_eff = true;
            cout << "--reml-est-fix" << endl;
        } else if (strcmp(argv[i], "--reml-alg") == 0) {
            reml_mtd = atoi(argv[++i]);
            cout << "--reml-alg " << reml_mtd << endl;
            if (reml_mtd < 0 || reml_mtd > 2) throw ("\nError: --reml-alg should be 0, 1 or 2.\n");
        } else if (strcmp(argv[i], "--reml-no-constrain") == 0) {
            reml_flag = true;
            no_constrain = true;
            cout << "--reml-no-constrain" << endl;
        } else if (strcmp(argv[i], "--reml-priors") == 0) {
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                reml_priors.push_back(atof(argv[i]));
            }
            i--;
            cout << "--reml-priors ";
            bool err_flag = false;
            for (j = 0; j < reml_priors.size(); j++) {
                cout << reml_priors[j] << " ";
                if (reml_priors[j] > 1.0 || reml_priors[j] < -10.0) err_flag = true;
            }
            cout << endl;
            if (err_flag || reml_priors.empty()) throw ("\nError: --reml-priors. Prior values should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--reml-priors-var") == 0) {
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                reml_priors_var.push_back(atof(argv[i]));
            }
            i--;
            cout << "--reml-priors-var ";
            bool err_flag = false;
            for (j = 0; j < reml_priors_var.size(); j++) {
                cout << reml_priors_var[j] << " ";
                if (reml_priors_var[j] < 0.0) err_flag = true;
            }
            cout << endl;
            if (err_flag || reml_priors_var.empty()) throw ("\nError: --reml-priors-var. Prior values should be positive.\n");
        } else if (strcmp(argv[i], "--reml-no-lrt") == 0) {
            no_lrt = true;
            cout << "--reml-no-lrt" << endl;
        } else if (strcmp(argv[i], "--reml-lrt") == 0) {
            no_lrt = false;
            reml_lrt_flag = true;
            reml_drop.clear();
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                reml_drop.push_back(atoi(argv[i]));
            }
            i--;
            cout << "--reml-lrt ";
            bool err_flag = false;
            for (j = 0; j < reml_drop.size(); j++) {
                cout << reml_drop[j] << " ";
                if (reml_drop[j] < 1) err_flag = true;
            }
            cout << endl;
            if (err_flag || reml_drop.empty()) throw ("\nError: invalid values specified after --reml-lrt.\n");
        } else if (strcmp(argv[i], "--reml-maxit") == 0) {
            MaxIter = atoi(argv[++i]);
            cout << "--reml-maxit " << MaxIter << endl;
            if (MaxIter < 1 || MaxIter > 10000) throw ("\nError: --reml-maxit should be within the range from 1 to 10000.\n");
        } else if (strcmp(argv[i], "--reml-bending") == 0) {
            reml_bending = true;
            cout << "--reml-bending " << endl;
        } else if (strcmp(argv[i], "--reml-diag-one") == 0) {
            reml_diag_one = true;
            cout << "--reml-diag-one " << endl;
        } else if (strcmp(argv[i], "--pheno") == 0) {
            phen_file = argv[++i];
            cout << "--pheno " << phen_file << endl;
            CommFunc::FileExist(phen_file);
        } else if (strcmp(argv[i], "--mpheno") == 0) {
            mphen = atoi(argv[++i]);
            cout << "--mpheno " << mphen << endl;
            if (mphen < 1) throw ("Error: --mpheno should be > 0.");
        } else if (strcmp(argv[i], "--qcovar") == 0) {
            qcovar_file = argv[++i];
            cout << "--qcovar " << qcovar_file << endl;
            CommFunc::FileExist(qcovar_file);
        } else if (strcmp(argv[i], "--covar") == 0) {
            covar_file = argv[++i];
            cout << "--covar " << covar_file << endl;
            CommFunc::FileExist(covar_file);
        } else if (strcmp(argv[i], "--gxqe") == 0) {
            qgxe_file = argv[++i];
            cout << "--gxqe " << qgxe_file << endl;
            CommFunc::FileExist(qgxe_file);
        } else if (strcmp(argv[i], "--gxe") == 0) {
            gxe_file = argv[++i];
            cout << "--gxe " << gxe_file << endl;
            CommFunc::FileExist(gxe_file);
        } else if (strcmp(argv[i], "--blup-snp") == 0) {
            blup_snp_flag = true;
            blup_indi_file = argv[++i];
            cout << "--blup-snp " << blup_indi_file << endl;
            CommFunc::FileExist(blup_indi_file);
        } else if (strcmp(argv[i], "--reml-wfam") == 0) {
            reml_flag = true;
            within_family = true;
            cout << "--reml-within-family " << endl;
        } else if (strcmp(argv[i], "--reml-bivar") == 0) {
            bivar_reml_flag = true;
            thread_flag = true;
            vector<int> mphen_buf;
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                mphen_buf.push_back(atoi(argv[i]));
            }
            i--;
            if (mphen_buf.size() < 2 && mphen_buf.size() > 0) throw ("\nError: --reml-bivar. Please specify two traits for the bivariate REML analysis.");
            if (mphen_buf.size() == 0) {
                mphen = 1;
                mphen2 = 2;
            } else {
                mphen = mphen_buf[0];
                mphen2 = mphen_buf[1];
            }
            if (mphen < 1 || mphen2 < 1 || mphen == mphen2) throw ("\nError: --reml-bivar. Invalid input parameters.");
            cout << "--reml-bivar " << mphen << " " << mphen2 << endl;
        } else if (strcmp(argv[i], "--reml-bivar-prevalence") == 0) {
            vector<double> K_buf;
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                K_buf.push_back(atof(argv[i]));
            }
            i--;
            if (K_buf.size() < 1 || K_buf.size() > 2) throw ("\nError: --reml-bivar-prevalence. Please specify the prevalences of the two diseases.");
            if (K_buf.size() == 2) {
                if (K_buf[0] < 0.0 || K_buf[0] > 1.0 || K_buf[1] < 0.0 || K_buf[1] > 1.0) throw ("\nError: --reml-bivar-prevalence. Disease prevalence should be betwen 0 and 1.");
                cout << "--reml-bivar-prevalence " << K_buf[0] << " " << K_buf[1] << endl;
                prevalence = K_buf[0];
                prevalence2 = K_buf[1];
            } else {
                if (K_buf[0] < 0.0 || K_buf[0] > 1.0) throw ("\nError: --reml-bivar-prevalence. Disease prevalence should be betwen 0 and 1.");
                cout << "--reml-bivar-prevalence " << K_buf[0] << endl;
                prevalence = prevalence2 = K_buf[0];
            }
        } else if (strcmp(argv[i], "--reml-bivar-nocove") == 0) {
            ignore_Ce = true;
            cout << "--reml-bivar-nocove" << endl;
        } else if (strcmp(argv[i], "--reml-bivar-lrt-rg") == 0) {
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                fixed_rg_val.push_back(atof(argv[i]));
            }
            i--;
            cout << "--reml-bivar-lrt-rg ";
            bool err_flag = false;
            for (j = 0; j < fixed_rg_val.size(); j++) {
                cout << fixed_rg_val[j] << " ";
                if (fixed_rg_val[j] > 1.0 || fixed_rg_val[j]<-1.0) err_flag = true;
            }
            cout << endl;
            if (err_flag || fixed_rg_val.empty()) throw ("\nError: --reml-bivar-lrt-rg. Any input paramter should be within the range from -1 to 1.\n");
            bool haveZero = false;
            if (CommFunc::FloatEqual(fixed_rg_val[0], 0.0)) haveZero = true;
            for (j = 1; j < fixed_rg_val.size(); j++) {
                if ((CommFunc::FloatNotEqual(fixed_rg_val[0], 0.0) && haveZero) || (CommFunc::FloatEqual(fixed_rg_val[0], 0.0) && !haveZero)) throw ("\nError: --reml-bivar-lrt-rg. Input paramters should be all zero or all non-zero values.\n");
            }
        } else if (strcmp(argv[i], "--reml-bivar-no-constrain") == 0) {
            bivar_no_constrain = true;
            cout << "--reml-bivar-no-constrain" << endl;
        } else if (strcmp(argv[i], "--cojo-file") == 0) {
            massoc_file = argv[++i];
            cout << "--cojo-file " << massoc_file << endl;
            CommFunc::FileExist(massoc_file);
        } else if (strcmp(argv[i], "--cojo-slct") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 0;
            cout << "--cojo-slct" << endl;
        } else if (strcmp(argv[i], "--cojo-stepwise") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 0;
            cout << "--cojo-stepwise" << endl;
        } else if (strcmp(argv[i], "--cojo-forward") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 1;
            cout << "--cojo-forward" << endl;
        } else if (strcmp(argv[i], "--cojo-backward") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 2;
            cout << "--cojo-backward" << endl;
        } else if (strcmp(argv[i], "--cojo-top-SNPs") == 0) {
            massoc_slct_flag = true;
            massoc_top_SNPs = atoi(argv[++i]);
            cout << "--cojo-top-SNPs " << massoc_top_SNPs << endl;
            if (massoc_top_SNPs < 1 || massoc_top_SNPs > 10000) throw ("\nError: --cojo-top-SNPs should be within the range from 1 to 10000.\n");
        } else if (strcmp(argv[i], "--cojo-actual-geno") == 0) {
            massoc_actual_geno_flag = true;
            cout << "--cojo-actual-geno" << endl;
        } else if (strcmp(argv[i], "--cojo-p") == 0) {
            massoc_p = atof(argv[++i]);
            cout << "--cojo-p " << massoc_p << endl;
            if (massoc_p > 0.05 || massoc_p <= 0) throw ("\nError: --cojo-p should be within the range from 0 to 0.05.\n");
        } else if (strcmp(argv[i], "--cojo-collinear") == 0) {
            massoc_collinear = atof(argv[++i]);
            cout << "--cojo-collinear " << massoc_collinear << endl;
            if (massoc_collinear > 0.99 || massoc_collinear < 0.01) throw ("\nError: --cojo-collinear should be within the ragne from 0.01 to 0.99.\n");
        } else if (strcmp(argv[i], "--cojo-wind") == 0) {
            massoc_wind = atoi(argv[++i]);
            cout << "--cojo-wind " << massoc_wind << endl;

            // debug
            if (massoc_wind > 100000) throw ("\nError: invalid value for --cojo-wind. Valid range: 100 ~ 100000\n");

            //if (massoc_wind < 100 || massoc_wind > 100000) throw ("\nError: invalid value for --cojo-wind. Valid range: 100 ~ 100000\n");
            massoc_wind *= 1000;
        } else if (strcmp(argv[i], "--cojo-joint") == 0) {
            massoc_joint_flag = true;
            cout << "--cojo-joint" << endl;
        } else if (strcmp(argv[i], "--cojo-cond") == 0) {
            massoc_cond_snplist = argv[++i];
            cout << "--cojo-cond " << massoc_cond_snplist << endl;
        } else if (strcmp(argv[i], "--cojo-gc") == 0) {
            massoc_gc_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                massoc_gc_val = -1;
                i--;
            } else {
                massoc_gc_val = atof(argv[i]);
                if (massoc_gc_val < 1 || massoc_gc_val > 10) throw ("\nError: invalid value specified after --cojo-gc.\n");
            }
            cout << "--cojo-gc " << ((massoc_gc_val < 0) ? "" : argv[i]) << endl;
        } else if (strcmp(argv[i], "--cojo-sblup") == 0) {
            massoc_sblup_flag = true;
            massoc_sblup_fac = atof(argv[++i]);
            cout << "--cojo-sblup " << massoc_sblup_fac << endl;
            if (massoc_sblup_fac < 0) throw ("\nError: invalid value for --cojo-sblup.\n");
        } else if (strcmp(argv[i], "--mlma") == 0) {
            reml_flag = false;
            mlma_flag = true;
            thread_flag = true;
            cout << "--mlm-assoc " << endl;
        } else if (strcmp(argv[i], "--mlma-loco") == 0) {
            reml_flag = false;
            mlma_loco_flag = true;
            thread_flag = true;
            cout << "--mlma-loco " << endl;
        } else if (strcmp(argv[i], "--mlma-no-adj-covar") == 0) {
            mlma_no_adj_covar = true;
            cout << "--mlma-no-adj-covar " << endl;
        } else if (strcmp(argv[i], "--fst") == 0) {
            subpopu_file = argv[++i];
            cout << "--fst " << subpopu_file << endl;
            CommFunc::FileExist(subpopu_file);
        } else if (strcmp(argv[i], "--sbat") == 0) {
            sbat_sAssoc_file = argv[++i];
            cout << "--sbat " << sbat_sAssoc_file << endl;
            CommFunc::FileExist(sbat_sAssoc_file);
        } else if (strcmp(argv[i], "--sbat-gene-list") == 0) {
            sbat_gAnno_file = argv[++i];
            cout << "--sbat-gene-list " << sbat_gAnno_file << endl;
            CommFunc::FileExist(sbat_gAnno_file);
        } else if (strcmp(argv[i], "--sbat-set-list") == 0) {
            sbat_snpset_file = argv[++i];
            cout << "--sbat-set-list " << sbat_snpset_file << endl;
            CommFunc::FileExist(sbat_snpset_file);
        } else if (strcmp(argv[i], "--sbat-wind") == 0) {
            sbat_wind = atoi(argv[++i]);
            cout << "--sbat-wind " << sbat_wind << endl;
            if (sbat_wind < 20 || sbat_wind > 1000) throw ("\nError: invalid value for --sbat-wind. Valid range: 20 ~ 1000\n");
            sbat_wind *= 1000;
        }
        else if (strcmp(argv[i], "--efile") == 0) {
            efile = argv[++i];
            efile_flag = true;
            cout << "--efile " << efile << endl;
            CommFunc::FileExist(efile);
        } 
        else if (strcmp(argv[i], "--e-cor") == 0) {
            eR_file = argv[++i];
            eR_file_flag = true;
            cout << "--e-cor " << eR_file << endl;
            CommFunc::FileExist(eR_file);
        } 
        else if (strcmp(argv[i], "--ecojo") == 0) {
            ecojo_ma_file = argv[++i];
            cout << "--ecojo " << ecojo_ma_file << endl;
            CommFunc::FileExist(ecojo_ma_file);
        } 
        else if (strcmp(argv[i], "--ecojo-slct") == 0) {
            ecojo_slct_flag = true;
            cout << "--ecojo-slct" << endl;
        } 
        else if (strcmp(argv[i], "--ecojo-p") == 0) {
            ecojo_p = atof(argv[++i]);
            cout << "--ecojo-p " << ecojo_p << endl;
            if (ecojo_p > 0.05 || ecojo_p <= 0) throw ("\nError: --ecojo-p should be within the range from 0 to 0.05.\n");
        } 
        else if (strcmp(argv[i], "--ecojo-collinear") == 0) {
            ecojo_collinear = atof(argv[++i]);
            cout << "--ecojo-collinear " << ecojo_collinear << endl;
            if (ecojo_collinear > 1 || ecojo_collinear < 0.01) throw ("\nError: --ecojo-collinear should be within the ragne from 0.01 to 0.99.\n");
        }
        else if (strcmp(argv[i], "--ecojo-blup") == 0) {
            ecojo_blup_flag = true;
            ecojo_lambda = atof(argv[++i]);
            cout << "--ecojo-blup " << ecojo_lambda << endl;
            if (ecojo_lambda < 0.01 || ecojo_lambda > 0.99) throw ("\nError: --ecojo-blup should be within the ragne from 0.01 to 0.99.\n");
        } 
        else if (strcmp(argv[i], "--make-erm") == 0) {
            make_erm_flag = true;
            thread_flag = true;
            cout << argv[i] << endl;
        }
        else if (strcmp(argv[i], "--make-erm-gz") == 0) {
            make_erm_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            cout << "--make-erm-gz" << endl;
        }
        else if (strcmp(argv[i], "--make-erm-alg") == 0) {
            make_erm_flag = true;
            make_erm_mtd = atoi(argv[++i]);
            thread_flag = true;
            cout << "--make-erm-alg " << make_erm_mtd << endl;
            if (make_erm_mtd < 1 || make_erm_mtd > 3) throw ("\nError: --make-erm-alg should be 1, 2 or 3.\n");
        } 
        else if (strcmp(argv[i], "gcta") == 0) break;
        else {
            stringstream errmsg;
            errmsg << "\nError: invalid option \"" << argv[i] << "\".\n";
            throw (errmsg.str());
        }
    }
    // conflicted options
    cout << endl;
    if (bfile2_flag && !bfile_flag) throw ("Error: the option --bfile2 should always go with the option --bfile.");
    if(bfile_flag && grm_cutoff>-1.0) throw ("Error: the --grm-cutoff option is invalid when used in combined with the --bfile option.");
    if (m_grm_flag) {
        if (grm_flag) {
            grm_flag = false;
            cout << "Warning: --grm option suppressed by the --mgrm option." << endl;
        }
        if (grm_cutoff>-1.0) {
            grm_cutoff = -2.0;
            cout << "Warning: --grm-cutoff option suppressed by the --mgrm option." << endl;
        }
    }
    if (pca_flag) {
        if (grm_adj_fac>-1.0) {
            grm_adj_fac = -2.0;
            cout << "Warning: --grm-adj option suppressed by the --pca option." << endl;
        } else if (dosage_compen>-1) {
            grm_adj_fac = -2;
            cout << "Warning: --dosage-compen option suppressed by the --pca option." << endl;
        }
    }
    if (!gxe_file.empty() && !grm_flag && !m_grm_flag) {
        cout << "Warning: --gxe option is ignored because there is no --grm or --mgrm option specified." << endl;
        gxe_file = "";
    }
    if (pred_rand_eff && !grm_flag && !m_grm_flag) {
        cout << "Warning: --reml-pred-rand option is ignored because there is no --grm or --mgrm option specified." << endl;
        pred_rand_eff = false;
    }
    if (dosage_compen>-1 && update_sex_file.empty()) throw ("Error: you need to specify the sex information for the individuals by the option --update-sex because of the option --dc.");
    if (bfile2_flag && update_freq_file.empty()) throw ("Error: you need to update the allele frequency by the option --update-freq because there are two datasets.");
    if ((dose_beagle_flag || dose_mach_flag || dose_mach_gz_flag) && dominance_flag) throw ("Error: unable to calculate the GRM for dominance effect using imputed dosage data.");
    if (make_grm_xchar_flag && dominance_flag) throw ("Error: unable to calculate the GRM for dominance effect for the X chromosome.");
    if (mlma_flag || mlma_loco_flag) {
        if (!gxe_file.empty()) cout << "Warning: the option --gxe option is disabled in this analysis." << endl;
        if (!update_sex_file.empty()) cout << "Warning: the option --update-sex option is disabled in this analysis." << endl;
        if (grm_adj_fac>-1.0) cout << "Warning: the option --grm-adj option is disabled in this analysis." << endl;
        if (dosage_compen>-1.0) cout << "Warning: the option --dc option is disabled in this analysis." << endl;
        if (est_fix_eff) cout << "Warning: the option --reml-est-fix option is disabled in this analysis." << endl;
        if (pred_rand_eff) cout << "Warning: the option --reml-pred-rand option is disabled in this analysis." << endl;
        if (reml_mtd != 0) cout << "Warning: the option --reml-alg option is disabled in this analysis. The default algorithm AI-REML is used." << endl;
        if (reml_lrt_flag) cout << "Warning: the option --reml-lrt option is disabled in this analysis." << endl;
    }
    if(bivar_reml_flag && prevalence_flag) throw("Error: --prevalence option is not compatible with --reml-bivar option. Please check the --reml-bivar-prevalence option!");

    // OpenMP
    stringstream ss;
    ss << thread_num;
    setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
    omp_set_num_threads(thread_num);
    if (thread_flag) {
        if (thread_num == 1) cout << "Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to speed up the computation if there are multiple processors in your machine." << endl;
        else cout << "Note: the program will be running on " << thread_num << " threads." << endl;
    }

    // set autosome
    if (autosome_flag) {
        extract_chr_start = 1;
        extract_chr_end = autosome_num;
    }
    if (make_grm_xchar_flag) extract_chr_start = extract_chr_end = (autosome_num + 1);

    // Implement
    cout << endl;
    gcta *pter_gcta = new gcta(autosome_num, rm_high_ld_cutoff, out); //, *pter_gcta2=new gcta(autosome_num, rm_high_ld_cutoff, out);
    if (grm_bin_flag || m_grm_bin_flag) pter_gcta->enable_grm_bin_flag();
    //if(simu_unlinked_flag) pter_gcta->simu_geno_unlinked(simu_unlinked_n, simu_unlinked_m, simu_unlinked_maf);
    if (!RG_fname_file.empty()) {
        if (RG_summary_file.empty()) throw ("Error: please input the summary information for the raw data files by the option --raw-summary.");
        pter_gcta->read_IRG_fnames(RG_summary_file, RG_fname_file, GC_cutoff);
    } 
    else if(efile_flag){
        pter_gcta->read_efile(efile);
        if (make_erm_flag) pter_gcta->make_erm(make_erm_mtd - 1, grm_out_bin_flag);
        else if(ecojo_slct_flag) pter_gcta->run_ecojo_slct(ecojo_ma_file, ecojo_p, ecojo_collinear);
        else if(ecojo_blup_flag) pter_gcta->run_ecojo_blup_efile(ecojo_ma_file, ecojo_lambda);
    }
    else if(eR_file_flag){
        pter_gcta->read_eR(eR_file);
        pter_gcta->run_ecojo_blup_eR(ecojo_ma_file, ecojo_lambda);
    }
    else if (bfile_flag) {
        if (hapmap_genet_dst) pter_gcta->genet_dst(bfile, hapmap_genet_dst_file);
        else {
            if (bfile2_flag) {
                cout << "There are two datasets specified (in PLINK binary PED format).\nReading dataset 1 ..." << endl;
                if (update_freq_file.empty()) throw ("Error: since there are two dataset, you should update the allele frequencies that are calculated in the combined dataset.");
            }
            pter_gcta->read_famfile(bfile + ".fam");
            if (!kp_indi_file.empty()) pter_gcta->keep_indi(kp_indi_file);
            if (!rm_indi_file.empty()) pter_gcta->remove_indi(rm_indi_file);
            if (!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
            if (!blup_indi_file.empty()) pter_gcta->read_indi_blup(blup_indi_file);
            pter_gcta->read_bimfile(bfile + ".bim");
            if (!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
            if (!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
            if (extract_chr_start > 0) pter_gcta->extract_chr(extract_chr_start, extract_chr_end);
            if (!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
            if (!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
            if (!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
            if (LD) pter_gcta->read_LD_target_SNPs(LD_file);
            pter_gcta->read_bedfile(bfile + ".bed");
            if (!update_impRsq_file.empty()) pter_gcta->update_impRsq(update_impRsq_file);
            if (!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
            if (dose_Rsq_cutoff > 0.0) pter_gcta->filter_impRsq(dose_Rsq_cutoff);
            if (maf > 0) pter_gcta->filter_snp_maf(maf);
            if (max_maf > 0.0) pter_gcta->filter_snp_max_maf(max_maf);
            if (out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
            else if (!paa_file.empty()) pter_gcta->paa(paa_file);
            else if (ibc) pter_gcta->ibc(ibc_all);
            else if (make_grm_flag){
                if(make_grm_ldwt_mtd == 2) pter_gcta->make_grm_pca(dominance_flag, make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, ldwt_wind, false);
                else pter_gcta->make_grm(dominance_flag, make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, false, make_grm_ldwt_mtd, i_ld_file, ldwt_wind, LD_rsq_cutoff, make_grm_f3_flag);
            }
            else if (recode || recode_nomiss) pter_gcta->save_XMat(recode_nomiss);
            else if (LD) pter_gcta->LD_Blocks(LD_step, LD_wind, LD_sig, LD_i, save_ram);
            else if (LD_prune_rsq>-1.0) pter_gcta->LD_pruning_mkl(LD_prune_rsq, LD_wind);
            else if (ld_mean_rsq_flag) pter_gcta->calcu_mean_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
            else if (ld_max_rsq_flag) pter_gcta ->calcu_max_ld_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
            else if (blup_snp_flag) pter_gcta->blup_snp_geno();
            else if (mlma_flag) pter_gcta->mlma(grm_file, phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, within_family, make_grm_inbred_flag, mlma_no_adj_covar);
            else if (mlma_loco_flag) pter_gcta->mlma_loco(phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
            else if (massoc_slct_flag | massoc_joint_flag) pter_gcta->run_massoc_slct(massoc_file, massoc_wind, massoc_p, massoc_collinear, massoc_top_SNPs, massoc_joint_flag, massoc_gc_flag, massoc_gc_val, massoc_actual_geno_flag, massoc_mld_slct_alg);
            else if (!massoc_cond_snplist.empty()) pter_gcta->run_massoc_cond(massoc_file, massoc_cond_snplist, massoc_wind, massoc_collinear, massoc_gc_flag, massoc_gc_val, massoc_actual_geno_flag);
            else if (massoc_sblup_flag) pter_gcta->run_massoc_sblup(massoc_file, massoc_wind, massoc_sblup_fac);
            else if (simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_seed, simu_output_causal, simu_emb_flag);
            else if (make_bed_flag) pter_gcta->save_plink();
            else if (!subpopu_file.empty()) pter_gcta->Fst(subpopu_file);
            else if (!sbat_sAssoc_file.empty()){
                if(!sbat_gAnno_file.empty()) pter_gcta->sbat_gene(sbat_sAssoc_file, sbat_gAnno_file, sbat_wind);
                if(!sbat_snpset_file.empty()) pter_gcta->sbat(sbat_sAssoc_file, sbat_snpset_file);
            }
        }
    } else if (dose_beagle_flag || dose_mach_flag || dose_mach_gz_flag) {
        if (massoc_slct_flag | massoc_joint_flag | !massoc_cond_snplist.empty()) throw ("Error: the --dosage option can't be used in combined with the --cojo options.");
        if (dose_mach_flag) pter_gcta->read_imp_info_mach(dose_info_file);
        else if (dose_mach_gz_flag) pter_gcta->read_imp_info_mach_gz(dose_info_file);
        else if (dose_beagle_flag) pter_gcta->read_imp_info_beagle(dose_info_file);
        if (!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
        if (!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
        if (!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
        if (!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
        if (extract_chr_start > 0) cout << "Warning: the option --chr, --autosome or --nonautosome is inactive for dosage data." << endl;
        if (!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
        if (dose_mach_flag) pter_gcta->read_imp_dose_mach(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
        else if (dose_mach_gz_flag) pter_gcta->read_imp_dose_mach_gz(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
        else if (dose_beagle_flag) pter_gcta->read_imp_dose_beagle(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
        if (!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
        if (!update_impRsq_file.empty()) pter_gcta->update_impRsq(update_impRsq_file);
        if (!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
        if (dose_Rsq_cutoff > 0.0) pter_gcta->filter_impRsq(dose_Rsq_cutoff);
        if (maf > 0.0) pter_gcta->filter_snp_maf(maf);
        if (max_maf > 0.0) pter_gcta->filter_snp_max_maf(max_maf);
        if (out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
        else if (make_grm_flag){
            if(make_grm_ldwt_mtd == 2) pter_gcta->make_grm_pca(dominance_flag, make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, ldwt_wind, false);
            else pter_gcta->make_grm(dominance_flag, make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, false, make_grm_ldwt_mtd, i_ld_file, ldwt_wind, LD_rsq_cutoff, make_grm_f3_flag);
        }
        else if (recode || recode_nomiss) pter_gcta->save_XMat(recode_nomiss);
        else if (LD_prune_rsq>-1.0) pter_gcta->LD_pruning_mkl(LD_prune_rsq, LD_wind);
        else if (ld_mean_rsq_flag) pter_gcta->calcu_mean_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
        else if (ld_max_rsq_flag) pter_gcta ->calcu_max_ld_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
        else if (blup_snp_flag) pter_gcta->blup_snp_dosage();
        else if (massoc_sblup_flag) pter_gcta->run_massoc_sblup(massoc_file, massoc_wind, massoc_sblup_fac);
        else if (simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_seed, simu_output_causal, simu_emb_flag);
        else if (make_bed_flag) pter_gcta->save_plink();
        else if (!subpopu_file.empty()) pter_gcta->Fst(subpopu_file);
        else if (mlma_flag) pter_gcta->mlma(grm_file, phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, within_family, make_grm_inbred_flag, mlma_no_adj_covar);
        else if (mlma_loco_flag) pter_gcta->mlma_loco(phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
    } else if (HE_reg_flag) pter_gcta->HE_reg(grm_file, phen_file, kp_indi_file, rm_indi_file, mphen);
    else if ((reml_flag || bivar_reml_flag) && phen_file.empty()) throw ("\nError: phenotype file is required for reml analysis.\n");
    else if (bivar_reml_flag) {
        pter_gcta->fit_bivar_reml(grm_file, phen_file, qcovar_file, covar_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, mphen2, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, prevalence2, no_constrain, ignore_Ce, fixed_rg_val, bivar_no_constrain);
    } else if (reml_flag) {
        pter_gcta->fit_reml(grm_file, phen_file, qcovar_file, covar_file, qgxe_file, gxe_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, no_constrain, mlma_flag, within_family, reml_bending, reml_diag_one);
    } else if (grm_flag || m_grm_flag) {
        if (pca_flag) pter_gcta->pca(grm_file, kp_indi_file, rm_indi_file, grm_cutoff, m_grm_flag, out_pc_num);
        else if (make_grm_flag) pter_gcta->save_grm(grm_file, kp_indi_file, rm_indi_file, update_sex_file, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, grm_out_bin_flag);
    } else throw ("Error: no analysis has been launched by the option(s).\n");

    delete pter_gcta;
}
