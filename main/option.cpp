/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * GCTA options
 *
 * 2010-2017 by Jian Yang <jian.yang@uq.edu.au> and others
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 *
 * Mocked by Zhili.
 *
 * This folder contains most majority functions of GCTA.
 * We are moving toward version 2 by huge mocks on original version.
 * 
 */

#include <stdlib.h>
#include "gcta.h"
#include "Logger.h"

void option(int option_num, char* option_str[]);

int main_v1(int argc, char* argv[])
{
    LOGGER << "*******************************************************************" << endl;
    LOGGER << "* Genome-wide Complex Trait Analysis (GCTA)" << endl;
    LOGGER << "* version 1.30 beta2" << endl;
    LOGGER << "* (C) 2010-2017, The University of Queensland" << endl;
    LOGGER << "* MIT License" << endl;
    LOGGER << "* Please report bugs to: Jian Yang <jian.yang@uq.edu.au>" << endl;
    LOGGER << "*******************************************************************" << endl;

    long int time_used = 0, start = time(NULL);
    time_t curr = time(0);
    LOGGER << "Analysis started: " << ctime(&curr) << endl;
    LOGGER << "Options:" << endl;
    try {
        option(argc, argv);
    } catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    } catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }
    curr = time(0);
    LOGGER << "\nAnalysis finished: " << ctime(&curr);
    time_used = time(NULL) - start;
    LOGGER << "Computational time: " << time_used / 3600 << ":" << (time_used % 3600) / 60 << ":" << time_used % 60 << endl;

    return 0;
}

void option(int option_num, char* option_str[])
{
    int i = 0, j = 0;

    // OpenMP
    bool thread_flag = false;
    int thread_num = omp_get_max_threads();

    // raw genotype data
    string RG_fname_file = "", RG_summary_file = "";
    double GC_cutoff = 0.7;

    // data management
    string bfile = "", bfile2 = "", bfile_list = "", update_sex_file = "", update_freq_file = "", update_refA_file = "", kp_indi_file = "", rm_indi_file = "", extract_snp_file = "", exclude_snp_file = "", extract_snp_name = "", exclude_snp_name = "", out = "gcta";
    bool SNP_major = false, make_bed_flag = false, dose_mach_flag = false, dose_mach_gz_flag = false, dose_beagle_flag = false, bfile2_flag = false, out_freq_flag = false, out_ssq_flag = false;
    bool ref_A = false, recode = false, recode_nomiss = false, recode_std = false, save_ram = false, autosome_flag = false;
    int bfile_flag = 0, autosome_num = 22, extract_chr_start = 0, extract_chr_end = 0, extract_region_chr = 0, extract_region_bp = 0, extract_region_wind = 0, exclude_region_chr = 0, exclude_region_bp = 0, exclude_region_wind = 0;
    string dose_file = "", dose_info_file = "", update_impRsq_file = "";
    double maf = 0.0, max_maf = 0.0, dose_Rsq_cutoff = 0.0;
    vector<string> multi_bfiles;

    // GRM
    bool ibc = false, ibc_all = false, grm_flag = false, grm_bin_flag = true, m_grm_flag = false, m_grm_bin_flag = true, make_grm_flag = false, make_grm_inbred_flag = false, dominance_flag = false, make_grm_xchar_flag = false, grm_out_bin_flag = true, make_grm_f3_flag = false;
    bool align_grm_flag = false;
    bool pca_flag = false, pcl_flag = false;
    bool project_flag = false;
    double grm_adj_fac = -2.0, grm_cutoff = -2.0, rm_high_ld_cutoff = -1.0, bK_threshold = -10.0;
    int dosage_compen = -2, out_pc_num = 20, make_grm_mtd = 0;
    string grm_file = "", paa_file = "", pc_file = "";
    //pca projection
    string project_file = "";
    int project_N = 0;


    // LD
    string LD_file = "", ld_score_multi_file = "";
    bool LD = false, LD_search = false, LD_i = false, ld_score_flag = false, ld_max_rsq_flag = false, ld_mean_rsq_seg_flag = false, ldscore_adj_flag = false;
    int LD_step = 10;
    double LD_wind = 1e7, LD_sig = 0.05, LD_prune_rsq = -1.0, LD_rsq_cutoff = 0.0, LD_seg = 1e5;

    // initialize paramters for simulation based on real genotype data
    bool simu_qt_flag = false, simu_cc = false, simu_emb_flag = false, simu_output_causal = false;
    int simu_rep = 1, simu_case_num = 0, simu_control_num = 0, simu_eff_mod = 0;
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
    bool prevalence_flag = false, reml_force_inv_fac_flag = false, reml_force_converge_flag = false, reml_no_converge_flag = false, reml_fixed_var_flag = false;
    int mphen = 1, mphen2 = 2, reml_mtd = 0, MaxIter = 100;
    double prevalence = -2.0, prevalence2 = -2.0;
    bool reml_flag = false, pred_rand_eff = false, est_fix_eff = false, blup_snp_flag = false, no_constrain = false, reml_lrt_flag = false, no_lrt = false, bivar_reml_flag = false, ignore_Ce = false, within_family = false, reml_bending = false, HE_reg_flag = false, reml_diag_one = false, bivar_no_constrain = false;
    bool cv_blup = false;
    bool HE_reg_bivar_flag = false;
    string phen_file = "", qcovar_file = "", covar_file = "", qgxe_file = "", gxe_file = "", blup_indi_file = "";
    vector<double> reml_priors, reml_priors_var, fixed_rg_val;
    vector<int> reml_drop;
    reml_drop.push_back(1);

    // Joint analysis of GWAS MA
    string massoc_file = "", massoc_init_snplist = "", massoc_cond_snplist = "";
    int massoc_wind = 1e7, massoc_top_SNPs = -1, massoc_mld_slct_alg = 0;
    double massoc_p = 5e-8, massoc_collinear = 0.9, massoc_sblup_fac = -1, massoc_gc_val = -1;
    bool massoc_slct_flag = false, massoc_joint_flag = false, massoc_sblup_flag = false, massoc_gc_flag = false, massoc_actual_geno_flag = false;
    double massoc_out_pC_thresh = -1;

    // mixed linear model association 
    bool mlma_flag = false, mlma_loco_flag = false, mlma_no_adj_covar = false;
    string subtract_grm_file = "";

    // Fst
    bool fst_flag = false;
    string subpopu_file = "";

    // gene-based association test
    bool sbat_seg_flag = false;
    double sbat_ld_cutoff = sqrt(0.9); //option to remove overly correlated snps in SBAT test
    bool sbat_write_snpset = false; //write snplist - used in conjunction with sbat_ld_cutoff
    string sbat_sAssoc_file = "", sbat_gAnno_file = "", sbat_snpset_file = "";
    int sbat_wind = 50000, sbat_seg_size = 1e5;

    // gene expression data
    string efile="", eR_file = "", ecojo_ma_file="";
    int make_erm_mtd = 1;
    double ecojo_p = 5e-6, ecojo_collinear = 0.9, ecojo_lambda = -1;
    bool efile_flag=false, eR_file_flag = false, ecojo_slct_flag = false, ecojo_blup_flag = false, make_erm_flag = false;

    // mtCOJO
    char chbuf = '\0';
    string mtcojolist_file="", mtcojo_bxy_file="", ref_ld_dirt="", w_ld_dirt="";
    int nsnp_gsmr=10;
    double freq_thresh = 0.2, gwas_thresh=5e-8, global_heidi_thresh = 0.0, indi_heidi_thresh = 0.01, ld_fdr_thresh=0.05, clump_wind_size=10000, clump_r2_thresh=0.05;
    bool mtcojo_flag=false, ref_ld_flag=false, w_ld_flag=false;

    // GSMR
    bool gsmr_flag = false, o_snp_instru_flag = false, gsmr_so_flag = false, gsmr_snp_update_flag = false;
    int gsmr_alg_flag = 0, gsmr_so_alg = -9;
    string expo_file_list = "", outcome_file_list = "";
    
    int argc = option_num;
    vector<char *> argv(option_num + 2);
    for (i = 0; i < option_num; i++) argv[i] = option_str[i];
    argv[option_num] = const_cast<char*>("gcta");
    argv[option_num + 1] = const_cast<char*>("gcta");
    LOGGER << "Accepted options:" << endl;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--thread-num") == 0) {
            thread_num = atoi(argv[++i]);
            LOGGER << "--thread-num " << thread_num << endl;
            if (thread_num < 1 || thread_num > 1000) LOGGER.e(0, "\n  --thread-num should be from 1 to 1000.\n");
        }
        else if (strcmp(argv[i], "--threads") == 0) {
            thread_num = atoi(argv[++i]);
            LOGGER << "--threads " << thread_num << endl;
            if (thread_num < 1 || thread_num > 1000) LOGGER.e(0, "\n  --threads should be from 1 to 1000.\n");
        }// raw genotype data
        else if (strcmp(argv[i], "--raw-files") == 0) {
            RG_fname_file = argv[++i];
            LOGGER << "--raw-files " << argv[i] << endl;
        } else if (strcmp(argv[i], "--raw-summary") == 0) {
            RG_summary_file = argv[++i];
            LOGGER << "--raw-summary " << argv[i] << endl;
        } else if (strcmp(argv[i], "--gencall") == 0) {
            GC_cutoff = atof(argv[++i]);
            LOGGER << "--gencall " << GC_cutoff << endl;
            if (GC_cutoff < 0.0 || GC_cutoff > 1.0) LOGGER.e(0, "\n  --gencall should be within the range from 0 to 1.\n");
        }            // data management
        else if (strcmp(argv[i], "--bfile") == 0) {
            bfile_flag = 1;
            bfile = argv[++i];
            LOGGER << "--bfile " << argv[i] << endl;
        } else if (strcmp(argv[i], "--mbfile") == 0) {
            bfile_flag = 2;
            bfile_list = argv[++i];
            LOGGER << "--mbfile " << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-bed") == 0) {
            make_bed_flag = true;
            LOGGER << "--make-bed " << endl;
        } else if (strcmp(argv[i], "--bfile2") == 0) {
            bfile2_flag = true;
            bfile2 = argv[++i];
            LOGGER << "--bfile2 " << argv[i] << endl;
        } else if (strcmp(argv[i], "--dosage-mach") == 0) {
            dose_mach_flag = true;
            dose_beagle_flag = false;
            dose_file = argv[++i];
            dose_info_file = argv[++i];
            LOGGER << "--dosage-mach " << dose_file << " " << dose_info_file << endl;
        } else if (strcmp(argv[i], "--dosage-mach-gz") == 0) {
            dose_mach_gz_flag = true;
            dose_beagle_flag = false;
            dose_file = argv[++i];
            dose_info_file = argv[++i];
            LOGGER << "--dosage-mach-gz " << dose_file << " " << dose_info_file << endl;
        } else if (strcmp(argv[i], "--dosage-beagle") == 0) {
            dose_beagle_flag = true;
            dose_mach_flag = false;
            dose_mach_gz_flag = false;
            dose_file = argv[++i];
            dose_info_file = argv[++i];
            LOGGER << "--dosage-beagle " << dose_file << " " << dose_info_file << endl;
        } else if (strcmp(argv[i], "--imput-rsq") == 0) {
            dose_Rsq_cutoff = atof(argv[++i]);
            LOGGER << "--imput-rsq " << dose_Rsq_cutoff << endl;
            if (dose_Rsq_cutoff < 0.0 || dose_Rsq_cutoff > 1.0) LOGGER.e(0, "\n  --imput-rsq should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--update-imput-rsq") == 0) {
            update_impRsq_file = argv[++i];
            LOGGER << "--update-imput-rsq " << update_impRsq_file << endl;
            CommFunc::FileExist(update_impRsq_file);
        } else if (strcmp(argv[i], "--update-freq") == 0) {
            update_freq_file = argv[++i];
            LOGGER << "--update-freq " << update_freq_file << endl;
            CommFunc::FileExist(update_freq_file);
        } else if (strcmp(argv[i], "--update-ref-allele") == 0) {
            update_refA_file = argv[++i];
            LOGGER << "--update-ref-allele " << update_refA_file << endl;
            CommFunc::FileExist(update_refA_file);
        } else if (strcmp(argv[i], "--keep") == 0) {
            kp_indi_file = argv[++i];
            LOGGER << "--keep " << kp_indi_file << endl;
            CommFunc::FileExist(kp_indi_file);
        } else if (strcmp(argv[i], "--remove") == 0) {
            rm_indi_file = argv[++i];
            LOGGER << "--remove " << rm_indi_file << endl;
            CommFunc::FileExist(rm_indi_file);
        } else if (strcmp(argv[i], "--update-sex") == 0) {
            update_sex_file = argv[++i];
            LOGGER << "--update-sex " << update_sex_file << endl;
            CommFunc::FileExist(update_sex_file);
        } else if (strcmp(argv[i], "--chr") == 0) {
            extract_chr_start = extract_chr_end = atoi(argv[++i]);
            LOGGER << "--chr " << extract_chr_start << endl;
            if (extract_chr_start < 1 || extract_chr_start > 100) LOGGER.e(0, "\n  --chr should be within the range from 1 to 100.\n");
        } else if (strcmp(argv[i], "--autosome-num") == 0) {
            autosome_num = atoi(argv[++i]);
            LOGGER << "--autosome-num " << autosome_num << endl;
            if (autosome_num < 1 || autosome_num > 100) LOGGER.e(0, "\n  invalid number specified after the option --autosome-num.\n");
        } else if (strcmp(argv[i], "--autosome") == 0) {
            autosome_flag = true;
            LOGGER << "--autosome" << endl;
        } else if (strcmp(argv[i], "--extract") == 0) {
            extract_snp_file = argv[++i];
            LOGGER << "--extract " << extract_snp_file << endl;
            CommFunc::FileExist(extract_snp_file);
        } else if (strcmp(argv[i], "--exclude") == 0) {
            exclude_snp_file = argv[++i];
            LOGGER << "--exclude " << exclude_snp_file << endl;
            CommFunc::FileExist(exclude_snp_file);
        } else if (strcmp(argv[i], "--extract-snp") == 0) {
            extract_snp_name = argv[++i];
            LOGGER << "--extract-snp " << extract_snp_name << endl;
        } else if (strcmp(argv[i], "--extract-region-snp") == 0) {
            extract_snp_name = argv[++i];
            extract_region_wind = atoi(argv[++i]);
            LOGGER << "--extract-region-snp " << extract_snp_name << " " << extract_region_wind << "Kb" << endl;
            extract_region_wind *= 1000;
            if(extract_region_wind < 1000 || extract_region_wind > 1e8) LOGGER.e(0, "\n  the second paramter of --extract-region is distance in Kb unit. It should take value between 1 and 1e5.");
        } else if (strcmp(argv[i], "--extract-region-bp") == 0) {
            extract_region_chr = atoi(argv[++i]);
            extract_region_bp = atoi(argv[++i]);
            extract_region_wind = atoi(argv[++i]);
            LOGGER << "--extract-region-bp " << extract_region_chr << " " << extract_region_bp << " " << extract_region_wind << "Kb" << endl;
            extract_region_wind *= 1000;
            if(extract_region_wind < 1000 || extract_region_wind > 1e8) LOGGER.e(0, "\n  the second paramter of --extract-region is distance in Kb unit. It should take value between 1 and 1e5.");
        } else if (strcmp(argv[i], "--exclude-snp") == 0) {
            exclude_snp_name = argv[++i];
            LOGGER << "--exclude-snp " << exclude_snp_name << endl;
        } else if (strcmp(argv[i], "--exclude-region-snp") == 0) {
            exclude_snp_name = argv[++i];
            exclude_region_wind = atoi(argv[++i]);
            LOGGER << "--exclude-region-snp " << exclude_snp_name << exclude_region_wind << "Kb" << endl;
            exclude_region_wind *= 1000;
            if(exclude_region_wind < 1000 || exclude_region_wind > 1e8) LOGGER.e(0, "\n  the second paramter of --exclude-region is distance in Kb unit. It should take value between 1 and 1e5.");
        } else if (strcmp(argv[i], "--exclude-region-bp") == 0) {
            exclude_region_chr = atoi(argv[++i]);
            exclude_region_bp = atoi(argv[++i]);
            exclude_region_wind = atoi(argv[++i]);
            LOGGER << "--exclude-region-bp " << exclude_region_chr << " " << exclude_region_bp << " " << exclude_region_wind << "Kb" << endl;
            exclude_region_wind *= 1000;
            if(exclude_region_wind < 1000 || exclude_region_wind > 1e8) LOGGER.e(0, "\n  the second paramter of --exclude-region is distance in Kb unit. It should take value between 1 and 1e5.");
        } else if (strcmp(argv[i], "--maf") == 0) {
            maf = atof(argv[++i]);
            LOGGER << "--maf " << maf << endl;
            if (maf < 0 || maf > 0.5) LOGGER.e(0, "\n  --maf should be within the range from 0 to 0.5.\n");
        } else if (strcmp(argv[i], "--max-maf") == 0) {
            max_maf = atof(argv[++i]);
            LOGGER << "--max-maf " << max_maf << endl;
            if (max_maf <= 0) LOGGER.e(0, "\n  --max-maf should be > 0.\n");
        } else if (strcmp(argv[i], "--out") == 0) {
            out = argv[++i];
            LOGGER << "--out " << out << endl;
        } else if (strcmp(argv[i], "--freq-v1") == 0) {
            out_freq_flag = true;
            thread_flag = true;
            LOGGER << "--freq-v1" << endl;
        } else if (strcmp(argv[i], "--freq") == 0) {
            out_freq_flag = true;
            thread_flag = true;
            LOGGER << "--freq" << endl;
        } else if (strcmp(argv[i], "--ssq") == 0) {
            out_ssq_flag = true;
            LOGGER << "--ssq" << endl;
        } else if (strcmp(argv[i], "--recode") == 0) {
            recode = true;
            thread_flag = true;
            LOGGER << "--recode" << endl;
        } else if (strcmp(argv[i], "--recode-nomiss") == 0) {
            recode_nomiss = true;
            thread_flag = true;
            LOGGER << "--recode-nomiss" << endl;
        } else if (strcmp(argv[i], "--recode-std") == 0) {
            recode_std = true;
            thread_flag = true;
            LOGGER << "--recode-std" << endl;
        } else if (strcmp(argv[i], "--save-ram") == 0) {
            save_ram = true;
            LOGGER << "--save-ram" << endl;
        }// GRM
        else if (strcmp(argv[i], "--paa") == 0) {
            paa_file = argv[++i];
            LOGGER << "--paa " << paa_file << endl;
            CommFunc::FileExist(paa_file);
        } else if (strcmp(argv[i], "--ibc") == 0) {
            ibc = true;
            LOGGER << "--ibc" << endl;
        } else if (strcmp(argv[i], "--ibc-all") == 0) {
            ibc = ibc_all = true;
            LOGGER << "--ibc-all" << endl;
        } else if (strcmp(argv[i], "--mgrm") == 0 || strcmp(argv[i], "--mgrm-bin") == 0) {
            m_grm_flag = true;
            grm_file = argv[++i];
            LOGGER << argv[i - 1] << " " << grm_file << endl;
        } else if (strcmp(argv[i], "--mgrm-gz") == 0) {
            m_grm_flag = true;
            m_grm_bin_flag = false;
            grm_bin_flag = false;
            grm_file = argv[++i];
            LOGGER << "--mgrm-gz " << grm_file << endl;
        } else if (strcmp(argv[i], "--grm") == 0 || strcmp(argv[i], "--grm-bin") == 0) {
            grm_flag = true;
            grm_file = argv[++i];
            LOGGER << argv[i - 1] << " " << grm_file << endl;
        } else if (strcmp(argv[i], "--grm-gz") == 0) {
            grm_flag = true;
            m_grm_bin_flag = false;
            grm_bin_flag = false;
            grm_file = argv[++i];
            LOGGER << "--grm-gz " << grm_file << endl;
        } else if (strcmp(argv[i], "--rm-high-ld") == 0) {
            rm_high_ld_cutoff = atof(argv[++i]);
            LOGGER << "--rm-high-ld " << rm_high_ld_cutoff << endl;
            if (rm_high_ld_cutoff <= 0 || rm_high_ld_cutoff >= 1) LOGGER.e(0, "\n  the value to be specified after --rm-high-ld should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--make-grm") == 0 || strcmp(argv[i], "--make-grm-v1") == 0 || strcmp(argv[i], "--make-grm-bin") == 0) {
            make_grm_flag = true;
            thread_flag = true;
            LOGGER << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-gz") == 0) {
            make_grm_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            LOGGER << "--make-grm-gz" << endl;
        } else if (strcmp(argv[i], "--make-grm-alg") == 0) {
            make_grm_flag = true;
            make_grm_mtd = atoi(argv[++i]);
            thread_flag = true;
            LOGGER << "--make-grm-alg " << make_grm_mtd << endl;
            if (make_grm_mtd < 0 || make_grm_mtd > 1) LOGGER.e(0, "\n  --make-grm-alg should be 0 or 1.\n");
        } else if (strcmp(argv[i], "--make-grm-f3") == 0) {
            make_grm_flag = true;
            make_grm_f3_flag = true;
            grm_out_bin_flag = true;
            thread_flag = true;
            LOGGER << "--make-grm-f3" << endl;
        } else if (strcmp(argv[i], "--make-grm-d-v1") == 0 || strcmp(argv[i], "--make-grm-d-bin") == 0) {
            make_grm_flag = true;
            dominance_flag = true;
            thread_flag = true;
            LOGGER << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-d-gz") == 0) {
            make_grm_flag = true;
            dominance_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            LOGGER << "--make-grm-d-gz" << endl;
        } else if (strcmp(argv[i], "--dominance") == 0) {
            dominance_flag = true;
            thread_flag = true;
            LOGGER <<"--dominance"<< endl;
        } else if (strcmp(argv[i], "--make-grm-xchr") == 0 || strcmp(argv[i], "--make-grm-xchr-bin") == 0) {
            make_grm_flag = true;
            make_grm_xchar_flag = true;
            thread_flag = true;
            LOGGER << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-xchr-gz") == 0) {
            make_grm_flag = true;
            make_grm_xchar_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            LOGGER << "--make-grm-xchr-gz" << endl;
        } else if (strcmp(argv[i], "--make-grm-inbred") == 0 || strcmp(argv[i], "--make-grm-inbred-bin") == 0) {
            make_grm_flag = true;
            make_grm_inbred_flag = true;
            thread_flag = true;
            LOGGER << argv[i] << endl;
        } else if (strcmp(argv[i], "--make-grm-inbred-gz") == 0) {
            make_grm_flag = true;
            grm_out_bin_flag = false;
            make_grm_inbred_flag = true;
            thread_flag = true;
            LOGGER << "--make-grm-inbred-gz" << endl;
        } else if (strcmp(argv[i], "--grm-adj") == 0) {
            grm_adj_fac = atof(argv[++i]);
            LOGGER << "--grm-adj " << grm_adj_fac << endl;
            if (grm_adj_fac < 0 || grm_adj_fac > 1) LOGGER.e(0, "\n  the value to be specified after --grm-adj should be within the range from 0 to 1.\n");
        } else if (strcmp(argv[i], "--dc") == 0) {
            dosage_compen = atoi(argv[++i]);
            LOGGER << "--dc " << dosage_compen << endl;
            if (dosage_compen != 0 && dosage_compen != 1) LOGGER.e(0, "\n  the value to be specified after --dc should be 0 or 1.\n");
        } else if (strcmp(argv[i], "--grm-cutoff") == 0 || strcmp(argv[i], "--grm-cutoff-v1") == 0) {
            grm_cutoff = atof(argv[++i]);
            if (grm_cutoff >= -1 && grm_cutoff <= 2) LOGGER << "--grm-cutoff" << grm_cutoff << endl;
            else grm_cutoff = -2;
        } else if (strcmp(argv[i], "--grm-align") == 0) {
            align_grm_flag = true;
            thread_flag = true;
        } else if (strcmp(argv[i], "--make-bK") == 0) {
            bK_threshold = atof(argv[++i]);
            if (bK_threshold < 0 || bK_threshold > 1) LOGGER.e(0, "\n  --make-bK threshold should be range from 0 to 1.\n");
            else LOGGER << "--make-bK " << bK_threshold << endl;
        } else if (strcmp(argv[i], "--pca") == 0) {
            pca_flag = true;
            thread_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                out_pc_num = 20;
                i--;
            } else out_pc_num = atoi(argv[i]);
            LOGGER << "--pca " << out_pc_num << endl;
            if (out_pc_num < 1) LOGGER.e(0, "\n  the value to be specified after --pca should be positive.\n");
        } else if (strcmp(argv[i], "--pc-loading") == 0) {
            pcl_flag = true;
            thread_flag = true;
            pc_file = argv[++i];
            //pcl_grm_N = atoi(argv[++i]);
            LOGGER << "--pc-loading " << pc_file <<  endl;
            //if(pcl_grm_N < 1 || pcl_grm_N > 1e20) LOGGER.e(0, "\n  invalid number of SNPs used to calculate PCs."); 
        }else if (strcmp(argv[i], "--project-loading") == 0 ){
            project_flag = true;
            thread_flag = true;
            project_file = argv[++i];
            project_N = atoi(argv[++i]);
            LOGGER << "--project-loading " << project_file << " " << project_N << endl;
            if(project_N < 1 || project_N > 1e3) LOGGER.e(0, "\n  invalid number of PCs to output");
        }
        // estimation of LD structure
        else if (strcmp(argv[i], "--ld") == 0) {
            LD = true;
            LD_file = argv[++i];
            LOGGER << "--ld " << LD_file << endl;
            CommFunc::FileExist(LD_file);
        } else if (strcmp(argv[i], "--ld-step") == 0) {
            LD_search = true;
            LD_step = atoi(argv[++i]);
            LOGGER << "--ld-step " << LD_step << endl;
            if (LD_step < 1 || LD_step > 20) LOGGER.e(0, "\n  --ld-step should be within the range from 1 to 20.\n");
        } else if (strcmp(argv[i], "--ld-wind") == 0 || strcmp(argv[i], "--ld-pruning-wind") == 0 || strcmp(argv[i], "--make-grm-wt-wind") == 0) {
            LD_wind = atof(argv[++i]);
            LOGGER << argv[i - 1] << " " << LD_wind << endl;
            LD_wind *= 1000;
            if (LD_wind < 1e3 || LD_wind > 2e7) {
                stringstream err_msg;
                err_msg << "\n  " << argv[i - 1] << " should be 1Kb or 20Mb.\n";
                LOGGER.e(0, err_msg.str());
            }
        } else if (strcmp(argv[i], "--ld-sig") == 0) {
            LD_sig = atof(argv[++i]);
            LOGGER << "--ld-sig " << LD_sig << endl;
            if (LD_sig <= 0) LOGGER.e(0, "\n  --ld-sig should be > 0.\n");
        } else if (strcmp(argv[i], "--ld-i") == 0) {
            LD_i = true;
            LOGGER << "--ld-i" << endl;
        } else if (strcmp(argv[i], "--ld-pruning") == 0) {
            thread_flag = true;
            LD_prune_rsq = atof(argv[++i]);
            LOGGER << "--ld-pruning " << LD_prune_rsq << endl;
            if (LD_prune_rsq < 0.0001 || LD_prune_rsq > 0.9999) LOGGER.e(0, "\n  --ld-pruning should be within the range from 0.0001 to 0.9999.\n");
        } else if (strcmp(argv[i], "--ld-score") == 0) {
            ld_score_flag = true;
            thread_flag = true;
            LOGGER << "--ld-score" << endl;
        } else if (strcmp(argv[i], "--ld-score-adj") == 0) {
            ldscore_adj_flag = true;
            LOGGER << "--ld-score-adj" << endl;
        } else if (strcmp(argv[i], "--ld-score-multi") == 0) {
            ld_score_flag = true;
            thread_flag = true;
            ld_score_multi_file = argv[++i];
            LOGGER << "--ld-score-multi " << ld_score_multi_file << endl;
            CommFunc::FileExist(ld_score_multi_file);
        } else if (strcmp(argv[i], "--ld-rsq-cutoff") == 0) {
            LD_rsq_cutoff = atof(argv[++i]);
            LOGGER << "--ld-rsq-cutoff " << LD_rsq_cutoff << endl;
            if (LD_rsq_cutoff < 0.0 || LD_rsq_cutoff > 1.0) {
                stringstream err_msg;
                err_msg << "\n  " << argv[i - 1] << " should be within the range from 0 to 1.\n";
                LOGGER.e(0, err_msg.str());
            }
        } else if (strcmp(argv[i], "--ld-max-rsq") == 0) {
            ld_max_rsq_flag = true;
            thread_flag = true;
            LOGGER << "--ld-max-rsq" << endl;
        } else if (strcmp(argv[i], "--ld-score-region") == 0) {
            ld_mean_rsq_seg_flag = true;
            thread_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                LD_seg = 200;
                i--;
            } else LD_seg = atoi(argv[i]);
            LOGGER << "--ld-score-region" << endl;
            if (LD_seg < 10) LOGGER.e(0, "\n  input value for --ld-score-region needs to be > 10.\n");
            LD_seg *= 1000;
        } else if (strcmp(argv[i], "--ld-file") == 0) {
            LD_file = argv[++i];
            LOGGER << "--ld-file " << LD_file << endl;
        }
        // simulation based on real genotype data
        else if (strcmp(argv[i], "--simu-qt") == 0) {
            simu_qt_flag = true;
            LOGGER << "--simu-qt" << endl;
        } else if (strcmp(argv[i], "--simu-cc") == 0) {
            simu_cc = true;
            simu_case_num = atoi(argv[++i]);
            simu_control_num = atoi(argv[++i]);
            LOGGER << "--simu-cc " << simu_case_num << " " << simu_control_num << endl;
            if (simu_case_num < 10) LOGGER.e(0, "--simu-cc, Invalid number of cases. Minimun number 10.");
            if (simu_control_num < 10) LOGGER.e(0, "--simu-cc, Invalid number of controls. Minimum number 10.");
        } else if (strcmp(argv[i], "--simu-rep") == 0) {
            simu_rep = atoi(argv[++i]);
            LOGGER << "--simu-rep " << simu_rep << endl;
            if (simu_rep < 1 || simu_rep > 10000) LOGGER.e(0, "--simu-rep should be within the range from 1 to 10000.");
        } else if (strcmp(argv[i], "--simu-hsq") == 0) {
            simu_h2 = atof(argv[++i]);
            LOGGER << "--simu-hsq " << simu_h2 << endl;
            if (simu_h2 > 1.0 || simu_h2 < 0.0) LOGGER.e(0, "--simu-h2 should be within the range from 0 to 1.");
        } else if (strcmp(argv[i], "--simu-k") == 0) {
            simu_K = atof(argv[++i]);
            LOGGER << "--simu-k " << simu_K << endl;
            if (simu_K > 0.5 || simu_K < 0.0001) LOGGER.e(0, "--simu-K should be within the range from 0.0001 to 0.5.");
        } else if (strcmp(argv[i], "--simu-causal-loci") == 0) {
            simu_causal = argv[++i];
            LOGGER << "--simu-causal-loci " << simu_causal << endl;
            CommFunc::FileExist(simu_causal);
        } else if (strcmp(argv[i], "--simu-embayesb") == 0) { // internal
            simu_emb_flag = true;
            LOGGER << "--simu-embayesb" << endl;
        } else if (strcmp(argv[i], "--simu-ouput-causal") == 0) { // internal
            simu_output_causal = true;
            LOGGER << "--simu-output-causal" << endl;
        } else if (strcmp(argv[i], "--simu-seed") == 0) {
            simu_seed = atof(argv[++i]);
            LOGGER << "--simu-seed " << simu_seed << endl;
            if (simu_seed <= 100) LOGGER.e(0, "--simu-seed should be >100.");
        } else if (strcmp(argv[i], "--simu-eff-mod") == 0) {
            simu_eff_mod = atoi(argv[++i]);
            LOGGER << "--simu-eff-mod " << simu_eff_mod << endl;
            if (simu_eff_mod != 0 && simu_eff_mod !=1) LOGGER.e(0, "--simu-eff-mod should be 0 or 1.");
        }
        else if (strcmp(argv[i], "--hapmap-genet-dst") == 0) { // calculate genetic dst based on HapMap data
            hapmap_genet_dst = true;
            hapmap_genet_dst_file = argv[++i];
            LOGGER << "--hapmap-genet-dst " << hapmap_genet_dst_file << endl;
        }// estimate variance explained by all SNPs
        else if (strcmp(argv[i], "--HEreg") == 0) {
            HE_reg_flag = true;
            thread_flag = true;
            LOGGER << "--HEreg" << endl;
        } else if (strcmp(argv[i], "--HEreg-bivar") == 0) {
            HE_reg_bivar_flag = true;
            thread_flag = true;
            vector<int> mphen_buf;
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                mphen_buf.push_back(atoi(argv[i]));
            }
            i--;
            if (mphen_buf.size() < 2 && mphen_buf.size() > 0) LOGGER.e(0, "\n  --HEreg-bivar. Please specify two traits for the HE regression for covariance analysis.");
            if (mphen_buf.size() == 0) {
                mphen = 1;
                mphen2 = 2;
            } else {
                mphen = mphen_buf[0];
                mphen2 = mphen_buf[1];
            }
            if (mphen < 1 || mphen2 < 1 || mphen == mphen2) LOGGER.e(0, "\n --HEreg-bivar. Invalid input parameters.");
            LOGGER << "--HEreg-bivar " << mphen << " " << mphen2 << endl;
        } else if (strcmp(argv[i], "--reml") == 0) {
            reml_flag = true;
            thread_flag = true;
            LOGGER << "--reml" << endl;
            if (m_grm_flag) no_lrt = true;
        } else if (strcmp(argv[i], "--prevalence") == 0) {
            prevalence_flag = true;
            prevalence = atof(argv[++i]);
            LOGGER << "--prevalence " << prevalence << endl;
            if (prevalence <= 0 || prevalence >= 1) LOGGER.e(0, "\n  --prevalence should be between 0 to 1.\n");
        } else if (strcmp(argv[i], "--reml-pred-rand") == 0) {
            pred_rand_eff = true;
            LOGGER << "--reml-pred-rand" << endl;
        } else if(strcmp(argv[i], "--cvblup") == 0){
            cv_blup = true;
            LOGGER << "--cvblup" << endl;
        } else if (strcmp(argv[i], "--reml-est-fix") == 0) {
            est_fix_eff = true;
            LOGGER << "--reml-est-fix" << endl;
        } else if (strcmp(argv[i], "--reml-alg") == 0) {
            reml_mtd = atoi(argv[++i]);
            LOGGER << "--reml-alg " << reml_mtd << endl;
            if (reml_mtd < 0 || reml_mtd > 2) LOGGER.e(0, "\n  --reml-alg should be 0, 1 or 2.\n");
        } else if (strcmp(argv[i], "--reml-no-constrain") == 0) {
            reml_flag = true;
            no_constrain = true;
            LOGGER << "--reml-no-constrain" << endl;
        } else if (strcmp(argv[i], "--reml-priors") == 0) {
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                reml_priors.push_back(atof(argv[i]));
            }
            i--;
            LOGGER << "--reml-priors ";
            bool err_flag = false;
            for (j = 0; j < reml_priors.size(); j++) {
                LOGGER << reml_priors[j] << " ";
                if (reml_priors[j] > 1.0 || reml_priors[j] < -10.0) err_flag = true;
            }
            LOGGER << endl;
            if (err_flag || reml_priors.empty()) LOGGER.e(0, "\n  --reml-priors. Prior values of variance explained should be between 0 and 1.\n");
        } else if (strcmp(argv[i], "--reml-priors-var") == 0  || strcmp(argv[i], "--reml-fixed-var") == 0) {
            string s_buf = argv[i];
            if(s_buf == "--reml-fixed-var") reml_fixed_var_flag = true;
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                reml_priors_var.push_back(atof(argv[i]));
            }
            i--;
            LOGGER << s_buf << " ";
            bool err_flag = false;
            for (j = 0; j < reml_priors_var.size(); j++) {
                LOGGER << reml_priors_var[j] << " ";
                if (reml_priors_var[j] < 0.0) err_flag = true;
            }
            LOGGER << endl;
            if (reml_priors_var.empty()) LOGGER.e(0, "\n  " + s_buf + ". Prior values of variance components are required.\n");
        } else if (strcmp(argv[i], "--reml-no-lrt") == 0) {
            no_lrt = true;
            LOGGER << "--reml-no-lrt" << endl;
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
            LOGGER << "--reml-lrt ";
            bool err_flag = false;
            for (j = 0; j < reml_drop.size(); j++) {
                LOGGER << reml_drop[j] << " ";
                if (reml_drop[j] < 1) err_flag = true;
            }
            LOGGER << endl;
            if (err_flag || reml_drop.empty()) LOGGER.e(0, "\n  invalid values specified after --reml-lrt.\n");
        } else if (strcmp(argv[i], "--reml-maxit") == 0) {
            MaxIter = atoi(argv[++i]);
            LOGGER << "--reml-maxit " << MaxIter << endl;
            if (MaxIter < 1 || MaxIter > 10000) LOGGER.e(0, "\n  --reml-maxit should be within the range from 1 to 10000.\n");
        } else if (strcmp(argv[i], "--reml-bendV") == 0) {
            reml_force_inv_fac_flag = true;
            LOGGER << "--reml-bendV " << endl;
        } else if (strcmp(argv[i], "--reml-force-converge") == 0) {
            reml_force_converge_flag = true;
            LOGGER << "--reml-force-converge " << endl;
        } else if (strcmp(argv[i], "--reml-allow-no-converge") == 0) {
            reml_no_converge_flag = true;
            LOGGER << "--reml-allow-no-converge " << endl;
        } else if (strcmp(argv[i], "--reml-bending") == 0) {
            reml_bending = true;
            LOGGER << "--reml-bending " << endl;
        } else if (strcmp(argv[i], "--reml-diag-one") == 0) {
            reml_diag_one = true;
            LOGGER << "--reml-diag-one " << endl;
        } else if (strcmp(argv[i], "--pheno") == 0) {
            phen_file = argv[++i];
            LOGGER << "--pheno " << phen_file << endl;
            CommFunc::FileExist(phen_file);
        } else if (strcmp(argv[i], "--mpheno") == 0) {
            mphen = atoi(argv[++i]);
            LOGGER << "--mpheno " << mphen << endl;
            if (mphen < 1) LOGGER.e(0, "--mpheno should be > 0.");
        } else if (strcmp(argv[i], "--qcovar") == 0) {
            qcovar_file = argv[++i];
            LOGGER << "--qcovar " << qcovar_file << endl;
            CommFunc::FileExist(qcovar_file);
        } else if (strcmp(argv[i], "--covar") == 0) {
            covar_file = argv[++i];
            LOGGER << "--covar " << covar_file << endl;
            CommFunc::FileExist(covar_file);
        } else if (strcmp(argv[i], "--gxqe") == 0) {
            qgxe_file = argv[++i];
            LOGGER << "--gxqe " << qgxe_file << endl;
            CommFunc::FileExist(qgxe_file);
        } else if (strcmp(argv[i], "--gxe") == 0) {
            gxe_file = argv[++i];
            LOGGER << "--gxe " << gxe_file << endl;
            CommFunc::FileExist(gxe_file);
        } else if (strcmp(argv[i], "--blup-snp") == 0) {
            blup_snp_flag = true;
            blup_indi_file = argv[++i];
            LOGGER << "--blup-snp " << blup_indi_file << endl;
            CommFunc::FileExist(blup_indi_file);
        } else if (strcmp(argv[i], "--reml-wfam") == 0) {
            reml_flag = true;
            within_family = true;
            LOGGER << "--reml-wfam " << endl;
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
            if (mphen_buf.size() < 2 && mphen_buf.size() > 0) LOGGER.e(0, "\n  --reml-bivar. Please specify two traits for the bivariate REML analysis.");
            if (mphen_buf.size() == 0) {
                mphen = 1;
                mphen2 = 2;
            } else {
                mphen = mphen_buf[0];
                mphen2 = mphen_buf[1];
            }
            if (mphen < 1 || mphen2 < 1 || mphen == mphen2) LOGGER.e(0, "\n  --reml-bivar. Invalid input parameters.");
            LOGGER << "--reml-bivar " << mphen << " " << mphen2 << endl;
        } else if (strcmp(argv[i], "--reml-bivar-prevalence") == 0) {
            vector<double> K_buf;
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                K_buf.push_back(atof(argv[i]));
            }
            i--;
            if (K_buf.size() < 1 || K_buf.size() > 2) LOGGER.e(0, "\n  --reml-bivar-prevalence. Please specify the prevalences of the two diseases.");
            if (K_buf.size() == 2) {
                if (K_buf[0] < 0.0 || K_buf[0] > 1.0 || K_buf[1] < 0.0 || K_buf[1] > 1.0) LOGGER.e(0, "\n  --reml-bivar-prevalence. Disease prevalence should be between 0 and 1.");
                LOGGER << "--reml-bivar-prevalence " << K_buf[0] << " " << K_buf[1] << endl;
                prevalence = K_buf[0];
                prevalence2 = K_buf[1];
            } else {
                if (K_buf[0] < 0.0 || K_buf[0] > 1.0) LOGGER.e(0, "\n  --reml-bivar-prevalence. Disease prevalence should be between 0 and 1.");
                LOGGER << "--reml-bivar-prevalence " << K_buf[0] << endl;
                prevalence = prevalence2 = K_buf[0];
            }
        } else if (strcmp(argv[i], "--reml-bivar-nocove") == 0) {
            ignore_Ce = true;
            LOGGER << "--reml-bivar-nocove" << endl;
        } else if (strcmp(argv[i], "--reml-bivar-lrt-rg") == 0) {
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                fixed_rg_val.push_back(atof(argv[i]));
            }
            i--;
            LOGGER << "--reml-bivar-lrt-rg ";
            bool err_flag = false;
            for (j = 0; j < fixed_rg_val.size(); j++) {
                LOGGER << fixed_rg_val[j] << " ";
                if (fixed_rg_val[j] > 1.0 || fixed_rg_val[j]<-1.0) err_flag = true;
            }
            LOGGER << endl;
            if (err_flag || fixed_rg_val.empty()) LOGGER.e(0, "\n  --reml-bivar-lrt-rg. Any input paramter should be within the range from -1 to 1.\n");
            bool haveZero = false;
            if (CommFunc::FloatEqual(fixed_rg_val[0], 0.0)) haveZero = true;
            for (j = 1; j < fixed_rg_val.size(); j++) {
                if ((CommFunc::FloatNotEqual(fixed_rg_val[0], 0.0) && haveZero) || (CommFunc::FloatEqual(fixed_rg_val[0], 0.0) && !haveZero)) LOGGER.e(0, "\n  --reml-bivar-lrt-rg. Input paramters should be all zero or all non-zero values.\n");
            }
        } else if (strcmp(argv[i], "--reml-bivar-no-constrain") == 0) {
            bivar_no_constrain = true;
            LOGGER << "--reml-bivar-no-constrain" << endl;
        } else if (strcmp(argv[i], "--cojo-file") == 0) {
            massoc_file = argv[++i];
            LOGGER << "--cojo-file " << massoc_file << endl;
            CommFunc::FileExist(massoc_file);
        } else if (strcmp(argv[i], "--cojo-slct") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 0;
            LOGGER << "--cojo-slct" << endl;
        } else if (strcmp(argv[i], "--cojo-stepwise") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 0;
            LOGGER << "--cojo-stepwise" << endl;
        } else if (strcmp(argv[i], "--cojo-forward") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 1;
            LOGGER << "--cojo-forward" << endl;
        } else if (strcmp(argv[i], "--cojo-backward") == 0) {
            massoc_slct_flag = true;
            massoc_mld_slct_alg = 2;
            LOGGER << "--cojo-backward" << endl;
        } else if (strcmp(argv[i], "--cojo-top-SNPs") == 0) {
            massoc_slct_flag = true;
            massoc_top_SNPs = atoi(argv[++i]);
            LOGGER << "--cojo-top-SNPs " << massoc_top_SNPs << endl;
            if (massoc_top_SNPs < 1 || massoc_top_SNPs > 10000) LOGGER.e(0, "\n  --cojo-top-SNPs should be within the range from 1 to 10000.\n");
        } else if (strcmp(argv[i], "--cojo-actual-geno") == 0) {
            massoc_actual_geno_flag = true;
            LOGGER << "--cojo-actual-geno" << endl;
        } else if (strcmp(argv[i], "--cojo-p") == 0) {
            massoc_p = atof(argv[++i]);
            LOGGER << "--cojo-p " << massoc_p << endl;
            if (massoc_p > 0.05 || massoc_p <= 0) LOGGER.e(0, "\n  --cojo-p should be within the range from 0 to 0.05.\n");
        } else if (strcmp(argv[i], "--restrict-output-pC") == 0) {
            massoc_out_pC_thresh = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--cojo-collinear") == 0) {
            massoc_collinear = atof(argv[++i]);
            LOGGER << "--cojo-collinear " << massoc_collinear << endl;
            if (massoc_collinear > 0.99 || massoc_collinear < 0.01) LOGGER.e(0, "\n  --cojo-collinear should be within the ragne from 0.01 to 0.99.\n");
        } else if (strcmp(argv[i], "--cojo-wind") == 0) {
            massoc_wind = atoi(argv[++i]);
            LOGGER << "--cojo-wind " << massoc_wind << endl;

            // debug
            if (massoc_wind > 100000) LOGGER.e(0, "\n  invalid value for --cojo-wind. Valid range: 100 ~ 100000\n");

            //if (massoc_wind < 100 || massoc_wind > 100000) LOGGER.e(0, "\n  invalid value for --cojo-wind. Valid range: 100 ~ 100000\n");
            massoc_wind *= 1000;
        } else if (strcmp(argv[i], "--cojo-joint") == 0) {
            massoc_joint_flag = true;
            LOGGER << "--cojo-joint" << endl;
        } else if (strcmp(argv[i], "--cojo-cond") == 0) {
            massoc_cond_snplist = argv[++i];
            LOGGER << "--cojo-cond " << massoc_cond_snplist << endl;
        } else if (strcmp(argv[i], "--cojo-gc") == 0) {
            massoc_gc_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                massoc_gc_val = -1;
                i--;
            } else {
                massoc_gc_val = atof(argv[i]);
                if (massoc_gc_val < 1 || massoc_gc_val > 10) LOGGER.e(0, "\n  invalid value specified after --cojo-gc.\n");
            }
            LOGGER << "--cojo-gc " << ((massoc_gc_val < 0) ? "" : argv[i]) << endl;
        } else if (strcmp(argv[i], "--cojo-sblup") == 0) {
            massoc_sblup_flag = true;
            massoc_sblup_fac = atof(argv[++i]);
            LOGGER << "--cojo-sblup " << massoc_sblup_fac << endl;
            if (massoc_sblup_fac < 0) LOGGER.e(0, "\n  invalid value for --cojo-sblup.\n");
        } else if (strcmp(argv[i], "--mlma") == 0) {
            reml_flag = false;
            mlma_flag = true;
            thread_flag = true;
            LOGGER << "--mlma " << endl;
        } else if (strcmp(argv[i], "--mlma-subtract-grm") == 0) {
            subtract_grm_file = argv[++i];
            LOGGER << "--mlma-subtract-grm " << subtract_grm_file << endl;
        } else if (strcmp(argv[i], "--mlma-loco") == 0) {
            reml_flag = false;
            mlma_loco_flag = true;
            thread_flag = true;
            LOGGER << "--mlma-loco " << endl;
        } else if (strcmp(argv[i], "--mlma-no-adj-covar") == 0) {
            mlma_no_adj_covar = true;
            LOGGER << "--mlma-no-adj-covar " << endl;
        } else if (strcmp(argv[i], "--fst") == 0) {
            fst_flag = true;
            LOGGER << "--fst " << endl;
        } else if (strcmp(argv[i], "--sub-popu") == 0) {
            subpopu_file = argv[++i];
            LOGGER << "--sub-popu " << subpopu_file << endl;
            CommFunc::FileExist(subpopu_file);
        } else if (strcmp(argv[i], "--fastBAT-ld-cutoff") == 0) {
            sbat_ld_cutoff = sqrt(atof(argv[++i]));
            LOGGER << "--fastBAT-ld-cutoff " << sbat_ld_cutoff << endl;
            if (sbat_ld_cutoff <= 0.1) LOGGER.e(0, "\n  --fastBAT_ld_cutoff should be > 0.1\n");
        } else if (strcmp(argv[i], "--fastBAT-write-snpset") == 0) {
            sbat_write_snpset = true;
            LOGGER << "--fastBAT-write-snpset" << endl;
        } else if (strcmp(argv[i], "--fastBAT") == 0) {
            sbat_sAssoc_file = argv[++i];
            LOGGER << "--fastBAT " << sbat_sAssoc_file << endl;
            CommFunc::FileExist(sbat_sAssoc_file);
        } else if (strcmp(argv[i], "--fastBAT-gene-list") == 0) {
            sbat_gAnno_file = argv[++i];
            LOGGER << "--fastBAT-gene-list " << sbat_gAnno_file << endl;
            CommFunc::FileExist(sbat_gAnno_file);
        } else if (strcmp(argv[i], "--fastBAT-set-list") == 0) {
            sbat_snpset_file = argv[++i];
            LOGGER << "--fastBAT-set-list " << sbat_snpset_file << endl;
            CommFunc::FileExist(sbat_snpset_file);
        } else if (strcmp(argv[i], "--fastBAT-wind") == 0) {
            sbat_wind = atoi(argv[++i]);
            LOGGER << "--fastBAT-wind " << sbat_wind << endl;
            if (sbat_wind < 0 || sbat_wind > 1000) LOGGER.e(0, "\n  invalid value for --fastBAT-wind. Valid range: 0 ~ 1000\n");
            sbat_wind *= 1000;
        } else if (strcmp(argv[i], "--fastBAT-seg") == 0) {
            sbat_seg_flag = true;
            thread_flag = true;
            i++;
            if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) {
                sbat_seg_size = 100;
                i--;
            } else sbat_seg_size = atoi(argv[i]);
            LOGGER << "--fastBAT-seg " << sbat_seg_size << endl;
            if (sbat_seg_size < 10 || sbat_seg_size > 10000) LOGGER.e(0, "\n  invalid value for --fastBAT-seg. Valid range: 10 ~ 10000\n");
            sbat_seg_size *= 1000;
        }
        else if (strcmp(argv[i], "--efile") == 0) {
            efile = argv[++i];
            efile_flag = true;
            LOGGER << "--efile " << efile << endl;
            CommFunc::FileExist(efile);
        } 
        else if (strcmp(argv[i], "--e-cor") == 0) {
            eR_file = argv[++i];
            eR_file_flag = true;
            LOGGER << "--e-cor " << eR_file << endl;
            CommFunc::FileExist(eR_file);
        } 
        else if (strcmp(argv[i], "--ecojo") == 0) {
            ecojo_ma_file = argv[++i];
            LOGGER << "--ecojo " << ecojo_ma_file << endl;
            CommFunc::FileExist(ecojo_ma_file);
        } 
        else if (strcmp(argv[i], "--ecojo-slct") == 0) {
            ecojo_slct_flag = true;
            LOGGER << "--ecojo-slct" << endl;
        } 
        else if (strcmp(argv[i], "--ecojo-p") == 0) {
            ecojo_p = atof(argv[++i]);
            LOGGER << "--ecojo-p " << ecojo_p << endl;
            if (ecojo_p > 0.05 || ecojo_p <= 0) LOGGER.e(0, "\n  --ecojo-p should be within the range from 0 to 0.05.\n");
        } 
        else if (strcmp(argv[i], "--ecojo-collinear") == 0) {
            ecojo_collinear = atof(argv[++i]);
            LOGGER << "--ecojo-collinear " << ecojo_collinear << endl;
            if (ecojo_collinear > 1 || ecojo_collinear < 0.01) LOGGER.e(0, "\n  --ecojo-collinear should be within the ragne from 0.01 to 0.99.\n");
        }
        else if (strcmp(argv[i], "--ecojo-blup") == 0) {
            ecojo_blup_flag = true;
            ecojo_lambda = atof(argv[++i]);
            LOGGER << "--ecojo-blup " << ecojo_lambda << endl;
            if (ecojo_lambda < 0.01 || ecojo_lambda > 0.99) LOGGER.e(0, "\n  --ecojo-blup should be within the ragne from 0.01 to 0.99.\n");
        } 
        else if (strcmp(argv[i], "--make-erm") == 0) {
            make_erm_flag = true;
            thread_flag = true;
            LOGGER << argv[i] << endl;
        }
        else if (strcmp(argv[i], "--make-erm-gz") == 0) {
            make_erm_flag = true;
            grm_out_bin_flag = false;
            thread_flag = true;
            LOGGER << "--make-erm-gz" << endl;
        }
        else if (strcmp(argv[i], "--make-erm-alg") == 0) {
            make_erm_flag = true;
            make_erm_mtd = atoi(argv[++i]);
            thread_flag = true;
            LOGGER << "--make-erm-alg " << make_erm_mtd << endl;
            if (make_erm_mtd < 1 || make_erm_mtd > 3) LOGGER.e(0, "\n  --make-erm-alg should be 1, 2 or 3.\n");
        } else if (strcmp(argv[i], "--gsmr-file") == 0 ) {
            gsmr_flag = true;

            vector<string> gsmr_file_list;
            while (1) {
                i++;
                if (strcmp(argv[i], "gcta") == 0 || strncmp(argv[i], "--", 2) == 0) break;
                gsmr_file_list.push_back(argv[i]);
            }
            i--;
            if (gsmr_file_list.size() < 1 || gsmr_file_list.size() > 2) 
                LOGGER.e(0, "--gsmr-file, please specify the GWAS summary data for the exposure(s) and the outcome(s).");

            expo_file_list = gsmr_file_list[0];
            outcome_file_list = gsmr_file_list[1];
            LOGGER << "--gsmr-file " << expo_file_list << " " << outcome_file_list << endl;
            CommFunc::FileExist(expo_file_list);
            CommFunc::FileExist(outcome_file_list);
        } else if (strcmp(argv[i], "--gsmr-direction") == 0) {
            gsmr_alg_flag = atoi(argv[++i]);
            if(gsmr_alg_flag < 0 || gsmr_alg_flag > 2) 
               LOGGER.e(0, "--gsmr-direction should be 0 (forward-GSMR), 1 (reverse-GSMR) or 2 (bi-GSMR).");
            LOGGER << "--gsmr-direction " << gsmr_alg_flag << endl;
        } else if (strcmp(argv[i], "--gsmr-alg") == 0) {
            LOGGER.e(0, "--gsmr-alg has been superceded by --gsmr-direction.");
        } else if (strcmp(argv[i], "--gsmr-so") == 0) {
            gsmr_so_flag = true;
            //gsmr_so_alg = atoi(argv[++i]);
            if(gsmr_so_alg < 0 || gsmr_so_alg > 1) 
                LOGGER.e(0, "--gsmr-so should be 0 (LD score regression) or 1 (correlation of SNP effects).");
            LOGGER << "--gsmr-so " << gsmr_so_alg << endl;
        } else if (strcmp(argv[i], "--effect-plot") == 0) {
            o_snp_instru_flag = true;
            LOGGER << "--effect-plot" << endl;
        } else if (strcmp(argv[i], "--mtcojo-file") == 0) {
            mtcojo_flag = true;
            mtcojolist_file = argv[++i];
            LOGGER << "--mtcojo-file " << mtcojolist_file << endl;
            CommFunc::FileExist(mtcojolist_file);
        } else if (strcmp(argv[i], "--mtcojo-bxy") == 0) {
            mtcojo_bxy_file = argv[++i];
            LOGGER << "--mtcojo-bxy " << mtcojo_bxy_file << endl;
            CommFunc::FileExist(mtcojo_bxy_file);
        } else if (strcmp(argv[i], "--ref-ld-chr") == 0) {
            ref_ld_flag = true;
            ref_ld_dirt = argv[++i];
            chbuf = ref_ld_dirt.back();

#ifdef _WIN32
	    if(chbuf != '\\') ref_ld_dirt = ref_ld_dirt + '\\';
#elif defined __linux__ || defined __APPLE__
	    if(chbuf != '/') ref_ld_dirt = ref_ld_dirt + '/';
#else
#error Only Windows, Mac and Linux are supported.
#endif
            LOGGER << "--ref-ld-chr " << ref_ld_dirt << endl;
        } else if (strcmp(argv[i], "--w-ld-chr") == 0) {
            w_ld_flag = true;
            w_ld_dirt = argv[++i];
            chbuf = w_ld_dirt.back();
#ifdef _WIN32
            if(chbuf != '\\') w_ld_dirt = w_ld_dirt + '\\';
#elif defined __linux__ || defined __APPLE__
            if(chbuf != '/') w_ld_dirt = w_ld_dirt + '/';
#else
#error Only Windows, Mac and Linux are supported.
#endif

            LOGGER << "--w-ld-chr " << w_ld_dirt << endl;
        } else if (strcmp(argv[i], "--diff-freq") == 0) {
            freq_thresh = atof(argv[++i]);
            if(freq_thresh <0 || freq_thresh >1)
                LOGGER.e(0, "--diff-freq, Invalid threshold for difference of allele frequencies.");
            LOGGER<<"--diff-freq "<<freq_thresh<<endl;
        } else if (strcmp(argv[i], "--gwas-thresh") == 0) {
            gwas_thresh = atof(argv[++i]);
            if(gwas_thresh <0 || gwas_thresh >1)
                LOGGER.e(0, "--gwas-thresh, Invalid p-value threshold for GWAS summary data.");
            LOGGER<<"--gwas-thresh "<<gwas_thresh<<endl;
        } else if (strcmp(argv[i], "--heidi-thresh") == 0) {
            indi_heidi_thresh = atof(argv[++i]);
            if(indi_heidi_thresh <0 || indi_heidi_thresh >1)
                LOGGER.e(0, "--heidi-thresh, Invalid p-value threshold for HEIDI test.");
            LOGGER<<"--heidi-thresh "<<indi_heidi_thresh<<endl;
        } else if (strcmp(argv[i], "--heidi-snp") == 0) {
            LOGGER.e(0, "--heidi-snp is discontinued. Please use --gsmr-snp-min to specify minimum number of SNP instruments for the HEIDI-outlier analysis.");
        } else if ((strcmp(argv[i], "--gsmr-snp") == 0) || (strcmp(argv[i], "--gsmr-snp-min") == 0)) {
            if(strcmp(argv[i], "--gsmr-snp") == 0) gsmr_snp_update_flag = true;
            nsnp_gsmr = atoi(argv[++i]);
            if(nsnp_gsmr < 0 || nsnp_gsmr > 1e6)
                LOGGER.e(0, "--gsmr-snp-min, Invalid SNP number threshold for GSMR analysis.");
            LOGGER<<"--gsmr-snp-min "<<nsnp_gsmr<<endl;
        } else if (strcmp(argv[i], "--gsmr-ld-fdr") == 0) {
            ld_fdr_thresh = atoi(argv[++i]);
            if(ld_fdr_thresh < 0 || ld_fdr_thresh > 1)
                LOGGER.e(0, "--gsmr-ld-fdr, Invalid FDR threshold for LD correlation matrix.");
            LOGGER<<"--gsmr-ld-fdr "<<ld_fdr_thresh<<endl;
        } else if (strcmp(argv[i], "--clump-p1") == 0) {
            LOGGER.e(0, "--clump-p1 is discontinued. Please use --gwas-thresh to specify p-value threshold for index SNPs.");   
        } else if (strcmp(argv[i], "--clump-kb") == 0) {
            clump_wind_size = atof(argv[++i]);
            if(clump_wind_size <0 || clump_wind_size >1e6)
                LOGGER.e(0, "--clump-kb, Invalid window size for clumping analysis.");
            LOGGER<<"--clump-kb   "<<clump_wind_size<<endl;
        } else if (strcmp(argv[i], "--clump-r2") == 0) {
            clump_r2_thresh = atof(argv[++i]);
            if(clump_r2_thresh <0 || clump_r2_thresh >1)
                LOGGER.e(0, "--clump-r2, Invalid LD r2 threshold for clumping analysis.");
            LOGGER<<"--clump-r2 "<<clump_r2_thresh<<endl;
        } else if (strcmp(argv[i], "gcta") == 0) break;
        else {
            stringstream errmsg;
            errmsg << "\n  invalid option \"" << argv[i] << "\".\n";
            LOGGER.e(0, errmsg.str());
        }
    }
    // conflicted options
    LOGGER << endl;
    if (bfile2_flag && !bfile_flag) LOGGER.e(0, "the option --bfile2 should always go with the option --bfile.");
    if(bfile_flag && grm_cutoff>-1.0) LOGGER.e(0, "the --grm-cutoff option is invalid when used in combined with the --bfile option.");
    if (m_grm_flag) {
        if (grm_flag) {
            grm_flag = false;
            LOGGER << "Warning: --grm option suppressed by the --mgrm option." << endl;
        }
        if (grm_cutoff>-1.0) {
            grm_cutoff = -2.0;
            LOGGER << "Warning: --grm-cutoff option suppressed by the --mgrm option." << endl;
        }
    }
    if (pca_flag) {
        if (grm_adj_fac>-1.0) {
            grm_adj_fac = -2.0;
            LOGGER << "Warning: --grm-adj option suppressed by the --pca option." << endl;
        } else if (dosage_compen>-1) {
            grm_adj_fac = -2;
            LOGGER << "Warning: --dosage-compen option suppressed by the --pca option." << endl;
        }
    }
    if (!gxe_file.empty() && !grm_flag && !m_grm_flag) {
        LOGGER << "Warning: --gxe option is ignored because there is no --grm or --mgrm option specified." << endl;
        gxe_file = "";
    }
    if (pred_rand_eff && !grm_flag && !m_grm_flag) {
        LOGGER << "Warning: --reml-pred-rand option is ignored because there is no --grm or --mgrm option specified." << endl;
        pred_rand_eff = false;
    }
    if (cv_blup && !grm_flag && !m_grm_flag) {
        LOGGER << "Warning: --cvblup option is ignored because there is no --grm or --mgrm option specified." << endl;
        cv_blup = false;
    }
    if(cv_blup && pred_rand_eff){
        LOGGER << "Warning: --reml-pred-rand options is ignored because --cvblup does more than this option" << endl;
        pred_rand_eff = false;
    }

    if (dosage_compen>-1 && update_sex_file.empty()) LOGGER.e(0, "you need to specify the sex information for the individuals by the option --update-sex because of the option --dc.");
    if (bfile2_flag && update_freq_file.empty()) LOGGER.e(0, "you need to update the allele frequency by the option --update-freq because there are two datasets.");
    if ((dose_beagle_flag || dose_mach_flag || dose_mach_gz_flag) && dominance_flag) LOGGER.e(0, "unable to calculate the GRM for dominance effect using imputed dosage data.");
    if (make_grm_xchar_flag && dominance_flag) LOGGER.e(0, "unable to calculate the GRM for dominance effect for the X chromosome.");
    if (mlma_flag || mlma_loco_flag) {
        if (!gxe_file.empty()) LOGGER << "Warning: the option --gxe option is disabled in this analysis." << endl;
        if (!update_sex_file.empty()) LOGGER << "Warning: the option --update-sex option is disabled in this analysis." << endl;
        if (grm_adj_fac>-1.0) LOGGER << "Warning: the option --grm-adj option is disabled in this analysis." << endl;
        if (dosage_compen>-1.0) LOGGER << "Warning: the option --dc option is disabled in this analysis." << endl;
        if (est_fix_eff) LOGGER << "Warning: the option --reml-est-fix option is disabled in this analysis." << endl;
        if (pred_rand_eff) LOGGER << "Warning: the option --reml-pred-rand option is disabled in this analysis." << endl;
        if(cv_blup) LOGGER << "Warning: the option --cvblup option is disabled in this analysis." << endl;
        if (reml_lrt_flag) LOGGER << "Warning: the option --reml-lrt option is disabled in this analysis." << endl; 
    }
    if(bivar_reml_flag && prevalence_flag) LOGGER.e(0, "--prevalence option is not compatible with --reml-bivar option. Please check the --reml-bivar-prevalence option!");
    if(gsmr_flag || mtcojo_flag){
        if(ref_ld_flag && !w_ld_flag) LOGGER.e(0, "--ref-ld-chr, please specify the directory of LD score files.");
        if(!ref_ld_flag && w_ld_flag) LOGGER.e(0, "--w-ld-chr, please specify the directory of LD scores for the regression weights.");
        if(gsmr_snp_update_flag) LOGGER.w(0, "--gsmr-snp has been superseded by --gsmr-snp-min.");
        if(nsnp_gsmr < 5) LOGGER.w(0, "The number of SNP instruments included in the analysis is too small. There might not be enough SNPs to perform the HEIDI-outlier analysis.");
        // if(!gsmr_so_flag && ref_ld_flag && w_ld_flag) { gsmr_so_alg = 0; LOGGER.w(0, "--gsmr-so is not specified. The default value is 0. GSMR analysis will perform LD score regression to estimate sample overlap."); }
        // if(gsmr_so_alg == 1 && ref_ld_flag && w_ld_flag) { gsmr_so_alg = 0; LOGGER.w(0, "The LD score regression instead of correlation method will be used to estimate sample overlap."); }
        // if(gsmr_so_alg == 0 && !ref_ld_flag && !w_ld_flag) LOGGER.e(0, "Please specify the directory of LD score files to perform LD score regression analysis.");
        // if(!gsmr_so_flag && !ref_ld_flag && !w_ld_flag) LOGGER.w(0, "The GSMR analysis will be performed assuming no sample overlap between the GWAS data for exposure and outcome.");
    }
    // OpenMP
    if (thread_flag) {
        if (thread_num == 1) LOGGER << "Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to speed up the computation if there are multiple processors in your machine." << endl;
        else LOGGER << "Note: the program will be running on " << thread_num << " threads." << endl;
    }

    // set autosome
    if (autosome_flag) {
        if(extract_chr_start == extract_chr_end && extract_chr_start != 0){
            LOGGER << "Warning: --autosome option omitted. You have specified the chromosome to analysis." << endl;
        }else{
            extract_chr_start = 1;
            extract_chr_end = autosome_num;
        }
    }
    if (make_grm_xchar_flag) extract_chr_start = extract_chr_end = (autosome_num + 1);

    // Implement
    LOGGER << endl;
    gcta *pter_gcta = new gcta(autosome_num, rm_high_ld_cutoff, out); //, *pter_gcta2=new gcta(autosome_num, rm_high_ld_cutoff, out);
    if(ldscore_adj_flag) pter_gcta->set_ldscore_adj_flag(ldscore_adj_flag);
    if(reml_force_inv_fac_flag) pter_gcta->set_reml_force_inv();
    if(reml_force_converge_flag) pter_gcta->set_reml_force_converge();
    if(reml_no_converge_flag) pter_gcta->set_reml_no_converge();
    if(reml_fixed_var_flag) pter_gcta->set_reml_fixed_var();
    if(reml_mtd != 0) pter_gcta->set_reml_mtd(reml_mtd);
    if (grm_bin_flag || m_grm_bin_flag) pter_gcta->enable_grm_bin_flag();
    //if(simu_unlinked_flag) pter_gcta->simu_geno_unlinked(simu_unlinked_n, simu_unlinked_m, simu_unlinked_maf);
    if (!RG_fname_file.empty()) {
        if (RG_summary_file.empty()) LOGGER.e(0, "please input the summary information for the raw data files by the option --raw-summary.");
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
                LOGGER << "There are two datasets specified (in PLINK binary PED format).\nReading dataset 1 ..." << endl;
                if (update_freq_file.empty()) LOGGER.e(0, "since there are two dataset, you should update the allele frequencies that are calculated in the combined dataset.");
            }
            // Read the list, if there are multiple bfiles
            if(bfile_flag==2) multi_bfiles = pter_gcta->read_bfile_list(bfile_list);
            // Start to read the genotypes
            if(bfile_flag==1) pter_gcta->read_famfile(bfile + ".fam");
            else pter_gcta->read_multi_famfiles(multi_bfiles);
            if (!kp_indi_file.empty()) pter_gcta->keep_indi(kp_indi_file);
            if (!rm_indi_file.empty()) pter_gcta->remove_indi(rm_indi_file);
            if (!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
            if (!blup_indi_file.empty()) pter_gcta->read_indi_blup(blup_indi_file);
            if(bfile_flag==1) pter_gcta->read_bimfile(bfile + ".bim");
            else pter_gcta->read_multi_bimfiles(multi_bfiles);
            if (!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
            if (extract_chr_start > 0) pter_gcta->extract_chr(extract_chr_start, extract_chr_end);
            if(extract_region_chr>0) pter_gcta->extract_region_bp(extract_region_chr, extract_region_bp, extract_region_wind);
            if (!extract_snp_name.empty()){
                if(extract_region_wind>0) pter_gcta->extract_region_snp(extract_snp_name, extract_region_wind);
                else pter_gcta->extract_single_snp(extract_snp_name);
            }
            if (!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
            if(exclude_region_chr>0) pter_gcta->exclude_region_bp(exclude_region_chr, exclude_region_bp, exclude_region_wind);
            if (!exclude_snp_name.empty()) {
                if(exclude_region_wind>0) pter_gcta->exclude_region_snp(exclude_snp_name, exclude_region_wind);
                else pter_gcta->exclude_single_snp(exclude_snp_name);
            }
            if (!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
            if (LD) pter_gcta->read_LD_target_SNPs(LD_file);
            if(gsmr_flag) pter_gcta->read_gsmrfile(expo_file_list, outcome_file_list, gwas_thresh, nsnp_gsmr, gsmr_so_alg);
            if(mtcojo_flag) pter_gcta->read_mtcojofile(mtcojolist_file, gwas_thresh, nsnp_gsmr);
            if(bfile_flag==1) pter_gcta->read_bedfile(bfile + ".bed");
            else pter_gcta->read_multi_bedfiles(multi_bfiles);
            if (!update_impRsq_file.empty()) pter_gcta->update_impRsq(update_impRsq_file);
            if (!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
            if (dose_Rsq_cutoff > 0.0) pter_gcta->filter_impRsq(dose_Rsq_cutoff);
            if (maf > 0) pter_gcta->filter_snp_maf(maf);
            if (max_maf > 0.0) pter_gcta->filter_snp_max_maf(max_maf);
            if (out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
            else if (!paa_file.empty()) pter_gcta->paa(paa_file);
            else if (ibc) pter_gcta->ibc(ibc_all);
            else if (make_grm_flag) pter_gcta->make_grm(dominance_flag, make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, false, make_grm_f3_flag, subpopu_file);
            else if (recode || recode_nomiss || recode_std) pter_gcta->save_XMat(recode_nomiss, recode_std);
            else if (LD) pter_gcta->LD_Blocks(LD_step, LD_wind, LD_sig, LD_i, save_ram);
            else if (LD_prune_rsq>-1.0) pter_gcta->LD_pruning_mkl(LD_prune_rsq, LD_wind);
            else if (ld_score_flag){
                if(ld_score_multi_file.empty()) pter_gcta->calcu_mean_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
                else pter_gcta->calcu_mean_rsq_multiSet(ld_score_multi_file, LD_wind, LD_rsq_cutoff, dominance_flag);
            }
            else if (ld_mean_rsq_seg_flag) pter_gcta->ld_seg(LD_file, LD_seg, LD_wind, LD_rsq_cutoff, dominance_flag);
            else if (ld_max_rsq_flag) pter_gcta ->calcu_max_ld_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
            else if (blup_snp_flag) pter_gcta->blup_snp_geno();
            else if (mlma_flag) pter_gcta->mlma(grm_file, m_grm_flag, subtract_grm_file, phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, within_family, make_grm_inbred_flag, mlma_no_adj_covar);
            else if (mlma_loco_flag) pter_gcta->mlma_loco(phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
            else if (massoc_slct_flag | massoc_joint_flag) {pter_gcta->set_massoc_pC_thresh(massoc_out_pC_thresh); pter_gcta->set_diff_freq(freq_thresh);pter_gcta->run_massoc_slct(massoc_file, massoc_wind, massoc_p, massoc_collinear, massoc_top_SNPs, massoc_joint_flag, massoc_gc_flag, massoc_gc_val, massoc_actual_geno_flag, massoc_mld_slct_alg);}
            else if (!massoc_cond_snplist.empty()) {pter_gcta->set_massoc_pC_thresh(massoc_out_pC_thresh); pter_gcta->set_diff_freq(freq_thresh);pter_gcta->run_massoc_cond(massoc_file, massoc_cond_snplist, massoc_wind, massoc_collinear, massoc_gc_flag, massoc_gc_val, massoc_actual_geno_flag);}
            else if (massoc_sblup_flag) {pter_gcta->set_diff_freq(freq_thresh);pter_gcta->run_massoc_sblup(massoc_file, massoc_wind, massoc_sblup_fac);}
            else if (gsmr_flag) pter_gcta->gsmr(gsmr_alg_flag, ref_ld_dirt, w_ld_dirt, freq_thresh, gwas_thresh, clump_wind_size, clump_r2_thresh, global_heidi_thresh, indi_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, o_snp_instru_flag, gsmr_so_alg);
            else if (mtcojo_flag) pter_gcta->mtcojo(mtcojo_bxy_file, ref_ld_dirt, w_ld_dirt, freq_thresh, gwas_thresh, clump_wind_size, clump_r2_thresh, global_heidi_thresh, indi_heidi_thresh, ld_fdr_thresh, nsnp_gsmr);
            else if (simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_seed, simu_output_causal, simu_emb_flag, simu_eff_mod);
            else if (make_bed_flag) pter_gcta->save_plink();
            else if (fst_flag) pter_gcta->Fst(subpopu_file);
            else if (!sbat_sAssoc_file.empty()){
                if(!sbat_gAnno_file.empty()) pter_gcta->sbat_gene(sbat_sAssoc_file, sbat_gAnno_file, sbat_wind, sbat_ld_cutoff, sbat_write_snpset);
                else if(!sbat_snpset_file.empty()) pter_gcta->sbat(sbat_sAssoc_file, sbat_snpset_file, sbat_ld_cutoff, sbat_write_snpset);
                else if(sbat_seg_flag) pter_gcta->sbat_seg(sbat_sAssoc_file, sbat_seg_size, sbat_ld_cutoff, sbat_write_snpset);
            }
            else if(pcl_flag) pter_gcta->snp_pc_loading(pc_file);
            else if(project_flag) pter_gcta->project_loading(project_file, project_N);
        }
    } else if (dose_beagle_flag || dose_mach_flag || dose_mach_gz_flag) {
        if (massoc_slct_flag | massoc_joint_flag | !massoc_cond_snplist.empty()) LOGGER.e(0, "the --dosage option can't be used in combined with the --cojo options.");
        if (dose_mach_flag) pter_gcta->read_imp_info_mach(dose_info_file);
        else if (dose_mach_gz_flag) pter_gcta->read_imp_info_mach_gz(dose_info_file);
        else if (dose_beagle_flag) pter_gcta->read_imp_info_beagle(dose_info_file);
        if (!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
        if (!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
        if (!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
        if (!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
        if (extract_chr_start > 0) LOGGER << "Warning: the option --chr, --autosome or --nonautosome is inactive for dosage data." << endl;
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
        else if (make_grm_flag) pter_gcta->make_grm(dominance_flag, make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, false, make_grm_f3_flag, subpopu_file);
        else if (recode || recode_nomiss || recode_std) pter_gcta->save_XMat(recode_nomiss, recode_std);
        else if (LD_prune_rsq>-1.0) pter_gcta->LD_pruning_mkl(LD_prune_rsq, LD_wind);
        else if (ld_score_flag){
                if(ld_score_multi_file.empty()) pter_gcta->calcu_mean_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
                else pter_gcta->calcu_mean_rsq_multiSet(ld_score_multi_file, LD_wind, LD_rsq_cutoff, dominance_flag);
        }
        else if (ld_max_rsq_flag) pter_gcta ->calcu_max_ld_rsq(LD_wind, LD_rsq_cutoff, dominance_flag);
        else if (blup_snp_flag) pter_gcta->blup_snp_dosage();
        else if (massoc_sblup_flag) {pter_gcta->set_diff_freq(freq_thresh);pter_gcta->run_massoc_sblup(massoc_file, massoc_wind, massoc_sblup_fac);}
        else if (simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_seed, simu_output_causal, simu_emb_flag, simu_eff_mod);
        else if (make_bed_flag) pter_gcta->save_plink();
        else if (fst_flag) pter_gcta->Fst(subpopu_file);
        else if (mlma_flag) pter_gcta->mlma(grm_file, m_grm_flag, subtract_grm_file, phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, within_family, make_grm_inbred_flag, mlma_no_adj_covar);
        else if (mlma_loco_flag) pter_gcta->mlma_loco(phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
    } else if (HE_reg_flag) pter_gcta->HE_reg(grm_file, m_grm_flag, phen_file, kp_indi_file, rm_indi_file, mphen);
    else if (HE_reg_bivar_flag) pter_gcta->HE_reg_bivar(grm_file, m_grm_flag, phen_file, kp_indi_file, rm_indi_file, mphen, mphen2);
    else if ((reml_flag || bivar_reml_flag) && phen_file.empty()) LOGGER.e(0, "\n  phenotype file is required for reml analysis.\n");
    else if (bivar_reml_flag) {
        pter_gcta->set_cv_blup(cv_blup);
        pter_gcta->fit_bivar_reml(grm_file, phen_file, qcovar_file, covar_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, mphen2, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, prevalence2, no_constrain, ignore_Ce, fixed_rg_val, bivar_no_constrain);
    } else if (reml_flag) {
        pter_gcta->set_cv_blup(cv_blup);
        pter_gcta->fit_reml(grm_file, phen_file, qcovar_file, covar_file, qgxe_file, gxe_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, no_constrain, mlma_flag, within_family, reml_bending, reml_diag_one);
    } else if (grm_flag || m_grm_flag) {
        if (pca_flag) pter_gcta->pca(grm_file, kp_indi_file, rm_indi_file, grm_cutoff, m_grm_flag, out_pc_num);
        else if (make_grm_flag) pter_gcta->save_grm(grm_file, kp_indi_file, rm_indi_file, update_sex_file, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, grm_out_bin_flag);
        else if (align_grm_flag) pter_gcta->align_grm(grm_file);
        else if (bK_threshold > -1) pter_gcta->grm_bK(grm_file, kp_indi_file, rm_indi_file, bK_threshold, grm_out_bin_flag);
    }
    else if (ld_mean_rsq_seg_flag) pter_gcta->ld_seg(LD_file, LD_seg, LD_wind, LD_rsq_cutoff, dominance_flag);
    else LOGGER.e(0, "no analysis has been launched by the option(s).\n");

    delete pter_gcta;
}
