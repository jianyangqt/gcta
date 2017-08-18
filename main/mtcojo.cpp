#include "gcta.h"
#include "Logger.h"

int read_snp_metafile(string metafile) {
    int success=1;

    ifstream meta_snp(metafile.c_str());
    if (!meta_snp) { throw ("Error: can not open the file [" + mtcojolist_file + "] to read."); success=0; }
    //
    if (StrFunc::split_string(str_buf, vs_buf) < 7) throw ("Error: format error in the input file [" + metafile + "]."); 
}

void read_mtcojofile(string mtcojolist_file, string target_pheno, vector<string> covar_pheno, int ncovar) {

    LOGGER.i(0, "\nReading GWAS summary-level statistics from [" + mtcojolist_file + "] ...");
    ifstream meta_list(mtcojolist_file.c_str());
    if (!meta_list) throw ("Error: can not open the file [" + mtcojolist_file + "] to read.");

    string strbuf="", target_pheno_file="";
    vector<string> vs_buf, covar_pheno_file, snplist;
    int success=1;

    // Retrieve the GWAS summary data file
    // The 1st row: the target trait
    getline(meta_list, str_buf);
    if (StrFunc::split_string(str_buf, vs_buf) != 2) { throw ("Error: format error in the input file [" + mtcojolist_file + "]."); success=0; }
    target_pheno=vs_buf[1]; target_pheno_file=vs_buf[2];
 
    // The rest rows: the covariate traits
    ncovar=0;
    while(getline(meta_list, str_buf)) {
         if (StrFunc::split_string(str_buf, vs_buf) < 2) { throw ("Error: format error in the input file [" + mtcojolist_file + "]."); success=0;}
         covar_pheno.push_back(vs_buf[1]);
         covar_pheno_file.puch_back(vs_buf[2]);
         ncovar++;
    }

    read_snp_metafile(target_pheno_file, snplist);
    read_multi_metafile(target_pheno_file);   

    cout << "Matching the GWAS meta-analysis results to the genotype data ..." << endl;
    gcta::update_id_map_kp(snplist, _snp_name_map, _include);

}


void gcta::mtcojo(string mtcojolist_file, string string ref_ld_dirt, string w_ld_dirt,
                  double gwas_thresh, double heidi_thresh, int nsnp_heidi, int nsnp_gsmr) {
    string target_pheno="";
    vector<string> covar_pheno;
    int ncovar=0, readFlag=0;

    readFlag=read_mtcojofile(mtcojolist_file, target_pheno, covar_pheno, ncovar);
//    read_multi_metafile(metafile_list);
//   qc_metafile();
//    gsmr();
//    mtcojo_ldsc();
//    mtcojo_cond();
//    ifstream meta_list(metafile.c_str());
}

