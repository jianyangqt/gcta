/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for reading the raw genotype data
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "gcta.h"

void gcta::read_IRG_fnames(string snp_info_file, string fname_file, double GC_cutoff)
{
    // read SNP summary information
    LOGGER<<"Reading summary information of the SNPs ..."<<endl;
    ifstream i_snp_info(snp_info_file.c_str());
    if(!i_snp_info) LOGGER.e(0, "cannot open the file ["+snp_info_file+"] to read.");
    int i=0,j=0;
    string str_buf;
    vector<string> vs_buf;
    getline(i_snp_info, str_buf);
    while(getline(i_snp_info, str_buf)){
        StrFunc::split_string(str_buf, vs_buf, "\t ,;\n");
        j=0;
        _snp_name.push_back(vs_buf[++j]);
        _chr.push_back(atoi(vs_buf[++j].c_str()));
        _bp.push_back(atoi(vs_buf[++j].c_str()));
    }
    i_snp_info.close();
    _snp_num=_snp_name.size();
    LOGGER<<_snp_num<<" SNPs specified in ["+snp_info_file+"]."<<endl;

    // save PLINK map file
    string map_file=_out+".map";
    ofstream omap(map_file.c_str());
    LOGGER<<"Saving PLINK MAP file ..."<<endl;
    if(!omap) LOGGER.e(0, "cannot open the file ["+map_file+"] to write.");
    for(i=0; i<_snp_num; i++) omap<<_chr[i]<<"\t"<<_snp_name[i]<<"\t0\t"<<_bp[i]<<endl;
    omap.close();
    LOGGER<<"PLINK MAP file has been saved in ["+map_file+"].\n"<<endl;

    // read the filenames of the raw genotype data files
    ifstream i_fnames(fname_file.c_str());
    if(!i_fnames) LOGGER.e(0, "cannot open the file ["+fname_file+"] to read.");
    vector<string> fnames;
    while(getline(i_fnames, str_buf)){
        if(StrFunc::split_string(str_buf, vs_buf)==1) fnames.push_back(vs_buf[0]);
    }
    i_fnames.close();
    _indi_num=fnames.size();
    LOGGER<<_indi_num<<" raw genotype data filenames specified in ["+fname_file+"]."<<endl;

    // read raw genotype file
    _snp_1.resize(_snp_num);
    _snp_2.resize(_snp_num);
    for(i=0; i<_snp_num; i++){
        _snp_1[i].reserve(_indi_num);
        _snp_2[i].reserve(_indi_num);
    }
    LOGGER<<"Reading the raw genotype files and saving the genotype data in PLINK PED format ..."<<endl;
    LOGGER<<"(SNP genotypes with GenCall rate < "<<GC_cutoff<<" are regarded as missing)"<<endl;
    string ped_file=_out+".ped";
    ofstream oped(ped_file.c_str());
    if(!oped) LOGGER.e(0, "cannot open the file ["+ped_file+"] to read.");
    for(i=0; i<_indi_num; i++){
        read_one_IRG(oped, i, fnames[i], GC_cutoff);
        LOGGER<<i+1<<" of "<<_indi_num<<" files.\r";
    }
    LOGGER<<"Genotype data for "<<_indi_num<<" individuals have been save in the file ["+ped_file+"]."<<endl;
    oped.close();
}

char gcta::flip_allele(char a)
{
    if(a=='A') return('T');
    else if(a=='C') return('G');
    else if(a=='G') return('C');
    else if(a=='T') return('A');
    else return('0');
}

void gcta::read_one_IRG(ofstream &oped, int ind, string IRG_fname, double GC_cutoff)
{
    ifstream i_IRG(IRG_fname.c_str());
    if(!i_IRG) LOGGER.e(0, "cannot open the file ["+IRG_fname+"] to read.");
    char a1='0', a2='0';
    int i=0, j=0, snp_num=0;
    double GC=0.0;
    string str_buf, fid, pid;
    vector<string> vs_buf;
    for(i=0; i<2; i++) getline(i_IRG, str_buf);
    StrFunc::split_string(str_buf, vs_buf);
    bool oldversion=false;
    if(vs_buf[2]=="1.1.9") oldversion=true;
    for(i=0; i<3; i++) getline(i_IRG, str_buf);
    StrFunc::split_string(str_buf, vs_buf);
    snp_num=atoi(vs_buf[2].c_str());
    if(snp_num!=_snp_num) LOGGER.e(0, "the number of SNPs specified in the summary data file does not match that in the raw genotype file ["+IRG_fname+"].");
    for(i=0; i<6; i++) getline(i_IRG, str_buf);
    for(i=0; i<snp_num; i++){
        getline(i_IRG, str_buf);
        StrFunc::split_string(str_buf, vs_buf, "\t ,;\n");
        if(vs_buf[0]!=_snp_name[i]) LOGGER.e(0, "the SNP ["+vs_buf[0]+"] specified in the summary data file does not match that in the raw genotype file ["+IRG_fname+"]. Has the order of the SNPs been changed?");
        pid=vs_buf[1];
        if(oldversion) fid=pid;
        else fid=vs_buf[2];
        if(oldversion){
            a1=vs_buf[2][0];
            a2=vs_buf[3][0];
            GC=atof(vs_buf[8].c_str());
        }
        else{
            GC=atof(vs_buf[3].c_str());
            a1=vs_buf[6][0];
            a2=vs_buf[7][0];
        }
        if(GC<GC_cutoff) a1=a2='0';
        if(i==0) oped<<fid<<" "<<pid<<" -9 -9 -9 -9 ";
        oped<<a1<<" "<<a2<<" ";
    }
    oped<<endl;
    i_IRG.close();
}
