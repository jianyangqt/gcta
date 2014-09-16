/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Read imputed data
 *
 * 2012 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::read_imp_info_mach(string zinfofile)
{
    _dosage_flag=true;

    int i=0;
    const int MAX_LINE_LENGTH = 1000;
    char buf[MAX_LINE_LENGTH];
    gzifstream zinf;
    zinf.open( zinfofile.c_str() );
    if(! zinf.is_open()) throw("Error: can not open the file ["+zinfofile+"] to read.");

    string str_buf, errmsg="Reading dosage data failed. Please check the format of the map file.";
    char c_buf;
    double f_buf=0.0;
    cout<<"Reading map file of the imputed dosage data from ["+zinfofile+"]."<<endl;
    zinf.getline(buf, MAX_LINE_LENGTH, '\n'); // skip the header
    vector<string> vs_buf;
    int col_num=StrFunc::split_string(buf, vs_buf);
    if(vs_buf[col_num-1]!="Rsq") throw(errmsg);
    _snp_name.clear();
    _allele1.clear();
    _allele2.clear();
    _dosage_Rsq.clear();
    while(1){
        zinf.getline(buf, MAX_LINE_LENGTH, '\n');
        if(zinf.fail() || !zinf.good()) break;
        stringstream ss(buf);
        string nerr=errmsg+"\nError line: "+ss.str();
        if(!(ss>>str_buf)) throw(nerr);
        _snp_name.push_back(str_buf);
        if(!(ss>>c_buf)) throw(nerr);
        _allele1.push_back(c_buf);
        if(!(ss>>c_buf)) throw(nerr);
        _allele2.push_back(c_buf);
        for(i=0; i<col_num-3; i++) if(!(ss>>f_buf)) throw(nerr);
        _dosage_Rsq.push_back(f_buf);
        if(ss>>str_buf) throw(nerr);
    }
    zinf.clear();
    zinf.close();
    _snp_num=_snp_name.size();
    cout<<_snp_num<<" SNPs to be included from ["+zinfofile+"]."<<endl;
    _chr.resize(_snp_num);
    _bp.resize(_snp_num);
    _genet_dst.resize(_snp_num);
    _ref_A=_allele1;

	// Initialize _include
	init_include();
}

void gcta::read_imp_prob_mach(string filenames)
{
    ifstream ifilenames(filenames.c_str());
    if(!ifilenames) throw("Error: can not open the file ["+filenames+"] to read.");
    string str_buf;
    vector<string> mlinfo_files, mlprob_files, vs_buf;
    while(getline(ifilenames, str_buf)){
        if(!str_buf.empty()){
            if(StrFunc::split_string(str_buf, vs_buf)!=2) continue;
            mlinfo_files.push_back(vs_buf[0]);
            mlprob_files.push_back(vs_buf[1]);
        }
    }
    cout<<"There are "<<mlinfo_files.size()<<" pairs of *.mlinfo and *.mlprob filenames specified in ["+filenames+"]."<<endl;
    if(mlinfo_files.size()<1) throw("Error: no *.mlinfo and *.mlprob filename is found in ["+filenames+"].");
    
    int i=0;
    const int MAX_LINE_LENGTH = 1000;
    char c_buf, buf[MAX_LINE_LENGTH];
    double d_buf;
    vector<string> snp;
    vector<char> A1, A2;
    vector<double> freq, maf, qual, rsq;
    for(i=0; i<mlinfo_files.size(); i++){
        // read *.mlinfo file
        vector<string> snp_buf;
        vector<char> A1_buf, A2_buf;
        vector<double> freq_buf, maf_buf, qual_buf, rsq_buf;
        string errmsg="Error: format error in the file ["+mlinfo_files[i]+"].";
        gzifstream zinf;
        zinf.open(mlinfo_files[i].c_str());
        if(! zinf.is_open()) throw("Error: can not open the file ["+mlinfo_files[i]+"] to read.");     

        cout<<"Reading SNP summary information of the imputed data from ["+mlinfo_files[i]+"]..."<<endl;
        zinf.getline(buf, MAX_LINE_LENGTH, '\n'); // skip the header
        vector<string> vs_buf;
        int col_num=StrFunc::split_string(buf, vs_buf);
        if(vs_buf[col_num-1]!="Rsq") throw(errmsg);
        while(1){
            zinf.getline(buf, MAX_LINE_LENGTH, '\n');
            if(zinf.fail() || !zinf.good()) break;
            stringstream ss(buf);
            string nerr=errmsg+"\nError line: "+ss.str();
            if(!(ss>>str_buf)) throw(nerr);
            snp_buf.push_back(str_buf);
            if(!(ss>>c_buf)) throw(nerr);
            A1_buf.push_back(c_buf);
            if(!(ss>>c_buf)) throw(nerr);
            A2_buf.push_back(c_buf);
            if(!(ss>>d_buf)) throw(nerr);
            freq_buf.push_back(d_buf);
            if(!(ss>>d_buf)) throw(nerr);
            maf_buf.push_back(d_buf);
            if(!(ss>>d_buf)) throw(nerr);
            qual_buf.push_back(d_buf);
            if(!(ss>>d_buf)) throw(nerr);
            rsq_buf.push_back(d_buf);
            if(ss>>str_buf) throw(nerr);
        }
        zinf.clear();
        zinf.close();
        snp.insert(snp.end(), snp_buf.begin(), snp_buf.end());
        A1.insert(A1.end(), A1_buf.begin(), A1_buf.end());
        A2.insert(A2.end(), A2_buf.begin(), A2_buf.end());
        freq.insert(freq.end(), freq_buf.begin(), freq_buf.end());
        maf.insert(maf.end(), maf_buf.begin(), maf_buf.end());
        qual.insert(qual.end(), qual_buf.begin(), qual_buf.end());
        rsq.insert(rsq.end(), rsq_buf.begin(), rsq_buf.end());
        
        // read *.mlprob file
        
        
    }
    
    ///////////////////////
    int i=0;
    
    string str_buf, errmsg="Reading dosage data failed. Please check the format of the map file.";

    double f_buf=0.0;
    
	// Initialize _include
	init_include();
    
    //////////////////////
    
    
    
    
    int i=0, j=0;

    const int MAX_LINE_LENGTH = 10000000;
    char buf[MAX_LINE_LENGTH];
    gzifstream zinf;
    zinf.open( zprobfile.c_str() );
    if(! zinf.is_open()) throw("Error: can not open the file ["+zprobfile+"] to read.");

    string str_buf, id_buf, err_msg="Error: reading imputed posteriror probability data failed. Are the *.mlinfo file and the *.mlprob file matched?";
    double f_buf=0.0;
    vector<string> kept_id, vs_buf;
    vector<float> vf_buf;
    cout<<"Reading imputed posteriror probability data data from ["+zprobfile+"]..."<<endl;
    _fid.clear();
    _pid.clear();
    _geno_dose.clear();
    while(1){
        bool kp_flag=true;
        zinf.getline(buf, MAX_LINE_LENGTH, '\n');
        if(zinf.fail() || !zinf.good()) break;
        stringstream ss(buf);
        ss>>str_buf;
        StrFunc::split_string(str_buf, vs_buf, ">");
        if(vs_buf[0].empty()) throw("Error: family ID of the individual ["+str_buf+"] is missing.");
        else vs_buf[0].erase(vs_buf[0].end()-1);
        _fid.push_back(vs_buf[0]);
        _pid.push_back(vs_buf[1]);
        id_buf=vs_buf[0]+":"+vs_buf[1];
        ss>>str_buf; // skip one string
        vf_buf.clear();
        vf_buf.resize(_include.size());
        for(i=0, j=0; i<_snp_num; i++){
            if(!(ss>>f_buf)) throw(err_msg);
            if(rsnp[i]){ vf_buf[j]=(f_buf); j++; }
        }
        if(ss>>str_buf) throw("Error: reading dosage data failed. Are the map file and the dosage file matched?");
        if(kp_indi_flag && kp_id_map.find(id_buf)==kp_id_map.end()) kp_flag=false;
        if(kp_flag && blup_indi_flag && blup_id_map.find(id_buf)==blup_id_map.end()) kp_flag=false;
        if(kp_flag && rm_indi_flag && rm_id_map.find(id_buf)!=rm_id_map.end()) kp_flag=false;
        if(kp_flag){
            _geno_dose.push_back(vf_buf);
            kept_id.push_back(id_buf);
        }
    }
    zinf.clear();
    zinf.close();
    _indi_num=_fid.size();

    if(kp_indi_flag) cout<<"Because of the --keep option."<<endl;
    if(blup_indi_flag) cout<<"Because of the --blup-snp option."<<endl;
    if(rm_indi_flag) cout<<"Because of the --remove option."<<endl;
    cout<<"Imputed dosage data for "<<kept_id.size()<<" individuals are included from ["<<zdosefile<<"]."<<endl;
    _fa_id.resize(_indi_num);
    _mo_id.resize(_indi_num);
    _sex.resize(_indi_num);
    _pheno.resize(_indi_num);

 	// initialize keep
    init_keep();
    update_id_map_kp(kept_id, _id_map, _keep);
    if(_keep.size()==0) throw("Error: No individual is retained for analysis.");

    if(blup_indi_flag) read_indi_blup(blup_indi_file);

	// update data
	update_bim(rsnp);
}
