/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for mtcojo
 *
 * 2018 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */


#include "gcta.h"
#include "Logger.h"
#include "StatFunc.h"
#include "zlib.h"
#include <limits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iterator>
#include <map>
#include <set>
#include <algorithm>

double eps = 1e-6;

bool file_exists(string filestr) {
    if(FILE *input_file = fopen(filestr.c_str(), "r")) {
        fclose(input_file);
        return true;
    }
    else return false;
}

vector<vector<int>> split_extract(vector<int> extract) {
    int i = 0, nBlock = 0, nExtract = extract.size();
    vector<int> indexToExtract(2);
    vector<vector<int>> blockToExtract;
  
    indexToExtract[0] = indexToExtract[1] = extract[0];
    blockToExtract.push_back(indexToExtract);
    for( i = 1; i < nExtract; i++ ) {
        if(extract[i]==indexToExtract[1]+1) blockToExtract[nBlock][1] = extract[i];
        else {
            indexToExtract[0] = indexToExtract[1] = extract[i];
            blockToExtract.push_back(indexToExtract);
            nBlock++;
        }
    }

    return(blockToExtract);
}

void extractRow(eigenMatrix &matrix, vector<int> rowToExtract)
{
    int i = 0;
    int numRows = rowToExtract.size(), numCols = matrix.cols();
    eigenMatrix matrix_new(numRows, numCols);

    vector<vector<int>> blockToExtract = split_extract(rowToExtract);
    int nBlock = blockToExtract.size(), rowIndex = 0;
    for( i = 0; i < nBlock; i++ ) {
        int numRowsBuf = blockToExtract[i][1] - blockToExtract[i][0] + 1;
        matrix_new.block(rowIndex, 0, numRowsBuf, numCols) = matrix.block(blockToExtract[i][0], 0, numRowsBuf, numCols);
        rowIndex += numRowsBuf;
    }

    matrix = matrix_new;
}

void extractCol(eigenMatrix &matrix, vector<int> colToExtract)
{
    int i = 0;
    int numRows = matrix.rows(), numCols = colToExtract.size();
    eigenMatrix matrix_new(numRows, numCols);

    vector<vector<int>> blockToExtract = split_extract(colToExtract);
    int nBlock = blockToExtract.size(), colIndex = 0;
    for( i = 0; i < nBlock; i++ ) {
        int numColsBuf = blockToExtract[i][1] - blockToExtract[i][0] + 1;
        matrix_new.block(0, colIndex, numRows, numColsBuf) = matrix.block(0, blockToExtract[i][0], numRows, numColsBuf);
        colIndex += numColsBuf;
    }

    matrix = matrix_new;
}

void removeElement(eigenVector &vector, int indexToRemove)
{
    int numElements = vector.size()-1;

    if( indexToRemove < numElements )
        vector.segment(indexToRemove, numElements-indexToRemove) = vector.segment(indexToRemove+1, numElements-indexToRemove);

    vector.conservativeResize(numElements);
}

void removeRow(eigenMatrix &matrix, int rowToRemove)
{
    int numRows = matrix.rows()-1;
    int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(eigenMatrix &matrix, int colToRemove)
{
    int numRows = matrix.rows();
    int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

vector<string> update_common_snps(const vector<string> snplist, vector<string> common_snps) {
    
    int i = 0, nsnpbuf = snplist.size();
    vector<string> snpbuf;
    std::vector<string>::iterator iter;
    
    for(i=0; i<nsnpbuf; i++) {
        iter = find(common_snps.begin(), common_snps.end(), snplist[i]);
        if( iter != common_snps.end()) snpbuf.push_back(snplist[i]);
    }
    stable_sort(snpbuf.begin(), snpbuf.end());
    
    return snpbuf;
}

void read_metafile_list(string mtcojolist_file, string &target_pheno, string &target_pheno_file, vector<string> &covar_pheno, vector<string> &covar_pheno_file, vector<double> &popu_prev, vector<double> &smpl_prev) {

    ifstream meta_list(mtcojolist_file.c_str());
    if (!meta_list)
        LOGGER.e(0, "cannot open the file [" + mtcojolist_file + "] to read.");
    
    string strbuf="", prevbuf1="", prevbuf2="";
    double d_prev1 = 0.0, d_prev2 = 0.0;
    // Retrieve the GWAS summary data file
    // The 1st row: the target trait
    int line_number = 1;
    std::getline(meta_list, strbuf);
    std::istringstream linebuf(strbuf);
    std::istream_iterator<string> begin_title(linebuf), end_title;
    vector<string> line_elements(begin_title, end_title);
    if(line_elements.size() != 2 && line_elements.size() != 4) {
        LOGGER.e(0, "the format of file [" + mtcojolist_file + "] is incorrect, line " + to_string(line_number) + ".");
    }
    target_pheno=line_elements[0]; target_pheno_file=line_elements[1];
    // prevelance
    d_prev1 = nan(""); d_prev2 = nan("");
    if(line_elements.size()==4) {
        prevbuf1 = line_elements[2]; prevbuf2 = line_elements[3];
        StrFunc::to_upper(prevbuf1); StrFunc::to_upper(prevbuf2);
        // available, 1 - sample prevelance, 2 - population prevelance
        if(prevbuf1 != "NA"  &&  prevbuf1!= "NAN" && prevbuf1!= ".") {
            d_prev1 = atof(prevbuf1.c_str());
            if(d_prev1 < 0 || d_prev1 > 1)
                LOGGER.e(0, "invalid sample prevalence for trait [" + target_pheno + "].");
        }
        if(prevbuf2 != "NA"  &&  prevbuf2!= "NAN" && prevbuf2 != ".") {
            d_prev2 = atof(prevbuf2.c_str());
            if(d_prev2 < 0 || d_prev2 > 1)
                LOGGER.e(0, "invalid population prevalence for trait [" + target_pheno + "].");
        }
    }
    smpl_prev.push_back(d_prev1); popu_prev.push_back(d_prev2);
    
    // The rest rows: the covariate traits
    while(std::getline(meta_list, strbuf)) {
        line_number++;
        //std::istringstream linebuf(strbuf);
        linebuf.clear();
        linebuf.str(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 2 && line_elements.size() != 4) {
            LOGGER.e(0, "the format of file [" + mtcojolist_file + "] is incorrect, line " + to_string(line_number) + ".");
        }
        covar_pheno.push_back(line_elements[0]);
        covar_pheno_file.push_back(line_elements[1]);
        // prevelance
        d_prev1 = nan(""); d_prev2 = nan("");
        if(line_elements.size()==4) {
            prevbuf1 = line_elements[2]; prevbuf2 = line_elements[3];
            StrFunc::to_upper(prevbuf1); StrFunc::to_upper(prevbuf2);
            // available, 1 - sample prevelance, 2 - population prevelance
            if(prevbuf1 != "NA"  &&  prevbuf1!= "NAN" && prevbuf1!= ".") {
                d_prev1 = atof(prevbuf1.c_str());
                if(d_prev1 < 0 || d_prev1 > 1)
                    LOGGER.e(0, "invalid sample prevalence for trait [" + line_elements[0] + "].");
            }
            if(prevbuf2 != "NA"  &&  prevbuf2!= "NAN" && prevbuf2 != ".") {
                d_prev2 = atof(prevbuf2.c_str());
                if(d_prev2 < 0 || d_prev2 > 1)
                    LOGGER.e(0, "invalid population prevalence for trait [" + line_elements[0] + "].");
            }
        }
        smpl_prev.push_back(d_prev1); popu_prev.push_back(d_prev2);
    }
    meta_list.close();
}

void gcta::init_meta_snp_map(vector<string> snplist, map<string, int> &snp_name_map, vector<string> &snp_name, vector<int> &remain_snp) {
    int i=0, size=0, nsnp = snplist.size();

    snp_name_map.clear(); remain_snp.clear();
    remain_snp.resize(nsnp); 
    for(i=0; i<nsnp; i++) {
        snp_name_map.insert( pair<string,int>(snplist[i], i) );
        if(size==snp_name_map.size()) LOGGER.e(0, "Duplicated SNP ID found: " + snplist[i] + ".");
        size = snp_name_map.size();
        remain_snp[i] = i;
    }

    snp_name.clear(); snp_name = snplist;
}

void gcta::init_gwas_variable(vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, eigenMatrix &snp_freq, eigenMatrix &snp_b, eigenMatrix &snp_se, eigenMatrix &snp_pval, eigenMatrix &n, int npheno, int nsnp) {
    int i = 0;
    // two alleles
    snp_a1.clear(); snp_a2.clear();
    snp_a1.resize(npheno); snp_a2.resize(npheno);
    for( i=0; i<npheno; i++) {
        snp_a1[i].resize(nsnp); snp_a2[i].resize(nsnp);
    }
    // allele frequency
    snp_freq.resize(nsnp, npheno);
    // effect size, se, p-value and sample size
    snp_b.resize(nsnp, npheno); snp_se.resize(nsnp, npheno);
    snp_pval.resize(nsnp, npheno); n.resize(nsnp, npheno);
}

void gcta::update_meta_snp_list(vector<string> &snplist, map<string, int> snp_id_map) {
    int i = 0, nsnpbuf = 0;
    vector<string> snpbuf(snplist);
    nsnpbuf = snplist.size();

    snplist.clear();
    for(i=0; i<nsnpbuf; i++) {
        if(snp_id_map.find(snpbuf[i]) == snp_id_map.end()) continue;
        snplist.push_back(snpbuf[i]);
    }
}

void gcta::update_meta_snp_map(vector<string> snplist, map<string, int> &snp_id_map, vector<string> &snp_id, vector<int> &snp_indx, bool indx_flag) {
    int i = 0, j = 0, size = 0, nsnp = snplist.size();
    map<string,int> snp_add_map;
    
    for(i=0, j=snp_id_map.size(); i<nsnp; i++) {
        // Check duplicated SNP ID
        snp_add_map.insert(pair<string,int>(snplist[i], i));
        if(size==snp_add_map.size()) LOGGER.e(0, "Duplicated SNP ID found: " + snplist[i] + ".");
        size = snp_add_map.size();

        if(snp_id_map.find(snplist[i]) != snp_id_map.end()) continue;
        snp_id_map.insert(pair<string, int>(snplist[i], j));
        snp_id.push_back(snplist[i]);
        if(indx_flag) snp_indx.push_back(j);
        j++;
    }
}

void gcta::update_meta_snp(map<string,int> &snp_name_map, vector<string> &snp_name, vector<int> &snp_remain) {
    int i = 0, nsnp = snp_remain.size();
    vector<int> snp_remain_buf;
    vector<string> snp_name_buf;

    snp_name_buf = snp_name;
    snp_remain_buf = snp_remain;
    // reset snp_name and snp_remain
    snp_name.clear(); snp_name.resize(nsnp);
    snp_remain.clear(); snp_remain.resize(nsnp);
    for(i = 0; i<nsnp; i++) {
        snp_name[i] = snp_name_buf[snp_remain_buf[i]];
        snp_remain[i] = i;
    }
    // reset snp_name_map
    snp_name_map.clear();
    for(i=0; i<nsnp; i++) 
        snp_name_map.insert(pair<string, int>(snp_name[i], i));
}

void update_mtcojo_snp_kp(const vector<string> adjsnps, map<string,int> &snp_id_map, vector<int> &remain_snp_indx) {
    
    int i=0, nsnpbuf = adjsnps.size();
    map<string,int> snp_id_map_buf(snp_id_map);
    for(i=0; i<nsnpbuf; i++) snp_id_map_buf.erase(adjsnps[i]);
    
    std::map<string,int>::iterator iter;
    for(iter=snp_id_map_buf.begin(); iter!=snp_id_map_buf.end(); iter++)
        snp_id_map.erase(iter->first);
    
    nsnpbuf = snp_id_map.size();
    remain_snp_indx.clear(); remain_snp_indx.resize(nsnpbuf);
    for(iter=snp_id_map.begin(), i=0; iter!=snp_id_map.end(); iter++, i++) remain_snp_indx[i] = iter->second;
    
    stable_sort(remain_snp_indx.begin(), remain_snp_indx.end());
}

void gcta::update_mtcojo_snp_rm(vector<string> adjsnps, map<string,int> &snp_id_map, vector<int> &remain_snp_indx) {
    
    int i=0, nsnpbuf=adjsnps.size();
    std::map<string,int>::iterator iter;
    
    for(i=0; i<nsnpbuf; i++) snp_id_map.erase(adjsnps[i]);
    
    nsnpbuf = snp_id_map.size();
    remain_snp_indx.clear(); remain_snp_indx.resize(nsnpbuf);
    for(iter=snp_id_map.begin(), i=0; iter!=snp_id_map.end(); iter++, i++) remain_snp_indx[i] = iter->second;

    stable_sort(remain_snp_indx.begin(), remain_snp_indx.end());
}

vector<string> gcta::read_snp_metafile_txt(string metafile, map<string,int> &gws_snp_name_map, double thresh) {
    ifstream meta_snp(metafile.c_str());
    if (!meta_snp)
         LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    
    string strbuf="", valbuf = "";
    vector<string> snplist;
    int line_number=0;

    // Read the summary data
    while(std::getline(meta_snp, strbuf)) {
        line_number++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8) {
            LOGGER.e(0, "the GWAS summary data file [" + metafile + "] should be in GCTA-COJO format, line " + to_string(line_number) + ".");
        }
        if(line_number==1) continue;
        // keep significant SNPs
        snplist.push_back(line_elements[0]);
        valbuf = line_elements[6];
        double pval_buf = 1.0;
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") pval_buf = atof(valbuf.c_str());
        if( pval_buf < thresh && (gws_snp_name_map.find(line_elements[0]) == gws_snp_name_map.end()) ) {
            int size = gws_snp_name_map.size();
            gws_snp_name_map.insert(pair<string,int>(line_elements[0], size));
        }
    }
    meta_snp.close();
    return snplist;
}

vector<string> gcta::read_snp_metafile_gz(string metafile, map<string,int> &gws_snp_name_map, double thresh) {
    string strbuf="", snpbuf="";
    vector<string> snplist;
    int line_number=0;
    const int MAX_LINE_LENGTH = 1024;
    char buf[MAX_LINE_LENGTH];
    
    gzifstream meta_snp(metafile.c_str());
    if (!meta_snp)
         LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    
    string err_msg = "Failed to read [" + metafile + "]. An error occurs in line ";

    // Read the summary data
    while(1) {
        meta_snp.getline(buf, MAX_LINE_LENGTH, '\n');
        if(meta_snp.fail() || !meta_snp.good()) break;
        line_number ++;
        if(line_number==1) continue;
        stringstream ss(buf);
        if (!(ss >> snpbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        // Save the SNPs
        snplist.push_back(snpbuf);
        int i = 0;
        for( i=0; i<6; i++)
            if (!(ss >> strbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        // Save the p-value
        double pval_buf = 1.0;
        if(strbuf!="." && strbuf!="NA" && strbuf!="NAN") pval_buf = atof(strbuf.c_str());
        if( pval_buf < thresh && (gws_snp_name_map.find(snpbuf) == gws_snp_name_map.end()) ) {
            int size = gws_snp_name_map.size();
            gws_snp_name_map.insert(pair<string,int>(snpbuf, size));
        }
    }

    meta_snp.close();
    return snplist;
}

double gcta::read_single_metafile_txt(string metafile, map<string, int> id_map,
                         vector<string> &snp_a1, vector<string> &snp_a2,
                         eigenVector &snp_freq, eigenVector &snp_b,
                         eigenVector &snp_se, eigenVector &snp_pval,
                         eigenVector &snp_n, vector<bool> &snp_flag) {
   
    ifstream meta_raw(metafile.c_str());
    if (!meta_raw)
        LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    string strbuf="", valbuf="";
    int line_number=0, snp_indx=0;
    map<string, int>::iterator iter;
    // Read the summary data
    double pval_thresh = 0.5, h_buf = 0.0, vp_buf = 0.0, median_vp = 0.0;
    bool missing_flag = false;
    vector<double> vec_vp_buf;
     
    while(std::getline(meta_raw, strbuf)) {
        missing_flag = false;
        line_number++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8) {
            LOGGER.e(0, "the GWAS summary data file [" + metafile + "] should be in GCTA-COJO format, line " + to_string(line_number) + ".");
        }        
        // Read the summary data
        if(line_number==1) continue;
        iter = id_map.find(line_elements[0]);      
        if(iter == id_map.end()) continue;     

        snp_indx = iter->second;         
        snp_a1[snp_indx] = line_elements[1]; snp_a2[snp_indx] = line_elements[2];
        StrFunc::to_upper(snp_a1[snp_indx]); StrFunc::to_upper(snp_a2[snp_indx]);
        valbuf =line_elements[3];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_freq(snp_indx) = atof(valbuf.c_str());
        else { snp_freq(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[4];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_b(snp_indx) = atof(valbuf.c_str());
        else { snp_b(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[5];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_se(snp_indx) = atof(valbuf.c_str());
        else { snp_se(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[6];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_pval(snp_indx) = atof(valbuf.c_str());
        else { snp_pval(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[7];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_n(snp_indx) = atof(valbuf.c_str());
        else { snp_n(snp_indx) = nan(""); missing_flag=true; }
        snp_flag[snp_indx] = true;

        if(!missing_flag) {
            h_buf = 2 * snp_freq(snp_indx) * ( 1- snp_freq(snp_indx) );
            vp_buf = h_buf * snp_b(snp_indx) * snp_b(snp_indx) + h_buf * snp_n(snp_indx) * snp_se(snp_indx) * snp_se(snp_indx);
            vec_vp_buf.push_back(vp_buf);
        }        
    }
    meta_raw.close();
    if(vec_vp_buf.size()>0) median_vp = CommFunc::median(vec_vp_buf);

    return median_vp;
}

double gcta::read_single_metafile_gz(string metafile, map<string, int> id_map,
                         vector<string> &snp_a1, vector<string> &snp_a2,
                         eigenVector &snp_freq, eigenVector &snp_b,
                         eigenVector &snp_se, eigenVector &snp_pval,
                         eigenVector &snp_n, vector<bool> &snp_flag) {
   
    const int MAX_LINE_LENGTH = 1024;
    char buf[MAX_LINE_LENGTH];

    gzifstream meta_raw(metafile.c_str());
    if (!meta_raw)
        LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    string err_msg = "Failed to read [" + metafile + "]. An error occurs in line ";

    int line_number=0, snp_indx=0;
    map<string, int>::iterator iter;
    // Read the summary data
    double pval_thresh = 0.5, h_buf = 0.0, vp_buf = 0.0, median_vp = 0.0;
    bool missing_flag = false;
    string strbuf = "";
    vector<double> vec_vp_buf;
    while(1) {
        meta_raw.getline(buf, MAX_LINE_LENGTH, '\n');
        if(meta_raw.fail() || !meta_raw.good()) break;
        missing_flag = false;
        line_number++;
        if(line_number==1) continue;
         
        string snpbuf = "", strbuf = "", valbuf = "";
        stringstream ss(buf);
        // SNP
        if (!(ss >> snpbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        iter = id_map.find(snpbuf);
        if(iter == id_map.end()) continue;     
        snp_indx = iter->second;
        // a1 
        if (!(ss >> strbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        StrFunc::to_upper(strbuf);
        snp_a1[snp_indx] = strbuf;
        // a2
        if (!(ss >> strbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        StrFunc::to_upper(strbuf);
        snp_a2[snp_indx] = strbuf;
        // freq
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_freq(snp_indx) = atof(valbuf.c_str());
        else { snp_freq(snp_indx) = nan(""); missing_flag=true; }
        // b
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_b(snp_indx) = atof(valbuf.c_str());
        else { snp_b(snp_indx) = nan(""); missing_flag=true; }
        // se
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_se(snp_indx) = atof(valbuf.c_str());
        else { snp_se(snp_indx) = nan(""); missing_flag=true; }
        // p-value
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_pval(snp_indx) = atof(valbuf.c_str());
        else { snp_pval(snp_indx) = nan(""); missing_flag=true; }
        // n
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_n(snp_indx) = atof(valbuf.c_str());
        else { snp_n(snp_indx) = nan(""); missing_flag=true; }
        snp_flag[snp_indx] = true;

        if(!missing_flag) {
            h_buf = 2 * snp_freq(snp_indx) * ( 1- snp_freq(snp_indx) );
            vp_buf = h_buf * snp_b(snp_indx) * snp_b(snp_indx) + h_buf * snp_n(snp_indx) * snp_se(snp_indx) * snp_se(snp_indx);
            vec_vp_buf.push_back(vp_buf);
        }        
    }
    meta_raw.close();
    if(vec_vp_buf.size()>0) median_vp = CommFunc::median(vec_vp_buf);

    return median_vp;
}

vector<string> gcta::remove_bad_snps(vector<string> snp_name, vector<int> snp_remain, vector<vector<bool>> snp_flag, vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, eigenMatrix &snp_freq,  
                                    eigenMatrix &snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, map<string,int> plink_snp_name_map, vector<string> snp_ref_a1, vector<string> snp_ref_a2, 
                                    int ntarget, int ncovar, string outfile_name) {
    int i=0, j=0, nsnp = snp_remain.size(), npheno = ntarget+ncovar;
    double snp_freq_bak = 0.0;
    string snp_a1_bak = "", snp_a2_bak = ""; 
    vector<string> badsnps;
    vector<int> bad_indx;
    map<string,int>::iterator iter;
    vector<vector<bool>> flip_flag(npheno);

    for(i=0; i<npheno; i++) {
        flip_flag[i].resize(nsnp);
        for(j=0; j<nsnp; j++) flip_flag[i][j] = false;
    }

    for( i=0; i<nsnp; i++) {
        bool iterFlag=true;
        int snpindx = -9;
        vector<string> allelebuf;
        // Choose a SNP
        iter = plink_snp_name_map.find(snp_name[snp_remain[i]]);
        if(iter!=plink_snp_name_map.end()) snpindx = iter->second;
        // Alleles
        if(snpindx > 0) {
            allelebuf.push_back(snp_ref_a1[snpindx]); allelebuf.push_back(snp_ref_a2[snpindx]);
        } else {
            for(j=0; j<npheno; j++) {
                if(snp_flag[j][snp_remain[i]]) {
                    allelebuf.push_back(snp_a1[j][snp_remain[i]]); 
                    allelebuf.push_back(snp_a2[j][snp_remain[i]]);
                    break;
                }
            }
        }

        snp_a1_bak = allelebuf[0]; snp_a2_bak = allelebuf[1];
        for(j=0; j<npheno; j++) {
            if(!snp_flag[j][snp_remain[i]]) continue;
            // Alleles
            allelebuf.push_back(snp_a1[j][snp_remain[i]]);
            allelebuf.push_back(snp_a2[j][snp_remain[i]]);
            // Removing SNPs with missing value
            if( std::isnan(snp_b(snp_remain[i],j)) || std::isnan(snp_se(snp_remain[i],j)) || std::isnan(snp_pval(snp_remain[i],j)) || std::isnan(snp_n(snp_remain[i],j)) ) {
                iterFlag=false; break;
            }
            // Removing SNPs with extremely small SE
            if(snp_se(snp_remain[i],j) < eps) {
                iterFlag=false; break;
            }
            // Alignment of b
            if(iterFlag && (allelebuf[0] != snp_a1[j][snp_remain[i]])) {
                snp_b(snp_remain[i],j) = -1*snp_b(snp_remain[i],j);
                if(!std::isnan(snp_freq(snp_remain[i],j))) 
                    snp_freq(snp_remain[i],j)  =  1 - snp_freq(snp_remain[i],j);
                std::swap(snp_a1[j][snp_remain[i]], snp_a2[j][snp_remain[i]]);
            }
            snp_freq_bak = snp_freq(snp_remain[i],j);
        }

        // Remove SNPs with multiple alleles
        stable_sort(allelebuf.begin(), allelebuf.end());      
        allelebuf.erase(unique(allelebuf.begin(), allelebuf.end()), allelebuf.end());
        if(allelebuf.size()!=2) iterFlag=false;
        // Collect bad SNPs
        if(iterFlag==false) {
            badsnps.push_back(snp_name[snp_remain[i]]);
        }
        
        // update the estimates for the first phenotype
        if(!snp_flag[0][snp_remain[i]]) {
            snp_a1[0][snp_remain[i]] = snp_a1_bak;
            snp_a2[0][snp_remain[i]] = snp_a2_bak;
            snp_freq(snp_remain[i],0) = snp_freq_bak;
        }
    }

    if (!badsnps.empty()) {
        string badsnpfile = outfile_name + ".badsnps", strbuf="";
        ofstream obadsnp(badsnpfile.c_str());
        if(!obadsnp) LOGGER.e(0, "cannot open file [" + badsnpfile + "] to write bad SNPs.");
        int nbadsnps = badsnps.size();
        for (i = 0; i < nbadsnps; i++) 
            obadsnp << badsnps[i] << endl;
        obadsnp.close();
        LOGGER.i(0,  to_string(nbadsnps) + " SNPs have missing value or mismatched alleles. These SNPs have been saved in [" + badsnpfile + "].");
    }
   
    stable_sort(badsnps.begin(), badsnps.end());
    return badsnps;
}

vector<string> gcta::remove_freq_diff_snps(vector<string> meta_snp_name, vector<int> meta_snp_remain, map<string,int> snp_name_map, vector<double> ref_freq, eigenMatrix meta_freq, vector<vector<bool>> snp_flag, int ntrait, double freq_thresh, string outfile_name) {
    int i = 0, nsnp = meta_snp_remain.size(), nsnp_ttl = meta_snp_remain.size();
    string snpbuf="";
    vector<string> afsnps;
    map<string,int>::iterator iter_ref;

    for( i=0; i<nsnp; i++ ) {
        int refsnp_index = 0;
        vector<double> snpfreq;

        // AF from the reference sample
        snpbuf = meta_snp_name[meta_snp_remain[i]];
        iter_ref = snp_name_map.find(snpbuf);
        if( iter_ref != snp_name_map.end()) {
            refsnp_index = iter_ref -> second;
            snpfreq.push_back(ref_freq[refsnp_index]/2);
        }

        // AF from the GWAS summary
        int j = 0;
        for( j=0; j<ntrait; j++ ) {
            if(!snp_flag[j][meta_snp_remain[i]]) continue;
            if(std::isnan(meta_freq(meta_snp_remain[i],j))) continue;
            snpfreq.push_back(meta_freq(meta_snp_remain[i],j));
        }

        // Test the difference
        int t1 = 0, t2 = 0, nsnp_test = snpfreq.size();
        double max_freq_dev = 0.0, freq_dev = 0.0;
        for( t1=0; t1<nsnp_test-1; t1++) {
            for(t2=t1+1; t2<nsnp_test; t2++) {
                freq_dev = fabs(snpfreq[t1]-snpfreq[t2]);
                if(max_freq_dev < freq_dev) max_freq_dev = freq_dev;
            }
        }
        if(max_freq_dev > freq_thresh) afsnps.push_back(snpbuf);     
    }  

    if (!afsnps.empty()) {
        string afsnpfile = outfile_name + ".freq.badsnps", strbuf="";
        ofstream oafsnp(afsnpfile.c_str());
        if(!oafsnp) LOGGER.e(0, "cannot open file [" + afsnpfile + "] to write bad SNPs.");
        int nafsnps = afsnps.size();
        for (i = 0; i < nafsnps; i++) oafsnp << afsnps[i] << endl;
        oafsnp.close();
        LOGGER.i(0,  to_string(nafsnps) + " SNP(s) have large difference of allele frequency between the GWAS summary data and the reference sample. These SNPs have been saved in [" + afsnpfile + "].");
        if(nafsnps > nsnp_ttl*0.05) 
            LOGGER.e(0, "there are too many SNPs that have large difference in allele frequency. Please check the GWAS summary data.");
    }
    return(afsnps);
}

vector<string> gcta::remove_mono_snps(map<string,int> snp_name_map, vector<double> ref_snpfreq, string outfile_name) {
    int i = 0, n_raresnp = 0;
    vector<string> afsnps;
    map<string,int>::iterator iter;

    for(iter=snp_name_map.begin(); iter!=snp_name_map.end(); iter++) {
        double af = ref_snpfreq[iter->second]/2;
        if(af > 0.5) af = 1.0 - af;
        if(CommFunc::FloatEqual(af, 0.0)) afsnps.push_back(iter->first);
        if(af < 0.01) n_raresnp++;
    }

    if(n_raresnp > 0) LOGGER.w(0,  "There are " + to_string(n_raresnp) + " SNPs with MAF < 0.01 in the reference sample.");

    if (!afsnps.empty()) {
        string afsnpfile = outfile_name + ".mono.badsnps", strbuf="";
        ofstream oafsnp(afsnpfile.c_str());
        if(!oafsnp) LOGGER.e(0, "cannot open file [" + afsnpfile + "] to write bad SNPs.");
        int nafsnps = afsnps.size();
        for (i = 0; i < nafsnps; i++) oafsnp << afsnps[i] << endl;
        oafsnp.close();
        LOGGER.i(0,  to_string(nafsnps) + " monomorphic SNP(s) have been saved in [" + afsnpfile + "].");
    }

    return(afsnps);
}

vector<string> gcta::filter_meta_snp_pval(vector<string> snp_name, vector<int> remain_snp_indx,  eigenMatrix snp_pval, int start_indx, int end_indx, vector<vector<bool>> snp_flag, double pval_thresh) {
    
    int i=0, j=0, nsnpbuf = remain_snp_indx.size();
    vector<string> kpsnps;
   
    for(i=0; i<nsnpbuf; i++) {
        for(j=start_indx; j<end_indx; j++) {
            if((snp_pval(remain_snp_indx[i],j) <  pval_thresh) && (snp_flag[j][remain_snp_indx[i]])) {
                kpsnps.push_back(snp_name[remain_snp_indx[i]]);
                break;
            }
        }
    }
    return kpsnps;
}

int gcta::read_mtcojofile(string mtcojolist_file, double gwas_thresh, int nsnp_gsmr) {

    int ncovar=0, i=0, j=0;
    string target_pheno_file="";
    vector<string> covar_pheno_file, snplist;
    
    std::map<string,int>::iterator iter;

    // Read the summary data
    LOGGER.i(0, "\nReading GWAS summary data from [" + mtcojolist_file + "] ...");
    // Read the file list for mtCOJO analysis
    _covar_pheno_name.clear(); _meta_smpl_prev.clear(); _meta_popu_prev.clear(); 
    read_metafile_list(mtcojolist_file, _target_pheno_name, target_pheno_file, _covar_pheno_name, covar_pheno_file, _meta_popu_prev, _meta_smpl_prev);

    // Read the SNPs
    // Covariates
    map<string,int> gws_snp_name_map;
    ncovar = _covar_pheno_name.size();
    for( i=0; i<ncovar; i++) {
        if(covar_pheno_file[i].substr(covar_pheno_file[i].length()-3,3)!=".gz")
            snplist=read_snp_metafile_txt(covar_pheno_file[i], gws_snp_name_map, -9);
        else 
            snplist=read_snp_metafile_gz(covar_pheno_file[i], gws_snp_name_map, -9);
        if( i == 0 ) init_meta_snp_map(snplist, _meta_snp_name_map, _meta_snp_name, _meta_remain_snp);
        else update_meta_snp_map(snplist, _meta_snp_name_map, _meta_snp_name, _meta_remain_snp, true);
    }

    // Target trait
    if(target_pheno_file.substr(target_pheno_file.length()-3,3)!=".gz")
        snplist=read_snp_metafile_txt(target_pheno_file, gws_snp_name_map, -9);
    else 
        snplist=read_snp_metafile_gz(target_pheno_file, gws_snp_name_map, -9);
    update_id_map_kp(snplist, _meta_snp_name_map, _meta_remain_snp);

    // Initialization of variables
    int nsnp = _meta_snp_name_map.size();
    eigenMatrix snp_freq;
    vector<vector<string>> snp_a1, snp_a2;

    init_gwas_variable(snp_a1, snp_a2, snp_freq, _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, ncovar+1, nsnp); 

    // reset SNP variables
    update_meta_snp(_meta_snp_name_map, _meta_snp_name, _meta_remain_snp);
//    _meta_snp_name.clear(); _meta_remain_snp.clear();
//    _meta_snp_name.resize(nsnp); _meta_remain_snp.resize(nsnp);
//    for( i=0, iter=_meta_snp_name_map.begin(); iter!=_meta_snp_name_map.end(); i++, iter++) {
//        _meta_snp_name[i] = iter->first;
//        iter->second = i;
//        _meta_remain_snp[i] = i;
//    }
    LOGGER.i(0, to_string(nsnp) + " SNPs in common between the target trait and the covariate trait(s).");
    // Read the meta analysis
    eigenVector snp_freq_buf(nsnp), snp_b_buf(nsnp), snp_se_buf(nsnp), snp_pval_buf(nsnp), snp_n_buf(nsnp);

    _meta_vp_trait.resize(ncovar+1);
    _snp_val_flag.clear(); _snp_val_flag.resize(ncovar+1);
    for(i=0; i<ncovar+1; i++) {
        _snp_val_flag[i].resize(nsnp);
        for(j=0; j<nsnp; j++) _snp_val_flag[i][j] = false;
    }

    // Target trait
    if(target_pheno_file.substr(target_pheno_file.length()-3,3)!=".gz")
        _meta_vp_trait(0) = read_single_metafile_txt(target_pheno_file, _meta_snp_name_map,  snp_a1[0], snp_a2[0], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[0]);
    else _meta_vp_trait(0) = read_single_metafile_gz(target_pheno_file, _meta_snp_name_map,  snp_a1[0], snp_a2[0], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[0]);
    if(_meta_vp_trait(0) < 0) LOGGER.e(0, "negative phenotypic variance of the target trait, " + _covar_pheno_name[0] + ".");
    
    snp_freq.col(0) = snp_freq_buf;
    _meta_snp_b.col(0) = snp_b_buf;
    _meta_snp_se.col(0) = snp_se_buf;
    _meta_snp_pval.col(0) = snp_pval_buf;
    _meta_snp_n_o.col(0) = snp_n_buf;
   
    // Covariates
    for(i=0; i<ncovar; i++) {
        if(covar_pheno_file[i].substr(covar_pheno_file[i].length()-3,3)!=".gz")
            _meta_vp_trait(i+1) = read_single_metafile_txt(covar_pheno_file[i], _meta_snp_name_map,  snp_a1[i+1], snp_a2[i+1], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i+1]);
        else _meta_vp_trait(i+1) = read_single_metafile_gz(covar_pheno_file[i], _meta_snp_name_map,  snp_a1[i+1], snp_a2[i+1], snp_freq_buf, snp_b_buf, snp_se_buf, snp_pval_buf, snp_n_buf, _snp_val_flag[i+1]);
        if(_meta_vp_trait(i+1) < 0) LOGGER.e(0, "negative phenotypic variance of the covariate #" + to_string(i+1) + ", " + _covar_pheno_name[i+1] + ".");
        snp_freq.col(i+1) = snp_freq_buf;
        _meta_snp_b.col(i+1) = snp_b_buf;
        _meta_snp_se.col(i+1) = snp_se_buf;
        _meta_snp_pval.col(i+1) = snp_pval_buf;
        _meta_snp_n_o.col(i+1) = snp_n_buf;
    }

    // QC of SNPs
    LOGGER.i(0, "Filtering out SNPs with multiple alleles or missing value ...");
   
    vector<string> badsnps;
    badsnps = remove_bad_snps(_meta_snp_name, _meta_remain_snp, _snp_val_flag, snp_a1, snp_a2, snp_freq,  _meta_snp_b, _meta_snp_se, _meta_snp_pval, _meta_snp_n_o, 
                             _snp_name_map, _allele1, _allele2, 1, ncovar, _out);
                                
    if(badsnps.size()>0) {
        update_id_map_rm(badsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(badsnps, _meta_snp_name_map, _meta_remain_snp);
    }  
  
    // For output
    _meta_snp_a1 = snp_a1[0]; _meta_snp_a2 = snp_a2[0];
    _meta_snp_freq = snp_freq;  
    nsnp = _meta_remain_snp.size();
    if(nsnp<1) LOGGER.e(0, "no SNP is retained after filtering.");
    else LOGGER.i(0, to_string(nsnp) + " SNPs are retained after filtering.");
  
    // Only keep SNPs with p-value < threshold
    double pval_thresh =  gwas_thresh;
    vector<string> keptsnps; 
    
    keptsnps = filter_meta_snp_pval(_meta_snp_name, _meta_remain_snp, _meta_snp_pval, 1, 1+ncovar, _snp_val_flag, pval_thresh);
    update_id_map_kp(keptsnps, _snp_name_map, _include);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(1) << gwas_thresh;
    LOGGER.i(0, "There are " + to_string(_include.size()) + " genome-wide significant SNPs with p < " + ss.str() + ".\n");
    if(keptsnps.size() < nsnp_gsmr) LOGGER.w(0, "We will use LD score regression to estimate bxy, because there are not enough significant SNPs for GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required.\n");

    return(_include.size());
}

eigenVector read_external_bxy(string filestr, vector<string> covar_pheno_name) {
    // Read the external bxy
    int i = 0, ncovar = covar_pheno_name.size();
    string strbuf="";
    eigenVector bxy_est(ncovar);
    map<string,int> covar_pheno_map;
    bxy_est.setConstant(-999999.0);

    ifstream extern_bxy(filestr.c_str());
    if (!extern_bxy) LOGGER.e(0, "cannot open the file [" + filestr + "] to read.");

    for(i=0; i<ncovar; i++) covar_pheno_map.insert(pair<string,int>(covar_pheno_name[i], i));
    
    while(std::getline(extern_bxy, strbuf)) {
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 2) LOGGER.e(0, "the format of file [" + filestr + "] is not correct.");
        string phenobuf = line_elements[0].c_str();
        map<string,int>::iterator iter = covar_pheno_map.find(phenobuf);
        if(iter==covar_pheno_map.end()) continue;
        bxy_est(iter->second) = atof(line_elements[1].c_str());
    }

    return(bxy_est);
}

vector<string> gcta::clumping_meta(eigenVector snp_chival, vector<bool> snp_flag, double pval_thresh, int wind_size, double r2_thresh) {   
    wind_size = wind_size*1e3;
    // Only select SNPs that are matched with plink binary file
    vector<pair<double, int>> snp_pvalbuf;
    int i = 0, indx = 0, nsnp_plink = _include.size(), nsnp_meta = _meta_remain_snp.size();
    double pvalbuf = 0.0;
    string snpbuf = "", snpbuf_left="", snpbuf_right="";
    map<string,int>::iterator iter;

    // Sort the p-value
    for(i=0; i<nsnp_meta; i++) {
        if(!snp_flag[i]) continue;
        iter = _snp_name_map.find(_meta_snp_name[_meta_remain_snp[i]]);
        if(iter==_snp_name_map.end()) continue;
        snp_pvalbuf.push_back(make_pair(snp_chival(_meta_remain_snp[i]), _meta_remain_snp[i]));
    }

    std::stable_sort(snp_pvalbuf.begin(), snp_pvalbuf.end());
    std::reverse(snp_pvalbuf.begin(), snp_pvalbuf.end());

    int nsnp_clumped=snp_pvalbuf.size();
    pval_thresh = StatFunc::qchisq(pval_thresh, 1);
    // Start to clump
    map<string, bool> clumped_snp;
    for(i=0; i<nsnp_clumped; i++) {
        pvalbuf = snp_pvalbuf[i].first;
        if( pvalbuf <= pval_thresh) continue;
        indx = snp_pvalbuf[i].second;
        clumped_snp.insert(pair<string,bool>(_meta_snp_name[indx],false));
    }

    map<string, bool>::iterator iter_clump;
    int geno_indx=0, geno_indx_j = 0, geno_indx_buf = 0, geno_indx_center = 0, nindi=_keep.size();
    int left_indx = 0, right_indx = 0;
    double r2_left=0.0, r2_right=0.0;
    vector<string> indices_snp;
    eigenVector x(nindi), x_j(nindi);
    nsnp_clumped=clumped_snp.size();
    for(i=0; i<nsnp_clumped; i++) {
        indx = snp_pvalbuf[i].second;
        snpbuf = _meta_snp_name[indx];
        if( clumped_snp[snpbuf]) continue;
        iter = _snp_name_map.find(snpbuf);
        geno_indx = iter -> second;
        geno_indx_center = std::find(_include.begin(), _include.end(), geno_indx) - _include.begin();
        makex_eigenVector(geno_indx_center, x, false, true);
        geno_indx_j = geno_indx_center;
        // Left side
        while(1) {
            r2_left=-1;  geno_indx_j--;
            if(geno_indx_j<0) break;
            snpbuf_left = _snp_name[_include[geno_indx_j]];
            if(_chr[geno_indx]==_chr[_include[geno_indx_j]] && abs(_bp[geno_indx] - _bp[_include[geno_indx_j]]) < wind_size) {
                iter_clump = clumped_snp.find(snpbuf_left);
                if(iter_clump==clumped_snp.end()) continue;
                // r2
                makex_eigenVector(geno_indx_j, x_j, false, true);
                r2_left= x.dot(x_j) / sqrt(x.dot(x) * x_j.dot(x_j));
                r2_left = r2_left*r2_left;
                // save the SNP
                if(r2_left >= r2_thresh) iter_clump->second=true;
            } else{
                break;
            }
        }
        // Right side
        geno_indx_j = geno_indx_center;
        while (1) {
            r2_right=-1; geno_indx_j++;
            if(geno_indx_j >= nsnp_plink) break;
            snpbuf_right = _snp_name[_include[geno_indx_j]];
            if(_chr[geno_indx] == _chr[_include[geno_indx_j]] && abs(_bp[geno_indx] - _bp[_include[geno_indx_j]]) < wind_size ) {
                 iter_clump = clumped_snp.find(snpbuf_right);
                if(iter_clump==clumped_snp.end()) continue;
                // r2
                makex_eigenVector(geno_indx_j, x_j, false, true);
                r2_right= x.dot(x_j) / sqrt(x.dot(x) * x_j.dot(x_j));
                r2_right = r2_right*r2_right;
                // Save the SNP
                if(r2_right >= r2_thresh) iter_clump->second=true;
            } else{
                break;
            }
        }
        indices_snp.push_back(snpbuf);
    }
  
    return indices_snp;
}

vector<int> rm_cor_elements(eigenMatrix r_mat, double r2_thresh, bool r2_flag) {
    int i = 0, j = 0, n = r_mat.cols();
    vector<int> kept_ID;

    if(r2_flag) r2_thresh = sqrt(r2_thresh);
    // identify the positions where you see a value > than the threshold
    vector<int> rm_mat_ID1, rm_mat_ID2;
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            if (abs(r_mat(i, j)) > r2_thresh) {
                rm_mat_ID1.push_back(i);
                rm_mat_ID2.push_back(j);
            }
        }
    }

    if(rm_mat_ID1.size() > 0) {
        int i_buf = 0;
        // count the number of appearance of each "position" in the vector, which involves a few steps
        vector<int> rm_uni_ID(rm_mat_ID1);
        rm_uni_ID.insert(rm_uni_ID.end(), rm_mat_ID2.begin(), rm_mat_ID2.end());
        stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
        rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
        map<int, int> rm_uni_ID_count;
        for (i = 0; i < rm_uni_ID.size(); i++) {
            i_buf = count(rm_mat_ID1.begin(), rm_mat_ID1.end(), rm_uni_ID[i]) + count(rm_mat_ID2.begin(), rm_mat_ID2.end(), rm_uni_ID[i]);
            rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
        }

        // swapping
        map<int, int>::iterator iter1, iter2;
        int n_buf = rm_mat_ID1.size();
        for (i = 0; i < n_buf; i++) {
            iter1 = rm_uni_ID_count.find(rm_mat_ID1[i]);
            iter2 = rm_uni_ID_count.find(rm_mat_ID2[i]);
            if (iter1->second < iter2->second) {
                i_buf = rm_mat_ID1[i];
                rm_mat_ID1[i] = rm_mat_ID2[i];
                rm_mat_ID2[i] = i_buf;
            }
        }
    
        stable_sort(rm_mat_ID1.begin(), rm_mat_ID1.end());
        rm_mat_ID1.erase(unique(rm_mat_ID1.begin(), rm_mat_ID1.end()), rm_mat_ID1.end());
        
        int k = 0;
        kept_ID.resize(n-rm_mat_ID1.size());
        for(i=0, j=0, k=0; i<n; i++) {
            if(i==rm_mat_ID1[j]) { j++; continue; }
            kept_ID[k] = i; k++;
        }
    } else {
        kept_ID.resize(n);
        for(i=0; i<n; i++) kept_ID[i] = i;
    }

    return kept_ID;
}

void adjust_ld_r_fdr(eigenMatrix &ld_r_mat, vector<int> kept_ID, vector<pair<double, int>> ld_pval, int m, double thresh) {
    int i = 0, nproc = ld_pval.size(), row_indx=0, col_indx=0;
    vector<double> pval_buf(nproc);

    for(i=0; i<nproc; i++) pval_buf[i] = ld_pval[i].first;

    pval_buf = StatFunc::ControlFDR_BH(pval_buf);

    for(i=0; i<nproc; i++) {
        if(pval_buf[i] < thresh) break;
        row_indx = ld_pval[i].second/m; 
        col_indx = ld_pval[i].second%m;
        ld_r_mat(kept_ID[row_indx], kept_ID[col_indx]) = ld_r_mat(kept_ID[col_indx], kept_ID[row_indx]) = 0;
    }
}

vector<int> init_interval_bxy(eigenVector bxy_sort, eigenVector bxy, vector<int> kept_ID, double lower_prob, double upper_prob) {
    int i = 0, nsnp = kept_ID.size();
    vector<int> ci_index;

    double lower_bounder = CommFunc::quantile(bxy_sort, lower_prob);
    double upper_bounder = CommFunc::quantile(bxy_sort, upper_prob);
    for(i=0; i<nsnp; i++) {
        if(bxy(kept_ID[i]) >= lower_bounder && bxy(kept_ID[i]) <= upper_bounder)
            ci_index.push_back(kept_ID[i]);
    }

    return(ci_index);
}

int topsnp_bxy(eigenVector bxy, eigenVector bzx, eigenVector bzx_se, map<string, int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID, int n_indices_snp) {
    // find the top SNP
    int i = 0, topindex = -9;
    double lower_bounder=0.0, upper_bounder = 0.0, max_bzx_chival=0.0;
    eigenVector bxy_sort(n_indices_snp);
    map<string, int>::iterator iter;
    for(i=0; i<n_indices_snp; i++) bxy_sort(i) = bxy(kept_ID[i]);
    std::stable_sort(bxy_sort.data(), bxy_sort.data()+bxy_sort.size());
    lower_bounder = CommFunc::quantile(bxy_sort, 0.4);
    upper_bounder = CommFunc::quantile(bxy_sort, 0.6);
    for(i=0; i<n_indices_snp; i++) {
        if(bxy(kept_ID[i]) >= lower_bounder && bxy(kept_ID[i]) <= upper_bounder) {
            iter = meta_snp_name_map.find(indices_snp[kept_ID[i]]);
            int i_buf = iter->second;
            double chivalbuf = bzx(i_buf)*bzx(i_buf) / (bzx_se(i_buf)*bzx_se(i_buf));
            if( chivalbuf > max_bzx_chival) {
                max_bzx_chival = chivalbuf;
                topindex = kept_ID[i];
            }
        }
    }
    return(topindex);
}

void est_cov_bxy(eigenMatrix &cov_bxy_p1, eigenMatrix &cov_bxy_p2, eigenVector bzx, eigenVector bzx_se, eigenVector bzy_se, eigenMatrix ld_r_mat,
            map<string, int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID) {

    int n_indices_snp = kept_ID.size();
    eigenVector zscore_inv1(n_indices_snp), zscore_inv2(n_indices_snp);

    #pragma omp parallel for
    for(int i=0; i<n_indices_snp; i++) {
        map<string, int>::iterator iter;
        iter = meta_snp_name_map.find(indices_snp[kept_ID[i]]);
        zscore_inv1(i) = bzx_se(iter->second)/bzx(iter->second);
        zscore_inv2(i) = bzy_se(iter->second)/bzx(iter->second);
    }

    #pragma omp parallel for
    for(int i=0; i<n_indices_snp; i++) {
        for(int j=i; j<n_indices_snp; j++) {
            cov_bxy_p1(i,j) = cov_bxy_p1(j,i) = ld_r_mat(kept_ID[i],kept_ID[j])*zscore_inv2(i)*zscore_inv2(j);
            cov_bxy_p2(i,j) = cov_bxy_p2(j,i) = ld_r_mat(kept_ID[i],kept_ID[j])*zscore_inv1(i)*zscore_inv1(j);
        }
    }
}

void est_cov_bxy_gsmr(eigenMatrix &cov_bxy, double bxy_hat, eigenVector bzx, eigenVector bzx_se, eigenVector bzy_se, eigenMatrix ld_r_mat,
                            map<string, int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID) {
    int n_indices_snp = kept_ID.size();
    eigenMatrix cov_bxy_p1(n_indices_snp, n_indices_snp), cov_bxy_p2(n_indices_snp, n_indices_snp);
    est_cov_bxy(cov_bxy_p1, cov_bxy_p2, bzx, bzx_se, bzy_se, ld_r_mat,
                   meta_snp_name_map, indices_snp, kept_ID);
    cov_bxy = cov_bxy_p1 + cov_bxy_p2*(bxy_hat*bxy_hat);
}

void est_cov_bxy_heidi(eigenMatrix &cov_bxy, double bxy_hat, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, eigenMatrix ld_r_mat,
                            map<string, int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID) {
    int n_indices_snp = kept_ID.size();
    eigenMatrix cov_bxy_p1(n_indices_snp, n_indices_snp), cov_bxy_p2(n_indices_snp, n_indices_snp);

    est_cov_bxy(cov_bxy_p1, cov_bxy_p2, bzx, bzx_se, bzy_se, ld_r_mat,
                   meta_snp_name_map, indices_snp, kept_ID);
    
    eigenMatrix cov_bxy_buf(cov_bxy_p2);
    #pragma omp parallel for
    for(int i=0; i<n_indices_snp; i++) {
        map<string,int>::iterator iter = meta_snp_name_map.find(indices_snp[kept_ID[i]]);
        double bxy_i = bzy(iter->second)/bzx(iter->second);
        cov_bxy_p2.col(i) = cov_bxy_buf.col(i)*(bxy_i*bxy_hat);
    }
    cov_bxy = cov_bxy_p1 + cov_bxy_p2.transpose();
}

vector<double> est_bxy_gsmr(eigenVector bxy_sort, eigenVector bxy, vector<string> indices_snp, vector<int> kept_ID, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, eigenMatrix ld_r_mat, map<string, int> meta_snp_name_map, eigenVector &vec_1t_v) {
    double bxy_median = CommFunc::quantile(bxy_sort, 0.50);

    // cov(bxy_i, bxy_j)
    int i = 0, n_indices_snp = kept_ID.size();
    eigenMatrix cov_bxy(n_indices_snp, n_indices_snp);
    est_cov_bxy_gsmr(cov_bxy, bxy_median, bzx, bzx_se, bzy_se, ld_r_mat,
                    meta_snp_name_map, indices_snp, kept_ID);

    eigenVector bxy_kept(n_indices_snp);
    for(i = 0; i < n_indices_snp; i++) bxy_kept(i) = bxy(kept_ID[i]);

    // gsmr
    VectorXd vec_1= VectorXd::Ones(n_indices_snp);
    cov_bxy = cov_bxy + eps*eigenMatrix::Identity(n_indices_snp, n_indices_snp);
    LDLT<eigenMatrix> ldlt_cov_bxy(cov_bxy);

    if( ldlt_cov_bxy.vectorD().minCoeff() <= 0 )
        LOGGER.e(0, "the variance-covariance matrix of bxy is not invertible.");
    eigenMatrix cov_bxy_inv = eigenMatrix::Identity(n_indices_snp, n_indices_snp);
    ldlt_cov_bxy.solveInPlace(cov_bxy_inv);
    vec_1t_v = vec_1.transpose()*cov_bxy_inv;
    double mat_buf = vec_1t_v.dot(vec_1);
    double bxy_gsmr_se  = 1/mat_buf;
    mat_buf = vec_1t_v.dot(bxy_kept);
    double bxy_gsmr = bxy_gsmr_se*mat_buf;
    double bxy_gsmr_pval = StatFunc::pchisq( bxy_gsmr*bxy_gsmr/bxy_gsmr_se, 1);
    bxy_gsmr_se = sqrt(bxy_gsmr_se);
    // save the result
    vector<double> rst(4);
    rst[0] = bxy_gsmr; rst[1] = bxy_gsmr_se; rst[2] = bxy_gsmr_pval; rst[3] = (double)n_indices_snp;
    return rst;
}

void est_var_bxy(eigenVector &var_bxy, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se,
                        vector<string> indices_snp, vector<int> kept_ID, map<string,int> meta_snp_name_map) {
    int n_indices_snp = kept_ID.size();
    
    #pragma omp parallel for
    for(int i = 0; i < n_indices_snp; i++) {
        map<string, int>::iterator iter;
        iter = meta_snp_name_map.find(indices_snp[kept_ID[i]]);
        double d_buf = pow(bzx_se(iter->second),2.0)*pow(bzy(iter->second), 2.0)/pow(bzx(iter->second), 4.0) + pow(bzy_se(iter->second), 2.0)/pow(bzx(iter->second), 2.0);
        var_bxy(i) = d_buf;
    }                        
}

void est_diff_bxy(eigenVector &bxy_diff, double bxy_est, vector<string> indices_snp, vector<int> kept_ID, 
                             eigenVector bzx, eigenVector bzy, map<string,int> meta_snp_name_map) {
    int n_indices_snp = kept_ID.size(); 

    #pragma omp parallel for
    for(int i=0; i<n_indices_snp; i++ ) {
        map<string,int>::iterator iter = meta_snp_name_map.find(indices_snp[kept_ID[i]]);
        double d_buf = bzy(iter->second)/bzx(iter->second) - bxy_est;
        bxy_diff(i) = d_buf;
    }
}

void est_var_diff_bxy(eigenMatrix &var_d, eigenMatrix cov_bxy, eigenVector cov_bxy_bgsmr, double bxy_hat_se, int n_indices_snp) {
    var_d.setZero(n_indices_snp, n_indices_snp);

    #pragma omp parallel for
    for(int i=0; i<n_indices_snp; i++ ) {
        for(int j=i; j<n_indices_snp; j++) {
            double d_buf = cov_bxy(i,j) + bxy_hat_se*bxy_hat_se - cov_bxy_bgsmr(i) - cov_bxy_bgsmr(j);
            var_d(i,j) = var_d(j,i) = d_buf;
        }
    }
}

vector<double> indi_heidi_pvalue(double bxy_hat, double bxy_hat_se, eigenVector vec_1t_v, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, 
                              eigenMatrix ld_r_mat, map<string,int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID, vector<int> remain_index) {
    int i = 0, j = 0, n_indices_snp = kept_ID.size();

    // var(bxy)
    eigenVector var_bxy(n_indices_snp);
    est_var_bxy(var_bxy, bzx, bzx_se, bzy, bzy_se, indices_snp, kept_ID, meta_snp_name_map);
    // cov(bxy_i, bxy_j)
    eigenMatrix cov_bxy(n_indices_snp, n_indices_snp);
    est_cov_bxy_heidi(cov_bxy, bxy_hat, bzx, bzx_se, bzy, bzy_se, ld_r_mat, meta_snp_name_map, indices_snp, kept_ID);

    vector<int> slct_index(remain_index.size());
    for(i=0, j=0; i<n_indices_snp; i++) {
        if(kept_ID[i]!=remain_index[j]) continue;
        slct_index[j++] = i;
    }

    extractCol(cov_bxy, slct_index);

    eigenVector cov_bxy_bgsmr(n_indices_snp);
    #pragma omp parallel for
    for(int k=0; k<n_indices_snp; k++)
        cov_bxy_bgsmr(k) = bxy_hat_se*bxy_hat_se*vec_1t_v.dot(cov_bxy.row(k));
    
    // d = bxy_i - bxy_gsmr
    eigenVector d(n_indices_snp);
    est_diff_bxy(d, bxy_hat, indices_snp, kept_ID, bzx, bzy, meta_snp_name_map);

    eigenVector var_d = var_bxy + bxy_hat_se*bxy_hat_se*eigenVector::Ones(n_indices_snp) - 2*cov_bxy_bgsmr;
    vector<double> indi_het_pval(n_indices_snp);
    for(i = 0; i < n_indices_snp; i++) 
        indi_het_pval[i] = StatFunc::pchisq(d(i)*d(i)/var_d(i), 1);

    return(indi_het_pval);
}

double global_heidi_pvalue(double bxy_hat, double bxy_hat_se, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, 
                           eigenMatrix ld_r_mat, eigenVector vec_1t_v, vector<string> indices_snp, vector<int> kept_ID, map<string, int> meta_snp_name_map) {
    int n_indices_snp = kept_ID.size();

    // d = bxy_i - bxy_gsmr
    eigenVector d(n_indices_snp);
    est_diff_bxy(d, bxy_hat, indices_snp, kept_ID, bzx, bzy, meta_snp_name_map);
    // V matrix
    eigenMatrix cov_bxy(n_indices_snp, n_indices_snp);
    est_cov_bxy_gsmr(cov_bxy, bxy_hat, bzx, bzx_se, bzy_se, ld_r_mat, meta_snp_name_map, indices_snp, kept_ID);
    eigenVector cov_bxy_bgsmr(n_indices_snp);
    #pragma omp parallel for
    for(int i=0; i<n_indices_snp; i++) cov_bxy_bgsmr(i) = bxy_hat_se*bxy_hat_se*vec_1t_v.dot(cov_bxy.col(i));
    eigenMatrix var_d(n_indices_snp, n_indices_snp);
    est_var_diff_bxy(var_d, cov_bxy, cov_bxy_bgsmr, bxy_hat_se, n_indices_snp);

    // inverse Vd matrix
    double global_heidi_pval;
    var_d = var_d + eps*eigenMatrix::Identity(n_indices_snp, n_indices_snp);

    // pchisqsum method
    eigenMatrix corr_d(n_indices_snp, n_indices_snp);
    corr_d.setIdentity(n_indices_snp, n_indices_snp);
    #pragma omp parallel for
    for(int i=0; i<(n_indices_snp-1); i++) 
        for(int j=i+1; j<n_indices_snp; j++) 
            corr_d(i,j) = corr_d(j,i) = var_d(i,j) / sqrt(var_d(i,i)*var_d(j,j));
    double chival = (d.array()*d.array()/var_d.diagonal().array()).sum();
    SelfAdjointEigenSolver<eigenMatrix> saes(corr_d);
    global_heidi_pval = StatFunc::pchisqsum(chival, saes.eigenvalues().cast<double>());

    return(global_heidi_pval);  
}

// Excluding SNPs based on HEIDI-outlier p-value at individual SNP
vector<int> indi_heidi_outlier(eigenVector bxy, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, 
                               eigenMatrix ld_r_mat, map<string,int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID) {
    
    int i = 0, n_indices_snp = kept_ID.size();
    eigenVector bxy_sort(n_indices_snp);
    for(i=0; i<n_indices_snp; i++) bxy_sort(i) = bxy(kept_ID[i]);
    std::stable_sort(bxy_sort.data(), bxy_sort.data()+n_indices_snp);
    vector<int> ci_index = init_interval_bxy(bxy_sort, bxy, kept_ID, 0.1, 0.9);
    if(ci_index.size() == 0) return(ci_index);

    eigenVector vec_1t_v;
    vector<double> gsmr_rst = est_bxy_gsmr(bxy_sort, bxy, indices_snp, ci_index, bzx, bzx_se, bzy, bzy_se, ld_r_mat, meta_snp_name_map, vec_1t_v);    
    double bxy_hat = gsmr_rst[0], bxy_hat_se = gsmr_rst[1];

    vector<double> indi_het_pval;
    indi_het_pval = indi_heidi_pvalue(bxy_hat, bxy_hat_se, vec_1t_v, bzx, bzx_se, bzy, bzy_se, 
                                    ld_r_mat, meta_snp_name_map, indices_snp, kept_ID, ci_index);

    vector<int> kept_ID_buf(kept_ID);
    kept_ID.clear();
    double indi_heidi_thresh = min(0.01, 0.05/(double)n_indices_snp);
    for(i=0; i<n_indices_snp; i++) {
        if(indi_het_pval[i] < indi_heidi_thresh) continue;
        kept_ID.push_back(kept_ID_buf[i]); 
    }

    return(kept_ID);
}

vector<int> indi_heidi_outlier_iter(eigenVector bxy, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, 
                               eigenMatrix ld_r_mat, map<string,int> meta_snp_name_map, vector<string> indices_snp, vector<int> kept_ID) {
    int i = 0, n_indices_snp = kept_ID.size();
    eigenVector bxy_sort(n_indices_snp);
    for(i=0; i<n_indices_snp; i++) bxy_sort(i) = bxy(kept_ID[i]);
    std::stable_sort(bxy_sort.data(), bxy_sort.data()+n_indices_snp);
    vector<int> ci_index = init_interval_bxy(bxy_sort, bxy, kept_ID, 0.1, 0.9);

    if(ci_index.size() == 0) return(ci_index);

    eigenVector vec_1t_v;
    vector<double> gsmr_rst = est_bxy_gsmr(bxy_sort, bxy, indices_snp, ci_index, bzx, bzx_se, bzy, bzy_se, ld_r_mat, meta_snp_name_map, vec_1t_v);    

    double bxy_hat = gsmr_rst[0], bxy_hat_se = gsmr_rst[1];

    vector<double> indi_het_pval;
    indi_het_pval = indi_heidi_pvalue(bxy_hat, bxy_hat_se, vec_1t_v, bzx, bzx_se, bzy, bzy_se, 
                                    ld_r_mat, meta_snp_name_map, indices_snp, kept_ID, ci_index);

    double indi_heidi_thresh = min(0.01, 0.05/(double)n_indices_snp);
//    double indi_heidi_thresh = 0.01;
    double min_pval = 1.0;
    int min_index = 0;
    for(i=0; i<n_indices_snp; i++) {
        if(indi_het_pval[i] >= min_pval) continue;
        min_pval = indi_het_pval[i];
        min_index = i;
    }

    if(min_pval < indi_heidi_thresh) kept_ID.erase(kept_ID.begin()+min_index);

    return(kept_ID);
}

int topsnp_heidi_outlier(eigenVector bxy, vector<int> kept_ID, vector<int> ci_index, eigenVector bzx, eigenVector bzx_se, vector<string> indices_snp, map<string,int> meta_snp_name_map) {
    int nbins = 5, n_candidate_snp = ci_index.size();
    // 1. 5 bins of bxy
    double hist_prob[6] = {0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double lower_bounder[5], upper_bounder[5];

    eigenVector bxy_sort(n_candidate_snp);
    for(int i = 0; i < n_candidate_snp; i++) bxy_sort(i) = bxy(kept_ID[ci_index[i]]);
    std::stable_sort(bxy_sort.data(), bxy_sort.data()+n_candidate_snp);

    for(int i = 0; i < nbins; i++) {
        lower_bounder[i] = CommFunc::quantile(bxy_sort, hist_prob[i]);
        upper_bounder[i] = CommFunc::quantile(bxy_sort, hist_prob[i+1]);
    }
    
    // 2. calculate mean(chi^2) in each bin
    int topsnp_bin[5];
    double chisq_bin[5], maximum_chisq[5], nsnps_bin[5];
    for(int i = 0; i < nbins; i++) {
        topsnp_bin[i] = 0; 
        chisq_bin[i] = maximum_chisq[i] = nsnps_bin[i] = 0.0;
    }
    for(int i = 0; i < n_candidate_snp; i++) {
        for(int j = 0; j < nbins; j++) {
            if(bxy(kept_ID[ci_index[i]]) < lower_bounder[j] || bxy(kept_ID[ci_index[i]]) >= upper_bounder[j]) continue;
            map<string,int>::iterator iter = meta_snp_name_map.find(indices_snp[kept_ID[ci_index[i]]]);
            double chisqbuf = (bzx(iter->second)*bzx(iter->second))/(bzx_se(iter->second)*bzx_se(iter->second));
            chisq_bin[j] += chisqbuf; nsnps_bin[j]++;
            if(maximum_chisq[j] < chisqbuf) {
                maximum_chisq[j] = chisqbuf;
                topsnp_bin[j] = ci_index[i];
            }
            break;
        }
    }
    
    // 3. select the top SNP
    int topsnp_index = 0;
    double max_chisq_bin = 0.0;
    for(int i = 0; i < nbins; i++) {
        if(nsnps_bin[i] < 1.0) continue;
        chisq_bin[i] /= nsnps_bin[i];        
        if(max_chisq_bin < chisq_bin[i]) {
            max_chisq_bin = chisq_bin[i];
            topsnp_index = topsnp_bin[i];
        }
    }

    return(topsnp_index);
}

double heidi_outlier_topsnp(eigenVector &indi_het_pval, eigenVector bxy, eigenMatrix cov_bxy, vector<int> kept_ID, eigenVector bzx, eigenVector bzy, int topsnp_index, vector<string> indices_snp, map<string,int> meta_snp_name_map, bool multi_snp_heidi_flag) {
    int n_indices_snp = kept_ID.size();

    // 1. d = bxy_i - bxy_top
    eigenVector d(n_indices_snp);
    est_diff_bxy(d, bxy(kept_ID[topsnp_index]), indices_snp, kept_ID, bzx, bzy, meta_snp_name_map);

    // 2. var(d)  
    eigenMatrix var_d = cov_bxy + cov_bxy(topsnp_index, topsnp_index)*eigenMatrix::Ones(n_indices_snp, n_indices_snp);
    for(int i=0; i<n_indices_snp; i++) var_d.col(i) -= cov_bxy.col(topsnp_index);
    for(int i=0; i<n_indices_snp; i++) var_d.row(i) -= cov_bxy.col(topsnp_index);
    var_d(topsnp_index, topsnp_index) = 1;

    // 3. d^2/var(d)
    indi_het_pval.resize(n_indices_snp);
    for(int i = 0; i < n_indices_snp; i++) {
        indi_het_pval[i] = StatFunc::pchisq(d(i)*d(i)/var_d(i,i), 1);
    }
    if(!multi_snp_heidi_flag) {
        // single-SNP HEIDI-outlier
        return 0.0;
    } else {
        // multi-SNP HEIDI-outlier
        // exclude the top SNP
        removeElement(d, topsnp_index);
        removeRow(var_d, topsnp_index); removeColumn(var_d, topsnp_index);
        // Saddle point method
        eigenMatrix corr_d(n_indices_snp-1, n_indices_snp-1);
        corr_d.setIdentity(n_indices_snp-1, n_indices_snp-1);
        #pragma omp parallel for
        for(int i=0; i<(n_indices_snp-2); i++) 
            for(int j=i+1; j<(n_indices_snp-1); j++) 
                corr_d(i,j) = corr_d(j,i) = var_d(i,j) / sqrt(var_d(i,i)*var_d(j,j));
        double chival = (d.array()*d.array()/var_d.diagonal().array()).sum();            
        SelfAdjointEigenSolver<eigenMatrix> saes(corr_d);
        double global_heidi_pval = StatFunc::pchisqsum(chival, saes.eigenvalues().cast<double>());
        return(global_heidi_pval);
    }
}

vector<int> keep_non_associate_snp(eigenVector bzy_pval, vector<string> indices_snp, vector<int> kept_ID, map<string,int> meta_snp_name_map, double pval_thresh) {
    int n_indices_snp = kept_ID.size();
    vector<int> ci_index;

    for(int i=0; i<n_indices_snp; i++) {
        map<string,int>::iterator iter = meta_snp_name_map.find(indices_snp[kept_ID[i]]);
        if(bzy_pval(iter->second) < pval_thresh) continue;
        ci_index.push_back(i);
    }

    return(ci_index);
}

int indi_heidi_outlier_topsnp_iter(eigenVector bxy, eigenMatrix &cov_bxy, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, eigenVector bzy_pval, 
                               eigenMatrix ld_r_mat, map<string,int> meta_snp_name_map, vector<string> indices_snp, vector<int> &kept_ID, double heidi_thresh) {
    int n_indices_snp = kept_ID.size();
    double pval_thresh = 5e-8;
    vector<int> ci_index;
    // 1. removing SNPs with significant bzy
    //ci_index = keep_non_associate_snp(bzy_pval, indices_snp, kept_ID, meta_snp_name_map, pval_thresh);
    ci_index.resize(n_indices_snp);
    for(int i=0; i<n_indices_snp; i++) ci_index[i] = i;

    // 2. removing SNPs with extreme bxy, [0.05, 0.95]
    int n_candidate_snp = ci_index.size();
    eigenVector bxy_sort(n_candidate_snp), bxy_tmp(n_indices_snp);
    for(int i=0; i<n_candidate_snp; i++) bxy_sort(i) = bxy(kept_ID[ci_index[i]]);
    for(int i=0; i<n_indices_snp; i++) bxy_tmp(i) = bxy(kept_ID[i]);
    std::stable_sort(bxy_sort.data(), bxy_sort.data()+n_candidate_snp);
    ci_index = init_interval_bxy(bxy_sort, bxy_tmp, ci_index, 0.05, 0.95);
    if(ci_index.size() == 0) { kept_ID.clear(); return 0; }

    // 3. choose the top SNP
    int topsnp_index = topsnp_heidi_outlier(bxy, kept_ID, ci_index, bzx, bzx_se, indices_snp, meta_snp_name_map);

    // 4. single-SNP HEIDI-outlier analysis
    eigenVector indi_het_pval(n_indices_snp);
    heidi_outlier_topsnp(indi_het_pval, bxy, cov_bxy, kept_ID, bzx, bzy, topsnp_index, indices_snp, meta_snp_name_map, false);

    // 5. remove the SNP with the smallest HEIDI p-value
    std::ptrdiff_t filtered_snp_index = 0;
    double minimum_pvalue = 1.0;
    minimum_pvalue = indi_het_pval.minCoeff(&filtered_snp_index);

    if(minimum_pvalue < heidi_thresh) {
        kept_ID.erase(kept_ID.begin()+filtered_snp_index);
        removeRow(cov_bxy, filtered_snp_index);
        removeColumn(cov_bxy, filtered_snp_index);
    }

    return 1;
}

double global_heidi_outlier_topsnp_iter(eigenVector bxy, eigenMatrix &cov_bxy, eigenVector bzx, eigenVector bzx_se, eigenVector bzy, eigenVector bzy_se, eigenVector bzy_pval, 
                               eigenMatrix ld_r_mat, map<string,int> meta_snp_name_map, vector<string> indices_snp, vector<int> &kept_ID, double heidi_thresh) {
    int n_indices_snp = kept_ID.size();
    double pval_thresh = 5e-8;
    vector<int> ci_index;
    // 1. removing SNPs with significant bzy
    ci_index = keep_non_associate_snp(bzy_pval, indices_snp, kept_ID, meta_snp_name_map, pval_thresh);

    // 2. removing SNPs with extreme bxy, [0.05, 0.95]
    int n_candidate_snp = ci_index.size();
    eigenVector bxy_sort(n_candidate_snp), bxy_tmp(n_indices_snp);
    for(int i=0; i<n_candidate_snp; i++) bxy_sort(i) = bxy(kept_ID[ci_index[i]]);
    for(int i=0; i<n_indices_snp; i++) bxy_tmp(i) = bxy(kept_ID[i]);
    std::stable_sort(bxy_sort.data(), bxy_sort.data()+n_candidate_snp);
    ci_index = init_interval_bxy(bxy_sort, bxy_tmp, ci_index, 0.05, 0.95);

    if(ci_index.size() == 0) { kept_ID.clear(); return(-9); }

    // 3. choose the top SNP
   int topsnp_index = topsnp_heidi_outlier(bxy, kept_ID, ci_index, bzx, bzx_se, indices_snp, meta_snp_name_map);

    // 4. single-SNP and multi-SNPs HEIDI-outlier analyses
    double global_heidi_pvalue = 1.0;
    eigenVector indi_het_pval(n_indices_snp);
    global_heidi_pvalue = heidi_outlier_topsnp(indi_het_pval, bxy, cov_bxy, kept_ID, bzx, bzy, topsnp_index, indices_snp, meta_snp_name_map, true);

    // 5. return the index of the SNP with the smallest HEIDI p-value, if global_heidi_pval < threshold
    if(global_heidi_pvalue < heidi_thresh) {
        std::ptrdiff_t filtered_snp_index = 0;
        double minimum_pvalue = 1.0;
        minimum_pvalue = indi_het_pval.minCoeff(&filtered_snp_index);
        kept_ID.erase(kept_ID.begin()+filtered_snp_index);
        removeRow(cov_bxy, filtered_snp_index);
        removeColumn(cov_bxy, filtered_snp_index);
    }

    return(global_heidi_pvalue);
}

vector<double> gcta::gsmr_meta(vector<string> &snp_instru, eigenVector bzx, eigenVector bzx_se, eigenVector bzx_pval, eigenVector bzy, eigenVector bzy_se, eigenVector bzy_pval, double rho_pheno, vector<bool> snp_flag, double gwas_thresh, int wind_size, double r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, string &pleio_snps, string &err_msg) {

    int i=0, j=0, nsnp = _include.size(), nindi=_keep.size();    
    vector<string> indices_snp;
    vector<double> rst(_n_gsmr_rst_item);
    for(i=0; i<_n_gsmr_rst_item; i++) rst[i] = nan("");

    if(nsnp < nsnp_gsmr) {
        err_msg = "Not enough SNPs to perform the GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required."; 
        return rst;
    }

    // clumping analysis
    eigenVector bzx_chival=bzx.array()*bzx.array()/(bzx_se.array()*bzx_se.array());
    indices_snp = clumping_meta(bzx_chival, snp_flag, gwas_thresh, wind_size, r2_thresh);
    int n_indices_snp = indices_snp.size();

    //LOGGER.i(0, to_string(n_indices_snp) + " index SNPs are obtained from the clumping analysis of the " + nsnp + " genome-wide significant SNPs.");

    std::stringstream ss1, ss2;
    ss1 << std::scientific << std::setprecision(1) << gwas_thresh;
    ss2 << std::fixed << std::setprecision(2) << r2_thresh;
    if(n_indices_snp < nsnp_gsmr) {
        LOGGER.i(0, to_string(n_indices_snp) + " index SNPs are obtained from the clumping analysis with p < " + ss1.str() + " and LD r2 < " + ss2.str() + ".");
        err_msg = "Not enough SNPs to perform the GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required."; 
        return rst;
    }

    // LD r
    MatrixXf x_sub(nindi, n_indices_snp);
    vector<int> snp_sn(n_indices_snp);
    map<string, int>::iterator iter;
    int i_buf = 0;
    for( i = 0; i < n_indices_snp; i++ ) {
        iter = _snp_name_map.find(indices_snp[i]);
        i_buf = iter->second;
        snp_sn[i] = find(_include.begin(), _include.end(), i_buf) - _include.begin();
    }

    // construct x coefficient
    eigenMatrix ld_r_mat(n_indices_snp, n_indices_snp);
    make_XMat_subset(x_sub, snp_sn, true);

    double x_sd1 = 0.0, x_sd2 = 0.0, x_cov = 0.0;
    ld_r_mat = MatrixXd::Identity(n_indices_snp, n_indices_snp);
    for(i=0; i<(n_indices_snp-1); i++) {
        x_sd1 = x_sub.col(i).norm();
        for(j=(i+1); j<n_indices_snp; j++) {
            x_cov = x_sub.col(i).dot(x_sub.col(j));
            x_sd2 = x_sub.col(j).norm();
            ld_r_mat(i,j) = ld_r_mat(j,i) = x_cov/(x_sd1*x_sd2);
        }
    }
   
    // LD pruning
    vector<int> kept_ID;
    kept_ID = rm_cor_elements(ld_r_mat, r2_thresh, true);
    n_indices_snp = kept_ID.size();

    LOGGER.i(0, to_string(n_indices_snp) + " index SNPs are obtained from the clumping analysis with p < " + ss1.str() + " and LD r2 < " + ss2.str() + ".");
  
    if(n_indices_snp < nsnp_gsmr) {
        err_msg = "Not enough SNPs to perform the GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required."; 
        return rst;
    }

    // Adjust LD
    int k = 0;
    eigenMatrix ld_pval_mat(n_indices_snp, n_indices_snp);
    vector<pair<double,int>> ld_pval(n_indices_snp*(n_indices_snp-1)/2);

    for(i=0, k=0; i<(n_indices_snp-1); i++) {
        for(j=(i+1); j<n_indices_snp; j++) {
            ld_pval[k++] = make_pair(StatFunc::chi_prob(1, ld_r_mat(kept_ID[i],kept_ID[j])*ld_r_mat(kept_ID[i],kept_ID[j])*(double)nindi), i*n_indices_snp+j); 
        }
    }

    stable_sort(ld_pval.begin(), ld_pval.end(), [](const pair<double,int> a, const pair<double,int> b) {return a.first > b.first;});
    adjust_ld_r_fdr(ld_r_mat, kept_ID, ld_pval, n_indices_snp, ld_fdr_thresh);

    // estimate bxy
    eigenVector bxy(indices_snp.size());
    bxy.setZero();
    for(i=0; i<n_indices_snp; i++) {
        iter = _meta_snp_name_map.find(indices_snp[kept_ID[i]]);
        bxy(kept_ID[i]) = bzy(iter->second)/bzx(iter->second);
    }

    vector<int> remain_index(kept_ID);
    bool heidi_flag = CommFunc::FloatNotEqual(abs(global_heidi_thresh) + abs(std_heidi_thresh), 0) ? true : false;

    // estimate cov(bxy) matrix for HEIDI-outlier
    eigenMatrix cov_bxy_p1(n_indices_snp, n_indices_snp), cov_bxy_p2(n_indices_snp, n_indices_snp);
    est_cov_bxy(cov_bxy_p1, cov_bxy_p2, bzx, bzx_se, bzy_se, ld_r_mat, _meta_snp_name_map, indices_snp, kept_ID);
    eigenVector bxy_kept(n_indices_snp);
    eigenMatrix cov_bxy(n_indices_snp, n_indices_snp);
    eigenMatrix cov_bxy_buf(cov_bxy_p2);
    for(i = 0; i < n_indices_snp; i++) bxy_kept(i) = bxy(kept_ID[i]);
    for(i = 0; i < n_indices_snp; i++) {
        cov_bxy_p2.col(i) = cov_bxy_buf.col(i).array()*(bxy_kept*bxy_kept(i)).array();
    }
    cov_bxy = cov_bxy_p1 + cov_bxy_p2;

    if(_gsmr_beta_version & heidi_flag) {
        // HEIDI-outlier step 1
        int n_kept_snp = 0;
        while(1) {          
            n_kept_snp = kept_ID.size();
            int process = indi_heidi_outlier_topsnp_iter(bxy, cov_bxy, bzx, bzx_se, bzy, bzy_se, bzy_pval, 
                        ld_r_mat, _meta_snp_name_map, indices_snp, kept_ID, std_heidi_thresh);
            n_indices_snp = kept_ID.size();
            if(n_kept_snp - n_indices_snp == 0 ) break;
            // check
            if(n_indices_snp < nsnp_gsmr || process == 0) {
                err_msg = "Not enough SNPs to perform GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required.";
                return rst;
            }
        }

        // HEIDI-outlier step 2
        double global_heidi_pval = 0.0;
        while(1) {
            // bxy by GSMR
            n_kept_snp = kept_ID.size();
            global_heidi_pval = global_heidi_outlier_topsnp_iter(bxy, cov_bxy, bzx, bzx_se, bzy, bzy_se, bzy_pval, 
                               ld_r_mat, _meta_snp_name_map, indices_snp, kept_ID, global_heidi_thresh);
            n_indices_snp = kept_ID.size();         
            if(n_kept_snp - n_indices_snp == 0 ) break;
            // check
            if(n_indices_snp < nsnp_gsmr || global_heidi_pval < 0) {
                err_msg = "Not enough SNPs to perform GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required.";
                return rst;
            }         
        }
        rst[4] = global_heidi_pval;
    } else if(!_gsmr_beta_version & heidi_flag) {
        // HEIDI-outlier using top SNP
        // Step 1, find the top SNP
        int topsnp_index = -9;
        topsnp_index = topsnp_bxy(bxy, bzx, bzx_se, _meta_snp_name_map, indices_snp, kept_ID, n_indices_snp);
        if(topsnp_index < 0)  {
            err_msg = "Not enough SNPs to perform the GSMR analysis. At least " + to_string(nsnp_gsmr) + " SNPs are required.";
            return rst;
        }

        // Step 2, estimate HEIDI-outlier p-value
        // d = bxy_i - bxy_gsmr
        eigenVector indi_het_pval(n_indices_snp);
        heidi_outlier_topsnp(indi_het_pval, bxy, cov_bxy, kept_ID, bzx, bzy, topsnp_index, indices_snp, _meta_snp_name_map, false);        
        vector<int> kept_ID_tmp;
        for(i = 0; i < n_indices_snp; i++) {
            if(i==topsnp_index) {
                kept_ID_tmp.push_back(kept_ID[i]);
                continue;
            }
            if(indi_het_pval(i) >= std_heidi_thresh) kept_ID_tmp.push_back(kept_ID[i]);
        }
        kept_ID = kept_ID_tmp;
        n_indices_snp = kept_ID.size();
    }
    // estimate bxy(gsmr)
    bxy_kept.resize(n_indices_snp);   
    for(i = 0; i < n_indices_snp; i++) bxy_kept(i) = bxy(kept_ID[i]);
    std::stable_sort(bxy_kept.data(), bxy_kept.data()+n_indices_snp);

    eigenVector vec_1t_v;
    vector<double> gsmr_rst = est_bxy_gsmr(bxy_kept, bxy, indices_snp, kept_ID, bzx, bzx_se, bzy, bzy_se, ld_r_mat, _meta_snp_name_map, vec_1t_v);
    rst[0] = gsmr_rst[0]; rst[1] = gsmr_rst[1]; 
    rst[2] = StatFunc::pchisq(rst[0]*rst[0]/(rst[1]*rst[1]), 1); rst[3] = n_indices_snp;

    // save pleiotropic SNPs
    int npleio = remain_index.size() - n_indices_snp;
    if(npleio > 0) {
        int nsnp_kept = remain_index.size();
        for(i=0, j=0; i<nsnp_kept && j<n_indices_snp; i++) {
            if(remain_index[i] == kept_ID[j]) { j++; continue; }
            pleio_snps += indices_snp[remain_index[i]] + ",";
        }
        if(i < nsnp_kept) {
            for(; i < nsnp_kept; i++) pleio_snps += indices_snp[remain_index[i]] + ",";
        }
        pleio_snps = pleio_snps.substr(0, pleio_snps.size()-1);
        LOGGER.i(0, to_string(npleio) + " pleiotropic SNPs are filtered by HEIDI-outlier analysis.");	
    }
    // save the SNPs after HEIDI-outlier
    snp_instru.clear(); snp_instru.resize(n_indices_snp);
    for(i=0; i<n_indices_snp; i++)
        snp_instru[i] = indices_snp[kept_ID[i]];
    LOGGER.i(0, to_string(n_indices_snp) + " index SNPs are retained after HEIDI-outlier analysis.");    
    return rst;
}

int read_ld_marker(string ref_ld_dirt) {
    // Read the number of markers
    int i = 0, i_buf = 0, ttl_mk_num = 0, chr_num = 22;
    string filestr = "", strbuf = "";
    for(i=0; i<chr_num; i++) {
        filestr = ref_ld_dirt + to_string(i+1)+ ".l2.M_5_50";
        ifstream ref_marker(filestr.c_str());
        if (!ref_marker) LOGGER.e(0, "cannot open the file [" + filestr + "] to read.");
        std::getline(ref_marker, strbuf);
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 1) LOGGER.e(0, "the format of file [" + filestr + "] is incorrect.");
        i_buf = atoi(line_elements[0].c_str());
        ttl_mk_num += i_buf;
        ref_marker.close();
    }
    return ttl_mk_num;
}

vector<string> read_ld_score_txt(string filestr, map<string,int> snplist_map, vector<double> &ld_score) {
    
    int line_number = 0, indxbuf = 0, ldsc_index = 0;
    double ldscbuf = 0.0;
    string strbuf = "", snpbuf = "";
    vector<string> ld_score_snps;
    map<string,int>::iterator snp_iter;

    ifstream ldsc_marker(filestr.c_str());
    while(std::getline(ldsc_marker, strbuf)) {
        line_number ++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() == 4) {
            ldsc_index = 3;
        } else if(line_elements.size() == 6) {
            ldsc_index = 5;
        } else {
            LOGGER.e(0, "the format of file [" + filestr + "] is incorrect, line " + to_string(line_number) + ".");
        }
        if(line_number==1) continue;
        snpbuf = line_elements[1];
        ldscbuf  = atof(line_elements[ldsc_index].c_str());
        // save the data
        snp_iter = snplist_map.find(snpbuf);
        if(snp_iter!=snplist_map.end()) {
            indxbuf = snp_iter->second;
            ld_score[indxbuf] = ldscbuf;
            ld_score_snps.push_back(snpbuf);
        }
    }
    ldsc_marker.close();
    return ld_score_snps;
}

vector<string> read_ld_score_gz(string filestr, map<string,int> snplist_map, vector<double> &ld_score) {
    int indxbuf = 0, line_number = 0, ldsc_index = 0;
    double ldscbuf = 0.0;
    string snpbuf = "";
    vector<string> ld_score_snps;
    map<string,int>::iterator snp_iter;

    string err_msg = "Failed to read [" + filestr + "]. An error occurs in line ";

    const int MAX_LINE_LENGTH = 1024;
    char strbuf[MAX_LINE_LENGTH];
    gzifstream ldsc_marker(filestr.c_str());
    while(1) {
        ldsc_marker.getline(strbuf, MAX_LINE_LENGTH, '\n');
        if (ldsc_marker.fail() || !ldsc_marker.good()) break;
        line_number ++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() == 4) {
            ldsc_index = 3;
        } else if(line_elements.size() == 6) {
            ldsc_index = 5;
        } else {
            LOGGER.e(0, "the format of file [" + filestr + "] is not correct, line " + to_string(line_number) + ".");
        }
        if(line_number==1) continue;
        snpbuf = line_elements[1];
        ldscbuf  = atof(line_elements[ldsc_index].c_str());

        // save the data
        snp_iter = snplist_map.find(snpbuf);
        if(snp_iter!=snplist_map.end()) {
            indxbuf = snp_iter->second;
            ld_score[indxbuf] = ldscbuf;
            ld_score_snps.push_back(snpbuf);
        }
    }
    ldsc_marker.close();
    return ld_score_snps;
}

vector<string> read_ld_score(string ld_dirt, map<string,int> snplist_map, int nsnp, vector<double> &ld_score) {

    // Read the reference / weighted LD score
    int i = 0, chr_num = 22;
    string filestr_t1 = "", filestr_t2 = "";
    vector<string> ld_score_snps, snpbuf;
    
    ld_score.clear(); ld_score.resize(nsnp);
    for(i=0; i<nsnp; i++) ld_score[i] = -9.0;
    for(i=0; i<chr_num; i++) {
        filestr_t1 = ld_dirt  + to_string(i+1)+ ".l2.ldscore";
        filestr_t2 = ld_dirt  + to_string(i+1)+ ".l2.ldscore.gz";
        if(file_exists(filestr_t1)) {
            snpbuf = read_ld_score_txt(filestr_t1, snplist_map, ld_score);
            ld_score_snps.insert(ld_score_snps.end(), snpbuf.begin(), snpbuf.end());
        } else if(file_exists(filestr_t2)) {
            snpbuf = read_ld_score_gz(filestr_t2, snplist_map, ld_score);
            ld_score_snps.insert(ld_score_snps.end(), snpbuf.begin(), snpbuf.end());
        }
        else LOGGER.e(0, "cannot open the file [" + filestr_t1 + "] or [" + filestr_t2 + "] to read.");
    }
    
    return ld_score_snps;
}

eigenVector update_weights_hsq(double intercept, double h, int ttl_mk_num, int n_ld_snp, eigenVector ref_ld, eigenVector w_ld, eigenVector n) {
    int i=0;
    if(h < 0.0) h = 0.0;
    if(h > 1.0) h = 1.0;
    for(i=0; i<n_ld_snp; i++) {
        if(ref_ld(i) < 1.0) ref_ld(i) = 1.0;
        if(w_ld(i) < 1.0) w_ld(i) = 1.0;
    }
    eigenVector intercept_vec = intercept * eigenVector::Ones(n_ld_snp);
    eigenVector denominator = (intercept_vec + h/(double)ttl_mk_num*n.cwiseProduct(ref_ld));
    denominator = 2*w_ld.cwiseProduct(denominator.cwiseProduct(denominator));
    eigenVector numerator = eigenVector::Ones(n_ld_snp);
    numerator  = numerator.cwiseQuotient(denominator);
    return numerator;
}

eigenVector compute_irls(double &intercept, double &hsq, eigenMatrix x, eigenVector y, eigenVector wt, eigenVector ref_ld, eigenVector w_ld, eigenVector n, int ttl_mk_num, int n_ld_snp, bool intercept_flag, bool x_flag) {
    eigenVector w;
    eigenVector wx(x.col(0));
    // x = x*sqrt(w); y = y*sqrt(w);
    wt = wt.cwiseSqrt();
    wt = wt/wt.sum();
    x.col(0) = x.col(0).cwiseProduct(wt);
    x.col(1) = x.col(1).cwiseProduct(wt);
    y = y.cwiseProduct(wt);
    // b = (x'x)^-1x'y
    if(intercept_flag) {
        eigenMatrix xt_x;
        eigenVector xt_y, b_coeff;  
        xt_x = x.transpose()*x;    
        xt_y = x.transpose()*y;
        b_coeff = xt_x.ldlt().solve(xt_y);     
        // update h2
        hsq=b_coeff(0)*(double)ttl_mk_num/n.mean();      
        intercept = b_coeff(1);     
    } else {      
        double xt_x = 0.0, xt_y = 0.0, b_coeff = 0.0;
        xt_x = x.col(0).dot(x.col(0));
        xt_y = x.col(0).dot(y);
        b_coeff = xt_y/xt_x;
        // update h2
        hsq = b_coeff*(double)ttl_mk_num/n.mean();
    }
  
    // weights
    if(x_flag) {       
        w = update_weights_hsq(intercept, hsq, ttl_mk_num, n_ld_snp, wx, w_ld, n);
    } else {     
        w = update_weights_hsq(intercept, hsq, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n);
    }
   
    return w;
}

vector<double> est_hsq_trait_1_step(eigenVector chival, eigenVector n, eigenVector ref_ld, eigenVector w_ld, int n_ld_snp, int ttl_mk_num) {

    // Estimate prior of the weights
    eigenVector denominator = ref_ld.cwiseProduct(n);
    double h_prior = (chival.mean() - 1)*(double)ttl_mk_num/denominator.mean();
    double intercept_prior = 1.0;
    eigenVector wt_ttl = update_weights_hsq(intercept_prior, h_prior, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n);
    eigenMatrix x(n_ld_snp, 2);
    x.col(1).setOnes(); x.col(0) = ref_ld.cwiseProduct(n)/n.mean();
    
    // Iteratively reweighted least squares
    wt_ttl = compute_irls(intercept_prior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, true, false);
    wt_ttl = compute_irls(intercept_prior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, true, false);
    // Estimate intercept and slope
    double intercept_posterior = 0.0, h_posterior = 0.0;
    wt_ttl = compute_irls(intercept_posterior, h_posterior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, true, false);
    
    vector<double> ldsc_est(2);
    ldsc_est[0]  = intercept_posterior; ldsc_est[1] = h_posterior;
    return ldsc_est;
}

vector<double> est_hsq_trait_2_steps(eigenVector chival, eigenVector n, eigenVector ref_ld, eigenVector w_ld, int n_ld_snp, int ttl_mk_num) {
    int i = 0, n_subset_snp = 0;
    double thresh = 30;
    vector<int> subset_indx;
    eigenVector subset_chi, subset_n, subset_ref_ld, subset_w_ld;

    // Estimate prior of the weights
    eigenVector denominator = ref_ld.cwiseProduct(n);
    double h_prior = (chival.mean() - 1)*(double)ttl_mk_num/denominator.mean();
    double intercept_prior = 1.0;
    eigenVector wt_ttl = update_weights_hsq(intercept_prior, h_prior, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n);
    eigenVector subset_wt;
    eigenMatrix x(n_ld_snp,2), subset_x;
    x.col(1).setOnes(); x.col(0) = ref_ld.cwiseProduct(n)/n.mean();

    // SNPs with chi^2 < 30
    for(i=0; i<n_ld_snp; i++) {
        if(chival(i) < thresh) subset_indx.push_back(i);
    }
    n_subset_snp = subset_indx.size();
    
    subset_chi.resize(n_subset_snp); subset_n.resize(n_subset_snp);
    subset_ref_ld.resize(n_subset_snp);  subset_w_ld.resize(n_subset_snp);
    subset_wt.resize(n_subset_snp); subset_x.resize(n_subset_snp,2);  
    for(i=0; i<n_subset_snp; i++) {
        subset_chi(i) = chival(subset_indx[i]);
        subset_n(i) = n(subset_indx[i]);
        subset_ref_ld(i) = ref_ld(subset_indx[i]);
        subset_w_ld(i) = w_ld(subset_indx[i]);
        subset_wt(i) = wt_ttl(subset_indx[i]);
        subset_x.row(i) = x.row(subset_indx[i]);
    }
 
    // Iteratively reweighted least squares
    subset_wt = compute_irls(intercept_prior, h_prior, subset_x, subset_chi, subset_wt, subset_ref_ld, subset_w_ld, subset_n, ttl_mk_num, n_subset_snp, true, true);
    subset_wt = compute_irls(intercept_prior, h_prior, subset_x, subset_chi, subset_wt, subset_ref_ld, subset_w_ld, subset_n, ttl_mk_num, n_subset_snp, true, true);
 
    // Estimate intercept
    double intercept_posterior = 0.0;
    subset_wt = compute_irls(intercept_posterior, h_prior, subset_x, subset_chi, subset_wt, subset_ref_ld, subset_w_ld, subset_n, ttl_mk_num, n_subset_snp, true, true);
  
    // Iteratively reweighted least squares
    chival = chival - eigenVector::Ones(n_ld_snp)*intercept_posterior;
    wt_ttl = compute_irls(intercept_posterior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, false, false);
    wt_ttl = compute_irls(intercept_posterior, h_prior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, false, false);
   
    // Estimate slope
    double h_posterior = 0.0;
    wt_ttl = compute_irls(intercept_posterior, h_posterior, x, chival, wt_ttl, ref_ld, w_ld, n, ttl_mk_num, n_ld_snp, false, false);
   
    vector<double> ldsc_est(2);
    ldsc_est[0]  = intercept_posterior; ldsc_est[1] = h_posterior;

    return ldsc_est;
}

eigenVector update_weights_gcov(double intercept1, double h1, double intercept2, double h2, double intercept_gcov, double gcov, int ttl_mk_num, int n_ld_snp, eigenVector ref_ld, eigenVector w_ld, eigenVector n1, eigenVector n2, eigenVector n_gcov ) {
    int i=0;
    if(h1 < 0.0) h1 = 0.0; if(h1 > 1.0) h1 = 1.0;
    if(h2 < 0.0) h2 = 0.0; if(h2 > 1.0) h2 = 1.0;
    if(gcov < -1.0) gcov = -1.0; if(gcov > 1.0) gcov = 1.0;
    for(i=0; i<n_ld_snp; i++) {
        if(ref_ld(i) < 1.0) ref_ld(i) = 1.0;
        if(w_ld(i) < 1.0) w_ld(i) = 1.0;
    }
    
    double d1=0.0, d2=0.0, d3=0.0;
    eigenVector w(n_ld_snp);

    for(i=0; i<n_ld_snp; i++) {
        d1 = n1(i)*h1*ref_ld(i)/(double)ttl_mk_num + intercept1;
        d2 = n2(i)*h2*ref_ld(i)/(double)ttl_mk_num + intercept2;
        d3 = n_gcov(i)*gcov*ref_ld(i)/(double)ttl_mk_num + intercept_gcov;
        w(i) = 1/( w_ld(i) * (d1*d2 + d3*d3) );
    }
    return w;
}

eigenVector compute_irls_gcov(double &intercept_gcov, double &gcov, eigenMatrix x, eigenVector y, eigenVector wt, eigenVector ref_ld, eigenVector w_ld, eigenVector n_gcov, double intercept1,  double hsq1, eigenVector n1, double intercept2,  double hsq2, eigenVector n2, int ttl_mk_num, int n_ld_snp, bool intercept_flag, bool x_flag) {

    eigenVector w;
    // x = x*sqrt(w); y = y*sqrt(w);
    wt = wt.cwiseSqrt();
    wt = wt/wt.sum();
    x.col(0) = x.col(0).cwiseProduct(wt);
    x.col(1) = x.col(1).cwiseProduct(wt);
    y = y.cwiseProduct(wt);
    // b = (x'x)^-1x'y
    if(intercept_flag) {
        eigenMatrix xt_x;
        eigenVector xt_y, b_coeff;
        xt_x = x.transpose()*x;
        xt_y = x.transpose()*y;
        b_coeff = xt_x.ldlt().solve(xt_y);
        // update h2
        gcov=b_coeff(0)*(double)ttl_mk_num/n_gcov.mean();
        intercept_gcov = b_coeff(1);
    } else {
        double xt_x = 0.0, xt_y = 0.0, b_coeff = 0.0;
        xt_x = x.col(0).transpose()*x.col(0);
        xt_y = x.col(0).transpose()*y;
        b_coeff = xt_y/xt_x;
        // update h2
        gcov = b_coeff*(double)ttl_mk_num/n_gcov.mean();
    }
    
    // weights
    if(x_flag) {
        w = update_weights_gcov(intercept1, hsq1, intercept2, hsq2, intercept_gcov, gcov, ttl_mk_num, n_ld_snp, x.col(0), w_ld, n1,  n2, n_gcov);
    } else {
        w = update_weights_gcov(intercept1, hsq1, intercept2, hsq2, intercept_gcov, gcov, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n1,  n2, n_gcov);
    }
    return w;
}

vector<double> est_gcov_trait_1_step(eigenVector zscore, eigenVector n_gcov, eigenVector ref_ld, eigenVector w_ld, double intercept1, double hsq1, eigenVector n1, double intercept2, double hsq2, eigenVector n2, int n_ld_snp, int ttl_mk_num) {
    
    // Estimate prior of the weights
    eigenVector denominator = ref_ld.cwiseProduct(n_gcov);
    double gcov_prior = zscore.mean()*(double)ttl_mk_num/denominator.mean();
    double intercept_gcov_prior = 0.0;

    eigenVector wt_ttl = update_weights_gcov(intercept1, hsq1, intercept2, hsq2, intercept_gcov_prior, gcov_prior, ttl_mk_num, n_ld_snp, ref_ld, w_ld, n1, n2, n_gcov);
    eigenMatrix x(n_ld_snp, 2);
    x.col(1).setOnes(); x.col(0) = ref_ld.cwiseProduct(n_gcov)/n_gcov.mean();
    
    // Estimate intercept and slope by Iteratively reweighted least squares
    double intercept_gcov_posterior = 0.0, gcov_posterior = 0.0;

    wt_ttl = compute_irls_gcov(intercept_gcov_posterior, gcov_posterior, x, zscore, wt_ttl, ref_ld, w_ld, n_gcov, intercept1,  hsq1, n1, intercept2,  hsq2, n2, ttl_mk_num, n_ld_snp, true, false);
    wt_ttl = compute_irls_gcov(intercept_gcov_posterior, gcov_posterior, x, zscore, wt_ttl, ref_ld, w_ld, n_gcov, intercept1,  hsq1, n1, intercept2,  hsq2, n2, ttl_mk_num, n_ld_snp, true, false);
    wt_ttl = compute_irls_gcov(intercept_gcov_posterior, gcov_posterior, x, zscore, wt_ttl, ref_ld, w_ld, n_gcov, intercept1,  hsq1, n1, intercept2,  hsq2, n2, ttl_mk_num, n_ld_snp, true, false);
    
    vector<double> ldsc_est(2);
    ldsc_est[0]  = intercept_gcov_posterior; ldsc_est[1] = gcov_posterior/sqrt(hsq1*hsq2);
    return ldsc_est;
}

double transform_hsq_L(double P, double K, double hsq) {
    double t = StatFunc::qnorm(1.0 - K);
    double z = StatFunc::dnorm(t);
    double C = (K * (1 - K) / (z * z))*(K * (1 - K) / (P * (1 - P)));
    return (hsq * C);
}

vector<string> gcta::read_snp_ldsc(map<string,int> ldsc_snp_name_map, vector<string> snp_name, vector<int> snp_remain, int &ttl_mk_num, 
                                   string ref_ld_dirt, string w_ld_dirt, vector<double> &ref_ld_vec, vector<double> &w_ld_vec) {
    int i=0, nsnp = snp_remain.size();
    vector<string> ref_ld_snps, w_ld_snps;

    ref_ld_vec.clear(); w_ld_vec.clear();
    // Read the total number of markers
    ttl_mk_num = read_ld_marker(ref_ld_dirt);
    // Read the reference LD scores
    ref_ld_snps = read_ld_score(ref_ld_dirt, ldsc_snp_name_map, nsnp, ref_ld_vec);
    // Read the weighted LD scores
    w_ld_snps = read_ld_score(w_ld_dirt, ldsc_snp_name_map, nsnp, w_ld_vec);
    // SNPs in common
    map<string,int> w_ld_snp_map;
    vector<string> cm_ld_snps;
    int n_ref_ld_snps = ref_ld_snps.size(), n_w_ld_snps = w_ld_snps.size();
    for(i=0; i<n_w_ld_snps; i++)
        w_ld_snp_map.insert(pair<string, int>(w_ld_snps[i], i));
    for(i=0; i<n_ref_ld_snps; i++) {
        if(w_ld_snp_map.find(ref_ld_snps[i])!=w_ld_snp_map.end()) {
            cm_ld_snps.push_back(ref_ld_snps[i]);
        }
    }

    return cm_ld_snps;
}

void gcta::reorder_snp_effect(vector<int> snp_remain, eigenMatrix &bhat_z, eigenMatrix &bhat_n, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, 
                              vector<vector<bool>> &snp_flag, vector<vector<bool>> snp_val_flag, vector<int> &nsnp_cm_trait,
                              vector<string> cm_ld_snps, map<string,int> ldsc_snp_name_map,
                              eigenVector &ref_ld, eigenVector &w_ld, vector<double> ref_ld_vec, vector<double> w_ld_vec, int ntrait) {
    // Re-order the variables
    int i = 0, j = 0, n_cm_ld_snps = cm_ld_snps.size(), indxbuf = 0;
    map<string,int>::iterator iter;

    ref_ld.resize(n_cm_ld_snps); ref_ld.setZero(n_cm_ld_snps);
    w_ld.resize(n_cm_ld_snps); w_ld.setZero(n_cm_ld_snps);
    bhat_z.resize(n_cm_ld_snps, ntrait); bhat_z.setZero(n_cm_ld_snps, ntrait);
    bhat_n.resize(n_cm_ld_snps, ntrait); bhat_n.setZero(n_cm_ld_snps, ntrait);
    snp_flag.clear(); snp_flag.resize(ntrait);
    for(i=0; i<ntrait; i++) snp_flag[i].resize(n_cm_ld_snps);
    for(i=0; i<n_cm_ld_snps; i++) {
        iter = ldsc_snp_name_map.find(cm_ld_snps[i]);
        indxbuf = iter->second;
        ref_ld(i) = ref_ld_vec[indxbuf]; w_ld(i) = w_ld_vec[indxbuf];
        for(j=0; j<ntrait; j++) {
            snp_flag[j][i] = snp_val_flag[j][snp_remain[indxbuf]];
            if(!snp_flag[j][i]) continue;
            nsnp_cm_trait[j]++;
            bhat_z(i,j) = snp_b(snp_remain[indxbuf], j)/snp_se(snp_remain[indxbuf], j);
            bhat_n(i,j) = snp_n(snp_remain[indxbuf], j);
        }
    }
}

eigenMatrix gcta::ldsc_snp_h2(eigenMatrix bhat_z, eigenMatrix bhat_n, eigenVector ref_ld, eigenVector w_ld, vector<vector<bool>> snp_flag, vector<int> nsnp_cm_trait, int n_cm_ld_snps, int ttl_mk_num, vector<string> trait_name, int ntrait) {
    // Estimate SNP h2
    int i = 0, j = 0;
    vector<double> rst_ldsc(2);
    vector<vector<vector<double>>> ldsc_trait(3);
    eigenMatrix ldsc_var(ntrait, 2);
    int nsnp_ldsc_thresh = 5e5;

    for(i=0; i<ntrait; i++) {
        // Remove missing value
        if(nsnp_cm_trait[i] < nsnp_ldsc_thresh) 
            LOGGER.w(0, "Only " + to_string(nsnp_cm_trait[i]) + " are retained in the univariate LD score regression analysis for " 
                    + trait_name[i] + ". The estimate may not be accurate.");                    
        int k = 0, nsnp_buf = 0;
        eigenVector chi_val_buf(nsnp_cm_trait[i]), n_buf(nsnp_cm_trait[i]), ref_ld_buf(nsnp_cm_trait[i]), w_ld_buf(nsnp_cm_trait[i]);
        for(j = 0, k = 0; j < n_cm_ld_snps; j++) {
            if(!snp_flag[i][j]) continue;
            chi_val_buf(k) = bhat_z(j,i)*bhat_z(j,i);
            n_buf(k) = bhat_n(j,i);
            ref_ld_buf(k) = ref_ld(j);
            w_ld_buf(k) = w_ld(j); 
            k++;
        }
        // h2
        rst_ldsc = est_hsq_trait_2_steps(chi_val_buf, n_buf, ref_ld_buf, w_ld_buf, nsnp_cm_trait[i], ttl_mk_num);
        ldsc_var(i, 0) = rst_ldsc[0];
        if(rst_ldsc[1]>0) {
            ldsc_var(i,1) = rst_ldsc[1];
            LOGGER.i(0, trait_name[i] + ": " + to_string(rst_ldsc[0]) + " " + to_string(rst_ldsc[1]));  
        } else {
            LOGGER.e(0, "negative SNP heritability estimate for " + trait_name[i] + ". Exiting ...");
        }
    }  
    return ldsc_var;
}

eigenMatrix gcta::ldsc_snp_rg(eigenMatrix ldsc_var_h2, eigenMatrix bhat_z, eigenMatrix bhat_n, eigenVector ref_ld, eigenVector w_ld, vector<vector<bool>> snp_flag, vector<int> trait_indx1, vector<int> trait_indx2, int n_cm_ld_snps, int ttl_mk_num, vector<string> trait_name) {
    // Estimate SNP rg
    int i = 0, nproc = trait_indx1.size();
    eigenMatrix ldsc_var_rg(nproc, 2);
    eigenVector gcov_n1n2(n_cm_ld_snps);
    vector<double> rst_ldsc(2);
    int nsnp_ldsc_thresh = 5e5;
    for( i = 0; i < nproc; i++ ) {
        int k = 0, n_cm_snps_buf = 0;
        vector<bool> snp_pair_flag(n_cm_ld_snps);
        for(k=0; k<n_cm_ld_snps; k++) {
            snp_pair_flag[k] = snp_flag[trait_indx1[i]][k] && snp_flag[trait_indx2[i]][k];
            n_cm_snps_buf += snp_pair_flag[k];
        }
        if(n_cm_snps_buf < nsnp_ldsc_thresh) 
            LOGGER.w(0, "Only " + to_string(n_cm_snps_buf) + " are retained in the bivariate LD score regression analysis for " 
                    + trait_name[trait_indx1[i]] + " and " + trait_name[trait_indx2[i]] + ". The estimate may not be accurate.");
        eigenVector gcov_z1z2(n_cm_snps_buf), gcov_n1n2(n_cm_snps_buf), ref_ld_buf(n_cm_snps_buf), w_ld_buf(n_cm_snps_buf), n_buf_i(n_cm_snps_buf), n_buf_j(n_cm_snps_buf);
        int t = 0;
        for(k=0, t=0; k<n_cm_ld_snps; k++) {
            if(!snp_pair_flag[k]) continue;
            gcov_z1z2(t) = bhat_z(k,trait_indx1[i])*bhat_z(k,trait_indx2[i]);
            gcov_n1n2(t) = sqrt(bhat_n(k,trait_indx1[i])*bhat_n(k,trait_indx2[i]));
            ref_ld_buf(t) = ref_ld(k);
            w_ld_buf(t) = w_ld(k);
            n_buf_i(t) = bhat_n(k,trait_indx1[i]);
            n_buf_j(t) = bhat_n(k,trait_indx2[i]);
            t++;
        }
        rst_ldsc = est_gcov_trait_1_step(gcov_z1z2, gcov_n1n2, ref_ld_buf, w_ld_buf, 
                                        ldsc_var_h2(trait_indx1[i],0), ldsc_var_h2(trait_indx1[i],1), n_buf_i, 
                                        ldsc_var_h2(trait_indx2[i],0), ldsc_var_h2(trait_indx2[i],1), n_buf_j, n_cm_snps_buf, ttl_mk_num);
        ldsc_var_rg(i,0) = rst_ldsc[0];
        ldsc_var_rg(i,1) = rst_ldsc[1];
    } 
    return ldsc_var_rg;
}

bool gcta::mtcojo_ldsc(vector<vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, int ntrait, 
                vector<string> snp_name, vector<int> snp_remain, string ref_ld_dirt, string w_ld_dirt, vector<string> trait_name,
                eigenMatrix &ldsc_intercept, eigenMatrix &ldsc_slope) {

    int i = 0, j = 0, nsnp = snp_remain.size();
    int ttl_mk_num = 0;
    vector<int> nsnp_cm_trait(ntrait);    
    vector<double> ref_ld_vec, w_ld_vec;
    vector<string> cm_ld_snps;
    vector<vector<bool>> snp_flag;
    map<string,int> ldsc_snp_name_map;
    eigenVector ref_ld, w_ld;
    eigenMatrix bhat_z, bhat_n;
    bool proc_flag = true;

    // Initialize the parameters
    for(i=0; i<nsnp; i++) {
        ldsc_snp_name_map.insert(pair<string,int>(snp_name[snp_remain[i]], i));
    }
    for(i=0; i<ntrait; i++) nsnp_cm_trait[i] = 0;
    // Read the LD scores
    cm_ld_snps = read_snp_ldsc(ldsc_snp_name_map, snp_name, snp_remain, ttl_mk_num, ref_ld_dirt, w_ld_dirt, ref_ld_vec, w_ld_vec);
    // Re-order the SNP effects and LD scores
    reorder_snp_effect(snp_remain, bhat_z, bhat_n, snp_b, snp_se, snp_n, snp_flag, snp_val_flag, nsnp_cm_trait,
                       cm_ld_snps, ldsc_snp_name_map, ref_ld, w_ld, ref_ld_vec, w_ld_vec, ntrait);
    for(i=0; i<ntrait; i++) {
        if(nsnp_cm_trait[i]==0)
            LOGGER.e(0, "\nno SNP in common between the summary data and the LD score files. Please check.");
    }
    // Check the p-value
    for(i=0; i<ntrait; i++) {
        double chi_val_buf = bhat_z.col(i).squaredNorm()/nsnp_cm_trait[i];
        if(chi_val_buf > 10) {
            proc_flag = false; return(proc_flag);
        }
    }
    // LDSC analysis
    int n_cm_ld_snps = cm_ld_snps.size();
    ldsc_intercept.resize(ntrait, ntrait);
    ldsc_slope.resize(ntrait, ntrait);
    // Estimate SNP h2
    LOGGER.i(0, "\nUnivariate LD score regression analysis to estimate SNP-based heritability ..."); 
    eigenMatrix ldsc_var_h2;
    ldsc_var_h2 = ldsc_snp_h2(bhat_z, bhat_n, ref_ld, w_ld, snp_flag, nsnp_cm_trait, n_cm_ld_snps, ttl_mk_num, trait_name, ntrait);
    for( i = 0; i < ntrait; i++ ) {
        ldsc_intercept(i,i) = ldsc_var_h2(i,0);
        ldsc_slope(i,i) = ldsc_var_h2(i,1);
    }
    // Estimate SNP rg
    LOGGER.i(0, "Bivariate LD score regression analysis to estimate genetic correlation between each pair of traits ...");
    int k = 0, nproc = ntrait*(ntrait-1)/2;
    eigenMatrix ldsc_var_rg;
    vector<int> trait_indx1(nproc), trait_indx2(nproc);
    for( i=0, k=0; i<(ntrait-1); i++ ) {
        for( j=i+1; j<ntrait; j++, k++) {
            trait_indx1[k] = i; trait_indx2[k] = j;
        }
    }
    ldsc_var_rg = ldsc_snp_rg(ldsc_var_h2, bhat_z, bhat_n, ref_ld, w_ld, snp_flag, trait_indx1, trait_indx2, n_cm_ld_snps, ttl_mk_num, trait_name);
    for( i=0; i<nproc; i++) {
        ldsc_intercept(trait_indx1[i], trait_indx2[i]) = ldsc_intercept(trait_indx2[i], trait_indx1[i]) = ldsc_var_rg(i,0);
        ldsc_slope(trait_indx1[i], trait_indx2[i]) = ldsc_slope(trait_indx2[i], trait_indx1[i]) = ldsc_var_rg(i,1);
    }
    // Print the intercept of rg
    stringstream ss;
    LOGGER.i(0, "Intercept:");
    for(i=0; i<ntrait; i++) {
        ss.str("");
        for(j=0; j<ntrait; j++)
            ss << ldsc_intercept(i,j) << " ";
        LOGGER.i(0, ss.str());
    }
    LOGGER.i(0, "rg:");
    for(i=0; i<ntrait; i++) {
        ss.str("");
        for(j=0; j<ntrait; j++) 
            ss << ldsc_slope(i,j) << " ";
        LOGGER.i(0, ss.str());
    }
    LOGGER.i(0, "The LD score regression analyses completed.");

    return(proc_flag);
}

eigenMatrix mtcojo_cond_single_covar(eigenVector bzy, eigenVector bzy_se,  eigenMatrix bzx, eigenMatrix bzx_se, double bxy, eigenMatrix ldsc_intercept, eigenMatrix ldsc_slope, int nsnp) {
    int i=0;
    double var_bzx_buf = 0.0, cov_bzx_bzy = 0.0;
    eigenMatrix mtcojo_est(nsnp, 3);

    for(i=0; i<nsnp; i++) {
        mtcojo_est(i,0) = bzy(i) - bzx(i,0)*bxy;
        var_bzx_buf = bxy*bxy*bzx_se(i,0)*bzx_se(i,0);
        cov_bzx_bzy = bxy*ldsc_intercept(0,1)*bzx_se(i,0)*bzy_se(i);
        mtcojo_est(i,1) = sqrt(bzy_se(i)*bzy_se(i) + var_bzx_buf - 2*cov_bzx_bzy);
        mtcojo_est(i,2) = StatFunc::pchisq(mtcojo_est(i,0)*mtcojo_est(i,0)/mtcojo_est(i,1)/mtcojo_est(i,1), 1);
     }
    
    return mtcojo_est;
}

eigenMatrix mtcojo_cond_multiple_covars(eigenVector bzy, eigenVector bzy_se,  eigenMatrix bzx, eigenMatrix bzx_se, eigenVector bxy, eigenMatrix ldsc_intercept, eigenMatrix ldsc_slope, eigenVector vp_trait, int nsnp, int ncovar) {
    int i=0, j=0;

    double var_bzx_buf = 0.0, cov_bzx_bzy = 0.0;
    eigenVector bjxy(ncovar);
    MatrixXd d_mat = MatrixXd::Zero(ncovar,ncovar), r_mat = MatrixXd::Identity(ncovar, ncovar);
    
    for(i=0; i<ncovar; i++) {
        d_mat(i,i) = sqrt(ldsc_slope(i+1, i+1)*vp_trait(i+1));
        for(j=1; j<ncovar; j++) {
            if(i==j) continue;
            r_mat(i,j) = r_mat(j,i) = ldsc_slope(i+1, j+1);
        }
    }
    
    bjxy = d_mat.ldlt().solve(r_mat.ldlt().solve(d_mat*bxy));
    
    double d_buf = 0.0;
    eigenVector bjxy_buf(ncovar), se_bzx_bzy(ncovar);
    eigenMatrix bzx_se_buf=MatrixXd::Zero(ncovar,ncovar), bzx_intercept=MatrixXd::Identity(ncovar,ncovar);
    for(i=0; i<(ncovar-1); i++) {
        for(j=i+1; j<ncovar; j++) {
            bzx_intercept(i,j) = bzx_intercept(j,i) = ldsc_intercept(i+1, j+1);
        }
    }

    eigenMatrix mtcojo_est(nsnp, 3);
    for(i=0; i<nsnp; i++) {
        d_buf = 0.0; cov_bzx_bzy = 0.0;
        for(j=0; j<ncovar; j++) {
            d_buf += bzx(i,j) * bjxy(j);
            bzx_se_buf(j,j) = bzx_se(i,j);
            cov_bzx_bzy += bzx_se(i,j)*bzy_se(i)*bjxy(j)*ldsc_intercept(0, j+1);
        }
        var_bzx_buf = (bjxy.transpose()*bzx_se_buf*bzx_intercept*bzx_se_buf*bjxy)(0,0);
        mtcojo_est(i,0) = bzy(i) - d_buf;
        mtcojo_est(i,1) = sqrt(bzy_se(i)*bzy_se(i) + var_bzx_buf - 2*cov_bzx_bzy);
        mtcojo_est(i,2) = StatFunc::pchisq(mtcojo_est(i,0)*mtcojo_est(i,0)/mtcojo_est(i,1)/mtcojo_est(i,1), 1);
    }
    
    return mtcojo_est;
}

void mtcojo_cond_output(string output_file, vector<string> snp_name, vector<int> snp_remain, vector<string> snp_a1, vector<string> snp_a2, eigenVector snp_freq, eigenVector snp_b, eigenVector snp_se, eigenVector snp_pval, eigenMatrix mtcojo_est, eigenVector snp_n, int nsnp) {
    
    int i=0, snpindx=0;
    ofstream ofile(output_file.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + output_file + "] to write.");
    
    ofile << "SNP\tA1\tA2\tfreq\tb\tse\tp\tN\tbC\tbC_se\tbC_pval" <<endl;
    for (i = 0; i < nsnp; i++) {
        snpindx = snp_remain[i];
        ofile << snp_name[snpindx] << "\t" <<snp_a1[snpindx] << "\t" << snp_a2[snpindx] << "\t" << snp_freq[snpindx] 
              << "\t" << snp_b[snpindx] << "\t" << snp_se[snpindx] << "\t"  << snp_pval[snpindx] << "\t" 
              << snp_n[snpindx] << "\t" << mtcojo_est(i,0) << "\t" << mtcojo_est(i,1) <<"\t"<< mtcojo_est(i,2) << "\t" <<endl;
    }
    ofile.close();
}

void gcta::mtcojo(string mtcojo_bxy_file, string ref_ld_dirt, string w_ld_dirt, double freq_thresh, double gwas_thresh, int clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, int gsmr_beta_version) {
    if(gsmr_beta_version) { _n_gsmr_rst_item = 5; _gsmr_beta_version = 1; }
    else { _n_gsmr_rst_item = 4; _gsmr_beta_version = 0; }

    int i=0, j=0, nsnp=_meta_remain_snp.size(), ncovar=_covar_pheno_name.size();
    vector<bool> snp_pair_flag(nsnp);
    vector<string> afsnps;
    eigenVector bxy_est(ncovar);
    bxy_est.setConstant(-999999.0);

    // Read the bxy
    if (!mtcojo_bxy_file.empty()) bxy_est = read_external_bxy(mtcojo_bxy_file, _covar_pheno_name);

    // Calculate allele frequency
    if (_mu.empty()) calcu_mu();
    LOGGER.i(0, "Checking the difference in allele frequency between the GWAS summary datasets and the LD reference sample...");
    afsnps = remove_freq_diff_snps(_meta_snp_name, _meta_remain_snp, _snp_name_map, _mu, _meta_snp_freq, _snp_val_flag, ncovar+1, freq_thresh, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(afsnps, _meta_snp_name_map, _meta_remain_snp);
    }
    // Remove monomorphic SNPs
    afsnps = remove_mono_snps(_snp_name_map, _mu, _out);
    // Update SNPs set
    if( afsnps.size()>0 ) {
        update_id_map_rm(afsnps, _snp_name_map, _include);
        update_mtcojo_snp_rm(afsnps, _meta_snp_name_map, _meta_remain_snp);
    }

    // Only keep the AF for the target trait
    int nsnp_freq=_meta_snp_freq.size();
    _meta_snp_freq.conservativeResize(nsnp_freq,1);

    // LDSC analysis to estimate h2 and rg
    eigenMatrix ldsc_intercept, ldsc_slope;
    vector<string> trait_name;
    trait_name.push_back(_target_pheno_name);
    trait_name.insert(trait_name.end(), _covar_pheno_name.begin(), _covar_pheno_name.end());
   
    bool proc_flag = true;
    proc_flag = mtcojo_ldsc(_snp_val_flag, _meta_snp_b, _meta_snp_se, _meta_snp_n_o, ncovar+1, _meta_snp_name, _meta_remain_snp, ref_ld_dirt, w_ld_dirt, trait_name, ldsc_intercept, ldsc_slope);
    if(!proc_flag) {
        if(ncovar==1) {
            ldsc_intercept.setZero(ncovar+1, ncovar+1);
            ldsc_slope.setZero(ncovar+1, ncovar+1);
            LOGGER.e(0, "not enough SNPs to perform the univariate LD score regression analysis. The mtCOJO analysis will be conducted assuming no sample overlap between the GWAS data for the target and covariate traits.");
        } else {
            LOGGER.e(0, "not enough SNPs to perform the univariate and bivariate LD score regression analyses. The mtCOJO analysis requires SNP-based heritability from univariate LD score regression analysis and genetic correlation from bivariate LD score regression analysis.");
        }
    } 

    // GSMR analysis
    vector<double> gsmr_rst(_n_gsmr_rst_item);
    std::stringstream ss_pleio;
    bool pleio_flag = false;
    for(i=1; i<=ncovar; i++) {
        if(fabs(bxy_est(i-1))<100) continue;
        string pleio_snps = "";
        vector<string> snp_instru;
        string err_msg;
        LOGGER.i(0, "\nGSMR analysis for covariate #" + to_string(i) + " (" + trait_name[i] + ")" + " ...");
        for(j=0; j<nsnp; j++) snp_pair_flag[j] = _snp_val_flag[0][_meta_remain_snp[j]] && _snp_val_flag[i][_meta_remain_snp[j]];
        gsmr_rst =  gsmr_meta(snp_instru, _meta_snp_b.col(i), _meta_snp_se.col(i), _meta_snp_pval.col(i), _meta_snp_b.col(0), _meta_snp_se.col(0), _meta_snp_pval.col(0), 0, snp_pair_flag, gwas_thresh, clump_wind_size, clump_r2_thresh, std_heidi_thresh, global_heidi_thresh, ld_fdr_thresh, nsnp_gsmr, pleio_snps, err_msg);
        if(std::isnan(gsmr_rst[3])) {
            LOGGER.w(0, err_msg);
            LOGGER.i(0, "bxy is estimated from rg.");
            bxy_est(i-1) = ldsc_slope(i,0)*sqrt(ldsc_slope(0,0)/ldsc_slope(i,i));
            LOGGER.i(0, "bxy " + to_string(gsmr_rst[0]) + " " + to_string(gsmr_rst[1]));
        }
        // Saving pleiotropic SNPs
        if(pleio_snps.size() > 0) {
            ss_pleio << _covar_pheno_name[i-1] << " " << _target_pheno_name << " " << pleio_snps << endl;
            pleio_flag = true;
        }
        bxy_est(i-1) = gsmr_rst[0];
        LOGGER.i(0, "bxy " + to_string(gsmr_rst[0]) + " " + to_string(gsmr_rst[1]));
        LOGGER.i(0, "GSMR analysis for covariate #" + to_string(i) + " (" + trait_name[i] + ") completed.");
    }
    if(pleio_flag) {
        string pleio_snpfile = _out + ".pleio_snps";	
        ofstream o_pleio_snp(pleio_snpfile.c_str());	
        if(!o_pleio_snp) LOGGER.e(0, "cannot open file [" + pleio_snpfile + "] to write pleiotropic SNPs.");
        o_pleio_snp << ss_pleio.str();            
        o_pleio_snp.close();	
        LOGGER.i(0,  "The pleiotropic SNPs filtered by HEIDI-outlier analysis have been saved in [" + pleio_snpfile + "].");
    }

    // mtcojo analysis
    // SNPs in common across the whole traits
    int n_buf = 0;
    vector<string> snp_buf;
    for(i=0; i<nsnp; i++) {
        n_buf = 0;
        for(j=0; j<ncovar+1; j++) n_buf += _snp_val_flag[j][_meta_remain_snp[i]];
        if(n_buf < ncovar+1) continue;
        snp_buf.push_back(_meta_snp_name[_meta_remain_snp[i]]);
    }
    update_id_map_kp(snp_buf, _meta_snp_name_map, _meta_remain_snp);

    nsnp = _meta_remain_snp.size();
    eigenVector snp_bzy(nsnp), snp_bzy_se(nsnp);
    eigenMatrix mtcojo_est, snp_bzx(nsnp, ncovar), snp_bzx_se(nsnp, ncovar);
    for(i=0; i<nsnp; i++) {
        snp_bzy(i) = _meta_snp_b(_meta_remain_snp[i], 0);
        snp_bzy_se(i) = _meta_snp_se(_meta_remain_snp[i], 0);
        for(j = 0; j<ncovar; j++) {
            snp_bzx(i,j) = _meta_snp_b(_meta_remain_snp[i], j+1);
            snp_bzx_se(i,j) = _meta_snp_se(_meta_remain_snp[i], j+1);
        }
    }

    LOGGER.i(0, "\nmtCOJO analysis ...");
    LOGGER.i(0, "There are " + to_string(nsnp) + " SNPs in common between the target trait and all the covariate trait(s).");
    if(ncovar==1) {
        mtcojo_est = mtcojo_cond_single_covar(snp_bzy, snp_bzy_se, snp_bzx, snp_bzx_se, bxy_est(0), ldsc_intercept, ldsc_slope, nsnp);
    } else {
        mtcojo_est = mtcojo_cond_multiple_covars(snp_bzy, snp_bzy_se, snp_bzx, snp_bzx_se, bxy_est, ldsc_intercept, ldsc_slope, _meta_vp_trait, nsnp, ncovar);
    }
    
    // Output
    string output_filename = _out + ".mtcojo.cma";
    LOGGER.i(0, "Saving the mtCOJO analysis results of " + to_string(nsnp) + " remaining SNPs to [" + output_filename + "] ...");
    LOGGER.i(0, "mtCOJO analysis completed.");
    mtcojo_cond_output(output_filename, _meta_snp_name, _meta_remain_snp, _meta_snp_a1, _meta_snp_a2, _meta_snp_freq, _meta_snp_b.col(0), _meta_snp_se.col(0), _meta_snp_pval.col(0), mtcojo_est, _meta_snp_n_o.col(0), nsnp);
}

