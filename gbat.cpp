/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementation of gene-based association test (GBAT) in GCTA
 *
 * 2013 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::gbat_read_snpAssoc(string snpAssoc_file, vector<string> &snp_name, vector<int> &snp_chr, vector<int> &snp_bp, vector<double> &snp_pval)
{
    ifstream in_snpAssoc(snpAssoc_file.c_str());
    if (!in_snpAssoc) throw ("Error: can not open the file [" + snpAssoc_file + "] to read.");
    cout << "\nReading SNP association results from [" + snpAssoc_file + "]." << endl;
    string str_buf;
    vector<string> vs_buf;
    map<string, int>::iterator iter;
    cout << "Reading association p-values from [" << snpAssoc_file << "]." << endl;
    while (getline(in_snpAssoc, str_buf)) {
        if (StrFunc::split_string(str_buf, vs_buf) != 2) throw ("Error: in line \"" + str_buf + "\".");
        iter = _snp_name_map.find(vs_buf[0]);
        if (iter == _snp_name_map.end()) continue;
        snp_name.push_back(vs_buf[0]);
        snp_pval.push_back(atof(vs_buf[1].c_str()));
    }
    in_snpAssoc.close();
    cout << "Association p-values of " << snp_name.size() << " SNPs have been included." << endl;

    update_id_map_kp(snp_name, _snp_name_map, _include);
    vector<string> snp_name_buf(snp_name);
    vector<double> snp_pval_buf(snp_pval);
    snp_name.clear();
    snp_pval.clear();
    snp_name.resize(_include.size());
    snp_pval.resize(_include.size());
    int i = 0;
    map<string, int> snp_name_buf_map;
    for (i = 0; i < snp_name_buf.size(); i++) snp_name_buf_map.insert(pair<string,int>(snp_name_buf[i], i));
    #pragma omp parallel for
    for (i = 0; i < _include.size(); i++) {
        map<string, int>::iterator iter = snp_name_buf_map.find(_snp_name[_include[i]]);
        snp_name[i] = snp_name_buf[iter->second];
        snp_pval[i] = snp_pval_buf[iter->second];
    }
    snp_chr.resize(_include.size());
    snp_bp.resize(_include.size());
    #pragma omp parallel for
    for (i = 0; i < _include.size(); i++) {
        snp_chr[i] = _chr[_include[i]];
        snp_bp[i] = _bp[_include[i]];
    }
    if (_include.size() < 1) throw ("Error: no SNP is included in the analysis.");
    else if (_chr[_include[0]] < 1) throw ("Error: chromosome information is missing.");
    else if (_bp[_include[0]] < 1) throw ("Error: bp information is missing.");
}

void gcta::gbat_read_geneAnno(string gAnno_file, vector<string> &gene_name, vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2) {
    ifstream in_gAnno(gAnno_file.c_str());
    if (!in_gAnno) throw ("Error: can not open the file [" + gAnno_file + "] to read.");
    cout << "Reading physical positions of the genes from [" + gAnno_file + "]." << endl;
    string str_buf;
    vector<string> vs_buf;
    while (getline(in_gAnno, str_buf)) {
        if (StrFunc::split_string(str_buf, vs_buf) != 4) throw ("Error: in line \"" + str_buf + "\".");
        gene_chr.push_back(atoi(vs_buf[0].c_str()));
        gene_bp1.push_back(atoi(vs_buf[1].c_str()));
        gene_bp2.push_back(atoi(vs_buf[2].c_str()));
        gene_name.push_back(vs_buf[3]);
    }
    in_gAnno.close();
    cout << "Physical positions of " << gene_name.size() << " genes have been include." << endl;
}

void gcta::gbat_calcu_ld(MatrixXf &X, eigenVector &sumsq_x, int snp1_indx, int snp2_indx, MatrixXf &C)
{
    int i = 0, j = 0;
    int size = snp2_indx - snp1_indx + 1;

    C.resize(0, 0);
    if (size == 1) {
        C.resize(1, 1);
        C(1, 1) = 1.0;
        return;
    }

    //MatrixXf X_sub = X.block(0,snp1_indx,_keep.size(),size);
    C = X.block(0,snp1_indx,_keep.size(),size).transpose() * X.block(0,snp1_indx,_keep.size(),size);
    #pragma omp parallel for private(j)
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            double d_buf = sqrt(sumsq_x[snp1_indx + i] * sumsq_x[snp1_indx + j]);
            if(d_buf>0.0) C(i, j) /= d_buf;
            else C(i, j) = 0.0;
        }
    }
}

void gcta::gbat(string sAssoc_file, string gAnno_file, int wind)
{
    int i = 0, j = 0;

    // read SNP association results
    vector<string> snp_name;
    vector<int> snp_chr, snp_bp;
    vector<double> snp_pval;
    gbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
    vector<double> snp_chisq(snp_pval.size());
    for (i = 0; i < snp_pval.size(); i++) snp_chisq[i] = StatFunc::qchisq(snp_pval[i], 1);

    // get start and end of chr
    int snp_num = snp_name.size();
    map<int, string> chr_begin_snp, chr_end_snp;
    chr_begin_snp.insert(pair<int, string>(snp_chr[0], snp_name[0]));
    for (i = 1; i < snp_num; i++) {
        if (snp_chr[i] != snp_chr[i - 1]) {
            chr_begin_snp.insert(pair<int, string>(snp_chr[i], snp_name[i]));
            chr_end_snp.insert(pair<int, string>(snp_chr[i - 1], snp_name[i - 1]));
        }
    }
    chr_end_snp.insert(pair<int, string>(snp_chr[snp_num - 1], snp_name[snp_num - 1]));
    
    // read gene list
    vector<string> gene_name;
    vector<int> gene_chr, gene_bp1, gene_bp2;
    gbat_read_geneAnno(gAnno_file, gene_name, gene_chr, gene_bp1, gene_bp2);

    // map genes to SNPs
    cout << "Mapping the physical positions of genes to SNP data (gene bounaries: " << wind / 1000 << "Kb away from UTRs) ..." << endl;

    int gene_num = gene_name.size();
    vector<string> gene2snp_1(gene_num), gene2snp_2(gene_num);
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chr_iter;
    vector<locus_bp> snp_vec;
    for (i = 0; i < snp_num; i++) snp_vec.push_back(locus_bp(snp_name[i], snp_chr[i], snp_bp[i]));
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp1[i] - wind));
        if (iter != snp_vec.end()) gene2snp_1[i] = iter->locus_name;
        else gene2snp_1[i] = "NA";
    }
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] == "NA") {
            gene2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp2[i] + wind));
        if (iter != snp_vec.end()){
            if (iter->bp ==  gene_bp2[i] + wind) gene2snp_2[i] = iter->locus_name;
            else {
                if(iter!=snp_vec.begin()){
                    iter--;
                    gene2snp_2[i] = iter->locus_name;
                }
                else gene2snp_2[i] = "NA";
            }
        }
        else {
            chr_iter = chr_end_snp.find(gene_chr[i]);
            if (chr_iter == chr_end_snp.end()) gene2snp_2[i] = "NA";
            else gene2snp_2[i] = chr_iter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] != "NA" && gene2snp_2[i] != "NA") mapped++;
    }
    if (mapped < 1) throw ("Error: no gene can be mapped to the SNP data. Please check the input data regarding chr and bp.");
    else cout << mapped << " genes have been mapped to SNP data." << endl;

    // recoding genotype
    MatrixXf X;
    eigenVector sumsq_x(_include.size());
    make_XMat(X);
    #pragma omp parallel for private(j)
    for(i = 0; i < _keep.size(); i++){
        for(j = 0; j < _include.size(); j++){
            if(X(i,j) < 1e5) X(i,j) -= _mu[_include[j]];
            else X(i,j) = 0.0;
        }
    }
    for (i = 0; i < _include.size(); i++) sumsq_x[i] = X.col(i).dot(X.col(i));

    // run gene-based test
    cout << "\nRunning gene-based association test (GBAT)..." << endl;
    vector<double> gene_pval(gene_num), chisq_o(gene_num);
    vector<int> snp_num_in_gene(gene_num);
    map<string, int>::iterator iter1, iter2;
    map<string, int> snp_name_map;
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(pair<string,int>(snp_name[i], i));
    for (i = 0; i < gene_num; i++) {
        iter1 = snp_name_map.find(gene2snp_1[i]);
        iter2 = snp_name_map.find(gene2snp_2[i]);
        bool skip = false;
        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) skip = true;
        snp_num_in_gene[i] = iter2->second - iter1->second + 1;
        if(!skip && snp_num_in_gene[i] > 10000){
            cout<<"Warning: Too many SNPs in the gene region ["<<gene_name[i]<<"]. Maximum limit is 10000. This gene is ignored in the analysis."<<endl;
            skip = true;  
        } 
        if(skip){
            gene_pval[i] = 2.0;
            snp_num_in_gene[i] = 0;
            continue;
        }
        chisq_o[i] = 0;
        for (j = iter1->second; j <= iter2->second; j++) chisq_o[i] += snp_chisq[j];
        MatrixXf C;
        gbat_calcu_ld(X, sumsq_x, iter1->second, iter2->second, C); // iter2->second-1 because iter2 is one step further
        if(snp_num_in_gene[i] == 1) gene_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            SelfAdjointEigenSolver<MatrixXf> saes(C);
            gene_pval[i] = StatFunc::pchisqsum(chisq_o[i], saes.eigenvalues().cast<double>());
        }

        if((i + 1) % 100 == 0 || (i + 1) == gene_num) cout << i + 1 << " of " << gene_num << " genes.\r";
    }

    string filename = _out + ".gbat";
    cout << "\nSaving the results of the gene-based association analysese to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Gene\tChr\tStart\tEnd\tNo.SNPs\tSNP_start\tSNP_end\tChisq(Obs)\tPvalue" << endl;
    for (i = 0; i < gene_num; i++) {
        if(gene_pval[i]>1.5) continue;
        ofile << gene_name[i] << "\t" << gene_chr[i] << "\t" << gene_bp1[i] << "\t" << gene_bp2[i] << "\t";
        ofile << snp_num_in_gene[i] << "\t" << gene2snp_1[i] << "\t" << gene2snp_2[i] << "\t" << chisq_o[i] << "\t" << gene_pval[i] << endl;
        //else ofile << "0\tNA\tNA\tNA\tNA" << endl;
    }
    ofile.close();
}

void gcta::sbat(string sAssoc_file, string snpset_file)
{
    int i = 0, j = 0;

    // read SNP association results
    vector<string> snp_name;
    vector<int> snp_chr, snp_bp;
    vector<double> snp_pval;
    gbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
    vector<double> snp_chisq(snp_pval.size());
    for (i = 0; i < snp_pval.size(); i++) snp_chisq[i] = StatFunc::qchisq(snp_pval[i], 1);

    // get start and end of chr
    int snp_num = snp_name.size();
    map<int, string> chr_begin_snp, chr_end_snp;
    chr_begin_snp.insert(pair<int, string>(snp_chr[0], snp_name[0]));
    for (i = 1; i < snp_num; i++) {
        if (snp_chr[i] != snp_chr[i - 1]) {
            chr_begin_snp.insert(pair<int, string>(snp_chr[i], snp_name[i]));
            chr_end_snp.insert(pair<int, string>(snp_chr[i - 1], snp_name[i - 1]));
        }
    }
    chr_end_snp.insert(pair<int, string>(snp_chr[snp_num - 1], snp_name[snp_num - 1]));
    
    // read gene list
    vector<string> gene_name;
    vector<int> gene_chr, gene_bp1, gene_bp2;
    gbat_read_geneAnno(gAnno_file, gene_name, gene_chr, gene_bp1, gene_bp2);

    // map genes to SNPs
    cout << "Mapping the physical positions of genes to SNP data (gene bounaries: " << wind / 1000 << "Kb away from UTRs) ..." << endl;

    int gene_num = gene_name.size();
    vector<string> gene2snp_1(gene_num), gene2snp_2(gene_num);
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chr_iter;
    vector<locus_bp> snp_vec;
    for (i = 0; i < snp_num; i++) snp_vec.push_back(locus_bp(snp_name[i], snp_chr[i], snp_bp[i]));
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp1[i] - wind));
        if (iter != snp_vec.end()) gene2snp_1[i] = iter->locus_name;
        else gene2snp_1[i] = "NA";
    }
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] == "NA") {
            gene2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp2[i] + wind));
        if (iter != snp_vec.end()){
            if (iter->bp ==  gene_bp2[i] + wind) gene2snp_2[i] = iter->locus_name;
            else {
                if(iter!=snp_vec.begin()){
                    iter--;
                    gene2snp_2[i] = iter->locus_name;
                }
                else gene2snp_2[i] = "NA";
            }
        }
        else {
            chr_iter = chr_end_snp.find(gene_chr[i]);
            if (chr_iter == chr_end_snp.end()) gene2snp_2[i] = "NA";
            else gene2snp_2[i] = chr_iter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] != "NA" && gene2snp_2[i] != "NA") mapped++;
    }
    if (mapped < 1) throw ("Error: no gene can be mapped to the SNP data. Please check the input data regarding chr and bp.");
    else cout << mapped << " genes have been mapped to SNP data." << endl;

    // recoding genotype
    MatrixXf X;
    eigenVector sumsq_x(_include.size());
    make_XMat(X);
    #pragma omp parallel for private(j)
    for(i = 0; i < _keep.size(); i++){
        for(j = 0; j < _include.size(); j++){
            if(X(i,j) < 1e5) X(i,j) -= _mu[_include[j]];
            else X(i,j) = 0.0;
        }
    }
    for (i = 0; i < _include.size(); i++) sumsq_x[i] = X.col(i).dot(X.col(i));

    // run gene-based test
    cout << "\nRunning gene-based association test (GBAT)..." << endl;
    vector<double> gene_pval(gene_num), chisq_o(gene_num);
    vector<int> snp_num_in_gene(gene_num);
    map<string, int>::iterator iter1, iter2;
    map<string, int> snp_name_map;
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(pair<string,int>(snp_name[i], i));
    for (i = 0; i < gene_num; i++) {
        iter1 = snp_name_map.find(gene2snp_1[i]);
        iter2 = snp_name_map.find(gene2snp_2[i]);
        bool skip = false;
        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) skip = true;
        snp_num_in_gene[i] = iter2->second - iter1->second + 1;
        if(!skip && snp_num_in_gene[i] > 10000){
            cout<<"Warning: Too many SNPs in the gene region ["<<gene_name[i]<<"]. Maximum limit is 10000. This gene is ignored in the analysis."<<endl;
            skip = true;  
        } 
        if(skip){
            gene_pval[i] = 2.0;
            snp_num_in_gene[i] = 0;
            continue;
        }
        chisq_o[i] = 0;
        for (j = iter1->second; j <= iter2->second; j++) chisq_o[i] += snp_chisq[j];
        MatrixXf C;
        gbat_calcu_ld(X, sumsq_x, iter1->second, iter2->second, C); // iter2->second-1 because iter2 is one step further
        if(snp_num_in_gene[i] == 1) gene_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            SelfAdjointEigenSolver<MatrixXf> saes(C);
            gene_pval[i] = StatFunc::pchisqsum(chisq_o[i], saes.eigenvalues().cast<double>());
        }

        if((i + 1) % 100 == 0 || (i + 1) == gene_num) cout << i + 1 << " of " << gene_num << " genes.\r";
    }

    string filename = _out + ".gbat";
    cout << "\nSaving the results of the gene-based association analysese to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Gene\tChr\tStart\tEnd\tNo.SNPs\tSNP_start\tSNP_end\tChisq(Obs)\tPvalue" << endl;
    for (i = 0; i < gene_num; i++) {
        if(gene_pval[i]>1.5) continue;
        ofile << gene_name[i] << "\t" << gene_chr[i] << "\t" << gene_bp1[i] << "\t" << gene_bp2[i] << "\t";
        ofile << snp_num_in_gene[i] << "\t" << gene2snp_1[i] << "\t" << gene2snp_2[i] << "\t" << chisq_o[i] << "\t" << gene_pval[i] << endl;
        //else ofile << "0\tNA\tNA\tNA\tNA" << endl;
    }
    ofile.close();
}

