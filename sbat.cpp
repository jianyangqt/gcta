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

void gcta::sbat_read_snpAssoc(string snpAssoc_file, vector<string> &snp_name, vector<int> &snp_chr, vector<int> &snp_bp, vector<double> &snp_pval)
{
    ifstream in_snpAssoc(snpAssoc_file.c_str());
    if (!in_snpAssoc) throw ("Error: can not open the file [" + snpAssoc_file + "] to read.");
    cout << "\nReading SNP association results from [" + snpAssoc_file + "]." << endl;
    string str_buf;
    vector<string> vs_buf;
    map<string, int>::iterator iter;
    while (getline(in_snpAssoc, str_buf)) {
        if (StrFunc::split_string(str_buf, vs_buf) != 2) throw ("Error: in line \"" + str_buf + "\".");
        iter = _snp_name_map.find(vs_buf[0]);
        if (iter == _snp_name_map.end()) continue;
        snp_name.push_back(vs_buf[0]);
        snp_pval.push_back(atof(vs_buf[1].c_str()));
    }
    in_snpAssoc.close();
    snp_name.erase(unique(snp_name.begin(), snp_name.end()), snp_name.end());
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

void gcta::sbat_read_geneAnno(string gAnno_file, vector<string> &gene_name, vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2) {
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

void gcta::sbat_gene(string sAssoc_file, string gAnno_file, int wind, bool reduce_cor, bool write_snpset)
{
    int i = 0, j = 0;
    int snp_count;

    // read SNP association results
    vector<string> snp_name;
    vector<int> snp_chr, snp_bp;
    vector<double> snp_pval;
    sbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
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
    sbat_read_geneAnno(gAnno_file, gene_name, gene_chr, gene_bp1, gene_bp2);

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

    // run gene-based test
    if (_mu.empty()) calcu_mu();
    cout << "\nRunning set-based association test (SBAT) for genes ..." << endl;
    vector<double> gene_pval(gene_num), chisq_o(gene_num), min_snp_pval(gene_num);
    vector<string> min_snp_name(gene_num);
    vector<int> snp_num_in_gene(gene_num);
    map<string, int>::iterator iter1, iter2;
    map<string, int> snp_name_map;
    string rgoodsnpfile = _out + ".gene.snpset";
    ofstream rogoodsnp;
    if (write_snpset) rogoodsnp.open(rgoodsnpfile.c_str());
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
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = iter1->second; j <= iter2->second; j++) {
            if (min_snp_pval[i] > snp_pval[j]) { 
                min_snp_pval[i] = snp_pval[j];
                min_snp_name[i] = snp_name[j];
            }
            chisq_o[i] += snp_chisq[j];
        }
        if(snp_num_in_gene[i] == 1) gene_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            vector<int> snp_indx;
            for (j = iter1->second; j <= iter2->second; j++) snp_indx.push_back(j);            
            snp_count=snp_num_in_gene[i];
            VectorXd eigenval;
            vector<int> sub_indx;
            sbat_calcu_lambda(snp_indx, eigenval, snp_count, reduce_cor, sub_indx);
            //recalculate chisq value from low correlation snp subset
            if (reduce_cor) {
                chisq_o[i] = 0;
                for (j = 0; j < sub_indx.size(); j++) chisq_o[i] += snp_chisq[snp_indx[sub_indx[j]]];
            } 
            gene_pval[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);
            snp_num_in_gene[i] = snp_count;

            if (write_snpset) {
                rogoodsnp << gene_name[i] << endl;
                for (int k = 0; k < sub_indx.size(); k++) rogoodsnp << snp_name[(iter1->second)+sub_indx[k]] << endl;
                rogoodsnp << "END" << endl << endl;
            }

        }

        if((i + 1) % 100 == 0 || (i + 1) == gene_num) cout << i + 1 << " of " << gene_num << " genes.\r";
    }

    string filename = _out + ".gene.sbat";
    cout << "\nSaving the results of the SBAT analyses to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Gene\tChr\tStart\tEnd\tNo.SNPs\tSNP_start\tSNP_end\tChisq(Obs)\tPvalue\tTopSNP_Pvalue\tTopSNP" << endl;
    for (i = 0; i < gene_num; i++) {
        if(gene_pval[i]>1.5) continue;
        ofile << gene_name[i] << "\t" << gene_chr[i] << "\t" << gene_bp1[i] << "\t" << gene_bp2[i] << "\t";
        ofile << snp_num_in_gene[i] << "\t" << gene2snp_1[i] << "\t" << gene2snp_2[i] << "\t" << chisq_o[i];
        ofile << "\t" << gene_pval[i] << "\t" << min_snp_pval[i] << "\t" << min_snp_name[i] << endl;
        //else ofile << "0\tNA\tNA\tNA\tNA" << endl;
    }
    ofile.close();
    if (write_snpset) {
        cout << "Writing snpset ..." << endl;
        rogoodsnp.close();
    }
}

void gcta::sbat_read_snpset(string snpset_file, vector<string> &set_name, vector< vector<string> > &snpset)
{
    ifstream in_snpset(snpset_file.c_str());
    if (!in_snpset) throw ("Error: can not open the file [" + snpset_file + "] to read.");
    cout << "\nReading SNP set from [" + snpset_file + "]." << endl;
    string str_buf;
    vector<string> vs_buf, snpset_buf, snp_name;
    int i = 0;
    while (in_snpset>>str_buf) {
        if(str_buf!="END" && str_buf!="end") vs_buf.push_back(str_buf);
        else{
            if(vs_buf.empty()) continue;
            set_name.push_back(vs_buf[0]);
            for(i = 1; i < vs_buf.size(); i++){
                if (_snp_name_map.find(vs_buf[i]) != _snp_name_map.end()){
                    snpset_buf.push_back(vs_buf[i]);
                    snp_name.push_back(vs_buf[i]);
                }
            }
            vs_buf.clear();
            snpset.push_back(snpset_buf);
            snpset_buf.clear();
        }
    }
    in_snpset.close();
    snp_name.erase(unique(snp_name.begin(), snp_name.end()), snp_name.end());
    update_id_map_kp(snp_name, _snp_name_map, _include);
    cout << snp_name.size() << " SNPs in " << snpset.size() << " sets have been included." << endl;
}

void gcta::sbat(string sAssoc_file, string snpset_file, bool reduce_cor, bool write_snpset)
{
    int i = 0, j = 0;
    int snp_count;

    // read SNP set file
    vector<string> set_name;
    vector< vector<string> > snpset;
    sbat_read_snpset(snpset_file, set_name, snpset);
    int set_num = set_name.size();

    // read SNP association results
    vector<string> snp_name;
    vector<int> snp_chr, snp_bp;
    vector<double> snp_pval;
    sbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
    vector<double> snp_chisq(snp_pval.size());
    for (i = 0; i < snp_pval.size(); i++) snp_chisq[i] = StatFunc::qchisq(snp_pval[i], 1);   

    // run gene-based test
    if (_mu.empty()) calcu_mu();
    cout << "\nRunning set-based association test (SBAT)..." << endl;
    vector<double> set_pval(set_num), chisq_o(set_num), min_snp_pval(set_num);
    vector<string> min_snp_name(set_num);
    vector<int> snp_num_in_set(set_num);
    map<string, int>::iterator iter;
    map<string, int> snp_name_map;

    string rgoodsnpfile = _out + ".snpset";
    ofstream rogoodsnp;
    if (write_snpset) rogoodsnp.open(rgoodsnpfile.c_str());
 
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(pair<string,int>(snp_name[i], i));
    for (i = 0; i < set_num; i++) {
        bool skip = false;
        if(snpset[i].size() < 1) skip = true;
        vector<int> snp_indx;
        for(j = 0; j < snpset[i].size(); j++){
            iter = snp_name_map.find(snpset[i][j]);
            if(iter!=snp_name_map.end()) snp_indx.push_back(iter->second);
        }
        snp_num_in_set[i] = snp_indx.size();
        if(!skip && snp_num_in_set[i] > 20000){
            cout<<"Warning: Too many SNPs in the set ["<<set_name[i]<<"]. Maximum limit is 20000. This gene is ignored in the analysis."<<endl;
            skip = true;  
        } 
        if(skip){
            set_pval[i] = 2.0;
            snp_num_in_set[i] = 0;
            continue;
        }
        chisq_o[i] = 0;
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = 0; j < snp_indx.size(); j++) {
            if (min_snp_pval[i] > snp_pval[snp_indx[j]]) {
                min_snp_pval[i] = snp_pval[snp_indx[j]];
                min_snp_name[i] = snp_name[snp_indx[j]];
            }
            chisq_o[i] += snp_chisq[snp_indx[j]];
        }
        if(snp_num_in_set[i] == 1) set_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            snp_count=snp_num_in_set[i];
            VectorXd eigenval;
            vector<int> sub_indx;
            sbat_calcu_lambda(snp_indx, eigenval, snp_count, reduce_cor, sub_indx);

            //recalculate chisq value from low correlation snp subset
            if (reduce_cor) {
                chisq_o[i] = 0;
                for (j = 0; j < sub_indx.size(); j++) chisq_o[i] += snp_chisq[snp_indx[sub_indx[j]]];
            }
            set_pval[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);
            snp_num_in_set[i]=snp_count;

            if (write_snpset) {
                rogoodsnp << set_name[i] << endl;
                for (int k = 0; k < sub_indx.size(); k++) rogoodsnp << snp_name[snp_indx[sub_indx[k]]] << endl;
                rogoodsnp << "END" << endl << endl;
            }

        }

        if((i + 1) % 100 == 0 || (i + 1) == set_num) cout << i + 1 << " of " << set_num << " sets.\r";
    }

    string filename = _out + ".sbat";
    cout << "\nSaving the results of the SBAT analyses to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Set\tNo.SNPs\tChisq(Obs)\tPvalue\tTopSNP_Pvalue\tTopSNP" << endl;
    for (i = 0; i < set_num; i++) {
        if(set_pval[i]>1.5) continue;
        ofile << set_name[i] << "\t" << snp_num_in_set[i] << "\t" << chisq_o[i] << "\t";
        ofile << set_pval[i] << "\t" << min_snp_pval[i] << "\t" << min_snp_name[i] << endl;
    }
    ofile.close();
    if (write_snpset) {
        cout << "Writing snpset ..." << endl;
        rogoodsnp.close();
    }
}

void gcta::sbat_seg(string sAssoc_file, int seg_size, bool reduce_cor, bool write_snpset)
{
    int i = 0, j = 0;
    int snp_count;

       // read SNP association results
    vector<string> snp_name;
    vector<int> snp_chr, snp_bp;
    vector<double> snp_pval;
    sbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
    vector<double> snp_chisq(snp_pval.size());
    for (i = 0; i < snp_pval.size(); i++) snp_chisq[i] = StatFunc::qchisq(snp_pval[i], 1);   

    // run gene-based test
    if (_mu.empty()) calcu_mu();
    cout << "\nRunning set-based association test (SBAT) at genomic segments with a length of " << seg_size/1000 << "Kb ..." << endl;
    vector< vector<int> > snp_set_indx;
    vector<int> set_chr, set_start_bp, set_end_bp;
    get_sbat_seg_blk(seg_size, snp_set_indx, set_chr, set_start_bp, set_end_bp);
    int set_num = snp_set_indx.size();
    vector<double> set_pval(set_num), chisq_o(set_num), min_snp_pval(set_num);
    vector<string> min_snp_name(set_num);
    vector<int> snp_num_in_set(set_num);

    string rgoodsnpfile = _out + ".seg.snpset";
    ofstream rogoodsnp;
    if (write_snpset) rogoodsnp.open(rgoodsnpfile.c_str());
 
    for (i = 0; i < set_num; i++) {
        bool skip = false;
        vector<int> snp_indx = snp_set_indx[i];
        if(snp_indx.size() < 1) skip = true;
        snp_num_in_set[i] = snp_indx.size();
        if(!skip && snp_num_in_set[i] > 20000){
            cout<<"Warning: Too many SNPs in the set on [chr" << set_chr[i] << ":" << set_start_bp[i] << "-" << set_end_bp[i] << "]. Maximum limit is 20000. This gene is ignored in the analysis."<<endl;
            skip = true;  
        } 
        if(skip){
            set_pval[i] = 2.0;
            snp_num_in_set[i] = 0;
            continue;
        }
        chisq_o[i] = 0; 
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = 0; j < snp_indx.size(); j++) {
            if (min_snp_pval[i] > snp_pval[snp_indx[j]]) {
                min_snp_pval[i] = snp_pval[snp_indx[j]];
                min_snp_name[i] = snp_name[snp_indx[j]];
            }
            chisq_o[i] += snp_chisq[snp_indx[j]];
        }
        if(snp_num_in_set[i] == 1) set_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            snp_count=snp_num_in_set[i];
            VectorXd eigenval;
            vector<int> sub_indx;
            sbat_calcu_lambda(snp_indx, eigenval, snp_count, reduce_cor, sub_indx);
            //recalculate chisq value from low correlation snp subset
            if (reduce_cor) {
                chisq_o[i] = 0;
                for (j = 0; j < sub_indx.size(); j++) chisq_o[i] += snp_chisq[snp_indx[sub_indx[j]]]; 
            }
            set_pval[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);
            snp_num_in_set[i] = snp_count;

           if (write_snpset) {
                rogoodsnp << "SEG" << i << ":" << set_start_bp[i] << "-" << set_end_bp[i] << endl;
                for (int k = 0; k < sub_indx.size(); k++) rogoodsnp << snp_name[snp_indx[sub_indx[k]]] << endl;
                rogoodsnp << "END" << endl << endl;
            }

        }

        if((i + 1) % 100 == 0 || (i + 1) == set_num) cout << i + 1 << " of " << set_num << " sets.\r";
    }

    string filename = _out + ".seg.sbat";
    cout << "\nSaving the results of the segment-based SBAT analyses to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Chr\tStart\tEnd\tNo.SNPs\tChisq(Obs)\tPvalue\tTopSNP_Pvalue\tTopSNP" << endl;
    for (i = 0; i < set_num; i++) {
        if(set_pval[i]>1.5) continue;
        ofile << set_chr[i] << "\t" << set_start_bp[i] << "\t"<< set_end_bp[i] << "\t";
        ofile << snp_num_in_set[i] << "\t" << chisq_o[i] << "\t" << set_pval[i] << "\t";
        ofile << min_snp_pval[i] << "\t" << min_snp_name[i] << endl;
 
    }
    ofile.close();
    if (write_snpset) {
        cout << "Writing snpset ..." << endl;
        rogoodsnp.close();
    }
}

void gcta::get_sbat_seg_blk(int seg_size, vector< vector<int> > &snp_set_indx, vector<int> &set_chr, vector<int> &set_start_bp, vector<int> &set_end_bp)
{
    int i = 0, j = 0, k = 0, m = _include.size();

    vector<int> brk_pnt;
    brk_pnt.push_back(0);
    for (i = 1, j = 0; i < m; i++) {
        if (i == (m - 1)) brk_pnt.push_back(m - 1);
        else if (_chr[_include[i]] != _chr[_include[brk_pnt[j]]]) {
            brk_pnt.push_back(i - 1);
            j++;
            brk_pnt.push_back(i);
            j++;
        }
        else if (_bp[_include[i]] - _bp[_include[brk_pnt[j]]] > seg_size) {
            brk_pnt.push_back(i - 1);
            j++;
            brk_pnt.push_back(i);
            j++;
        }
    }

    snp_set_indx.clear();
    set_chr.clear();
    set_start_bp.clear();
    set_end_bp.clear();
    for (i = 0; i < brk_pnt.size() - 1; i++) {
        int size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if(size < 3 && (i%2 != 0)) continue;
        vector<int> snp_indx(size);
        for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++){
            snp_indx[k] = j;  
        } 
        snp_set_indx.push_back(snp_indx);
        set_chr.push_back(_chr[_include[brk_pnt[i]]]);
        set_start_bp.push_back(_bp[_include[brk_pnt[i]]]);
        set_end_bp.push_back(_bp[_include[brk_pnt[i + 1]]]);
    }
}

void gcta::sbat_calcu_lambda(vector<int> &snp_indx, VectorXd &eigenval, int &snp_count, bool reduce_cor, vector<int> &sub_indx)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = snp_indx.size();

    MatrixXf X;
    make_XMat_subset(X, snp_indx, false);
    vector<int> rm_ID1;
    double R_cutoff = 0.9486833; //sqrt(0.9) 
    int qi = 0; //alternate index

    VectorXd sumsq_x(m);
    for (j = 0; j < m; j++) sumsq_x[j] = X.col(j).dot(X.col(j));

    MatrixXf C = X.transpose() * X;
    X.resize(0,0);
    #pragma omp parallel for private(j)
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            double d_buf = sqrt(sumsq_x[i] * sumsq_x[j]);
            if(d_buf>0.0) C(i,j) /= d_buf;
            else C(i,j) = 0.0;
        }
    }

    if (reduce_cor) rm_cor_sbat(C, R_cutoff, m, rm_ID1);
        //Create new index
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sub_indx.push_back(i);
            else {
                if (rm_ID1[qi] == i) qi++; //Skip removed snp
                else sub_indx.push_back(i);
            }
        }
        snp_count = sub_indx.size();
        if (sub_indx.size() < C.size()) { //Build new matrix
            MatrixXf D(sub_indx.size(),sub_indx.size());
            for (i = 0 ; i < sub_indx.size() ; i++) {
               for (j = 0 ; j < sub_indx.size() ; j++) {
                   D(i,j) = C(sub_indx[i],sub_indx[j]);
               }
            }
            C = D; 
        }
    
    SelfAdjointEigenSolver<MatrixXf> saes(C);
    eigenval = saes.eigenvalues().cast<double>();
}

void gcta::rm_cor_sbat(MatrixXf &R, double R_cutoff, int m, vector<int> &rm_ID1) {
    //Modified version of rm_cor_indi from grm.cpp
    
    int i = 0, j = 0, i_buf = 0;
    vector<int> rm_ID2;

    //float tmpr = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < i; j++) {
            //tmpr = R(i,j);
            //if (fabs(tmpr) > R_cutoff ) { 
            if (fabs(R(i,j)) > R_cutoff ) { 
                rm_ID1.push_back(i);
                rm_ID2.push_back(j);
            }
        }
    }

    // count the number of appearance of each "position" in the vector, which involves a few steps
    vector<int> rm_uni_ID(rm_ID1);
    rm_uni_ID.insert(rm_uni_ID.end(), rm_ID2.begin(), rm_ID2.end());
    stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    map<int, int> rm_uni_ID_count;
    for (i = 0; i < rm_uni_ID.size(); i++) {
        i_buf = count(rm_ID1.begin(), rm_ID1.end(), rm_uni_ID[i]) + count(rm_ID2.begin(), rm_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
    }

    // swapping
    map<int, int>::iterator iter1, iter2;
    for (i = 0; i < rm_ID1.size(); i++) {
        iter1 = rm_uni_ID_count.find(rm_ID1[i]);
        iter2 = rm_uni_ID_count.find(rm_ID2[i]);
        if (iter1->second < iter2->second) {
            i_buf = rm_ID1[i];
            rm_ID1[i] = rm_ID2[i];
            rm_ID2[i] = i_buf;
        }
    }
    stable_sort(rm_ID1.begin(), rm_ID1.end());
    rm_ID1.erase(unique(rm_ID1.begin(), rm_ID1.end()), rm_ID1.end());
}

