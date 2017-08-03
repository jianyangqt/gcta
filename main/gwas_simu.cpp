/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for GWAS simulation
 *
 * 2014 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

int gcta::read_QTL_file(string qtl_file, vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, vector<int> &have_eff)
{
    qtl_name.clear();
    qtl_pos.clear();
    qtl_eff.clear();
    have_eff.clear();

    ifstream i_qtl(qtl_file.c_str());
    if (!i_qtl) throw ("Error: can not open the file [" + qtl_file + "] to read.");
    string qtl_buf, str_buf;
    double qtl_eff_buf = 0.0;
    cout << "Reading a list of SNPs (as causal variants) from [" + qtl_file + "]." << endl;
    map<string, int>::iterator iter, End = _snp_name_map.end();
    vector<string> vs_buf;
    vector<int> confirm(_snp_num);
    int icount = 0;
    while (i_qtl) {
        i_qtl >> qtl_buf;
        if (i_qtl.eof()) break;
        iter = _snp_name_map.find(qtl_buf);
        if (getline(i_qtl, str_buf) && StrFunc::split_string(str_buf, vs_buf, " \t\n") > 0) {
            have_eff.push_back(1);
            qtl_eff_buf = atof(vs_buf[0].c_str());
            if (fabs(qtl_eff_buf) > 1e5) throw ("Error: invalid effect size specified for the causal variant [" + str_buf + "].");
        } else {
            have_eff.push_back(0);
            qtl_eff_buf = 0.0;
        }
        if (iter != End) {
            qtl_name.push_back(qtl_buf);
            qtl_pos.push_back(iter->second);
            qtl_eff.push_back(qtl_eff_buf);
        }
    }
    vector<string> qtl_name_buf=qtl_name;
    stable_sort(qtl_name_buf.begin(), qtl_name_buf.end());
    qtl_name_buf.erase(unique(qtl_name_buf.begin(), qtl_name_buf.end()), qtl_name_buf.end());
    i_qtl.close();
    
    if(qtl_name_buf.size() < qtl_name.size()) throw("Error: there are duplicated SNP IDs.");
    cout << qtl_pos.size() << " SNPs (as causal variants) to be included from [" + qtl_file + "]." << endl;
    return (qtl_pos.size());
}

void gcta::output_simu_par(vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, double Vp)
{
    int i = 0;
    string out_parfile = _out + ".par";
    ofstream out_par(out_parfile.c_str());
    if (!out_par) throw ("Error: can not open par file [" + out_parfile + "] to write!");
    out_par << "QTL\tRefAllele\tFrequency\tEffect" << endl;
    for (i = 0; i < qtl_eff.size(); i++) out_par << qtl_name[i] << "\t" << _ref_A[qtl_pos[i]] << "\t" << 0.5 * _mu[qtl_pos[i]] << "\t" << qtl_eff[i] << endl;
    out_par.close();
    cout << "Simulated QTL effect(s) have been saved in [" + out_parfile + "]." << endl;
}

void gcta::save_phenfile(vector< vector<double> > &y)
{
    string phenfile = _out + ".phen";
    ofstream phen(phenfile.c_str());
    if (!phen) throw ("Error: can not open the file [" + phenfile + "] to write.");
    int i = 0, j = 0;
    for (i = 0; i < _keep.size(); i++) {
        phen << _fid[_keep[i]] << " " << _pid[_keep[i]] << " ";
        for (j = 0; j < y.size(); j++) phen << y[j][i] << " ";
        phen << endl;
    }
    phen.close();
}

void gcta::GWAS_simu(string bfile, int simu_num, string qtl_file, int case_num, int control_num, double hsq, double K, int seed, bool output_causal, bool simu_emb_flag, int eff_mod)
{
    int i = 0, j = 0;
    bool cc_flag = false;
    if (case_num > 0 || control_num > 0) cc_flag = true;

    cout << "Simulation parameters:" << endl;
    cout << "Number of simulation replicate(s) = " << simu_num << " (Default = 1)" << endl;
    cout << "Heritability " << (cc_flag ? "of liability = " : " = ") << hsq << " (Default = 0.1)" << endl;
    if (cc_flag) {
        cout << "Disease prevalence = " << K << " (Default = 0.1)" << endl;
        cout << "Number of cases = " << case_num << endl;
        cout << "Number of controls = " << control_num << endl;
    }
    cout << endl;

    // Read QTL file
    vector<string> qtl_name;
    vector<int> qtl_pos, have_eff;
    vector<double> qtl_eff;
    int qtl_num = read_QTL_file(qtl_file, qtl_name, qtl_pos, qtl_eff, have_eff);
    update_id_map_kp(qtl_name, _snp_name_map, _include);
    
    // Generate QTL effects
    int Seed = -CommFunc::rand_seed();
    if (hsq > 0.0) {
        int num_gener_qtl_eff = 0;
        for (i = 0; i < qtl_num; i++) {
            if (have_eff[i] == 0) {
                qtl_eff[i] = StatFunc::gasdev(Seed);
                num_gener_qtl_eff++;
            }
        }
        if (qtl_num - num_gener_qtl_eff > 0) cout << qtl_num - num_gener_qtl_eff << " user-specified QTL effects." << endl;
        if (num_gener_qtl_eff > 0) cout << num_gener_qtl_eff << " unspecified QTL effects are generated from standard normal distribution." << endl;

        vector<string> vs_buf(qtl_num);
        for (i = 0; i < qtl_num; i++) vs_buf[i] = _snp_name[_include[i]];
        vector<int> indx;
        StrFunc::match(vs_buf, qtl_name, indx);
        qtl_name = vs_buf;
        vector<double> qtl_eff_buf(qtl_eff);
        vector<int> qtl_pos_buf(qtl_pos);
        for (i = 0; i < qtl_num; i++){
            qtl_eff[i] = qtl_eff_buf[indx[i]];
            qtl_pos[i] = qtl_pos_buf[indx[i]]; 
        }
    } else {
        qtl_eff.clear();
        qtl_eff.resize(qtl_num, 0.0);
    }

    // Calculate allele frequency
    MatrixXf X;
    make_XMat(X);
    if(eff_mod == 0){
        eigenVector sd_SNP;
        std_XMat(X, sd_SNP, false, true, true);
    } else {
        eigenVector sd_SNP;
        std_XMat(X, sd_SNP, false, true, false);
    }

    // Calculate Ve and threhold
    double var_g = 0.0, var_e = 1.0;
    vector<double> g(_keep.size());
    if (hsq > 0.0) {
        for (i = 0; i < _keep.size(); i++) {
            for (j = 0; j < qtl_num; j++) g[i] += X(i,j) * qtl_eff[j];
        }
        var_g = CommFunc::var(g);
        var_e = var_g * (1.0 / hsq - 1.0);
    }
    double sd_e = sqrt(var_e);

    // Output par file
    output_simu_par(qtl_name, qtl_pos, qtl_eff, var_e + var_g);

    // Output phenotype file
    cout << "Simulating GWAS based on the real genotyped data with " << simu_num << " replicate(s) ..." << endl;
    vector< vector<double> > y(simu_num);
    int case_num_buf = 0, control_num_buf = 0;
    for (i = 0; i < simu_num; i++) {
        y[i].resize(_keep.size());
        for (j = 0; j < _keep.size(); j++) {
            if (hsq < 1.0) y[i][j] = g[j] + sd_e * StatFunc::gasdev(Seed);
            else y[i][j] = g[j];
        }
        if (cc_flag) {
            case_num_buf = 0;
            control_num_buf = 0;
            vector<double> y_buf(y[i]);
            stable_sort(y_buf.begin(), y_buf.end());
            int n = (int) (_indi_num * (1.0 - K));
            double Th = 0.5 * (y_buf[n] + y_buf[n - 1]);
            for (j = 0; j < _keep.size(); j++) {
                if (y[i][j] > Th) {
                    if (case_num_buf < case_num) {
                        y[i][j] = 2;
                        case_num_buf++;
                    } else y[i][j] = -9;
                } else {
                    if (control_num_buf < control_num) {
                        y[i][j] = 1;
                        control_num_buf++;
                    } else y[i][j] = -9;
                }
            }
        }
    }

    if (!simu_emb_flag) {
        save_phenfile(y);
        if (cc_flag) cout << "Simulated " << case_num_buf << " cases and " << control_num << " controls have been saved in [" + _out + ".phen" + "]." << endl;
        else cout << "Simulated phenotypes of " << _keep.size() << " individuals have been saved in [" + _out + ".phen" + "]." << endl;
    } else {
        // emBayesB format
        if (!output_causal) update_id_map_rm(qtl_name, _snp_name_map, _include);
        string out_rstfile = _out + ".emb";
        ofstream out_emBayesB(out_rstfile.c_str());
        if (!out_emBayesB) throw ("Error: can not open the file [" + out_rstfile + "] to write.");
        cout << "Saving the simulated data to the file [" + out_rstfile + "] (in emBayesB format)." << endl;
        for (i = 0; i < _keep.size(); i++) {
            if (y[0][i] == -9) continue;
            out_emBayesB << _pid[_keep[i]] << " " << g[i] << " " << y[0][i] << endl;
            for (j = 0; j < _include.size(); j++) {
                if (_snp_1[_include[j]][_keep[i]] && !_snp_2[_include[j]][_keep[i]]) out_emBayesB << _mu[_include[j]] << " ";
                else out_emBayesB << (double) (_snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]]) << " ";
            }
            out_emBayesB << endl;
        }
        out_emBayesB.close();
        cout << "Simulated data (" << _keep.size() << " individuals and " << _include.size() << " SNPs) has been saved in [" + out_rstfile + "]." << endl;
    }
}

/////////////////////////////////////////
/* These codes are not used any more */
/////////////////////////////////////////
/*
void gcta::GenerCases(string bfile, string qtl_file, int case_num, int control_num, double hsq, double K, bool curr_popu, double gnrt)
{
        int i=0, j=0, k=0;

        if(gnrt<1.0 || gnrt>1e5) throw("Error: --simu-gener should be within the range from 1 to 100000.");
        if(hsq>1.0 || hsq<0.0) throw("Error: --simu-h2 should be within the range from 0 to 1.");
        if(K>0.5 || K<0.0001) throw("Error: --simu-K should be within the range from 0.0001 to 0.5.");
        if(!curr_popu && (case_num>1e5 || case_num<1)) throw("Error: --simu-cc, Invalid number of cases.");
        if(!curr_popu && (control_num>1e6 || control_num<1)) throw("Error: --simu-cc, Invalid number of controls.");

        // Read bim file: recombination rate is defined between SNP i and SNP i-1
        read_bimfile(bfile+".bim");
        kosambi();

        // calcualte the accumulative recombination probability over xxx generation
        for(i=0; i<_snp_num; i++) _rc_rate[i]=1.0-pow((double)(1.0-_rc_rate[i]), (double)gnrt);

    // Read QTL file
    vector<string> qtl_name;
        read_snplist(qtl_file, qtl_name);
    vector<int> qtl_pos;
        map<string, int>::iterator iter;
        for(int i=0; i<qtl_name.size(); i++){
        iter=_snp_name_map.find(qtl_name[i]);
        if(iter!=_snp_name_map.end()) qtl_pos.push_back(iter->second);
    }
    stable_sort(qtl_pos.begin(), qtl_pos.end());
        int qtl_num=qtl_pos.size();
        cout<<qtl_num<<" SNPs (as QTLs) to be included from ["+qtl_file+"]."<<endl;

        // Read fam and bed files
        read_famfile(bfile+".fam");
        read_bedfile(bfile+".bed");

        // Generate QTL effects
        int Seed=-CommFunc::rand_seed();
        cout<<"Generating QTL effects ("<<"Random seed = "<<abs(Seed)<<")."<<endl;
        vector<double> qtl_eff(qtl_num);
        for(i=0; i<qtl_num && hsq>0.0; i++) qtl_eff[i]=StatFunc::gasdev(Seed);

        // Calculate Ve and threhold
    int n=0, x=0;
    double var_g=0.0, var_e=1.0;
    vector<double> AF(qtl_num), y(_indi_num), y_backup;
    if(hsq>0.0){
        for(n=0; n<_indi_num; n++){
            for(i=0; i<qtl_num; i++){
                j=qtl_pos[i];
                if(_snp_a[n][j] && !_snp_b[n][j]) continue;
                x=_snp_a[n][j]+_snp_b[n][j];
                AF[i]+=x;
                y[n]+=x*qtl_eff[i];
            }
        }
        var_g=CommFunc::var(y);
            var_e=var_g*(1.0/hsq-1.0);
    }
        double sd_e=sqrt(var_e);
        for(n=0; n<_indi_num; n++) y[n]+=sd_e*StatFunc::gasdev(Seed);
        if(curr_popu) y_backup=y;
        stable_sort(y.begin(), y.end());
        n=(int)(_indi_num*(1.0-K));
        double Th=0.5*(y[n]+y[n-1]);

        // ouput par file
        vector<double> qsq(qtl_num);
        for(i=0; i<qtl_num; i++){
        AF[i]/=_indi_num;
        AF[i]*=0.5;
        qsq[i]=2.0*AF[i]*(1.0-AF[i])*qtl_eff[i]*qtl_eff[i];
        qsq[i]/=(var_e+var_g);
    }
        string out_parfile=_out+".par";
        ofstream out_par(out_parfile.c_str());
        if(!out_par) throw("Error: can not open par file ["+out_parfile+"] to write!");
        cout<<"Writing simulation parameters to file ["<<out_parfile<<"]."<<endl;
        for(i=0; i<qtl_num; i++) out_par<<qtl_name[i]<<"\t"<<AF[i]<<"\t"<<qtl_eff[i]<<"\t"<<qsq[i]<<endl;
    out_par.close();

        // If generate cases and controls based on current population
        if(curr_popu){
            cout<<"Generating cases and control based on the current population."<<endl;
        for(i=0; i<_indi_num; i++){
            if(y_backup[i]>Th) _pheno[i]=2;
            else _pheno[i]=1;
        }
        save_famfile();
        return;
        }

        // Generate population
    int fa_hap_buf=0, mo_hap_buf, fa=0, mo=0;
    double score=0.0;
        vector< vector<int> > case_fa, case_mo, case_fa_hap, case_mo_hap, control_fa, control_mo, control_fa_hap, control_mo_hap;
        cout<<"Generating population derived from the current population."<<endl;
        while(case_fa.size()<case_num || control_fa.size()<control_num){
        vector<int> fa_vbuf, mo_vbuf, fa_hap_vbuf, mo_hap_vbuf;
            for(i=0; i<_snp_num; i++){
            if(StatFunc::ran1(Seed)<_rc_rate[i]){
                fa_hap_vbuf.push_back(i);
                fa_vbuf.push_back((int)(StatFunc::ran1(Seed)*_indi_num*2));
            }
            if(StatFunc::ran1(Seed)<_rc_rate[i]){
                mo_hap_vbuf.push_back(i);
                mo_vbuf.push_back((int)(StatFunc::ran1(Seed)*_indi_num*2));
            }
        }
        score=0.0;
        for(i=0; i<qtl_num; i++){
            j=qtl_pos[i];
            fa_hap_buf=upper_bound(fa_hap_vbuf.begin(), fa_hap_vbuf.end(), j)-fa_hap_vbuf.begin()-1;
            mo_hap_buf=upper_bound(mo_hap_vbuf.begin(), mo_hap_vbuf.end(), j)-mo_hap_vbuf.begin()-1;
            fa=(int)(fa_vbuf[fa_hap_buf]*0.5);
            mo=(int)(mo_vbuf[mo_hap_buf]*0.5);
            if( (_snp_a[fa][j] && !_snp_b[fa][j]) || (_snp_a[mo][j] && !_snp_b[mo][j]) ) continue;
            if(fa_vbuf[fa_hap_buf]%2==0) x=_snp_a[fa][j];
            else x=_snp_b[fa][j];
            if(mo_vbuf[mo_hap_buf]%2==0) x+=_snp_a[mo][j];
            else x+=_snp_b[mo][j];
            score+=x*qtl_eff[i];
        }
        score+=sd_e*StatFunc::gasdev(Seed);
        if(score>Th){
            if(case_fa.size()<case_num){
                case_fa.push_back(fa_vbuf);
                case_mo.push_back(mo_vbuf);
                case_fa_hap.push_back(fa_hap_vbuf);
                case_mo_hap.push_back(mo_hap_vbuf);
            }
        }
        else{
            if(control_fa.size()<control_num){
                control_fa.push_back(fa_vbuf);
                control_mo.push_back(mo_vbuf);
                control_fa_hap.push_back(fa_hap_vbuf);
                control_mo_hap.push_back(mo_hap_vbuf);
            }
        }
    }
        cout<<case_num<<" cases and "<<control_num<<" controls have been generated."<<endl;
        _indi_num=case_num+control_num;

        // Output fam file
        string out_famfile=_out+".fam";
    ofstream out_fam(out_famfile.c_str());
    if(!out_fam) throw("Error: can not open fam file ["+out_famfile+"] to write!");
    cout<<"Writing IDs of cases and controls to file ["<<out_famfile<<"]."<<endl;
    for(i=0; i<case_num; i++) out_fam<<i+1<<"\t"<<i+1<<"\t0\t0\t0\t2"<<endl;
    for(i=case_num; i<_indi_num; i++) out_fam<<i+1<<"\t"<<i+1<<"\t0\t0\t0\t1"<<endl;
    out_fam.close();

        // Output bed file
        // Output first three bytes
    // Output genotype in binary format
    case_fa.insert(case_fa.end(), control_fa.begin(), control_fa.end());
    case_mo.insert(case_mo.end(), control_mo.begin(), control_mo.end());
    case_fa_hap.insert(case_fa_hap.end(), control_fa_hap.begin(), control_fa_hap.end());
    case_mo_hap.insert(case_mo_hap.end(), control_mo_hap.begin(), control_mo_hap.end());
    save_bedfile(case_fa, case_mo, case_fa_hap, case_mo_hap);
}


void gcta::kosambi()
{
    int i=0;
        double c=0.0;
        _rc_rate.clear();
        _rc_rate.resize(_snp_num);
    for(i=0; i<_snp_num; i++){
        if(i==0 || _chr[i]!=_chr[i-1]) _rc_rate[i]=0.5;
        else{
            c=exp(_genet_dst[i]*0.04);
            _rc_rate[i]=0.5*(c-1.0)/(c+1.0);
        }
    }
}
 */

/*
void gcta::save_bedfile(vector< vector<int> > &fa_indx, vector< vector<int> > &mo_indx, vector< vector<int> > &fa_hap, vector< vector<int> > &mo_hap, bool GENOME)
{
    bool fa_geno=false, mo_geno=false;
    int i=0, pos=0, n=0, fa=0, mo=0, fa_hap_buf=0, mo_hap_buf=0;
    string OutBedFile=_out+".bed";
        fstream OutBed(OutBedFile.c_str(), ios::out|ios::binary);
        if(!OutBed) throw("Error: can not open the file ["+OutBedFile+"] to write.");
        cout<<"Writing genotypes to PLINK BED file ["+OutBedFile+"]."<<endl;
        bitset<8> b;
        char ch[1];
    b.reset();
    b.set(2);  b.set(3);  b.set(5);  b.set(6);
    ch[0] = (char)b.to_ulong();
    OutBed.write(ch,1);
    b.reset();
    b.set(0);  b.set(1);  b.set(3);  b.set(4);
    ch[0] = (char)b.to_ulong();
    OutBed.write(ch,1);
    b.reset();
    b.set(0);
    ch[0] = (char)b.to_ulong();
    OutBed.write(ch,1);
    for(i=0; i<_snp_num; i++){
        pos=0;
        b.reset();
        for(n=0; n<_indi_num; n++){
            if(!GENOME){
                fa_hap_buf=upper_bound(fa_hap[n].begin(), fa_hap[n].end(), i)-fa_hap[n].begin()-1;
                mo_hap_buf=upper_bound(mo_hap[n].begin(), mo_hap[n].end(), i)-mo_hap[n].begin()-1;
                fa=(int)(fa_indx[n][fa_hap_buf]*0.5);
                mo=(int)(mo_indx[n][mo_hap_buf]*0.5);
            }
            else{
                fa=n;
                mo=n;
            }
                        if( (_snp_a[fa][i] && !_snp_b[fa][i]) || (_snp_a[mo][i] && !_snp_b[mo][i]) ){
                b[pos++]=1;
                b[pos++]=0;
            }
            else{
                if(!GENOME){
                    if(fa_indx[n][fa_hap_buf]%2==0) fa_geno=_snp_a[fa][i];
                    else fa_geno=_snp_b[fa][i];
                    if(mo_indx[n][mo_hap_buf]%2==0) mo_geno=_snp_a[mo][i];
                    else mo_geno=_snp_b[mo][i];
                    fa_geno=(!fa_geno);
                    mo_geno=(!mo_geno);
                }
                else{
                   fa_geno=_snp_a[fa][i];
                   mo_geno=_snp_b[mo][i];
                }
                b[pos++]=min(fa_geno, mo_geno);
                b[pos++]=max(fa_geno, mo_geno);
            }
            if(pos>7 || n==_indi_num-1){
                ch[0]=(char)b.to_ulong();
                OutBed.write(ch,1);
                pos=0;
                b.reset();
            }
        }
    }
    OutBed.close();
}

void gcta::save_famfile()
{
    string famfile=_out+".fam";
        ofstream Fam(famfile.c_str());
        if(!Fam) throw("Error: can not open the fam file "+famfile+" to save!");
        cout<<"Writing PLINK FAM file to ["+famfile+"]."<<endl;
        int i=0;
        for(i=0; i<_indi_num; i++){
                Fam<<_fid[i]<<"\t"<<_pid[i]<<"\t"<<_fa_id[i]<<"\t"<<_mo_id[i]<<"\t"<<_sex[i]<<"\t"<<_pheno[i]<<endl;
        }
        Fam.close();
        cout<<_indi_num<<" individuals to be saved to ["+famfile+"]."<<endl;
}

void gcta::save_bimfile()
{
        int i=0;
        string bimfile=_out+".bim";
        ofstream Bim(bimfile.c_str());
        if(!Bim) throw("Error: can not open the file ["+bimfile+"] to write.");
        cout<<"Writing PLINK bim file to ["<<bimfile<<"]."<<endl;
        for(i=0; i<_snp_num; i++){
                Bim<<_chr[i]<<"\t"<<_snp_name[i]<<"\t"<<_genet_dst[i]<<"\t"<<_bp[i]<<"\t"<<_allele1[i]<<"\t"<<_allele2[i]<<endl;
        }
        Bim.close();
        cout<<_snp_num<<" SNPs to be saved to ["<<bimfile<<"]."<<endl;
}
 */
void gcta::genet_dst(string bfile, string hapmap_genet_map)
{
    // Read bim file
    read_bimfile(bfile + ".bim");
    int snp_num = _snp_name.size();

    // Read HAPMAP genetic map files
    int i = 0, j = 0;
    string str_buf;
    vector<string> vs_buf;
    string genet_mapfile;
    vector< vector< vector<double> > > hap_genet(_autosome_num);
    for (i = 0; i < _autosome_num; i++) {
        stringstream str_strm;
        str_strm << hapmap_genet_map << i + 1 << "_CEU_b36.txt";
        genet_mapfile = str_strm.str();
        ifstream i_genet_map(genet_mapfile.c_str());
        if (!i_genet_map) throw ("Error: can not open HAPMAP genetic map file " + genet_mapfile + "!");
        hap_genet[i].resize(2);
        getline(i_genet_map, str_buf);
        while (getline(i_genet_map, str_buf)) {
            if (StrFunc::split_string(str_buf, vs_buf) < 3) continue;
            hap_genet[i][0].push_back(atof(vs_buf[0].c_str()));
            hap_genet[i][1].push_back(atof(vs_buf[1].c_str()));
        }
        i_genet_map.close();
    }

    // calculate genetic distance
    int pos1 = 0, pos2 = 0;
    double prev_bp = 0.0;
    vector<double> dst(snp_num);
    vector<double>::iterator iter1, iter2;
    for (i = 0; i < snp_num; i++) {
        if (i == 0 || _chr[i - 1] != _chr[i]) {
            dst[i] = 0.0;
            continue;
        }
        iter1 = upper_bound(hap_genet[_chr[i] - 1][0].begin(), hap_genet[_chr[i] - 1][0].end(), _bp[i - 1]);
        iter2 = upper_bound(hap_genet[_chr[i] - 1][0].begin(), hap_genet[_chr[i] - 1][0].end(), _bp[i]);
        if (iter1 == hap_genet[_chr[i] - 1][0].end() || iter2 == hap_genet[_chr[i] - 1][0].end()) {
            dst[i] = _bp[i] - _bp[i - 1];
            continue;
        }
        pos1 = iter1 - hap_genet[_chr[i] - 1][0].begin();
        pos2 = iter2 - hap_genet[_chr[i] - 1][0].begin();
        prev_bp = _bp[i - 1];
        while (pos1 < pos2) {
            dst[i] += (hap_genet[_chr[i] - 1][0][pos1] - prev_bp) * hap_genet[_chr[i] - 1][1][pos1 - 1];
            prev_bp = hap_genet[_chr[i] - 1][0][pos1];
            pos1++;
        }
        dst[i] += (_bp[i] - prev_bp) * hap_genet[_chr[i] - 1][1][pos1 - 1];
    }

    // Output fam file
    string out_bimfile = _out + ".genetdst";
    ofstream out_bim(out_bimfile.c_str());
    if (!out_bim) throw ("Error: can not open file " + out_bimfile + " to write!");
    for (i = 0; i < snp_num; i++) out_bim << _chr[i] << "\t" << _snp_name[i] << "\t" << dst[i]*1e-6 << "\t" << _bp[i] << "\t" << _allele1[i] << "\t" << _allele2[i] << endl;
    out_bim.close();
    cout << "Genetic distances have been created, and been saved in [" + out_bimfile + "]." << endl;
}

/*
void gcta::simu_genome(
          string popSize, 						// effective population size of each population
                  int nSubPOP, 							// number of populations
                  vector<int> & nSubSample, 			// number of samples (haplotypes) draw from each populations
                  int numPieces, 						// number of fragments for each sample (chromosome)
                  int pieceLen, 						// length in base pair of each fragment
                  int numIndepRegion, 					// number of independent regions (independent chromosome)
                  int s,								// fixed number of SNPs want to simulate, randomly place s SNPs on the genealogy
                  string rec,							// recombination rate between consecutive fragments per generation
                  double mut, 							// mutation rate per generation per base pair
                  double mig 							// migration rate per generation
                  )
{
    int seed=CommFunc::rand_seed();
    vector< vector<bool> > data; // the return chromosome by connecting all independent chromosome into one long chromosome

    cout<<"\n*********************************************"<<endl;
    cout<<"The following output is generated by the program GENOME\n"<<endl;
    genome(popSize, nSubPOP, nSubSample, numPieces, pieceLen, numIndepRegion, s, rec, mut, mig, data, seed, false, false);
    cout<<"\nEnd of output by GENOME"<<endl;
    cout<<"*********************************************\n"<<endl;
    _snp_a.clear();
    _snp_b.clear();
    _pheno.clear();
    int i=0, j=0, k=0, index=0;
    for(k=0; k<nSubPOP; k++){
        for(i=0; i<nSubSample[k]; i++){
            _snp_a.push_back(data[index++]);
            _snp_b.push_back(data[index++]);
            _pheno.push_back(k+1);
            i++;
            if(i==nSubSample[k]-2) break;
        }
    }
    _indi_num=_snp_a.size();
    _snp_num=_snp_a[0].size();

    // save famfile
        stringstream strstrm;
    _fid.clear();
    _fid.resize(_indi_num);
    _pid.clear();
    _pid.resize(_indi_num);
    _fa_id.clear();
    _fa_id.resize(_indi_num);
    _mo_id.clear();
    _mo_id.resize(_indi_num);
    _sex.clear();
    _sex.resize(_indi_num);
    for(i=0; i<_indi_num; i++){
        strstrm.str("");
        strstrm<<i+1;
        _fid[i]=strstrm.str();
                _pid[i]=strstrm.str();
                _fa_id[i]="0";
                _mo_id[i]="0";
        }
        save_famfile();

        // save bimfile
        _chr.clear();
        _chr.resize(_snp_num);
        _snp_name.clear();
        _snp_name.resize(_snp_num);
        _genet_dst.clear();
        _genet_dst.resize(_snp_num);
        _bp.clear();
        _bp.resize(_snp_num);
        _allele1.clear();
        _allele1.resize(_snp_num);
        _allele2.clear();
        _allele2.resize(_snp_num);
    for(i=0; i<_snp_num; i++){
        _chr[i]=1;
        strstrm.str("");
        strstrm<<"SNP_"<<i+1;
        _snp_name[i]=strstrm.str();
        _allele1[i]='A';
        _allele2[i]='C';
    }
    save_bimfile();

    // save bedfile
    vector< vector<int> > tmp;
    save_bedfile(tmp, tmp, tmp, tmp, true);
}
 */

/*
void gcta::simu_geno_unlinked(int N, int M, double maf)
{
    int i=0, j=0, x=0;
    
    // fam file
    _keep.resize(N);
    _fid.resize(N);
    _pid.resize(N);
    _fa_id.resize(N);
    _mo_id.resize(N);
    _sex.resize(N);
    _pheno.resize(N);
    for(i=0; i<N; i++){
        _keep[i]=i;
        stringstream ss;
        ss<<i+1;
        _fid[i]=ss.str();
        _pid[i]=_fid[i];
        _fa_id[i]="0";
        _mo_id[i]="0";
        _sex[i]=-9;
        _pid[i]=-9;
    }
    
    // bim file
    _include.resize(M);
    _snp_num=M;
    _chr.resize(M);
    _snp_name.resize(M);
        _bp.resize(M);
    _genet_dst.resize(M);
    _allele1.resize(M);
    _allele2.resize(M);
    _snp_1.resize(M);
    _snp_2.resize(M);

//double p=0.0;
    std::tr1::minstd_rand eng;
        eng.seed((unsigned int)time(NULL));
	
        cout<<"maf "<<maf<<endl;
	
    for(j=0; j<M; j++){
        _include[j]=j;
        stringstream ss;
        ss<<"SNP"<<j+1;
        _chr[j]=1;
        _snp_name[j]=ss.str();
        _bp[j]=j+1;
        _genet_dst[j]=0.0;
        _allele1[j]="A";
        _allele2[j]="G";
        _snp_1[j].resize(N);
        _snp_2[j].resize(N);
                std::tr1::uniform_real<double> runiform(maf,1-maf);
        double p = runiform(eng)/1.0e10;
		
                //debug
                cout<<"p = "<<p<<endl;
		
        for(i=0; i<N; i++){
            std::tr1::binomial_distribution<int, double> rbinom(2,0.5);
            x=rbinom(eng);
			
                        //debug
                        cout<<x<<"\t";
			
			
            if(x==2) _snp_1[j][i]=_snp_2[j][i]=true;
            else if(x==1){
                _snp_1[j][i]=false;
                _snp_2[j][i]=true;
            }
            else _snp_1[j][i]=_snp_2[j][i]=false;
        }  
		
                //debug
                cout<<endl;
    }
    
    save_plink();
}
 */

