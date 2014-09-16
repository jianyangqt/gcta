void gcta::grm_m(bool output_bin, int grm_mtd)
{
    check_autosome();
    if(_mu.empty()) calcu_mu();
	rm_high_ld();

    int i=0, j=0, k=0;
    double denom=0.0;

	// Calcuate A matrix
	cout<<"\nCalculating the genetic relationship matrix ...(NOTE: save memory mode with lower speed)"<<endl;
	vector<double> var_SNP(_include.size()), var_SNP_i(_include.size()); // variance of each SNP, 2pq
	vector<double> x(_include.size()), vbuf(_include.size());
	for(j=0; j<_include.size(); j++){
        var_SNP[j]=_mu[_include[j]]*(1.0-0.5*_mu[_include[j]]);
        if(fabs(var_SNP[j])<1.0e-50) var_SNP_i[j]=0.0;
        else var_SNP_i[j]=1.0/var_SNP[j];
	}
	vector< vector<double> > A(_keep.size());
	vector< vector<double> > A_N(_keep.size());
	for(i=0; i<_keep.size(); i++){
	    A[i].resize(_keep.size());
	    A_N[i].resize(_keep.size());
	}
	for(i=0; i<_keep.size(); i++){
 		for(k=0, denom=0.0; k<_include.size(); k++){
            if(!_snp_a[i][_include[k]] || _snp_b[i][_include[k]]){
                if(_allele1[_include[k]]==_ref_A[_include[k]]) x[k]=_snp_a[i][_include[k]]+_snp_b[i][_include[k]];
                else x[k]=2.0-(_snp_a[i][_include[k]]+_snp_b[i][_include[k]]);
                if(grm_mtd==0){
                    A[i][i]+=(x[k]*(x[k]-1.0-_mu[_include[k]])+0.5*_mu[_include[k]]*_mu[_include[k]])*var_SNP_i[k];
                    denom+=1.0;
                }
                else{
                    A[i][i]+=(x[k]*(x[k]-1.0-_mu[_include[k]])+0.5*_mu[_include[k]]*_mu[_include[k]]);
                    denom+=var_SNP[k];
                }
                x[k]-=_mu[_include[k]];
                vbuf[k]=x[k]*var_SNP_i[k];
                A_N[i][i]++;
            }
            else x[k]=1e6;
        }
        if(denom>0.0) A[i][i]=1.0+A[i][i]/denom;
        else A[i][i]=1.0;
		for(j=0; j<i; j++){
 			for(k=0, denom=0.0; k<_include.size(); k++){
                if(x[k]<1e5 && (!_snp_a[j][_include[k]] || _snp_b[j][_include[k]])){
                    if(grm_mtd==0){
                        A[i][j]+=vbuf[k]*(_snp_a[j][_include[k]]+_snp_b[j][_include[k]]-_mu[_include[k]]);
                        denom+=1.0;
                    }
                    else{
                        A[i][j]+=x[k]*(_snp_a[j][_include[k]]+_snp_b[j][_include[k]]-_mu[_include[k]]);
                        denom+=var_SNP[k];
                    }
                    A_N[i][j]++;
                }
            }
            if(denom>0.0) A[i][j]/=denom;
            else A[i][j]=0.0;
		}
		cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}
	cout<<endl<<endl;
	output_grm_vec(A, A_N, output_bin);
}

// old version 23 May 2011
void gcta::grm_s(bool grm_xchr_flag, bool output_bin, int grm_mtd)
{
    if(grm_xchr_flag) check_chrX();
    else check_autosome();
	rm_high_ld();

	int i=0, j=0, k=0;
	double denom=0.0;
    vector< vector<float> > X;
    make_XMat(X, false);
    vector<double> sd_SNP_i, var_SNP;
    if(grm_mtd==0) std_XMat(X, sd_SNP_i, grm_xchr_flag, true);
    else{
        std_XMat(X, sd_SNP_i, grm_xchr_flag, false);
        var_SNP.resize(_include.size());
        for(i=0; i<_include.size(); i++){
            if(sd_SNP_i[i]>0.0) var_SNP[i]=1.0/(sd_SNP_i[i]*sd_SNP_i[i]);
        }
    }

	// Calcuate A matrix
	if(_dosage_flag){
        if(grm_xchr_flag) cout<<"\nCalculating the genetic relationship matrix for the X chromosome using imputed dosage data ... "<<endl;
	    else cout<<"\nCalculating the genetic relationship matrix using imputed dosage data ... "<<endl;
	}
	else{
        if(grm_xchr_flag) cout<<"\nCalculating the genetic relationship matrix for the X chromosome ... ((NOTE: default speed-optimized mode, may use huge RAM))"<<endl;
        else cout<<"\nCalculating the genetic relationship matrix ... (NOTE: default speed-optimized mode, may use huge RAM)"<<endl;
	}
	vector< vector<double> > A(_keep.size());
	vector< vector<int> > A_N(_keep.size());
	for(i=0; i<_keep.size(); i++){
	    A[i].resize(_keep.size());
	    A_N[i].resize(_keep.size());
	}
	for(i=0; i<_keep.size(); i++){
		for(k=0, denom=0.0; k<_include.size(); k++){
            if(X[i][k]<1e5){
                if(grm_mtd==0) A[i][i]+=X[i][k]*X[i][k];//X[i][k]*X[i][k]+(_mu[_include[k]]-1.0)*X[i][k]*sd_SNP_i[k];
                else{
                    A[i][i]+=(X[i][k]+_mu[_include[k]])*(X[i][k]-1.0)+0.5*_mu[_include[k]]*_mu[_include[k]];
                    denom+=var_SNP[k];
                }
                A_N[i][i]++;
            }
        }
        if(A_N[i][i]>0){
            if(grm_mtd==0) A[i][i]/=(double)A_N[i][i];
            else A[i][i]=1.0+A[i][i]/denom;
        }
        else A[i][i]=1.0;
		for(j=0; j<i; j++){
			for(k=0, denom=0.0; k<_include.size(); k++){
                if(X[i][k]<1e5 && X[j][k]<1e5){
                    A[i][j]+=X[i][k]*X[j][k];
                    if(grm_mtd==0) denom+=1.0;
                    else denom+=var_SNP[k];
                    A_N[i][j]++;
                }
            }
            if(denom>0.0) A[i][j]/=denom;
            else A[i][j]=0.0;
		}
		cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}
	cout<<endl<<endl;
	output_grm_vec(A, A_N, output_bin);
}

void gcta::grm_between(gcta *plk, bool grm_xchr_flag, bool output_bin, int grm_mtd)
{
    if(grm_xchr_flag) check_chrX();
    else check_autosome();

    // SNPs in common
    int i=0, j=0, k=0;
    cout<<"\nMatching the SNPs between the two datasets ... "<<endl;
    string str_buf;
    vector<string> comm_SNPs;
    for(i=0; i<plk->_include.size(); i++){
        str_buf=plk->_snp_name[plk->_include[i]];
        if(_snp_name_map.find(str_buf)!=_snp_name_map.end()) comm_SNPs.push_back(str_buf);
    }
    int N=comm_SNPs.size();
    if(N==0) throw("Error: no SNP in common between the two datasets.");
    cout<<"There are "<<N<<" SNPs in common between the two datasets."<<endl;
    update_id_map_kp(comm_SNPs, _snp_name_map, _include);
    plk->update_id_map_kp(comm_SNPs, plk->_snp_name_map, plk->_include);

    // check SNP name and reference allele
    stringstream errmsg;
    for(i=0; i<_include.size(); i++){
        if(_snp_name[_include[i]]!=plk->_snp_name[plk->_include[i]]){
            errmsg<<"Error: the "<<i+1<<" th SNP included in the analysis is "<<_snp_name[_include[i]]<<" in the first dataset but is "<<plk->_snp_name[plk->_include[i]]<<" in the second dataset.";
            throw(errmsg.str());
        }
        if(_ref_A[_include[i]]!=plk->_ref_A[plk->_include[i]]){
            errmsg<<"Error: the reference allele for the SNP "<<_snp_name[_include[i]]<<" is "<<_ref_A[_include[i]]<<" in the first dataset but is "<<plk->_ref_A[plk->_include[i]]<<" in the second dataset.";
            throw(errmsg.str());
        }
    }

    vector< vector<float> > X1, X2;
    cout<<"Dataset 1: ";
    make_XMat(X1, false);
    cout<<"Dataset 2: ";
    plk->make_XMat(X2, false);
    vector<double> sd_SNP_i, var_SNP;
    if(grm_mtd==0){
        std_XMat(X1, sd_SNP_i, grm_xchr_flag, true);
        plk->std_XMat(X2, sd_SNP_i, grm_xchr_flag, true);
    }
    else{
        std_XMat(X1, sd_SNP_i, grm_xchr_flag, false);
        plk->std_XMat(X2, sd_SNP_i, grm_xchr_flag, false);
        var_SNP.resize(_include.size());
        for(i=0; i<_include.size(); i++){
            if(sd_SNP_i[i]>0.0) var_SNP[i]=1.0/(sd_SNP_i[i]*sd_SNP_i[i]);
        }
    }

	// Calcuate A matrix
	cout<<"\nCalculating the genetic relationshp matrix between the individuals in the two datasets ... "<<endl;
	double denom=0.0;
	vector< vector<double> > A(_keep.size());
	vector< vector<int> > A_N(_keep.size());
	for(i=0; i<_keep.size(); i++){
	    A[i].resize(plk->_keep.size());
	    A_N[i].resize(plk->_keep.size());
	}
	for(i=0; i<_keep.size(); i++){
	    denom=0.0;
		for(j=0; j<plk->_keep.size(); j++){
			for(k=0; k<N; k++){
                if(X1[i][k]<1e5 && X2[j][k]<1e5){
                    A[i][j]+=X1[i][k]*X2[j][k];
                    A_N[i][j]++;
                    if(denom>0.0) denom+=var_SNP[k];
                }
            }
            if(grm_mtd=0.0){
                if(A_N[i][j]>0) A[i][j]/=(double)A_N[i][j];
                else A[i][j]=0.0;
            }
            else{
                if(denom>0.0) A[i][j]/=denom;
                else A[i][j]=0.0;
            }
		}
		cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}
	cout<<endl<<endl;

    string grm_file;
    if(output_bin){
        // Save A matrix in binary file
        grm_file=_out+".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out|ios::binary);
        if(!A_Bin) throw("Error: can not open the file ["+grm_file+"] to write.");
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<plk->_keep.size(); j++) A_Bin.write((char *) &(A[i][j]), sizeof(double));
        }
        A_Bin.close();
        cout<<"Genetic relationship matrix betwee "<<_keep.size()<<" individuals in dataset 1 and "<<plk->_keep.size()<<" individuals in dataset 2 has been saved in the file ["+grm_file+"] (in binary format)."<<endl;
    }
    else{
        grm_file=_out+".grm.gz";
        gzofstream zoutf;
        zoutf.open( grm_file.c_str() );
        if(!zoutf.is_open()) throw("Error: can not open the file ["+grm_file+"] to write.");
        cout<<"Saving the genetic relationship matrix to the file ["+grm_file+"] (in compressed text format)."<<endl;
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<plk->_keep.size(); j++) zoutf<<setprecision(0)<<i+1<<'\t'<<j+1+_keep.size()<<'\t'<<A_N[i][j]<<'\t'<<setprecision(15)<<A[i][j]<<endl;
        }
        zoutf.close();
        cout<<"Genetic relationship matrix between "<<_keep.size()<<" individuals in dataset 1 and "<<plk->_keep.size()<<" individuals in dataset 2 has been saved in the file ["+grm_file+"] (in compressed text format)."<<endl;
    }
    string famfile=_out+".grm.id";
	ofstream Fam(famfile.c_str());
	if(!Fam) throw("Error: can not open the file ["+famfile+"] to write.");
	for(i=0; i<_keep.size(); i++) Fam<<_fid[_keep[i]]+"\t"+_pid[_keep[i]]<<endl;
	for(i=0; i<plk->_keep.size(); i++) Fam<<plk->_fid[plk->_keep[i]]+"\t"+plk->_pid[plk->_keep[i]]<<endl;
	Fam.close();
	cout<<"IDs for the GRM file ["+grm_file+"] have been saved in the file ["+famfile+"]."<<endl;
}
