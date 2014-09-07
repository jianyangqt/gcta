/*
 * Implementations of the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "CommFunc.h"

double CommFunc::Abs(const double &x)
{
	complex<double> cld(x);
	double ldAbs = abs(cld);
	return(ldAbs);
}

double CommFunc::sum(const vector<double> &x)
{
    int size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    return (double)d_buf;
}

double CommFunc::mean(const vector<double> &x)
{
    int size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    d_buf/=(double)size;
    return (double)d_buf;
}

double CommFunc::median(const vector<double> &x)
{
    vector<double> b(x);
    int size = b.size();
    if(size==1) return b[0];
    stable_sort(b.begin(), b.end());
    if(size%2==1) return b[(size-1)/2];
    else return (b[size/2]+b[size/2-1])/2;
}

double CommFunc::var(const vector<double> &x)
{
    int size = x.size();
    if(size<=1) return(0.0);
    int i=0;
    double mu=0.0, s2=0.0;
    for(i=0; i<size; i++) mu+=x[i];
    mu/=(double)size;
    for(i=0; i<size; i++) s2+=(x[i]-mu)*(x[i]-mu);
    s2/=(double)(size-1);
    return (double)s2;
}

double CommFunc::cov(const vector<double> &x, const vector<double> &y)
{
    int size = x.size();
    int i=0;
    double mu1=0.0, mu2=0.0, c=0.0;
    for(i=0; i<size; i++){
        mu1+=x[i];
        mu2+=y[i];
    }
    mu1/=(double)size;
    mu2/=(double)size;

    for(i=0; i<size; i++) c+=(x[i]-mu1)*(y[i]-mu2);
    c/=(double)(size-1);
    return c;
}

bool CommFunc::FloatEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) < FloatErr) return true;
	return false;
}

bool CommFunc::FloatNotEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) >= FloatErr) return true;
	return false;
}

const double CommFunc::Sqr(const double &a)
{
	return a*a;
}

const double CommFunc::Max(const double &a, const double &b)
{
	return b > a ? (b) : (a);
}

const double CommFunc::Min(const double &a, const double &b)
{
	return b < a ? (b) : (a);
}

const double CommFunc::Sign(const double &a, const double &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

int CommFunc::rand_seed()
{
  	stringstream str_strm;
	str_strm<<time(NULL);
	string seed_str=str_strm.str();
    reverse(seed_str.begin(), seed_str.end());
    seed_str.erase(seed_str.begin()+7, seed_str.end());
	return(abs(atoi(seed_str.c_str())));
}

void CommFunc::FileExist(string filename)
{
    ifstream ifile(filename.c_str());
    if(!ifile) throw("Error: can not open the file ["+filename+"] to read.");
}

// find the missing value
VectorXi CommFunc::find_missing_vec(eigenVector X) {
    int i=0, j=0;
    int length = X.size();
    VectorXi miss_vec(length);
    
    miss_vec.setZero();
    #pragma omp parallel for private(i)
    for( i = 0; i<length; i++) {
        if( X(i)<1e9 ) miss_vec(i)=1;
    }
    return miss_vec;
}

MatrixXi CommFunc::find_missing_X(eigenMatrix X) {
    int i=0;
    int nrow = X.rows(), ncol = X.cols();
    MatrixXi miss_mat(nrow, ncol);
    eigenVector xbuf;
    VectorXi miss_vec;
   
    #pragma omp parallel for private(i)
    for( i = 0; i<nrow; i++) {
        xbuf = X.row(i);
        miss_vec = find_missing_vec(xbuf);
        miss_mat.row(i) = miss_vec;
    }
    return(miss_mat);
}

int *CommFunc::find_missing_pos(MatrixXi miss_mat, int *miss_indx, vector<missValue> &miss_vale_pos, bool byrow) {
    int i=0, j=0, k=0, n=0, m=0, nonmiss=0;
    int nrow=miss_mat.rows(), ncol=miss_mat.cols();
    VectorXi miss_vec;
    missValue mv_buf={0,""};
    
    if(byrow) {
        n = nrow; m = ncol;
    } else {
        m = nrow; n = ncol;
    }
    
    // initialize the variables
    miss_indx = (int*) malloc (n*sizeof(int));
    for( i=0; i<n; i++) miss_indx[i]=0;
    miss_vale_pos.clear();
    miss_vale_pos.push_back(mv_buf);
    
    ostringstream cvrt;
    
    for( i =0 ; i<n; i++) {
        byrow ? (miss_vec = miss_mat.row(i)) : (miss_vec = miss_mat.col(i));
        nonmiss = miss_vec.sum();
        if( nonmiss == m ) continue;
        else {
            cvrt.str("");
            mv_buf.num = m - nonmiss;
            for( j=0; j<m; j++) {
                if(miss_vec(j)) continue;
                else cvrt <<j<<"," ;
            }
            mv_buf.pos = cvrt.str();
            miss_vale_pos.push_back(mv_buf);
            miss_indx[i] = ++k;
        }
    }
    return miss_indx;
}

int *CommFunc::split_miss_pos(missValue miss_vale_pos, int *miss_sn) {
    int i=0, nmiss=miss_vale_pos.num;
    vector<string> indxstrbuf;
    
    miss_sn = (int*) malloc (nmiss*sizeof(int));
    StrFunc::split_string(miss_vale_pos.pos, indxstrbuf);
    for(i=0; i<nmiss; i++) {
        istringstream (indxstrbuf[i]) >> miss_sn[i];
    }
    return miss_sn;
}

// standardise the vector or matrix to z-score
eigenVector CommFunc::assignVale(eigenVector vec, int *indxbuf, int nmiss, double dbuf) {
    int i = 0, length = vec.size();
    if( nmiss > length ) cerr <<"Incorrect size of vector. CommFunc::assignVale"<<endl;
    
    for(i =0; i<nmiss; i++) {
        vec(indxbuf[i]) = dbuf;
    }
    return vec;
}

double CommFunc::eigenMean(eigenVector vec, int missindx, missValue miss_vale_pos) {
    int i=0, *indxbuf=0, nmiss=miss_vale_pos.num;
    int nonmiss = vec.size() - nmiss;
    
    // set the missing value to zero
    if(missindx) {
        indxbuf = CommFunc::split_miss_pos(miss_vale_pos, indxbuf);
        vec = CommFunc::assignVale(vec, indxbuf, nmiss);
    }
    return vec.sum()/(double)nonmiss;
}

double CommFunc::eigenSd(eigenVector vec, int missindx, missValue miss_vale_pos) {
    int i=0, *indxbuf=0, nmiss=miss_vale_pos.num, length=vec.size();
    int nonmiss = length - nmiss;
    double sdbuf = 0.0;
    
    // set the missing value to zero
    if(missindx) {
        indxbuf = CommFunc::split_miss_pos(miss_vale_pos, indxbuf);
        vec = CommFunc::assignVale(vec, indxbuf, nmiss);
    }
    double meanbuf = vec.sum()/nonmiss;
    // estimate the sd
    if(missindx) vec = CommFunc::assignVale(vec, indxbuf, nmiss, meanbuf);
    for(i=0; i<length; i++) {
        sdbuf += pow((vec(i) - meanbuf), 2);
    }
    sdbuf = sqrt(sdbuf/(nonmiss-1));
    
    return sdbuf;
}

eigenVector CommFunc::scale_vec(eigenVector vec, int missindx, missValue miss_vale_pos, bool onlymean) {
    int i = 0, j =0, *indxbuf=0;
    int length = vec.size(), nmiss = miss_vale_pos.num;
    int nonmiss = length - nmiss;
    missValue missvalebuf;
    
    // set the missing value to zero
    if(missindx) {
        indxbuf = CommFunc::split_miss_pos(miss_vale_pos, indxbuf);
        vec = CommFunc::assignVale(vec, indxbuf, nmiss);
    }
    
    // estimate the average
    double meanbuf = vec.sum()/(double)nonmiss;
    // estimate the sd
    if(missindx) vec = CommFunc::assignVale(vec, indxbuf, nmiss, meanbuf);
    double sdbuf = onlymean ? 1.0 : eigenSd(vec, 0, miss_vale_pos);
    // transformed to z-score
    for( i=0; i<length; i++) vec(i) = (vec(i) - meanbuf)/sdbuf;
    // set missing value as 1e10
    if(missindx) {
        vec = CommFunc::assignVale(vec, indxbuf, nmiss, 1e10);
    }
    return vec;
}

eigenMatrix CommFunc::scale(eigenMatrix mat, int *missindx, vector<missValue> miss_vale_pos, bool byrow, bool onlymean) {
    int i = 0, missnum=0, nrow=mat.rows(), ncol=mat.cols();
    eigenVector vecbuf;
    
    if(byrow) {
        for( i = 0; i < nrow; i++) {
            vecbuf = mat.row(i);
            missnum = missindx[i];
            vecbuf = CommFunc::scale_vec(vecbuf, missnum, miss_vale_pos[missnum], onlymean);
            mat.row(i) = vecbuf;
        }
    } else {
        for( i = 0; i< ncol; i++) {
            vecbuf = mat.col(i);
            missnum = missindx[i];
            vecbuf = CommFunc::scale_vec(vecbuf, missnum, miss_vale_pos[missnum], onlymean);
            mat.col(i) = vecbuf;
        }
    }
    return mat;
}