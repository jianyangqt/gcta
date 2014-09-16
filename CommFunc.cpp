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

// standardise the vector or matrix to z-score
double CommFunc::eigenMean(eigenVector vec, int nonmiss) {
    int n = vec.size();
    
    // set the missing value to zero
    vec = (vec.array() < 1e9).select(vec,0);

    return vec.sum()/(double)nonmiss;
}

double CommFunc::eigenSd(eigenVector vec, double xmean, int nonmiss ) {
    int n = vec.size();
    
    // set the missing value to average
    vec.array() = (vec.array() < 1e9).select(vec,xmean);
    vec.array() = vec.array() - xmean;
    return sqrt(vec.dot(vec) / ((double)nonmiss-1));
}

eigenVector CommFunc::scale_vec(eigenVector vec, int nonmiss, bool onlymean) {
    double meanbuf = CommFunc::eigenMean(vec, nonmiss);
    double sdbuf = onlymean ? 1.0 : CommFunc::eigenSd(vec, meanbuf, nonmiss);
    
    vec.array() = ( vec.array() - meanbuf) / sdbuf;
    // set missing value to 1e10
    double crt=1e9/sdbuf;
    vec = (vec.array() < crt ).select(vec, 1e10);
    
    return vec;
}

eigenMatrix CommFunc::scale(eigenMatrix mat, VectorXi nonmiss, bool byrow, bool onlymean) {
    int i = 0, nrow=mat.rows(), ncol=mat.cols();
    eigenVector vecbuf;
    
    if(byrow) {
        for( i = 0; i < nrow; i++) {
            mat.row(i).array() = CommFunc::scale_vec(mat.row(i).array(), nonmiss(i), onlymean);
        }
    } else {
        for( i = 0; i< ncol; i++) {
            mat.col(i).array() = CommFunc::scale_vec(mat.col(i).array(), nonmiss(i), onlymean);
        }
    }
    return mat;
}