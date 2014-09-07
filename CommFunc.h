/*
 * Interface to the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _COMMFUNC_H
#define _COMMFUNC_H

#include <limits>
#include <complex>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include "GlobalVarb.h"
#include "StrFunc.h"

using namespace Eigen;
using namespace std;

namespace CommFunc
{
	const double FloatErr=numeric_limits<double>::epsilon();
	double Abs(const double &x);
	double sum(const vector<double> &x);
    double mean(const vector<double> &x);
    double median(const vector<double> &x);
    double var(const vector<double> &x);
    double cov(const vector<double> &x, const vector<double> &y);
	bool FloatEqual(double lhs, double rhs);
	bool FloatNotEqual(double lhs, double rhs);
	const double Sqr(const double &a);
	const double Max(const double &a, const double &b);
	const double Min(const double &a, const double &b);
	const double Sign(const double &a, const double &b);
	int rand_seed(); //positive value, return random seed using the system time
    void FileExist(string filename);
    
    // find the missing value
    VectorXi find_missing_vec(eigenVector vec);
    MatrixXi find_missing_X(eigenMatrix mat);
    int *find_missing_pos(MatrixXi miss_mat, int *miss_indx, vector<missValue> &miss_vale_pos, bool byrow=false);
    int *split_miss_pos(missValue miss_vale_pos, int *miss_sn);
    // standardise the vector or matrix to z-score
    eigenVector assignVale(eigenVector vec, int *indxbuf, int nmiss, double dbuf=0.0) ;
    double eigenMean(eigenVector vec, int missindx, missValue miss_vale_pos);
    double eigenSd(eigenVector vec, int missindx, missValue miss_vale_pos);
    eigenVector scale_vec(eigenVector vec, int missindx, missValue miss_vale_pos, bool onlymean=false);
    eigenMatrix scale(eigenMatrix mat, int *missindx, vector<missValue> miss_vale_pos, bool byrow=false, bool onlymean=false);
}

#endif
