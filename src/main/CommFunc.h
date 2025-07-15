/*
 * Interface to the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#ifndef _COMMFUNC_H
#define _COMMFUNC_H

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <complex>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <Eigen/StdVector>
using namespace std;

namespace CommFunc
{
	const double FloatErr=numeric_limits<double>::epsilon();
	double Abs(const double &x);
	double sum(const vector<double> &x);
    double mean(const vector<double> &x);
    double median(const vector<double> &x);
    double quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob);
    double quantile(const std::vector<double> &vals, double prob);
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
}

#endif
