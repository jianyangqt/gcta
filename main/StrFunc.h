/*
 * Interface to the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _STRFUNC_H
#define _STRFUNC_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
using namespace std;

namespace StrFunc
{
	bool str_within_quto(const string &str, string &str_buf);
	int split_string(const string &str, vector<string> &vec_str, string separator=" ,\t;\n");
	string first_string(const string &str, const char separator);
	string last_string(const string &str, const char separator);
	void to_upper(string &str);
	void to_lower(string &str);
	string get_sub_str(const string & rst, int pos);
	bool StrEqual(const string &StrA, const string &StrB, bool NoCaseSens=true);
	bool StrVecEqual(const vector<string> &VsBufA, const vector<string> &VsBufB, int Pos);

	// find a string in a string vector ignoring upper or lower case
	vector<string>::iterator find(vector<string> &target_vs, const string &target_str);

	// find a char in a string ignoring upper or lower case
	string::iterator find(string &target_str, const char target_ch);

	// go to the postion of a give string in a stream ignoring upper or lower case
	bool goto_str(std::istream &in_file, const string &str);

	// rewind a stream
	void rewind_if(std::istream &in_file);

	// match two vectors
	void match(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC);
        // compare two string ignore the case
        bool i_compare(string const& a, string const& b);
}

#endif
