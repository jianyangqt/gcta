/*
 * Implementations of the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "StrFunc.h"

int StrFunc::split_string(const string &str, vector<string> &vec_str, string separator)
{
	if(str.empty()) return 0;
	vec_str.clear();

	int i=0;
	bool look=false;
	string str_buf;
	string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
	string::size_type pos;

	for(i=0; i<separator.size(); i++){
		pos=symbol_pool.find(separator[i]);
		if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
	}

	for(i=0; i<str.size(); i++){
		if( symbol_pool.find(str[i])!=string::npos ){
			if(!look) look=true;
			str_buf += str[i];
		}
		else{
			if(look){
				look=false;
				vec_str.push_back(str_buf);
				str_buf.erase(str_buf.begin(), str_buf.end());
			}
		}
	}
	if(look) vec_str.push_back(str_buf);

	return vec_str.size();
}

string StrFunc::first_string(const string &str, const char separator)
{
	int pos=str.find(separator);
	if(pos!=-1) return string(str.begin(), str.begin()+pos);
	return string("");
}

string StrFunc::last_string(const string &str, const char separator)
{
	int pos=str.find_last_of(separator);
	if(pos!=-1) return string(str.begin()+pos+1, str.end());
	return string("");
}

void StrFunc::to_upper(string &str)
{
	int i=0;
	for(i=0; i<str.size(); i++){
		if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
	}
}

void StrFunc::to_lower(string &str)
{
	int i=0;
	for(i=0; i<str.size(); i++){
		if(str[i]>='A' && str[i]<='Z') str[i]-='A'-'a';
	}
}

string StrFunc::get_sub_str(const string & rst, int pos)
{
	vector<string> vs_buf;
	StrFunc::split_string(rst, vs_buf);
	return vs_buf[pos];
}

bool StrFunc::StrEqual(const string &StrA, const string &StrB, bool NoCaseSens)
{
	if(!NoCaseSens) return StrA==StrB;
	string StrBufA=StrA, StrBufB=StrB;
	to_upper(StrBufA);
	to_upper(StrBufB);
	return StrBufA==StrBufB;
}

bool StrFunc::StrVecEqual(const vector<string> &VsBufA, const vector<string> &VsBufB, int Pos)
{
	int SizeA=VsBufA.size(), SizeB=VsBufB.size();
	if(SizeA!=SizeB) return false;
	if(Pos>=SizeA) throw("Invalid Pos! StrFunc::StrVecEqual");

	int i=0;
	for(i=Pos; i<SizeA; i++){
		if(VsBufA[i]!=VsBufB[i]) return false;
	}

	return true;
}

bool StrFunc::str_within_quto(const string &str, string &str_buf)
{
	unsigned int begin=str.find_first_of("\"");
	unsigned int end=str.find_last_of("\"");
	if(begin==string::npos || end==string::npos || begin==end) return false;

	str_buf="";
	str_buf.insert(str_buf.begin(), str.begin()+begin+1, str.begin()+end);
	return true;
}

vector<string>::iterator StrFunc::find(vector<string> &target_vs, const string &target_str)
{
	string str_buf=target_str;
	vector<string> vs_buf=target_vs;

	int i=0;
	for(i=0; i<vs_buf.size(); i++) to_upper(vs_buf[i]);
	to_upper(str_buf);
	return target_vs.begin()+(std::find(vs_buf.begin(), vs_buf.end(), str_buf)-vs_buf.begin());
}

string::iterator StrFunc::find(string &target_str, const char target_ch)
{
	char ch_buf=target_ch;
	string str_buf=target_str;
	to_upper(str_buf);
	if(ch_buf>'a' && ch_buf<'z') ch_buf+='A'-'a';
	return target_str.begin()+(std::find(str_buf.begin(), str_buf.end(), ch_buf)-str_buf.begin());
}

bool StrFunc::goto_str(std::istream &in_file, const string &str)
{
	string str_buf;
	string query_str=str;
	vector<string> vs_buf;
	StrFunc::to_upper(query_str);
	while(in_file>>str_buf){
		if( StrFunc::split_string(str_buf, vs_buf)>0 ) str_buf=vs_buf[0];
		else continue;
		StrFunc::to_upper(str_buf);
		if(str_buf=="#") { getline(in_file, str_buf); continue; }
		if(str_buf==query_str) return true;
	}

	return false;
}

void StrFunc::rewind_if(std::istream &in_file)
{
	in_file.clear(ios::goodbit);
	in_file.seekg(ios::beg);
}

void StrFunc::match(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC)
{
    int i=0;
    map<string, int> id_map;
    map<string, int>::iterator iter;
    VecC.clear();
    for(i=0; i<VecB.size(); i++) id_map.insert(pair<string,int>(VecB[i], i));
    for(i=0; i<VecA.size(); i++){
        iter=id_map.find(VecA[i]);
        if(iter==id_map.end()) VecC.push_back(-9);
        else VecC.push_back(iter->second);
    }
}

