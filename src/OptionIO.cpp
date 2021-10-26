#include "OptionIO.h"
#include "Logger.h"
#include <iostream>
#include <sstream>
#include <set>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <string>

using std::to_string;

void addOneFileOption(string key_store, string append_string, string key_name, map<string, vector<string>> options_in, map<string,string>& options) {
    if(options_in.find(key_name) != options_in.end()){
        if(options_in[key_name].size() == 1){
            options[key_store] = options_in[key_name][0] + append_string;
        }else if(options_in[key_name].size() > 1){
            options[key_store] = options_in[key_name][0] + append_string;
            LOGGER.w(0, "There are multiple " + key_name + ". Only the first item will be used in the analysis." );
        }else{
            LOGGER.e(0, "no " + key_name + " parameter found");
        }
        if(!checkFileReadable(options[key_store])){
            LOGGER.e(0, key_name + " " + options[key_store] + " not found");
        }
    }
}

void addMFileListsOption(string key_store, string append_string, string key_name, map<string, vector<string>> options_in, map<string,string>& options){
    string list_filename;
    if(options_in.find(key_name) != options_in.end()){
        if(options_in[key_name].size() == 1){
            list_filename = options_in[key_name][0];
        }else if(options_in[key_name].size() > 1){
            list_filename = options_in[key_name][0];
            LOGGER.w(0, "There are multiple " + key_name + ". Only the first item will be used in the analysis." );
        }else{
            LOGGER.e(0, "no " + key_name + " parameter found");
        }

        vector<string> m_files;
        std::ifstream m_file(list_filename);
        string line;
        while(std::getline(m_file, line)){
            std::istringstream iss(line);
            string temp_fname;
            if (!(iss >> temp_fname)){
                continue;
            }
            boost::replace_all(temp_fname, "\r", "");
            m_files.push_back(temp_fname);
        }
        m_file.close();

        if(m_files.size() > 0){
            std::set<string> unique_mfiles;
            auto ends = std::remove_if(m_files.begin(), m_files.end(),
                    [&unique_mfiles](const string& item){
                        if(unique_mfiles.find(item) != unique_mfiles.end()){
                            return true;
                        }
                        unique_mfiles.insert(item);
                        return false;
                    });
            m_files.erase(ends, m_files.end());
            //LOGGER.i(0, to_string(m_files.size()) + " unique items in [" + list_filename + "].");

            vector<string> mc_files;
            for(auto &temp_fname : m_files){
                string t2_fname;
                if(temp_fname.size() < append_string.size()){
                    t2_fname = temp_fname + append_string;
                }else{
                    string ext = temp_fname.substr(temp_fname.size() - append_string.size());
                    if(ext == append_string){
                        t2_fname = temp_fname;
                    }else{
                        t2_fname = temp_fname + append_string;
                    }
                }
                mc_files.push_back(t2_fname);
            }

            options[key_store] = boost::algorithm::join(mc_files, "\t");
        }else{
            LOGGER.e(0, "no item in [" + list_filename + "].");
        }
    }
}

// get the header line
// return true:   success, head empty, means no head line; ifstream will be set to the first line of contents 
//        false:  fail to get header. ncol = 0: all comment header, or contains blank line in the header;
//                                    ncol != 0: head is inconsistent with contents
// 
bool getListHeaderNew(std::ifstream *input, bool &ioerr, int &ncol, vector<string> &head, int &nHeader){
    int startByte = 0;
    ioerr = false;
    ncol = 0;
    string line;
    uint64_t prePos = 0;
    string preLine;
    if(!input->good()){
        ioerr = true;
        return false;
    }
    std::size_t found;
    int curReadLine = 0;
    while(std::getline(*input, line)){
        curReadLine++;
        found = line.find_first_not_of("\r\t ");
        if(found != std::string::npos && line.at(found) == '#'){
            preLine = line;
            prePos = input->tellg();
            continue;
        }else{
            break;
        }
    }

    if(line.at(found) == '#' || found == std::string::npos){
        nHeader = curReadLine;
        return false;
    }

    vector<string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    int ncol_con = line_elements.size();

    ncol = ncol_con;
    startByte = prePos;
    nHeader = curReadLine - 1;

    if(!preLine.empty()){
        vector<string> line_elements_head;
        boost::split(line_elements_head, preLine, boost::is_any_of("\t "));
        int ncol_head = line_elements_head.size();
        boost::replace_all(line_elements_head[ncol_head - 1], "\r", "");
        if(ncol_head != ncol_con){
            return false;
        }else{
            head = line_elements_head;
        }
    }
    input->clear();
    input->seekg(startByte);
    return true;
}

bool readTxtList(string fileName, int minFields, vector<string> &head, int &nHeader, map<int, vector<string>> &lists){
    std::ifstream input(fileName.c_str());
    int ncol;
    bool ioerr;
    if(getListHeaderNew(&input, ioerr, ncol, head, nHeader)){
        if(ncol < minFields){
            LOGGER.e(0, "the number of fields in [" + fileName + "] is less than " + to_string(minFields) + ".");
        }
        //init the list;
        for(int i = 0; i < ncol; i++){
            lists[i] = {};
        }
        string line;
        int line_number = nHeader;
        while(std::getline(input, line)){
            line_number++;
            vector<string> line_elements;
            boost::split(line_elements, line, boost::is_any_of("\t "));
            boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
            if(line_elements.size() < ncol){
                LOGGER.e(0, "the file [" + fileName + "] contains different number of elements in line " + to_string(line_number) + ".");
            }else if(line_elements.size() > ncol){
                LOGGER.w(0, "the file [" + fileName + "] contains extra number of elements in line " + to_string(line_number) + ".");
            }
            for(int i = 0; i < ncol; i++){
                lists[i].emplace_back(line_elements[i]);
            }
        }
        input.close();
        return true;
    }else{
        if(ioerr){
            LOGGER.e(0, "can't read [" + fileName + "].");
        }
        if(ncol == 0){
            LOGGER.e(0, "blank header line or no valid text data.");
        }else{
            LOGGER.e(0, "data are inconsistent with the headers.");
        }
        input.close();
        return false;
    }
}

bool checkFileReadable(string filename){
    std::ifstream f(filename.c_str());
    return f.good();
}

uint64_t getFileSize(FILE * file){
   if(file == NULL){
       return 0;
   }
   fseek(file, 0, SEEK_END);
   uint64_t f_size = ftell(file);
   rewind(file);
   return f_size;
}




