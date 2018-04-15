#include "OptionIO.h"
#include "Logger.h"
#include <iostream>

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

bool checkFileReadable(string filename){
    std::ifstream f(filename.c_str());
    return f.good();
}



