/*
 *
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * 2010-2017 by Jian Yang <jian.yang@uq.edu.au> and others
 *
 * Mock layer of version 2 by Zhili Zheng <zhilizheng@outlook.com>
 *
 * Options parser, and bridge unsupported flags to version 1.

 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is attached along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include "Logger.h"
#include "Geno.h"
#include "Marker.h"
#include "Pheno.h"
#include "GRM.h"
#include "ThreadPool.h"
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <chrono>
#include "main/option.h"

using std::bind;
using std::map;
using std::function;
using std::placeholders::_1;
using std::placeholders::_2;

void out_ver(bool flag_outFile){
    function<void (int, const string&)> log;
    if(flag_outFile){
        log = bind(&Logger::l, LOGGER_P, _1, _2);
    }else{
        log = bind(&Logger::m, LOGGER_P, _1, _2);
    }

    log(0, "*******************************************************************");
    log(0, "* Genome-wide Complex Trait Analysis (GCTA)");
    log(0, "* version 1.90 beta1");
    log(0, "* (C) 2010-2017, The University of Queensland");
    log(0, "* Please report bugs to: Jian Yang <jian.yang@uq.edu.au>");
    log(0, "*******************************************************************");
}

int main(int argc, char *argv[]){
    out_ver(false);
    auto start = std::chrono::steady_clock::now();
    map<string, vector<string>> options;
    vector<string> keys;
    string last_key = "";
    string cur_string = "";
   if(argc == 1){
        LOGGER.m(0, "The GCTA has lots of options, you can refer to our online document at URL http://cnsgenomics.com/software/gcta/");
        LOGGER.e(0, "No analysis has been launched by the option(s)");
    }
    for(int index = 1; index < argc; index++){
        cur_string = argv[index];
        if(cur_string.substr(0, 2) == "--"){
            last_key = cur_string;
            if(options.find(last_key) == options.end()){
                options[last_key] = {};
                keys.push_back(last_key);
            }else{
                LOGGER.e(0, "Find double options: " + cur_string);
            }
        }else{
            if(last_key == ""){
                LOGGER.e(0, "the first options isn't start with \"--\": " + cur_string);
            }else{
                options[last_key].push_back(cur_string);
            }
        }
    }

    int num_parts = 1;
    int cur_part = 1;
    bool is_parted = false;

    if(options.find("--part") != options.end()){
        num_parts = std::stoi(options["--part"][0]);
        cur_part = std::stoi(options["--part"][1]);
        is_parted = true;
    }

    options["parts"].push_back(std::to_string(num_parts));
    options["parts"].push_back(std::to_string(cur_part));

    //Find the --out options
    string log_name;
    if(options.find("--out") != options.end()){
        vector<string> outs = options["--out"];
        if(outs.size() == 0){
            LOGGER.e(0, "you might forget to specify the output filename in the \"--out\" option");
        }else{
            log_name = outs[0];
            if(is_parted){
                //log_name = outs[0] + ".parted_" + std::to_string(num_parts) + "_" + std::to_string(cur_part);
                std::string s_parts = std::to_string(num_parts);
                std::string c_parts = std::to_string(cur_part);
                log_name = outs[0] + ".parted_" + s_parts + "_" + std::string(s_parts.length() - c_parts.length(), '0') + c_parts;
            }
            options["out"].push_back(log_name);
        }
    }else{
        LOGGER.e(0, "No \"--out\" options find, you might forget to specify the output option");
    }

    // multi thread mode
    int thread_num = 1;
    int is_threaded = false;
    if(options.find("--thread-num") != options.end()){
        vector<string> thread_nums = options["--thread-num"];
        if(thread_nums.size() == 1) {
            try{
                thread_num = std::stoi(thread_nums[0]);
                is_threaded = true;
            }catch(std::invalid_argument&){
                LOGGER.w(0, "can't get thread number from --thread-num option");
            }
        }
    }


    //start register the options
    // Please take care of the order, C++ has few reflation feature, I did in a ugly way.
    vector<string> module_names = {"phenotype", "marker", "genotype", "genetic relationship matrix"};
    vector<bool (*)(map<string, vector<string>>&)> registers = {
            Pheno::registerOption,
            Marker::registerOption,
            Geno::registerOption,
            GRM::registerOption
    };
    vector<void (*)()> processMains = {
            Pheno::processMain,
            Marker::processMain,
            Geno::processMain,
            GRM::processMain
    };

    vector<int> mains;
    string out = "";
    for(int index = 0; index != module_names.size(); index++){
        if(registers[index](options)){
            mains.push_back(index);
        }
        out += module_names[index] + ", ";
    }
    if(mains.size() != 0){
        LOGGER.open(log_name + ".log");
        out_ver(true);
        // List all of the options
        LOGGER.i(0, "");
        LOGGER.i(0, "Options: ");
        LOGGER.i(0, " ");
        for(auto &key : keys){
            LOGGER << key << " ";
            for(auto & element : options[key]){
                LOGGER << element << " ";
            }
            LOGGER << std::endl;
        }
        LOGGER.i(0, "");

        if(mains.size() > 1) LOGGER.e(0, "multiple main functions are not available currently");
        if(is_threaded) {
            LOGGER.i(0, "The analysis will run in " + std::to_string(thread_num) + " threads");
            LOGGER.i(0, "Note: Some analysis will ignore the multi-thread mode");
        }
        ThreadPool *threadPool = ThreadPool::GetPool(thread_num - 1);
        processMains[mains[0]]();
    }else{
        try{
            option(argc, argv);
        }catch (const string &err_msg) {
            LOGGER.e(0, err_msg);         
        }catch (const char *err_msg) {
            LOGGER.e(0, string(err_msg));
        }
    }

    auto end = std::chrono::steady_clock::now();
    float duration = std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count();
    LOGGER.i(0, "Finished in " + std::to_string(duration) + " seconds");
}
