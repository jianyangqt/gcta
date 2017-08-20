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
#include <algorithm>
#include "main/option.h"
#include "utils.hpp" 

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
    log(0, "* version 1.90 beta1.02");
    log(0, "* (C) 2010-2017, The University of Queensland");
    log(0, "* Please report bugs to: Jian Yang <jian.yang@uq.edu.au>");
    log(0, "*******************************************************************");
    log(0, "Analysis started: " + getLocalTime());
    log(0, "Hostname: " + getHostName());
    log(0, "");
}

int main(int argc, char *argv[]){
    out_ver(false);
    auto start = std::chrono::steady_clock::now();
    vector<string> supported_flagsV2 = {"--bfile", "--bim", "--fam", "--bed", "--keep", "--remove", 
        "--chr", "--autosome-num", "--autosome", "--extract", "--exclude", "--maf", "--max-maf", 
        "--freq", "--out", "--make-grm", "--make-grm-part", "--thread-num", "--grm",
        "--grm-cutoff", "--grm-no-relative", "--cutoff-detail", "--make-fam"};
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

     //Find the --out options
    if(options.find("--out") != options.end()){
        vector<string> outs = options["--out"];
        if(outs.size() == 0){
            LOGGER.e(0, "you might forget to specify the output filename in the \"--out\" option");
        }else{
            options["out"].push_back(outs[0]);
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
                LOGGER.e(0, "can't get thread number from --thread-num option");
            }
        }
    }

    map<string, vector<string>> options_total = options;

    //start register the options
    // Please take care of the order, C++ has few reflation feature, I did in a ugly way.
    vector<string> module_names = {"phenotype", "marker", "genotype", "genetic relationship matrix"};
    vector<int (*)(map<string, vector<string>>&)> registers = {
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
    for(int index = 0; index != module_names.size(); index++){
        int num_reg = registers[index](options);
        if(num_reg == 1){
            mains.push_back(index);
        }else if(num_reg > 1){
            LOGGER.e(0, "multiple main functions are invalid currently");
        }
    }
    bool unKnownFlag = false;
    for(string &key : keys){
        if(std::find(supported_flagsV2.begin(), supported_flagsV2.end(), key) == supported_flagsV2.end()){
            unKnownFlag = true;
            break;
        }
    }
    if(mains.size() != 0 && !unKnownFlag){
        LOGGER.open(options["out"][0] + ".log");
        out_ver(true);
        // List all of the options
        LOGGER.i(0, "Options: ");
        LOGGER.i(0, " ");
        for(auto &key : keys){
            LOGGER << key << " ";
            for(auto & element : options_total[key]){
                LOGGER << element << " ";
            }
            LOGGER << std::endl;
        }
        LOGGER.i(0, "");

        if(mains.size() > 1) LOGGER.e(0, "multiple main functions are invalid currently");
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

    LOGGER.i(0, "");

    auto end = std::chrono::steady_clock::now();
    float duration = std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count();
    LOGGER.i(0, "Analysis finished: " + getLocalTime());
    LOGGER.i(0, "Computational time: " + std::to_string(duration) + " seconds");
}
