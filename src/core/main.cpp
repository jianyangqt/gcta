/*
 *
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * 2010-present Jian Yang <jian.yang.qt@gmail.com> and others
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
#include "Covar.h"
#include "FastFAM.h"
#include "LD.h"
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "main/option.h"
#include "utils.hpp" 
#include <omp.h>
#include <cstdlib>
#include "mem.hpp"
#include "config.h"

#include "acat.hpp" //add by fanghl 20210303
#define ACAT_MODULE

using std::bind;
using std::map;
using std::to_string;
using std::function;
using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

void out_ver(bool flag_outFile){
    function<void (int, const string&, const string&)> log;
    if(flag_outFile){
        log = bind(&Logger::l, LOGGER_P, _1, _2, _3);
    }else{
        log = bind(&Logger::m, LOGGER_P, _1, _2, _3);
    }

    log(0, "*******************************************************************", "");
    log(0, "* Genome-wide Complex Trait Analysis (GCTA)", "");
    log(0, std::string("* version ") + std::string(GCTA_VERSION) + std::string(" ") + getOSName(), "");
#if (defined __linux__ || defined __linux) && defined __GNUC__
    std::string compile_info("");
    std::ostringstream outstring;
    outstring << "* Built at " << __DATE__ << " " << __TIME__ << ", "
        << "by GCC " << __GNUC__ << "." << __GNUC_MINOR__;
    compile_info = outstring.str();
    log(0, compile_info, "");
#endif
    log(0, "* (C) 2010-present, Yang Lab, Westlake University", "");
    log(0, "* Please report bugs to Jian Yang <jian.yang@westlake.edu.cn>", "");
    log(0, "*******************************************************************", "");
    log(0, string("at ") + getLocalTime() + ".", "Analysis started");
    log(0, "Hostname: " + getHostName(), "");
    log(0, "", "");
}

int main(int argc, char *argv[]){
    out_ver(false);
    LOGGER.ts("main");
    vector<string> supported_flagsV2 = {"--test-covar",
        "--bfile", "--bim", "--fam", "--bed", "--keep", "--remove", 
        "--chr", "--autosome-num", "--autosome", "--extract", "--exclude", "--maf", "--max-maf", 
        "--freq", "--out", "--make-grm", "--make-grm-part", "--thread-num", "--threads", "--grm",
        "--grm-cutoff", "--grm-singleton", "--cutoff-detail", "--make-bK-sparse", "--make-bK", "--pheno",
        "--mpheno", "--ge", "--fastGWA", "--fastGWA-mlm", "--fastGWA-mlm-exact", "--fastGWA-lr", "--save-fastGWA-mlm-residual", "--grm-sparse", "--qcovar", "--covar", "--rcovar", "--covar-maxlevel", "--make-grm-d", "--make-grm-d-part",
        "--cg", "--ldlt", "--llt", "--pardiso", "--tcg", "--lscg", "--save-inv", "--load-inv",
        "--update-ref-allele", "--update-freq", "--update-sex", "--mbfile", "--freqx", "--make-grm-xchr", "--make-grm-xchr-part", "--dc", "--make-grm-alg",
        "--make-bed", "--recodet", "--sum-geno-x", "--sample", "--bgen", "--mbgen", "--hard-call-thresh", "--dosage-call", "--dosage", "--mgrm", "--unify-grm", "--rel-only", 
        "--ld-matrix", "--r", "--ld-wind", "--r2", "--subtract-grm", "--save-pheno", "--save-bin", "--no-marker", "--joint-covar", "--sparse-cutoff", "--noblas", "--fastGWA-gram",
        "--inv-t1", "--est-vg", "--force-gwa", "--reml-detail", "--h2-limit", "--gwa-no-constrain", "--verbose", "--c-inf", "--c-inf-no-filter", "--geno", "--info", "--nofilter",
        "--set-list", "--burden",
        "--pfile", "--bpfile", "--mpfile", "--mbpfile", "--model-only", "--load-model", "--seed", "--fastGWA-mlm-binary", "--num-vec", "--trace-exact", "--cv-threshold", "--tao-start",
        "--acat", "--gene-list", "--snp-list", "--min-mac", "--max-maf", "--wind",
        "--envir", "--optimal-rho", "--noSandwich", "--grid-size",
    };
    map<string, vector<string>> options;
    vector<string> keys;
    string last_key = "";
    string cur_string = "";
   if(argc == 1){
        LOGGER.m(0, "Error: no analysis has been launched by the option(s)");
        LOGGER.m(0, "Please see online documentation at https://yanglab.westlake.edu.cn/software/gcta/");
        return 1;
    }
   #ifdef ACAT_MODULE //fanghl
    string opt_acat = "--acat";
    string opt_gene_list = "--gene-list";
    string arg_gene_list = "";
    string opt_snp_list = "--snp-list";
    string arg_snp_list = "";
    string opt_max_af = "--max-maf";
    double arg_max_af = 0.01;
    string opt_min_sample = "--min-mac";
    unsigned int arg_min_sample = 20;
    string opt_extend_len = "--wind";
    unsigned int arg_extend_len = 0;
    string opt_out = "--out";
    string arg_out = "acat_res.csv";
    int ACAT_call = 0;

    for (int j = 1; j < argc; j ++){
        if (argv[j] == opt_acat){
            ACAT_call = 1;
        } else if ((j + 1) < argc && string(argv[j + 1]).substr(0, 2) != "--"){
            if (argv[j] == opt_gene_list)
                arg_gene_list = string(argv[j + 1]);
            else if (argv[j] == opt_snp_list)
                arg_snp_list = string(argv[j + 1]);
            else if (argv[j] == opt_max_af)
                arg_max_af = std::stof(argv[j + 1]);
            else if (argv[j] == opt_min_sample)
                arg_min_sample = std::stoi(argv[j + 1]);
            else if (argv[j] == opt_extend_len)
                arg_extend_len = std::stoi(argv[j + 1]);
            else if (argv[j] == opt_out)
                arg_out = string(argv[j + 1]);
        }
    }

    if (ACAT_call && arg_gene_list != "" && arg_snp_list != "" ){
        const char * gene_list_file_name = arg_gene_list.c_str();
        const char * snp_list_file_name = arg_snp_list.c_str();
        const char * arg_out_file = arg_out.c_str();
        acat_func(gene_list_file_name, snp_list_file_name, arg_max_af, arg_min_sample, arg_extend_len, arg_out_file);
        exit(0);
        }
  #endif
    for(int index = 1; index < argc; index++){
        cur_string = argv[index];
        if(cur_string.substr(0, 2) == "--"){
            last_key = cur_string;
            if(options.find(last_key) == options.end()){
                options[last_key] = {};
                keys.push_back(last_key);
            }else{
                LOGGER.e(0, "invalid multiple options: " + cur_string);
            }
        }else{
            if(last_key == ""){
                LOGGER.e(0, "an option must start with \"--\": " + cur_string);
            }else{
                options[last_key].push_back(cur_string);
            }
        }
    }

     //Find the --out options
    if(options.find("--out") != options.end()){
        vector<string> outs = options["--out"];
        if(outs.size() == 0){
            LOGGER.e(0, "no output file name is specified in the \"--out\" option");
        }else{
            options["out"].push_back(outs[0]);
        }
    }else{
        LOGGER.e(0, "the --out option is missing.");
    }

    // multi thread mode
    int thread_num = 1;
    int is_threaded = false;
    char *thread_char = getenv("OMP_NUM_THREADS");
    try{
        if(thread_char){
            string thread_string(thread_char);
            thread_num = std::stoi(thread_string);
            is_threaded = true;
        }
    }catch(std::invalid_argument&){
        ;
    }

    bool is_thread_set = false;
    if(options.find("--thread-num") != options.end()){
        vector<string> thread_nums = options["--thread-num"];
        if(thread_nums.size() == 1) {
            try{
                thread_num = std::stoi(thread_nums[0]);
                is_threaded = true;
                is_thread_set = true;
            }catch(std::invalid_argument&){
                LOGGER.e(0, "can't get the thread number from the --thread-num option.");
            }
        }else{
            LOGGER.e(0, "can't set the number of threads.");
        }
    }

    if(options.find("--threads") != options.end()){
        vector<string> thread_nums = options["--threads"];
        if(is_thread_set){
            LOGGER.e(0,"can't set the number of threads using both --thread-num and --threads");
        }
        if(thread_nums.size() == 1) {
            try{
                thread_num = std::stoi(thread_nums[0]);
                is_threaded = true;
            }catch(std::invalid_argument&){
                LOGGER.e(0, "can't get the number of threads from --threads option.");
            }
        }else{
            LOGGER.e(0, "can't set the number of threads.");
        }
    }

    #ifdef _WIN32
       _putenv_s("OMP_NUM_THREADS", to_string(thread_num).c_str());
    #elif defined __linux__ || defined __APPLE__
        setenv("OMP_NUM_THREADS", to_string(thread_num).c_str(), 1);
    #else
        #error Only Windows, Mac and Linux are supported.
    #endif
    omp_set_num_threads(thread_num);


    // end thread;
    map<string, vector<string>> options_total = options;

    //start register the options
    // Please take care of the order, C++ has few reflation feature, I did in a ugly way.
    vector<string> module_names = {"phenotype", "marker", "genotype", "covar", "GRM", "fastFAM", "LD"};
    vector<int (*)(map<string, vector<string>>&)> registers = {
            Pheno::registerOption,
            Marker::registerOption,
            Geno::registerOption,
            Covar::registerOption,
            GRM::registerOption,
            FastFAM::registerOption,
            LD::registerOption
    };
    vector<void (*)()> processMains = {
            Pheno::processMain,
            Marker::processMain,
            Geno::processMain,
            Covar::processMain,
            GRM::processMain,
            FastFAM::processMain,
            LD::processMain
    };

    vector<int> mains;
    for(int index = 0; index != module_names.size(); index++){
        int num_reg = registers[index](options);
        if(num_reg == 1){
            mains.push_back(index);
        }else if(num_reg > 1){
            LOGGER.e(0, "multiple main functions are not supported currently.");
        }
    }
    bool unKnownFlag = false;
    for(string &key : keys){
        if(std::find(supported_flagsV2.begin(), supported_flagsV2.end(), key) == supported_flagsV2.end()){
            unKnownFlag = true;
            break;
        }
    }
    // other cases of not handle;
    // can't handle --mgrm --make-grm currently.
    if(std::find(keys.begin(), keys.end(), "--make-grm") != keys.end() && 
            std::find(keys.begin(), keys.end(), "--mgrm") != keys.end()){
        unKnownFlag = true;
    }
    bool bShowAll = false;
    if(options.find("--verbose") != options.end()){
        bShowAll = true;
        options.erase("--verbose");
    }
 

    LOGGER.open(options["out"][0] + ".log");
    out_ver(true);
    if(mains.size() != 0 && !unKnownFlag){
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

        if(mains.size() > 1) LOGGER.e(0, "multiple main functions are not supported currently.");
        if(is_threaded) {
            LOGGER.i(0, "The program will be running with up to " + std::to_string(thread_num) + " threads.");
        }
        //ThreadPool *threadPool = ThreadPool::GetPool(thread_num - 1);
        //avoid auto parallel

        processMains[mains[0]]();
    }else{
        try {
            option(argc, argv);
        } catch (const string &err_msg) {
            LOGGER.e(0, err_msg);
        } catch (const char *err_msg) {
            LOGGER.e(0, string(err_msg));
        }
    }

    LOGGER.i(0, "");

    LOGGER.i(0,  "at " + getLocalTime(), "Analysis finished");
    float duration = LOGGER.tp("main");
    int hours = (int) duration / 3600;
    string time_str = (hours == 0) ? "" : (to_string(hours) + " hour" + ((hours == 1) ? " ": "s "));
    int mins = (int) (duration - 3600 * hours) / 60;
    time_str += (mins == 0) ? "" : (to_string(mins) + " minute" + ((mins == 1) ? " ": "s "));
    float seconds = duration - 3600 * hours - 60 * mins;
    time_str = time_str + ((time_str == "") ? to_string_precision(seconds,2) : to_string((int)seconds)) + " sec"; 

    #ifdef __linux__
    if(bShowAll){
        float vmem = roundf(1000.0 * getVMPeakKB() / 1024/1024) / 1000.0; 
        float mem = roundf(1000.0 * getMemPeakKB() / 1024/1024) / 1000.0;     
        LOGGER.setprecision(3);
        LOGGER << "Peak memory: " << mem << " GB; Virtual memory: " << vmem << " GB." << std::endl;
        //system("cat /proc/$PPID/status");
    }
    #endif
 
    LOGGER << "Overall computational time: " << time_str  <<  "." << std::endl;
    return EXIT_SUCCESS;
}
