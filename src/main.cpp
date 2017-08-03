/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Options parser

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
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

using std::bind;
using std::map;
using std::function;
using std::placeholders::_1;
using std::placeholders::_2;

int main(int argc, char *argv[]){
    auto start = std::chrono::steady_clock::now();
    map<string, vector<string>> options;
    vector<string> keys;
    string last_key = "";
    string cur_string = "";
    if(argc == 1){
        LOGGER.i(0, "The GCTA has lots of options, you can refer to our online document at URL http://cnsgenomics.com/software/gcta/");
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
    if(options.find("--out") != options.end()){
        vector<string> outs = options["--out"];
        if(outs.size() == 0){
            LOGGER.e(0, "you might forget to specify the output filename in the \"--out\" option");
        }else{
            string log_name = outs[0];
            if(is_parted){
                //log_name = outs[0] + ".parted_" + std::to_string(num_parts) + "_" + std::to_string(cur_part);
                std::string s_parts = std::to_string(num_parts);
                std::string c_parts = std::to_string(cur_part);
                log_name = outs[0] + ".parted_" + s_parts + "_" + std::string(s_parts.length() - c_parts.length(), '0') + c_parts;
            }
            options["out"].push_back(log_name);
            LOGGER.open(log_name + ".log");
        }
    }else{
        LOGGER.e(0, "No \"--out\" options find, you might forget to specify the output option");
    }

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

    // multi thread mode
    int thread_num = 4;
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

    if(is_threaded) {
        LOGGER.i(0, "The analysis will run in " + std::to_string(thread_num) + " threads");
    }else{
        LOGGER.i(0, "The analysis will run in 4 threads to calculate in default; You can specify --thread-num number to "
                "change this behavior");
    }
    LOGGER.i(0, "Note: Some analysis will ignore the multi-thread mode");
    ThreadPool *threadPool = ThreadPool::GetPool(thread_num - 1);


    //start register the options
    // Please take care the order, C++ has few reflation feature, I did in a ugly way.
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
    if(mains.size() > 1){
        LOGGER.e(0, "multiple main functions are" + out + "not available currently");
    }else if(mains.size() == 1){
        processMains[mains[0]]();
    }else{
        LOGGER.e(0, "No main function specified");
    }

    auto end = std::chrono::steady_clock::now();
    float duration = std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count();
    LOGGER.i(0, "Finished in " + std::to_string(duration) + " seconds");

}
