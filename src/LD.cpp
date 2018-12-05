#include "LD.h"
#include "Logger.h"
#include "constants.hpp"
#include "utils.hpp"
#include "StatLib.h"
#include <numeric>
#include <string>
#include <functional>

map<string, string> LD::options;
map<string, int> LD::options_i;
vector<string> LD::processFunctions;
bool chr_ends;
unique_ptr<double[]> geno_buffer[2];
int cur_buffer = 0;
uint64_t cur_buffer_offset[2] = {0, 0};

LD::LD(Geno * geno){
    this->geno = geno;
    num_indi = geno->pheno->count_keep();
}

void LD::readGeno(uint64_t *buf, int num_marker){
    double * ptr = geno_buffer[cur_buffer].get();
    uint64_t cur_offset = cur_buffer_offset[cur_buffer];
    ptr += cur_offset;

    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < num_marker; i++){
        geno->makeMarkerX(buf, i, ptr + i * num_indi, true, true);
    }
    cur_buffer_offset[cur_buffer] += (uint64_t) num_marker * num_indi;
}


int LD::registerOption(map<string, vector<string>>& options_in){
    int ret_val = 0;
    options["out"] = options_in["out"][0] + ".zld";

    string curFlag = "--ld-matrix";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("ld_matrix");
        ret_val++;
        options_in.erase(curFlag);
    }

    // get the calculation methods, defalut LD r;
    options["method"] = "r";
    curFlag = "--r";
    if(options_in.find(curFlag) != options_in.end()){
        options["method"] = "r";
        options_in.erase(curFlag);
    }

    curFlag = "--r2";
    if(options_in.find(curFlag) != options_in.end()){
        options["method"] = "r2";
        options_in.erase(curFlag);
    }

    // set the LD window size Kb
    options_i["LD_window"] = 1000;
    curFlag = "--ld-wind";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            try{
                options_i["ld_window"] = std::stoi(options_in[curFlag][0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "LD window is not an integer.");
            }
        }
        options_in.erase(curFlag);
    }


    return ret_val;
}

void LD::processMain(){
    vector<function<void (uint64_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "ld_matrix"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            LD ld(&geno);

            LOGGER.i(0, "Generating LD matrix...");
            uint32_t window = options_i["ld_window"] * 1000;
            uint32_t total_num_marker = marker.count_extract(); 
            uint32_t cur_index_marker = 0;
            callBacks.push_back(std::bind(&Geno::freq64, &geno, _1, _2));
            callBacks.push_back(std::bind(&LD::readGeno, &ld, _1, _2));

            while(cur_index_marker < total_num_marker){
                cur_buffer = !cur_buffer;
                vector<uint32_t> indices1 = marker.getNextWindowIndex(cur_index_marker, window, chr_ends);
                geno_buffer[cur_buffer].reset(new double[indices1.size() * num_indi]);
                cur_buffer_offset[cur_buffer] = 0;
                geno.loop_64block(indices1, callBacks);
                cur_index_marker += indices1.size();

                if((!chr_ends) && cur_index_marker < total_num_marker){
                    cur_buffer = !cur_buffer;
                    vector<uint32_t> indices2 = marker.getNextWindowIndex(cur_index_marker, window, chr_ends);
                    geno_buffer[cur_buffer].reset(new double[indices2.size() * num_indi]);
                    cur_buffer_offset[cur_buffer] = 0;
                    geno.loop_64block(indices2, callBacks);
                    cur_index_marker += indices2.size();
                }else{
                    geno_buffer[!cur_buffer].reset(nullptr);
                }
                //ld.deduceLD();
            }
        }
    }
}

