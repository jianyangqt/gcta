#include "LD.h"
#include "Logger.h"
#include "constants.hpp"
#include "utils.hpp"
#include "StatLib.h"
#include <numeric>
#include <string>
#include <functional>
#include <omp.h>
#include <mkl.h>
#include <cstdio>
#include <algorithm>

map<string, string> LD::options;
map<string, int> LD::options_i;
vector<string> LD::processFunctions;
bool LD::chr_ends = false;
unique_ptr<double[]> LD::geno_buffer[2];
int LD::cur_buffer = 0;
uint64_t LD::cur_buffer_offset[2] = {0, 0};
uint32_t LD::num_indi = 0;

LD::LD(Geno * geno){
    this->geno = geno;
    num_indi = geno->pheno->count_keep();
    ld_window = options_i["LD_window"] * 1000;
    is_r2 = false;
    cur_process_marker_index = 0;
    if(options["method"] == "r2"){
        is_r2 = true;
    }

    h_ld = fopen(options["out"].c_str(), "wb");
    if(!h_ld){
        LOGGER.e(0, "can't open " + options["out"] + " for writing.");
    }
}

LD::~LD(){
    fclose(h_ld);
}

template <typename T>
struct static_cast_func
{
  template <typename T1> 
  T operator()(const T1& x) const { return static_cast<T>(x); }
};

void LD::calcLD(){
    // the current buffer
    int cacl_index_buffer;
    if(geno_buffer[!cur_buffer]){
        cacl_index_buffer = !cur_buffer;
    }else{
        cacl_index_buffer = cur_buffer;
    }
    const char *trans = "T", *notrans = "N", *uplo = "L";
    double zero = 0.0;
    int nr = num_indi;
    int nc1 = cur_buffer_offset[cacl_index_buffer] / nr;
    double alpha = 1.0 / (nr - 1);
    double *ptr1 = geno_buffer[cacl_index_buffer].get();
    double *res1 = new double[nc1 * nc1];
    dsyrk(uplo, trans, &nc1, &nr, &alpha, ptr1, &nr, &zero, res1, &nc1);

    double *res2 = nullptr;
    // is previous buffer active?
    int nc2 = 0;
    if(geno_buffer[!cacl_index_buffer]){
        //TODO: reduce half calculation
        nc2 = cur_buffer_offset[!cacl_index_buffer] / nr;
        double *ptr2 = geno_buffer[!cacl_index_buffer].get();
        res2 = new double[nc2 * nc1];
        dgemm(trans, notrans, &nc2, &nc1, &nr, &alpha, ptr2, &nr, ptr1, &nr, &zero, res2, &nc2);
    }
    
    for(int i = 0; i < nc1; i++){
        uint32_t cur_size;
        float *buffer;
        if(res2){ 
            cur_size = geno->marker->getNextWindowSize(cur_process_marker_index, ld_window);
            buffer = new float[cur_size];
            uint32_t cur_size1 = nc1 - i;
            uint32_t cur_size2 = cur_size - cur_size1;
            double* temp_ptr = res1 + i + i*nc1;
            double* temp_ptr2 = res2 + i * nc2;
            std::transform(temp_ptr, temp_ptr + cur_size1, buffer, static_cast_func<float>());
            std::transform(temp_ptr2, temp_ptr2 + cur_size2, buffer + cur_size1, static_cast_func<float>());
        }else{
            cur_size = nc1 - i;
            buffer = new float[cur_size];
            double* temp_ptr = res1 + i + i * nc1;
            std::transform(temp_ptr, temp_ptr + cur_size, buffer, static_cast_func<float>());
        }
        // write to file
        if(fwrite(&cur_size, sizeof(uint32_t), 1, h_ld) != 1){
            LOGGER.e(0, "can't write to " + options["out"] + ".");
        }
        if(fwrite(buffer, sizeof(float), cur_size, h_ld) != cur_size){
            LOGGER.e(0, "can't write to " + options["out"] + ".");
        }
        delete[] buffer;
        cur_process_marker_index++;
    }
    //clean 
    if(res1){
        delete[] res1;
    }
    if(res2){
        delete[] res2;
    }
    geno_buffer[cacl_index_buffer].reset(nullptr);
    fflush(h_ld);
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
                options_i["LD_window"] = std::stoi(options_in[curFlag][0]);
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
            uint32_t window = options_i["LD_window"] * 1000;
            uint32_t total_num_marker = marker.count_extract(); 
            uint32_t cur_index_marker = 0;
            callBacks.push_back(std::bind(&Geno::freq64, &geno, _1, _2));
            callBacks.push_back(std::bind(&LD::readGeno, &ld, _1, _2));

            while(cur_index_marker < total_num_marker){
                cur_buffer = !cur_buffer;
                vector<uint32_t> indices1 = marker.getNextWindowIndex(cur_index_marker, window, chr_ends);
                geno_buffer[cur_buffer].reset(new double[indices1.size() * num_indi]);
                cur_buffer_offset[cur_buffer] = 0;
                geno.loop_64block(indices1, callBacks, false);
                cur_index_marker += indices1.size();

                // if another buffer is NA,  not chr ends and still in range;
                if((!geno_buffer[!cur_buffer]) && (!chr_ends)){
                    continue;
                }
                //ld.deduceLD();
                ld.calcLD();
                LOGGER.p(0, to_string_precision(cur_index_marker * 100.0 / total_num_marker, 2) + "% finished.");
            }
            //last block if not chr ends;
            if(geno_buffer[cur_buffer]){
                ld.calcLD();
            }
        }
    }
}

