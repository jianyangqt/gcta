/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: read and process genotype of plink format in block way.

   Depends on the class of marker and phenotype

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#include "Geno.h"
#include "constants.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <chrono>
#include <ctime>
#include <iostream>
#include <iterator>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "utils.hpp"

using std::thread;
using std::to_string;

map<string, string> Geno::options;
map<string, double> Geno::options_d;
vector<string> Geno::processFunctions;

Geno::Geno(Pheno* pheno, Marker* marker) {
    if(options.find("geno_file") != options.end()){
        bed_file = options["geno_file"];
    }else{
        LOGGER.e(0, "No genotype file specified");
    }

    this->pheno = pheno;
    this->marker = marker;
    num_byte_per_marker = (pheno->count_raw() + 3) / 4;
    num_byte_buffer = num_byte_per_marker * Constants::NUM_MARKER_READ;
    last_byte_NA_sample = (4 - (pheno->count_raw() % 4)) % 4;
    check_bed();

    init_AF();
    init_AsyncBuffer();
    filter_MAF();
}

void Geno::filter_MAF(){
    if((options_d["min_maf"] != 0.0) || (options_d["max_maf"] != 0.5)){
        LOGGER.i(0, "Computing allele frequencies...");
        vector<function<void (uint8_t *, int)>> callBacks;
        callBacks.push_back(bind(&Geno::freq, this, _1, _2));
        loop_block(callBacks);
        // We adopt the EPSILON from plink, because the double value may have some precision issue;
        double min_maf = options_d["min_maf"] * (1 - Constants::SMALL_EPSILON);
        double max_maf = options_d["max_maf"] * (1 + Constants::SMALL_EPSILON);
        LOGGER.d(0, "min_maf: " + to_string(min_maf) + " max_maf: " + to_string(max_maf));
        vector<int> remove_extract_index;
        double cur_AF;
        for(int index = 0; index != AFA1.size(); index++){
            cur_AF = AFA1[index];
            if(cur_AF > 0.5) cur_AF = 1.0 - cur_AF;
            if((cur_AF < min_maf) || (cur_AF > max_maf)){
                remove_extract_index.push_back(index);
                LOGGER.d(0, to_string(index) + ": " + to_string(cur_AF));
            }
        }

        for(auto index = remove_extract_index.rbegin(); index != remove_extract_index.rend(); ++index){
            auto rindex = *index;
            LOGGER.d(5, to_string(rindex));
            AFA1.erase(AFA1.begin() + rindex);
            countA1A1.erase(countA1A1.begin() + rindex);
            countA1A2.erase(countA1A2.begin() + rindex);
            countA2A2.erase(countA2A2.begin() + rindex);
            RDev.erase(RDev.begin() + rindex);
        }

        marker->remove_extracted_index(remove_extract_index);
        init_AsyncBuffer();
        num_blocks = marker->count_extract() / Constants::NUM_MARKER_READ +
                     (marker->count_extract() % Constants::NUM_MARKER_READ != 0);
        LOGGER.i(0, to_string(remove_extract_index.size()) + " SNPs removed due to --maf or --max-maf,  " + to_string(marker->count_extract()) + " remained");
    }

}

void Geno::init_AF() {
    AFA1.clear();
    countA1A2.clear();
    countA1A1.clear();
    countA2A2.clear();
    RDev.clear();
    uint32_t num_marker = marker->count_extract();
    AFA1.resize(num_marker);
    countA1A1.resize(num_marker);
    countA1A2.resize(num_marker);
    countA2A2.resize(num_marker);
    RDev = vector<double>(num_marker, 0.0); 
    num_blocks = num_marker / Constants::NUM_MARKER_READ +
                 (num_marker % Constants::NUM_MARKER_READ != 0);
    LOGGER.d(0, "The program will run in " + to_string(num_blocks) + " blocks");
}

void Geno::init_AsyncBuffer(){
    if(asyncBuffer){
        delete asyncBuffer;
    }
    asyncBuffer = new AsyncBuffer<uint8_t>(num_byte_buffer);
}

Geno::~Geno(){
    delete asyncBuffer;
}

void Geno::out_freq(string filename){
    string name_frq = filename + ".frq";
    LOGGER.i(0, "Saving allele frequencies...");
    std::ofstream o_freq(name_frq.c_str());
    if (!o_freq) { LOGGER.e(0, "can not open the file [" + name_frq + "] to write"); }
    vector<string> out_contents;
    out_contents.reserve(AFA1.size() + 1);
    out_contents.push_back("CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS");
    for(int i = 0; i != AFA1.size(); i++){
        out_contents.push_back(marker->get_marker(marker->getExtractIndex(i)) + "\t" + to_string(AFA1[i])
                               + "\t" + to_string( 2*(countA1A1[i] + countA1A2[i] + countA2A2[i])));
    }
    std::copy(out_contents.begin(), out_contents.end(), std::ostream_iterator<string>(o_freq, "\n"));
    o_freq.close();
    LOGGER.i(0, "Allele frequencies of " + to_string(AFA1.size()) + " SNPs have been saved in the file [" + name_frq + "]");
}

bool Geno::check_bed(){
    FILE *pFile;
    uint64_t f_size;
    uint8_t * buffer = new uint8_t[3];

    pFile = fopen(bed_file.c_str(), "rb");
    if(pFile == NULL){
        LOGGER.e(0, "can't open [" + bed_file + "] to read");
    }
    fseek(pFile, 0, SEEK_END);
    f_size = ftell(pFile);
    rewind(pFile);

    if((f_size - 3) != num_byte_per_marker * marker->count_raw()){
        LOGGER.e(0, "invalid bed file [" + bed_file +
                "]. The sample and SNP number in bed file are different from bim and fam file");
    }

    size_t read_count = fread(buffer, 1, 3, pFile);
    fclose(pFile);
    if((read_count != 3) &&
            (*buffer != 0x6c) &&
            (*(buffer+1) != 0x1b) &&
            (*(buffer+2) != 0x01)){
        LOGGER.e(0, "invalid bed file [" + bed_file +
                    "], please convert it into new format (SNP major).");
        return false;
    }else{
        delete[] buffer;
        LOGGER.d(0, "BED file check OK");
        return true;
    }
}

void Geno::read_bed(){
    LOGGER.i(0, "Reading PLINK BED file from [" + bed_file + "] in SNP-major format...");
    //LOGGER.i(0, "Genotype data for " + to_string(pheno->count_keep()) + " individuals and " + 
    //            to_string(marker->count_extract()) + " SNPs to be included from [" + bed_file + "]");
    //clock_t begin = clock();
    FILE *pFile;
    pFile = fopen(bed_file.c_str(), "rb");
    if(pFile == NULL){
        LOGGER.e(0, "can't open [" + bed_file + "] to read.");
    }

    fseek(pFile, 3, SEEK_SET);

    uint8_t *w_buf = NULL, *iw_buf = NULL;
    bool isEOF = false;
    uint64_t index_marker_extracted = 0;
    uint64_t last_index = marker->getExtractIndex(index_marker_extracted);
    fseek(pFile, num_byte_per_marker * last_index, SEEK_CUR);

    int cur_block = 0;
    int last_num_marker = marker->count_extract() % Constants::NUM_MARKER_READ;

    LOGGER.ts("read_geno");

    while(!isEOF){
        w_buf = asyncBuffer->start_write();
        
        LOGGER.d(0, "read start");
        int num_marker = 0;
        size_t read_count = 0;
        iw_buf = w_buf;
        //begin = clock();
        int cur_num_block = (cur_block != (num_blocks - 1))? Constants::NUM_MARKER_READ : last_num_marker;
        do{
            uint64_t  cur_index = marker->getExtractIndex(index_marker_extracted);
            uint64_t lag_index = cur_index - last_index;

            //very arbitrary value to skip, maybe precised in the future.
            if(lag_index > 10){
                fseek(pFile, (lag_index - 1) * num_byte_per_marker, SEEK_CUR);
            }else{
                for(int64_t ab_index = 1; ab_index < lag_index; ab_index++){
                    fread(w_buf, 1, num_byte_per_marker, pFile);
                }
            }
            read_count = fread(w_buf, 1, num_byte_per_marker, pFile);
            w_buf += num_byte_per_marker;

            last_index = cur_index;
            index_marker_extracted++;
            num_marker++;

            if(read_count != num_byte_per_marker || index_marker_extracted == marker->count_extract()){
                asyncBuffer->setEOF();
                isEOF = true;
                break;
            }
        }while(num_marker != cur_num_block);

        pheno->mask_geno_keep(iw_buf, num_marker);
        asyncBuffer->end_write();
        LOGGER.d(2, "read block success");
        //std::cout << "   read bed" << "@T: " << 1.0 * (clock() - begin) / CLOCKS_PER_SEC << std::endl;
        cur_block++;
    }
    fclose(pFile);
    LOGGER.i(0, "Read bed time: " + to_string(LOGGER.tp("read_geno")));
}
/*
void Geno::freq(uint8_t *buf, int num_marker){
    if(num_marker_freq >= marker->count_extract()) return;
    for(int cur_marker_index = 0; cur_marker_index < num_marker; ++cur_marker_index){
        uint32_t curA1A1 = 0, curA1A2 = 0, curA2A2 = 0;
        uint64_t *pbuf = (uint64_t *) (buf + cur_marker_index * num_byte_per_marker);
        for(auto &index : pheno->keep_block_index){
            uint64_t geno_temp = *(pbuf + index);
            if(pheno->mask_add_items[index]){
                geno_temp = (geno_temp & pheno->mask_items[index]) + pheno->mask_add_items[index];
            }
            vector<uint16_t> genos = {(uint16_t)(geno_temp), (uint16_t)(geno_temp >> 16), 
                                       (uint16_t)(geno_temp >> 32), (uint16_t)(geno_temp >> 48)};
            g_table.set_count(genos, curA1A1, curA1A2, curA2A2); 

            int raw_index_marker = num_marker_freq + cur_marker_index;
            countA1A1[raw_index_marker] = curA1A1;
            countA1A2[raw_index_marker] = curA1A2;
            countA2A2[raw_index_marker] = curA2A2;
            AFA1[raw_index_marker] = (2.0 * curA1A1 + curA1A2) / (2.0 * (curA1A1 + curA1A2 + curA2A2));
        }
    }
    num_marker_freq += num_marker;
}
*/

void Geno::freq(uint8_t *buf, int num_marker) {
    if(num_marker_freq >= marker->count_extract()) return;
    //pheno->mask_geno_keep(buf, num_marker);
    int cur_num_marker_read = num_marker;
    static bool isLastTrunkSingle = (num_byte_per_marker % 2 != 0);
    static int num_trunk_per_marker = num_byte_per_marker / 2 + isLastTrunkSingle;
    uint16_t *p_buf;
    uint16_t *trunk_buf;      
    uint32_t curA1A1, curA1A2, curA2A2;
    int raw_index_marker;
    for(int cur_marker_index = 0; cur_marker_index < cur_num_marker_read; ++cur_marker_index){
        //It will cause problems when memory in critical stage.
        //However, other part of this program will also goes wrong in this situation.
        p_buf = (uint16_t *) (buf + cur_marker_index * num_byte_per_marker);
        trunk_buf = p_buf;
        curA1A1 = 0;
        curA1A2 = 0;
        curA2A2 = 0;
        for(int cur_trunk = 0; cur_trunk < num_trunk_per_marker - 1; ++cur_trunk){
            curA1A1 += g_table.get(*trunk_buf, 0);
            curA1A2 += g_table.get(*trunk_buf, 1);
            curA2A2 += g_table.get(*trunk_buf, 2);
            trunk_buf++;
        }

        uint16_t last_trunk = *trunk_buf;
        if(isLastTrunkSingle){
            last_trunk = last_trunk & 255;
            curA1A1 += (g_table.get(last_trunk, 0) - 4 - last_byte_NA_sample);
        }else{
            curA1A1 += (g_table.get(last_trunk, 0) - last_byte_NA_sample); 
        }

        curA1A2 += g_table.get(last_trunk, 1);
        curA2A2 += g_table.get(last_trunk, 2);

        raw_index_marker = num_marker_freq + cur_marker_index;
        

        countA1A1[raw_index_marker] = curA1A1;
        countA1A2[raw_index_marker] = curA1A2;
        countA2A2[raw_index_marker] = curA2A2;
        AFA1[raw_index_marker] = (2.0 * curA1A1 + curA1A2) / (2.0 * (curA1A1 + curA1A2 + curA2A2));
    }
    num_marker_freq += num_marker;
}

void Geno::freq2(uint8_t *buf, int num_marker) {
    const static uint64_t MASK = 6148914691236517205UL; 
    const static uint64_t FF = UINT64_MAX;
    static int pheno_count = pheno->count_keep();
    static int num_trunk = num_byte_per_marker / 8;
    static int remain_trunk = num_byte_per_marker % 8;
    static uint64_t last_mask = (FF >> remain_trunk) << remain_trunk ;

    if(num_marker_freq >= marker->count_extract()) return;

    int cur_num_marker_read = num_marker;

    for(int cur_marker_index = 0; cur_marker_index < cur_num_marker_read; ++cur_marker_index){
        uint32_t curA1A1, curA1A2, curA2A2;
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0;
        uint64_t *p_buf = (uint64_t *) (buf + cur_marker_index * num_byte_per_marker);
        int index = 0;
        for(; index < num_trunk ; index++){
            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            odd_ct += popcount(g_buf & MASK);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf & g_buf_h);
        }
        if(remain_trunk){
            uint64_t g_buf = p_buf[index + 1];
            g_buf &= last_mask;
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            odd_ct += popcount(g_buf & MASK);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf & g_buf_h);
        }

        curA1A1 = pheno_count + both_ct - even_ct - odd_ct;
        curA1A2 = even_ct - both_ct;
        curA2A2 = both_ct;

        int raw_index_marker = num_marker_freq + cur_marker_index;

        countA1A1[raw_index_marker] = curA1A1;
        countA1A2[raw_index_marker] = curA1A2;
        countA2A2[raw_index_marker] = curA2A2;
        AFA1[raw_index_marker] = (2.0 * curA1A1 + curA1A2) / (2.0 * (curA1A1 + curA1A2 + curA2A2));
    }
    num_marker_freq += num_marker;
}


void Geno::loop_block(vector<function<void (uint8_t *buf, int num_marker)>> callbacks) {
    num_finished_markers = 0;
    thread read_thread([this](){this->read_bed();});
    read_thread.detach();

    uint8_t *r_buf = NULL;
    bool isEOF = false;
    int cur_num_marker_read;

    LOGGER.ts("LOOP_GENO_TOT");
    LOGGER.ts("LOOP_GENO_PRE");

    for(int cur_block = 0; cur_block < num_blocks; ++cur_block){
        std::tie(r_buf, isEOF) = asyncBuffer->start_read();
        //LOGGER.i(0, "time get buffer: " + to_string(LOGGER.tp("LOOP_GENO_PRE")));

        LOGGER.d(0, "Process block " + std::to_string(cur_block));
        if(isEOF && cur_block != (num_blocks - 1)){
            LOGGER.e(0, "can't read to the end of file [" + bed_file + "].");
        }
        //correct the marker read;
        if(cur_block == (num_blocks - 1)){
            cur_num_marker_read = marker->count_extract() - Constants::NUM_MARKER_READ * cur_block;
        }else{
            cur_num_marker_read = Constants::NUM_MARKER_READ;
        }

        for(auto callback : callbacks){
         //   LOGGER.i(0, "time1: " + to_string(LOGGER.tp("LOOP_GENO_PRE")));
            callback(r_buf, cur_num_marker_read);
          //  LOGGER.i(0, "time2: " + to_string(LOGGER.tp("LOOP_GENO_PRE")));
        }

        asyncBuffer->end_read();

        num_finished_markers += cur_num_marker_read;
        if(cur_block % 100 == 0){
            float time_p = LOGGER.tp("LOOP_GENO_PRE");
            if(time_p > 300){
                LOGGER.ts("LOOP_GENO_PRE");
                float elapse_time = LOGGER.tp("LOOP_GENO_TOT");
                float finished_percent = (float) cur_block / num_blocks;
                float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                std::ostringstream ss;
                ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 
                
                LOGGER.i(1, ss.str());
            }
        }
    }
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(1) << "100% Finished in " << LOGGER.tp("LOOP_GENO_TOT") / 60 << " min";
    LOGGER.i(1, ss.str());
}

void Geno::makeMarkerX(uint8_t *buf, int cur_marker, double *w_buf, bool center, bool std){
    uint32_t cur_raw_marker = num_finished_markers + cur_marker;
    uint8_t *cur_buf = buf + cur_marker * num_byte_per_marker;
    double af = AFA1[cur_raw_marker];
    double mu = 2.0 * af;
    double center_value = 0.0;
    if(center){
        center_value = mu;
    }
    double rdev = 1.0;
    if(std){
        rdev = 1.0 / sqrt(mu * (1.0 - af));
        RDev[cur_raw_marker] = rdev; 
    }

    double g1_lookup[4];
    g1_lookup[0] = (2.0 - center_value) * rdev;
    g1_lookup[1] = (mu - center_value) * rdev;
    g1_lookup[2] = (1.0 - center_value) * rdev;
    g1_lookup[3] = (0.0 - center_value) * rdev;

    double g_lookup[256][4];
    for(uint16_t i = 0; i <= 255; i++){
        for(uint16_t j = 0; j < 4; j++){
            g_lookup[i][j] = g1_lookup[(i >> (2 * j)) & 3];
        }
    }

    uint32_t n_keep = pheno->count_keep();
    vector<uint32_t>& index_keep = pheno->get_index_keep();

    for(uint32_t index = 0; index != n_keep; index++){
        uint32_t raw_index = index_keep[index];
        w_buf[index] = g_lookup[*(cur_buf + (raw_index / 4))][raw_index % 4]; 
    }
}

void Geno::addOneFileOption(string key_store, string append_string, string key_name,
                                     map<string, vector<string>> options_in) {
    if(options_in.find(key_name) != options_in.end()){
        if(options_in[key_name].size() == 1){
            options[key_store] = options_in[key_name][0] + append_string;
        }else if(options_in[key_name].size() > 1){
            options[key_store] = options_in[key_name][0] + append_string;
            LOGGER.w(0, "Geno: multiple " + key_name + ", use the first one only" );
        }else{
            LOGGER.e(0, "no " + key_name + " parameter found");
        }
        std::ifstream f(options[key_store].c_str());
        if(!f.good()){
            LOGGER.e(0, key_name + " " + options[key_store] + " not found");
        }
        f.close();
    }
}

int Geno::registerOption(map<string, vector<string>>& options_in) {
    int return_value = 0;
    addOneFileOption("geno_file", ".bed","--bfile", options_in);
    options_in.erase("--bfile");
    options_d["min_maf"] = 0.0;
    options_d["max_maf"] = 0.5;

    if(options_in.find("--maf") != options_in.end()){
        auto option = options_in["--maf"];
        if(option.size() == 1){
            try{
                options_d["min_maf"] = std::stod(option[0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "illegal value in --maf");
            }
            if(options_d["min_maf"]<0.0 || options_d["max_maf"]>0.5){
                LOGGER.e(0, "--maf can't be smaller than 0 or larger than 0.5");
            }

        }else{
            LOGGER.e(0, "multiple value in --maf, not supported currently");
        }
        options_in.erase("--maf");
    }

     if(options_in.find("--max-maf") != options_in.end()){
        auto option = options_in["--max-maf"];
        if(option.size() == 1){
            try{
                options_d["max_maf"] = std::stod(option[0]);
           }catch(std::invalid_argument&){
                LOGGER.e(0, "illegal value in --maf");
           }
           if(options_d["max_maf"] < 0.0 || options_d["max_maf"] > 0.5){
               LOGGER.e(0, "--max-maf can't be smaller than 0 or larger than 0.5");
           }
        }else{
            LOGGER.e(0, "multiple value in --maf, not supported currently");
        }
        options_in.erase("--max-maf");
    }

    if(options_d["min_maf"] > options_d["max_maf"]){
        LOGGER.e(0, "--maf can't be larger than --max-maf value");
    }

    if(options_in.find("--freq") != options_in.end()){
        processFunctions.push_back("freq");
        if(options_in["--freq"].size() != 0){
            LOGGER.w(0, "--freq should not follow by other parameters, if you want to calculate in founders only, "
                    "please specify by --founders option");
        }
        options_in.erase("--freq");

        options["out"] = options_in["--out"][0];

        return_value++;
    }

    addOneFileOption("update_freq_file", "", "--update-freq", options_in);
    addOneFileOption("update_ref_allele_file", "", "--update-ref-allele", options_in);
    

    return return_value;
}

void Geno::processMain() {
    vector<function<void (uint8_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "freq"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            if(geno.num_marker_freq == 0 ){
                LOGGER.i(0, "Computing allele frequencies...");
                callBacks.push_back(bind(&Geno::freq2, &geno, _1, _2));
                geno.loop_block(callBacks);
            }
            geno.out_freq(options["out"]);
            callBacks.clear();
        }
    }

}
