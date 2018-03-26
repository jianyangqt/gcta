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
#include "omp.h"
#include "ThreadPool.h"
#include <cstring>

#ifdef _WIN64
  #include <intrin.h>
  uint32_t __inline CTZ64U(uint64_t value){
      unsigned long tz = 0;
      _BitScanForward64(&tz, value);
      return tz;
  }
  
  uint32_t __inline CLZ64U(uint64_t value){
      unsigned long lz = 0;
      _BitScanReverse64(&lz, value);
      return 63 - lz;
  }
#else
  //#define CTZU __builtin_ctz
  //#define CLZU __builtin_clz
  #ifdef __linux__
  #pragma message("multiple target")
  __attribute__((target_clones("avx2","default")))
  #endif
  uint32_t CTZ64U(uint64_t value){
      return __builtin_ctzll(value);
  }
 
#endif

typedef uint32_t halfword_t;
const uintptr_t k1LU = (uintptr_t)1;

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
    num_raw_sample = pheno->count_raw();
    num_byte_per_marker = (num_raw_sample + 3) / 4;
    num_byte_buffer = num_byte_per_marker * Constants::NUM_MARKER_READ;
    last_byte_NA_sample = (4 - (num_raw_sample % 4)) % 4;

    num_keep_sample = pheno->count_keep();
    num_item_1geno = (num_keep_sample + 31) / 32;
    num_item_geno_buffer = num_item_1geno * Constants::NUM_MARKER_READ;

    vector<uint32_t>& index_keep = pheno->get_index_keep();
    keep_mask = new uint64_t[(num_raw_sample + 63)/64]();
    for(auto keep_item : index_keep){
        uint32_t cur_qword = keep_item / 64;
        uint32_t cur_offset = keep_item % 64;
        keep_mask[cur_qword] |= k1LU << cur_offset;
    }

    check_bed();

    init_AF();
    init_AsyncBuffer();
    filter_MAF();
}

Geno::~Geno(){
    delete asyncBuffer;
    delete[] keep_mask;
}

void Geno::filter_MAF(){
    if((options_d["min_maf"] != 0.0) || (options_d["max_maf"] != 0.5)){
        LOGGER.i(0, "Computing allele frequencies...");
        vector<function<void (uint64_t *, int)>> callBacks;
        callBacks.push_back(bind(&Geno::freq64, this, _1, _2));
        loop_64block(callBacks);
        // We adopt the EPSILON from plink, because the double value may have some precision issue;
        double min_maf = options_d["min_maf"] * (1 - Constants::SMALL_EPSILON);
        double max_maf = options_d["max_maf"] * (1 + Constants::SMALL_EPSILON);
        LOGGER.d(0, "min_maf: " + to_string(min_maf) + " max_maf: " + to_string(max_maf));
        vector<uint32_t> extract_index;
        double cur_AF;

        for(int index = 0; index != AFA1.size(); index++){
            cur_AF = AFA1[index];
            if(cur_AF > 0.5) cur_AF = 1.0 - cur_AF;
            if((cur_AF > min_maf) && (cur_AF < max_maf)){
                extract_index.push_back(index);
                LOGGER.d(0, to_string(index) + ": " + to_string(cur_AF));
            }
        }

        vector<double> AFA1o = AFA1;
        vector<uint32_t> countA1A1o = countA1A1;
        vector<uint32_t> countA1A2o = countA1A2;
        vector<uint32_t> countA2A2o = countA2A2;
        vector<double> RDevo = RDev;

        AFA1.resize(extract_index.size());
        countA1A1.resize(extract_index.size());
        countA1A2.resize(extract_index.size());
        countA2A2.resize(extract_index.size());
        RDev.resize(extract_index.size());

        #pragma omp parallel for
        for(uint32_t index = 0; index < extract_index.size(); index++){
            uint32_t cur_index = extract_index[index];
            AFA1[index] = AFA1o[cur_index];
            countA1A1[index] = countA1A1[cur_index];
            countA1A2[index] = countA1A2[cur_index];
            countA2A2[index] = countA2A2[cur_index];
            RDev[index] = RDevo[cur_index];
        }

        marker->keep_extracted_index(extract_index);

        init_AsyncBuffer();
        num_blocks = marker->count_extract() / Constants::NUM_MARKER_READ +
                     (marker->count_extract() % Constants::NUM_MARKER_READ != 0);
        LOGGER.i(0, to_string(extract_index.size()) + " SNPs remained due to --maf or --max-maf,  " + to_string(marker->count_extract()) + " remained");
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
    //omp_set_num_threads(THREADS.getThreadCount()+1);
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

        asyncBuffer->end_write();
        LOGGER.d(2, "read block success");
        //std::cout << "   read bed" << "@T: " << 1.0 * (clock() - begin) / CLOCKS_PER_SEC << std::endl;
        cur_block++;
    }
    fclose(pFile);
    //LOGGER.i(0, "Read bed time: " + to_string(LOGGER.tp("read_geno")));
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
    //pheno->mask_geno_keep(buf, num_marker);
    const static uint64_t MASK = 6148914691236517205UL; 
    static int pheno_count = pheno->count_keep();
    static int num_trunk = num_byte_per_marker / 8;
    static int remain_bit = (num_byte_per_marker % 8) * 8;
    static int move_bit = 64 - remain_bit;

    if(num_marker_freq >= marker->count_extract()) return;

    int cur_num_marker_read = num_marker;
    
    #pragma omp parallel for schedule(dynamic) 
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
        if(remain_bit){
            uint64_t g_buf = p_buf[index];
            g_buf = (g_buf << move_bit) >> move_bit;
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

void Geno::freq64(uint64_t *buf, int num_marker) {
    //pheno->mask_geno_keep(buf, num_marker);
    const static uint64_t MASK = 6148914691236517205UL; 
    if(num_marker_freq >= marker->count_extract()) return;

    int cur_num_marker_read = num_marker;
    
    #pragma omp parallel for schedule(dynamic) 
    for(int cur_marker_index = 0; cur_marker_index < cur_num_marker_read; ++cur_marker_index){
        uint32_t curA1A1, curA1A2, curA2A2;
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0;
        uint64_t *p_buf = buf + cur_marker_index * num_item_1geno;
        for(int index = 0; index < num_item_1geno ; index++){
            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            odd_ct += popcount(g_buf & MASK);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf & g_buf_h);
        }

        curA1A1 = num_keep_sample + both_ct - even_ct - odd_ct;
        curA1A2 = even_ct - both_ct;
        curA2A2 = both_ct;

        int raw_index_marker = num_marker_freq + cur_marker_index;

        countA1A1[raw_index_marker] = curA1A1;
        countA1A2[raw_index_marker] = curA1A2;
        countA2A2[raw_index_marker] = curA2A2;
        double cur_af = (2.0 * curA1A1 + curA1A2) / (2.0 * (curA1A1 + curA1A2 + curA2A2));
        if(marker->isEffecRev(raw_index_marker)){
            cur_af = 1.0 - cur_af;
        }
        AFA1[raw_index_marker] = cur_af;
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

// extracted and revised from plink2.0
// GPL v3, license detailed on github
// https://github.com/chrchang/plink-ng
void copy_quaterarr_nonempty_subset(uint64_t* raw_quaterarr[], const uint64_t* subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uint64_t* output_quaterarr[], const int num_marker) {
    // in plink 2.0, we probably want (0-based) bit raw_quaterarr_entry_ct of
    // subset_mask to be always allocated and unset.  This removes a few special
    // cases re: iterating past the end of arrays.
    static const uint32_t kBitsPerWordD2 = 32;

    uint64_t cur_output_word[num_marker];
    memset(cur_output_word, 0, num_marker * 8);

    uint64_t* output_quaterarr_iter[num_marker];
    uint64_t* output_quaterarr_last[num_marker];
    for(int i = 0; i != num_marker; i++){
        output_quaterarr_iter[i] = output_quaterarr[i];
        output_quaterarr_last[i] = &(output_quaterarr[i][subset_entry_ct / kBitsPerWordD2]);
    }
    const uint32_t word_write_halfshift_end = subset_entry_ct % kBitsPerWordD2;
    uint32_t word_write_halfshift = 0;
    // if <= 2/3-filled, use sparse copy algorithm
    // (tried copy_bitarr_subset() approach, that actually worsened things)
    if (subset_entry_ct * (3 * k1LU) <= raw_quaterarr_entry_ct * (2 * k1LU)) {
        uint32_t subset_mask_widx = 0;
        while (1) {
            const uint64_t cur_include_word = subset_mask[subset_mask_widx];
            if (cur_include_word) {
                uint32_t wordhalf_idx = 0;
                uint32_t cur_include_halfword = (halfword_t)cur_include_word;
                while (1) {
                    if (cur_include_halfword) {
                        uint64_t raw_quaterarr_word[num_marker];
                        uint32_t temp_index = subset_mask_widx * 2 + wordhalf_idx;
                        for(int i = 0; i != num_marker; i++){
                            raw_quaterarr_word[i] = raw_quaterarr[i][temp_index];
                        }
                        do {
                            uint32_t rqa_idx_lowbits = CTZ64U(cur_include_halfword);
                            uint32_t lshift = word_write_halfshift * 2; 
                            uint32_t rshift = rqa_idx_lowbits * 2;
                            for(int i = 0; i != num_marker; i++){
                                cur_output_word[i] |= ((raw_quaterarr_word[i] >> rshift) & 3) << lshift;
                            }
                            if (++word_write_halfshift == kBitsPerWordD2) {
                                for(int i = 0; i != num_marker; i++){
                                    *(output_quaterarr_iter[i])++ = cur_output_word[i];
                                }
                                word_write_halfshift = 0;
                                //cur_output_word = 0;
                                memset(cur_output_word, 0, num_marker * 8);
                            }
                            cur_include_halfword &= cur_include_halfword - 1;
                        } while (cur_include_halfword);
                    }
                    if (wordhalf_idx) {
                        break;
                    }
                    ++wordhalf_idx;
                    cur_include_halfword = cur_include_word >> kBitsPerWordD2;
                }
                if (output_quaterarr_iter[0] == output_quaterarr_last[0]) {
                    if (word_write_halfshift == word_write_halfshift_end) {
                        if (word_write_halfshift_end) {
                            for(int i = 0; i != num_marker; i++){
                                *(output_quaterarr_last[i]) = cur_output_word[i];
                            }
                        }
                        return;
                    }
                }
            }
            ++subset_mask_widx;
        }
    }

    const uint64_t* raw_quaterarr_iter[num_marker];
    for(int i = 0; i != num_marker; i++){
        raw_quaterarr_iter[i] = raw_quaterarr[i];
    }
    //const uint64_t* raw_quaterarr_iter = raw_quaterarr;
    while (1) {
        const uint64_t cur_include_word = *subset_mask++;
        uint32_t wordhalf_idx = 0;
        uint64_t cur_include_halfword = (halfword_t)cur_include_word;
        while (1) {
           // uintptr_t raw_quaterarr_word = *raw_quaterarr_iter++;
            uint64_t raw_quaterarr_word[num_marker];
            for(int i = 0; i != num_marker; i++){
                raw_quaterarr_word[i] = *(raw_quaterarr_iter[i]++);
            }
            while (cur_include_halfword) {
                uint32_t rqa_idx_lowbits = CTZ64U(cur_include_halfword); // tailing zero
                uint64_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
                uint64_t raw_quaterarr_curblock_unmasked[num_marker];
                int m_bit = rqa_idx_lowbits * 2;
                for(int i = 0; i != num_marker; i++){
                    raw_quaterarr_curblock_unmasked[i] = raw_quaterarr_word[i] >> m_bit; 
                }
                //uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2); // remove mask bit tailing zero, not to keep
                uint32_t rqa_block_len = CTZ64U(halfword_invshifted);  // find another keep
                uint32_t block_len_limit = kBitsPerWordD2 - word_write_halfshift;
                m_bit = 2 * word_write_halfshift;
                for(int i = 0; i != num_marker; i++){
                    cur_output_word[i] |= raw_quaterarr_curblock_unmasked[i] << m_bit; // avoid overwrite current saved bits
                }
                if (rqa_block_len < block_len_limit) { //2  16
                    word_write_halfshift += rqa_block_len; // 0 2
                    m_bit = 2 * word_write_halfshift;
                    uint64_t temp_mask = (k1LU << m_bit) - k1LU;
                    for(int i = 0; i != num_marker; i++){
                        cur_output_word[i] &= temp_mask; // mask high end, and keep low needed bits
                    }
                } else {
                    // no need to mask, extra bits vanish off the high end
                    for(int i = 0; i != num_marker; i++){
                        *(output_quaterarr_iter[i]++) = cur_output_word[i];
                    }
                    word_write_halfshift = rqa_block_len - block_len_limit;
                    if (word_write_halfshift) {
                        uint64_t t_mask = ((k1LU << (2 * word_write_halfshift)) - k1LU), mi_bit = 2 * block_len_limit;
                        for(int i = 0; i != num_marker; i++){
                            cur_output_word[i] = (raw_quaterarr_curblock_unmasked[i] >> mi_bit) & t_mask;
                        }
                    } else {
                        // avoid potential right-shift-[word length]
                        //cur_output_word = 0;
                        memset(cur_output_word, 0, num_marker * 8);
                    }
                }
                cur_include_halfword &= (~(k1LU << (rqa_block_len + rqa_idx_lowbits))) + k1LU;
            }
            if (wordhalf_idx) {
                break;
            }
            ++wordhalf_idx;
            cur_include_halfword = cur_include_word >> kBitsPerWordD2;
        }
        if (output_quaterarr_iter[0] == output_quaterarr_last[0]) {
            if (word_write_halfshift == word_write_halfshift_end) {
                if (word_write_halfshift_end) {
                    for(int i = 0; i != num_marker; i++){
                        *(output_quaterarr_last[i]) = cur_output_word[i];
                    }
                }
                return;
            }
        }
    }
}

void Geno::move_geno(uint8_t *buf, uint64_t *keep_list, uint32_t num_raw_sample, uint32_t num_keep_sample, uint32_t num_marker, uint64_t *geno_buf){
    /*
    if(num_raw_sample == num_keep_sample){
        for(int index = 0; index < num_marker; index++){
            memcpy(

        }
    }
    */

    uint32_t num_byte_keep_geno = (num_keep_sample + 3) / 4;
    uint32_t num_byte_per_marker = (num_raw_sample + 3) / 4;
    uint32_t num_qword_per_marker = (num_byte_keep_geno + 7) / 8;

    static int remain_bit = (num_byte_per_marker % 8) * 8;
    static int move_bit = 64 - remain_bit;
    static uint64_t MASK = (UINT64_MAX << move_bit) >> move_bit;

    const int move_markers = 5;

    #pragma omp parallel for schedule(dynamic) 
    for(uint32_t index = 0; index < num_marker; index += move_markers){
        uint64_t *pbuf[move_markers], *gbuf[move_markers];
        for(uint32_t i = 0; i < move_markers; i++){
            pbuf[i] = (uint64_t *) (buf + (index + i) * num_byte_per_marker);
            gbuf[i] = geno_buf + (index + i) * num_qword_per_marker; 
        }
        copy_quaterarr_nonempty_subset(pbuf, keep_list, num_raw_sample, num_keep_sample, gbuf, move_markers);
    }
}

void Geno::loop_64block(vector<function<void (uint64_t *buf, int num_marker)>> callbacks) {
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

        uint64_t *geno_buf = new uint64_t[num_item_geno_buffer];
        move_geno(r_buf, keep_mask, num_raw_sample, num_keep_sample, cur_num_marker_read, geno_buf);
        asyncBuffer->end_read();

        for(auto callback : callbacks){
         //   LOGGER.i(0, "time1: " + to_string(LOGGER.tp("LOOP_GENO_PRE")));
            callback(geno_buf, cur_num_marker_read);
          //  LOGGER.i(0, "time2: " + to_string(LOGGER.tp("LOOP_GENO_PRE")));
        }

        delete[] geno_buf;

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


void Geno::makeMarkerX(uint64_t *buf, int cur_marker, double *w_buf, bool center, bool std){
    static uint32_t last_sample = num_keep_sample % 32;
    static uint32_t last_8block = last_sample / 4;
    static uint32_t last_2block = last_sample % 4; 

    uint32_t cur_raw_marker = num_finished_markers + cur_marker;
    uint64_t *cur_buf = buf + cur_marker * num_item_1geno;
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

    int sub_index = 0;
    uint32_t index = 0;
    for(; index < num_item_1geno - 1; index++){
        uint64_t geno = cur_buf[index];
        int move_bit = 0;
        for(int ri = 0; ri < 8; ri++){
            uint8_t geno_temp = (uint8_t) (geno >> move_bit);
            w_buf[sub_index] = g_lookup[geno_temp][0]; 
            w_buf[sub_index + 1] = g_lookup[geno_temp][1]; 
            w_buf[sub_index + 2] = g_lookup[geno_temp][2]; 
            w_buf[sub_index + 3] = g_lookup[geno_temp][3]; 
            move_bit += 8;
            sub_index += 4;
        }
    }
    //last geno
    uint64_t geno = cur_buf[index];
    int move_bit = 0;
    for(int ri = 0; ri < last_8block; ri++){
        uint8_t geno_temp = (uint8_t) (geno >> move_bit);
        w_buf[sub_index] = g_lookup[geno_temp][0]; 
        w_buf[sub_index + 1] = g_lookup[geno_temp][1]; 
        w_buf[sub_index + 2] = g_lookup[geno_temp][2]; 
        w_buf[sub_index + 3] = g_lookup[geno_temp][3]; 
        move_bit += 8;
        sub_index += 4;
    }
    //last 4
    uint8_t geno_temp = (uint8_t) (geno >> move_bit);
    for(int ri = 0; ri < last_2block; ri++){
        w_buf[sub_index + ri] = g_lookup[geno_temp][ri];
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
    //vector<function<void (uint8_t *, int)>> callBacks;
    vector<function<void (uint64_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "freq"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            if(geno.num_marker_freq == 0 ){
                LOGGER.i(0, "Computing allele frequencies...");
                callBacks.push_back(bind(&Geno::freq64, &geno, _1, _2));
                geno.loop_64block(callBacks);
            }
            geno.out_freq(options["out"]);
            callBacks.clear();
        }
    }

}
