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

#define NOMINMAX
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
#include <boost/algorithm/string.hpp>
#include "OptionIO.h"
#include "zlib.h"
#include "zstd.h"
#include <cstring>
#include "cpu.h"
#include <Eigen/Eigen>
#include <algorithm>
#include "PgenReader.h"
#include <numeric>

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
  #if defined(__linux__) && GCTA_CPU_x86
  __attribute__((target("default")))
  #endif
  uint32_t CTZ64U(uint64_t value){
      return __builtin_ctzll(value);
  }
  #if defined(__linux__) && GCTA_CPU_x86
  __attribute__((target("popcnt")))
  uint32_t CTZ64U(uint64_t value){
      return __builtin_ctzll(value);
  }
  #endif
 
#endif

#if defined(__linux__) && GCTA_CPU_x86
__attribute__((target("default")))
#endif
uint64_t fill_inter_zero(uint64_t x) {
   uint64_t t;
   t = (x ^ (x >> 16)) & 0x00000000FFFF0000;
   x = x ^ t ^ (t << 16);
   t = (x ^ (x >> 8)) & 0x0000FF000000FF00;
   x = x ^ t ^ (t << 8);
   t = (x ^ (x >> 4)) & 0x00F000F000F000F0;
   x = x ^ t ^ (t << 4);
   t = (x ^ (x >> 2)) & 0x0C0C0C0C0C0C0C0C;
   x = x ^ t ^ (t << 2);
   t = (x ^ (x >> 1)) & 0x2222222222222222;
   x = x ^ t ^ (t << 1);
   return x;
}
#if defined(__linux__) && GCTA_CPU_x86
#include <x86intrin.h>
__attribute__((target("bmi2")))
uint64_t fill_inter_zero(uint64_t x) {
    return _pdep_u64(x, 0x5555555555555555U);
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
    geno_files.clear();
    int num_geno = 0;
    if(options.find("geno_file") != options.end()){
        genoFormat = "BED";
        geno_files.push_back(options["geno_file"]);
        num_geno++;
        hasInfo = false;
    }

    if(options.find("m_file") != options.end()){
        genoFormat = "BED";
        boost::split(geno_files, options["m_file"], boost::is_any_of("\t "));
        //std::transform(geno_files.begin(), geno_files.end(), geno_files.begin(), [](string r){return r + ".bed";});
        num_geno++;
        hasInfo = false;
    }

    if(options.find("pgen_file") != options.end()){
        genoFormat = "PGEN";
        geno_files.push_back(options["pgen_file"]);
        boost::split(geno_files, options["pgen_file"], boost::is_any_of("\t "));
        num_geno++;
        hasInfo = false;
    }

    if(options.find("mpgen_file") != options.end()){
        genoFormat = "PGEN";
        boost::split(geno_files, options["mpgen_file"], boost::is_any_of("\t "));
        num_geno++;
        hasInfo = false;
    }

    if(options.find("bgen_file") != options.end()){
        genoFormat = "BGEN";
        geno_files.push_back(options["bgen_file"]);
        num_geno++;
        hasInfo = true;
    }

    if(options.find("mbgen_file") != options.end()){
        genoFormat = "BGEN";
        boost::split(geno_files, options["mbgen_file"], boost::is_any_of("\t "));
        num_geno++;
        hasInfo = true;
    }

    if(num_geno == 0){
        LOGGER.e(0, "no genotype file is specified");
    }

    //open and check genotype files
    openGFiles();

    this->pheno = pheno;
    this->marker = marker;

    //this->sampleKeepIndex = pheno->get_index_keep();

    //register format handlers
    preProcessFuncs["BED"] = &Geno::preProcess_bed;
    preProcessFuncs["BGEN"] = &Geno::preProcess_bgen;

    getGenoArrayFuncs["BED"] = &Geno::getGenoArray_bed;
    getGenoArrayFuncs["BGEN"] = &Geno::getGenoArray_bgen;

    endProcessFuncs["BED"] = &Geno::endProcess_bed;
    endProcessFuncs["BGEN"] = &Geno::endProcess_bgen;

    // register format handlers subset manner
    preGenoDoubleFuncs["BED"] = &Geno::preGenoDouble_bed;
    preGenoDoubleFuncs["PGEN"] = &Geno::preGenoDouble_bed;
    preGenoDoubleFuncs["BGEN"] = &Geno::preGenoDouble_bgen;

    getGenoDoubleFuncs["BED"] = &Geno::getGenoDouble_bed;
    getGenoDoubleFuncs["PGEN"] = &Geno::getGenoDouble_bed;
    getGenoDoubleFuncs["BGEN"] = &Geno::getGenoDouble_bgen;

    endGenoDoubleFuncs["BED"] = &Geno::endGenoDouble_bed;
    endGenoDoubleFuncs["PGEN"] = &Geno::endGenoDouble_bed;
    endGenoDoubleFuncs["BGEN"] = &Geno::endGenoDouble_bgen;

    readGenoFuncs["BED"] = &Geno::readGeno_bed;
    readGenoFuncs["PGEN"] = &Geno::readGeno_bed;
    readGenoFuncs["BGEN"] = &Geno::readGeno_bgen;
    
    //BED legacy codes
    num_raw_sample = pheno->count_raw();
    num_byte_per_marker = (num_raw_sample + 3) / 4;
    num_byte_buffer = num_byte_per_marker * Constants::NUM_MARKER_READ;
    last_byte_NA_sample = (4 - (num_raw_sample % 4)) % 4;

    /*
    if(options.find("bed_file") != options.end()
            || options.find("m_file") != options.end()){
        check_bed();
    }
    */

    string alleleFileName = "";
    if(options.find("update_freq_file") != options.end()){
        alleleFileName = options["update_freq_file"];
    }

    setMAF(options_d["min_maf"]);
    setMaxMAF(options_d["max_maf"]);
    setFilterInfo(options_d["info_score"]);
    setFilterMiss(1.0 - options_d["geno_rate"]);

    string filterprompt = "Threshold to filter variants:";
    bool outFilterPrompt = false;
    if(options_d["min_maf"] != 0.0){
        filterprompt += " MAF > " + to_string(options_d["min_maf"]);
        outFilterPrompt = true;
    }
    if(options_d["max_maf"] != 0.5){
        filterprompt += string(outFilterPrompt ? "," : "") + " MAF < " + to_string(options_d["max_maf"]);
        outFilterPrompt = true;
    }
    if(options_d["info_score"] != 0.0){
        filterprompt += string(outFilterPrompt ? "," : "") + " imputation INFO score > " + to_string(options_d["info_score"]);
        outFilterPrompt = true;
    }

    if(options_d["geno_rate"] != 1.0){
        filterprompt += string(outFilterPrompt ? "," : "") + " missingness rate < " + to_string(options_d["geno_rate"]);
        outFilterPrompt = true;
    }
    if(outFilterPrompt){
        LOGGER << filterprompt << "." << std::endl;
    }
    if(options_d["dos_dc"] == 1.0){
        iGRMdc = 1;
        iDC = 1;
        LOGGER << "Switch to the full dosage compensation mode." << std::endl;
    }else if(options_d["dos_dc"] == 0.0){
        iGRMdc = 0;
        iDC = 0;
        LOGGER << "Switch to the no dosage compensation mode." << std::endl;
    }else{
        // default equal variance mode for grm
        iGRMdc = -1;
        // compensation mode for x
        iDC = 1;
    }

    init_AF(alleleFileName);

    //olds
    //init_AsyncBuffer();

    //num_keep_sample = 0;
    //init_keep();
    //olds;
}

Geno::~Geno(){
    if(asyncBuffer)delete asyncBuffer;
    if(keep_mask)delete[] keep_mask;
    if(keep_male_mask) delete[] keep_male_mask;
}

void Geno::init_keep(){
    num_keep_sample = pheno->count_keep();
    total_markers = 2 * num_keep_sample;
    num_byte_keep_geno1 = (num_keep_sample + 3) / 4;
    num_item_1geno = (num_keep_sample + 31) / 32;
    num_item_geno_buffer = num_item_1geno * Constants::NUM_MARKER_READ;

    if(keep_mask) delete[] keep_mask;
    keep_mask = new uint64_t[(num_raw_sample + 63)/64]();
    pheno->getMaskBit(keep_mask);

    isX = false;
    if(keep_male_mask) delete[] keep_male_mask;
    if(options.find("sex") != options.end()){
        isX = true;
        num_male_keep_sample = pheno->count_male();
        total_markers -= pheno->count_male();
        keep_male_mask = new uint64_t[(num_keep_sample + 63)/64]();
        pheno->getMaskBitMale(keep_male_mask);
    }
}

uint32_t Geno::getTotalMarker(){
    return total_markers;
}

void Geno::setSexMode(){
    std::map<string, vector<string>> t_option;
    t_option["--chrx"] = {};
    t_option["--filter-sex"] = {}; 
    Pheno::registerOption(t_option);
    Marker::registerOption(t_option);
    Geno::registerOption(t_option);
}


//
//true:  filtered; flase: not neccesory to filter
bool Geno::filterMAF(){
    if((options_d["min_maf"] != 0.0) || (options_d["max_maf"] != 0.5)){
        LOGGER.i(0, "Computing allele frequencies...");
        vector<function<void (uint64_t *, int)>> callBacks;
        callBacks.push_back(bind(&Geno::freq64, this, _1, _2));
        loop_64block(this->marker->get_extract_index(), callBacks);
        // We adopt the EPSILON from plink, because the double value may have some precision issue;
        double min_maf = options_d["min_maf"] * (1.0 - Constants::SMALL_EPSILON);
        double max_maf = options_d["max_maf"] * (1.0 + Constants::SMALL_EPSILON);
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
        //vector<uint32_t> countA1A1o = countA1A1;
        //vector<uint32_t> countA1A2o = countA1A2;
        //vector<uint32_t> countA2A2o = countA2A2;
        vector<uint32_t> countMarkerso = countMarkers;
        //vector<double> RDevo = RDev;

        AFA1.resize(extract_index.size());
        //countA1A1.resize(extract_index.size());
        //countA1A2.resize(extract_index.size());
        //countA2A2.resize(extract_index.size());
        countMarkers.resize(extract_index.size());
        //RDev.resize(extract_index.size());

        #pragma omp parallel for
        for(uint32_t index = 0; index < extract_index.size(); index++){
            uint32_t cur_index = extract_index[index];
            AFA1[index] = AFA1o[cur_index];
            //countA1A1[index] = countA1A1[cur_index];
            //countA1A2[index] = countA1A2[cur_index];
            //countA2A2[index] = countA2A2[cur_index];
            countMarkers[index] = countMarkerso[cur_index];
            //RDev[index] = RDevo[cur_index];
        }

        marker->keep_extracted_index(extract_index);

        //init_AsyncBuffer();
        num_blocks = marker->count_extract() / Constants::NUM_MARKER_READ +
                     (marker->count_extract() % Constants::NUM_MARKER_READ != 0);
        //LOGGER << "num_blocks: " << num_blocks << ", count extract: " << marker->count_extract() << std::endl;
        LOGGER.i(0, to_string(extract_index.size()) + " SNPs remain from --maf or --max-maf,  ");
        num_marker_freq = extract_index.size();
        bFreqFiltered = true;
        num_finished_markers = 0;
        return true;
    }else{
        return false;
    }

}

void Geno::init_AF(string alleleFileName) {
    AFA1.clear();
    //countA1A2.clear();
    //countA1A1.clear();
    //countA2A2.clear();
    countMarkers.clear();
    //RDev.clear();
    if(!alleleFileName.empty()){
        LOGGER.i(0, "Reading frequencies from [" + alleleFileName + "]...");
        vector<int> field_return = {2};
        vector<string> fields;
        vector<bool> a_rev;
        marker->matchSNPListFile(alleleFileName, 3, field_return, fields, a_rev, false);

        AFA1.resize(a_rev.size());
        vector<uint32_t> extract_index;
        bool filterByMaf = false;
        for(int i = 0; i < a_rev.size(); i++){
            double af;
            try{
                af = stod(fields[i]);
            }catch(std::out_of_range &){
                LOGGER.e(0, "the third column should be numeric");
            }
            if(af < 0 || af > 1.0){
                LOGGER.e(0, "frequency values should range from 0 to 1");
            }
            double maf = std::min(af, 1.0 - af);
            if(maf > min_maf && maf < max_maf){
                extract_index.push_back(i);
            }else{
                filterByMaf = true;
            }

            if(a_rev[i]){
                AFA1[i] = 1.0 - af;
            }else{
                AFA1[i] = af;
            }
        }
        LOGGER.i(0, "Frequencies of " + to_string(AFA1.size()) + " SNPs are updated.");

        marker->keep_extracted_index(extract_index);
        vector<double> AFA1o = AFA1;
        //vector<uint32_t> countMarkerso = countMarkers;

        AFA1.resize(extract_index.size());
        countMarkers.resize(extract_index.size());

        #pragma omp parallel for
        for(uint32_t index = 0; index < extract_index.size(); index++){
            uint32_t cur_index = extract_index[index];
            AFA1[index] = AFA1o[cur_index];
           // countMarkers[index] = countMarkerso[cur_index];
        }
        if(filterByMaf)LOGGER << "  " << extract_index.size() << " SNPs remain after MAF filtering." << std::endl;
        bHasPreAF = true;
    }

    //init_AsyncBuffer();
    num_blocks = marker->count_extract() / Constants::NUM_MARKER_READ +
        (marker->count_extract() % Constants::NUM_MARKER_READ != 0);

    
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
    if (!o_freq) { LOGGER.e(0, "cannot open the file [" + name_frq + "] to write"); }
    vector<string> out_contents;
    out_contents.reserve(AFA1.size() + 1);
    out_contents.push_back("CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS");
    for(int i = 0; i != AFA1.size(); i++){
        out_contents.push_back(marker->get_marker(marker->getRawIndex(i)) + "\t" + to_string(AFA1[i])
                               + "\t" + to_string(countMarkers[i]));
    }
    std::copy(out_contents.begin(), out_contents.end(), std::ostream_iterator<string>(o_freq, "\n"));
    o_freq.close();
    LOGGER.i(0, "Allele frequencies of " + to_string(AFA1.size()) + " SNPs have been saved in the file [" + name_frq + "]");
}

bool Geno::check_bed(){
    bool has_error = false;
    FILE *pFile;
    uint64_t f_size;
    uint8_t buffer[3];
    string message;
    uint32_t previous_size  = 0;

    for(int i = 0; i < geno_files.size(); i++){
        string bed_file = geno_files[i];
        uint32_t cur_size =  marker->count_raw(i);

        pFile = fopen(bed_file.c_str(), "rb");
        if(pFile == NULL){
            has_error = true;
            message += "Can't open [" + bed_file + "] to read.\n";
            previous_size = cur_size;
            continue;
        }
        fseek(pFile, 0, SEEK_END);
        f_size = ftell(pFile);
        rewind(pFile);

        if((f_size - 3) != ((uint64_t)num_byte_per_marker) * (cur_size - previous_size)){
            has_error = true;
            message += "Invalid bed file [" + bed_file +
                "]. The sample and SNP number in bed file are different from bim and fam file.\n";
            previous_size = cur_size;
            continue;
        }

        size_t read_count = fread(buffer, 1, 3, pFile);
        fclose(pFile);
        if((read_count != 3) &&
                (*buffer != 0x6c) &&
                (*(buffer+1) != 0x1b) &&
                (*(buffer+2) != 0x01)){
            has_error = true;
            message += "Invalid bed file [" + bed_file +
                "], please convert it into new format (SNP major).\n";
        }
        previous_size = cur_size;
    }

    //delete[] buffer;
    if(has_error){
        LOGGER.e(0, message);
    }else{
        LOGGER.i(0, "BED file(s) check OK.");
    }
    return has_error;
}

void Geno::read_bed(const vector<uint32_t> &raw_marker_index){

    // init start index for each file
    vector<int32_t> pos;
    pos.push_back(-1);
    for(int i = 0; i != geno_files.size() - 1; i++){
        pos.push_back(marker->count_raw(i) - 1);
    }

    //init files handles;
    vector<FILE *> pFiles;
    for(auto & cur_bed_file : geno_files){
        FILE *pFile = fopen(cur_bed_file.c_str(), "rb");
        if(pFile == NULL){
            LOGGER.e(0, "can't open [" + cur_bed_file + "] to read.");
        }

        fseek(pFile, 3, SEEK_SET);
        pFiles.push_back(pFile);
    }

    uint8_t *w_buf = NULL;
    int num_marker_read = 0;
    bool bNewWrite = true;
    for(auto & cur_marker_index : raw_marker_index){
        if(bNewWrite){
            w_buf = asyncBuffer->start_write();
            bNewWrite = false;
        }
        int cur_file_index = marker->getMIndex(cur_marker_index);
        FILE * pFile = pFiles[cur_file_index];
        int32_t lag_index = cur_marker_index - pos[cur_file_index];
        //very arbitary number to skip
        if(lag_index > 10){
            fseek(pFile, ((uint64_t)lag_index - 1) * num_byte_per_marker, SEEK_CUR);
        }else{
            for(int32_t ab_index = 1; ab_index < lag_index; ab_index++){
                if(fread(w_buf, 1, num_byte_per_marker, pFile) != num_byte_per_marker){
                    LOGGER.e(0, "error in reading data from [" + geno_files[cur_file_index] + "].\nThere might be some problems with your storage, or have you changed the files?");
                }
            }
        }


        size_t read_count = fread(w_buf, 1, num_byte_per_marker, pFile);
        if(read_count != num_byte_per_marker){
            LOGGER.e(0, "error in reading data from [" + geno_files[cur_file_index] + "].\nThere might be some problems with your storage, or the files have been changed?");
        }
        w_buf += num_byte_per_marker;
        pos[cur_file_index] = cur_marker_index;

        num_marker_read += 1;
        if(num_marker_read == Constants::NUM_MARKER_READ){
            asyncBuffer->end_write();
            bNewWrite = true;
            num_marker_read = 0;
        }
    }

    if(!bNewWrite){
        asyncBuffer->end_write();
    }

    for(auto & pFile : pFiles){
        fclose(pFile);
    }

}

void Geno::preGenoDouble(int numMarkerBuf, bool bMakeGeno, bool bGenoCenter, bool bGenoStd, bool bMakeMiss){
    sampleKeepIndex = pheno->get_index_keep();
    keepSampleCT = sampleKeepIndex.size();
    rawSampleCT = pheno->count_raw();
    numMarkerBlock = numMarkerBuf;
    //sex mode
    keepSexIndex = pheno->getSexValidRawIndex();
    keepMaleIndex = pheno->getMaleRawIndex();
    keepMaleExtractIndex = pheno->getMaleExtractIndex();
    keepSexSampleCT = keepSexIndex.size();
    keepMaleSampleCT = keepMaleIndex.size();

    this->bMakeGeno = bMakeGeno;
    this->bGenoCenter = bGenoCenter;
    this->bGenoStd = bGenoStd;
    this->bMakeMiss = bMakeMiss;

    numMarkersReadBlocks.resize(3);
    isMarkersSexXYs.resize(3);
    fileIndexBuf.resize(3);
 
    (this->*preGenoDoubleFuncs[genoFormat])();
    
    //init base SNP each file for read
    baseIndexLookup.clear();
    baseIndexLookup.push_back(0);
    int32_t sumIndex = 0;
    for(int i = 0; i < geno_files.size() - 1; i++){
        sumIndex += rawCountSNPs[i];
        baseIndexLookup.push_back(sumIndex);
    }

    missPtrSize = PgenReader::GetSubsetMaskSize(keepSampleCT);
}

//void Geno::loopDouble(const vector<uint32_t> &extractIndex, )

void Geno::preGenoDouble_pgen(){
    hasInfo = false;
    compressFormats.clear();
    rawCountSamples.clear();
    rawCountSNPs.clear();
 
    for(int i = 0; i < geno_files.size(); i++){
        MarkerParam curParam = marker->getMarkerParams(i); 
        compressFormats.push_back(curParam.compressFormat);
        rawCountSamples.push_back(rawSampleCT);
        rawCountSNPs.push_back(curParam.rawCountSNP);
    }

    // raw genotype buffer size, alligned
    pgenGenoPtrSize = (PgenReader::GetGenoBufPtrSize(keepSampleCT) + 63) / 64 * 64;
    pgenDosageMainPtrSize = (PgenReader::GetDosageMainSize(keepSampleCT) + 63)/64 * 64; 
    pgenDosagePresentPtrSize = (PgenReader::GetDosagePresentSize(keepSampleCT) + 63)/64 * 64;

    pgenGenoBuf1PtrSize = (pgenGenoPtrSize + pgenDosageMainPtrSize + pgenDosagePresentPtrSize + 1 + 63) /64 * 64;
    asyncBuf64 = new AsyncBuffer<uintptr_t>(pgenGenoBuf1PtrSize * numMarkerBlock);
    if(!asyncBuf64->init_status()){
        LOGGER.e(0, "can't allocate enough memory to read genotype.");
    }
}

void Geno::preGenoDouble_bed(){
    hasInfo = false;
    numBytePerMarker = (rawSampleCT + 3) / 4; 
    //checking the bed files

    compressFormats.clear();
    rawCountSamples.clear();
    rawCountSNPs.clear();
 
    for(int i = 0; i < geno_files.size(); i++){
        MarkerParam curParam = marker->getMarkerParams(i); 
        compressFormats.push_back(curParam.compressFormat);
        rawCountSamples.push_back(rawSampleCT);
        rawCountSNPs.push_back(curParam.rawCountSNP);
    }

    // raw genotype buffer size
    uint32_t raw_sample_ct = rawSampleCT;
    bedRawGenoBuf1PtrSize = PgenReader::GetGenoBufPtrSize(raw_sample_ct);
    asyncBuf64 = new AsyncBuffer<uintptr_t>(bedRawGenoBuf1PtrSize * numMarkerBlock);
    if(!asyncBuf64->init_status()){
        LOGGER.e(0, "can't allocate enough memory to read genotype.");
    }
 

    maskPtrSize = PgenReader::GetSubsetMaskSize(raw_sample_ct);
    keepMaskPtr = new uintptr_t[maskPtrSize]; 
    keepMaskInterPtr = new uintptr_t[maskPtrSize];
    PgenReader::SetSampleSubsets(sampleKeepIndex, raw_sample_ct, keepMaskPtr, keepMaskInterPtr);

    sexMaskPtr = new uintptr_t[maskPtrSize];
    sexMaskInterPtr = new uintptr_t[maskPtrSize];
    PgenReader::SetSampleSubsets(keepSexIndex, raw_sample_ct, sexMaskPtr, sexMaskInterPtr);

    maleMaskPtr = new uintptr_t[maskPtrSize];
    maleMaskInterPtr = new uintptr_t[maskPtrSize];
    PgenReader::SetSampleSubsets(keepMaleIndex, raw_sample_ct, maleMaskPtr, maleMaskInterPtr);

    // for missing pointer size of 1 genotype
}

void Geno::preGenoDouble_bgen(){
    hasInfo = true;

    compressFormats.clear();
    rawCountSamples.clear();
    rawCountSNPs.clear();
    for(int i = 0; i < geno_files.size(); i++){
        MarkerParam curParam = marker->getMarkerParams(i); 
        if(curParam.rawCountSample != pheno->count_raw()){
            LOGGER.e(0, "inconsistent sample sizes between the .bgen file [" + geno_files[i] + "] and the .sample file (specified by --sample).");
        }
        compressFormats.push_back(curParam.compressFormat);
        rawCountSamples.push_back(curParam.rawCountSample);
        rawCountSNPs.push_back(curParam.rawCountSNP);
    }


    bgenRawGenoBuf1PtrSize = marker->getMaxGenoMarkerUptrSize();
    asyncBuf64 = new AsyncBuffer<uintptr_t>(bgenRawGenoBuf1PtrSize * numMarkerBlock);
    if(!asyncBuf64->init_status()){
        LOGGER.e(0, "can't allocate enough memory to read genotype.");
    }
 

}

int nextBufIndex(int curIndex){
    return (curIndex + 1) % 3;
}


void Geno::readGeno(const vector<uint32_t> &extractIndex){
    (this->*readGenoFuncs[genoFormat])(extractIndex);
}

void Geno::readGeno_bed(const vector<uint32_t> &extractIndex){
    const vector<uint32_t> raw_marker_index = marker->get_extract_index();
    vector<uint32_t> rawIndices(extractIndex.size());
    std::transform(extractIndex.begin(), extractIndex.end(), rawIndices.begin(), 
            [&raw_marker_index](size_t pos){return raw_marker_index[pos];});

    uintptr_t *g_buf = NULL;
    uint32_t numMarker = extractIndex.size();
    uint32_t finishedMarker = 0;
    uint32_t nextSize;
    int preFileIndex = -1;
    int fileIndex = 0;
    int base_index = baseIndexLookup[fileIndex];
    PgenReader reader;
    bool chr_ends;
    uint8_t isSexXY;
    int curWriteBufIndex = 0;
    //std::ofstream oidx("rawidx.snplist");
    while(finishedMarker != numMarker && (nextSize = marker->getNextSize(rawIndices, finishedMarker, numMarkerBlock,fileIndex, chr_ends, isSexXY)) != 0){
        g_buf = asyncBuf64->start_write();
        for(int i = 0; i < nextSize; i++){
            int processIndex = finishedMarker + i;
            int rawIndex = rawIndices[processIndex];
            //oidx << rawIndex << "\n";
            int curExtractIndex = extractIndex[processIndex];
            if(preFileIndex != fileIndex){
                //LOGGER << "reading " << fileIndex << ", sample: " << rawCountSamples[fileIndex] << ", marker: " << rawCountSNPs[fileIndex] << std::endl;
                reader.Load(geno_files[fileIndex], &rawCountSamples[fileIndex], &rawCountSNPs[fileIndex], sampleKeepIndex);
                base_index = baseIndexLookup[fileIndex];
                preFileIndex = fileIndex;
            }
            int lag_index = rawIndex - base_index;
            //int al_idx = marker->isEffecRevRaw(rawIndex) ? 0 : 1;
            reader.ReadRawFullHard(g_buf, lag_index);

            g_buf += bedRawGenoBuf1PtrSize;
        }

        finishedMarker += nextSize;
        numMarkersReadBlocks[curWriteBufIndex] = nextSize;
        isMarkersSexXYs[curWriteBufIndex] = isSexXY;
        fileIndexBuf[curWriteBufIndex] = fileIndex;
        asyncBuf64->end_write();
        curWriteBufIndex = nextBufIndex(curWriteBufIndex);
    }
    //oidx.close();
}

void Geno::readGeno_pgen(const vector<uint32_t> &extractIndex){
    const vector<uint32_t> raw_marker_index = marker->get_extract_index();
    vector<uint32_t> rawIndices(extractIndex.size());
    std::transform(extractIndex.begin(), extractIndex.end(), rawIndices.begin(), 
            [&raw_marker_index](size_t pos){return raw_marker_index[pos];});

    uintptr_t *g_buf = NULL;
    uint32_t numMarker = extractIndex.size();
    uint32_t finishedMarker = 0;
    uint32_t nextSize;
    int preFileIndex = -1;
    int fileIndex = 0;
    int base_index = baseIndexLookup[fileIndex];
    PgenReader reader;
    bool chr_ends;
    uint8_t isSexXY;
    int curWriteBufIndex = 0;
    while((finishedMarker != numMarker) && (nextSize = marker->getNextSize(rawIndices, finishedMarker, numMarkerBlock,fileIndex, chr_ends, isSexXY)) != 0){
        g_buf = asyncBuf64->start_write();
        //LOGGER << "before: " << (void*)g_buf << std::endl;
        for(int i = 0; i < nextSize; i++){
            int processIndex = finishedMarker + i;
            int rawIndex = rawIndices[processIndex];
            int curExtractIndex = extractIndex[processIndex];
            if(preFileIndex != fileIndex){
                reader.Load(geno_files[fileIndex], &rawCountSamples[fileIndex], &rawCountSNPs[fileIndex], sampleKeepIndex);
                base_index = baseIndexLookup[fileIndex];
                preFileIndex = fileIndex;
            }
            int lag_index = rawIndex - base_index;
            int al_idx = marker->isEffecRevRaw(rawIndex) ? 0 : 1;
            reader.ReadDosage(g_buf, lag_index, al_idx);
            g_buf += pgenGenoBuf1PtrSize;
        }

        finishedMarker += nextSize;
        numMarkersReadBlocks[curWriteBufIndex] = nextSize;
        isMarkersSexXYs[curWriteBufIndex] = isSexXY;
        fileIndexBuf[curWriteBufIndex] = fileIndex;
        asyncBuf64->end_write();
        curWriteBufIndex = nextBufIndex(curWriteBufIndex);
    }
}

void Geno::getGenoDouble(uintptr_t *buf, int bufIndex, GenoBufItem* gbuf){
    (this->*getGenoDoubleFuncs[genoFormat])(buf, bufIndex, gbuf);
}

void Geno::setGenoItemSize(uint32_t &genoSize, uint32_t &missSize){
    genoSize = keepSampleCT;
    missSize = missPtrSize;
}

void Geno::getGenoDouble_pgen(uintptr_t *buf, int idx, GenoBufItem* gbuf){
    SNPInfo snpinfo;
    uintptr_t *cur_buf = buf + idx * pgenGenoBuf1PtrSize;
    uintptr_t *geno_buf = cur_buf;
    uintptr_t *dosage_present = cur_buf + pgenGenoPtrSize;
    uint16_t *dosage_main = reinterpret_cast<uint16_t*>(cur_buf + pgenGenoPtrSize + pgenDosagePresentPtrSize);
    uint32_t dosage_ct = cur_buf[pgenGenoBuf1PtrSize - 1];
    uint8_t isSexXY = isMarkersSexXYs[curBufferIndex];

    const vector<uint32_t> *curMaleIndex = NULL;
    if(isSexXY == 1){
        curMaleIndex = &keepMaleExtractIndex;
    }

    string err;
    if(!PgenReader::CountHardDosage(cur_buf, dosage_main, dosage_present, curMaleIndex, keepSampleCT, dosage_ct, &snpinfo, err)){
        LOGGER.e(0, err);
    }

    double af = snpinfo.af;
    double std = snpinfo.std;
    if(bHasPreAF){
        af = AFA1[gbuf->extractedMarkerIndex];
        std = 2.0 * af * (1.0 - af);
    }
    double maf = std::min(af, 1.0 - af);
    if(maf >= min_maf && maf <= max_maf && snpinfo.nMissRate >= dFilterMiss){
        gbuf->valid = true;
        gbuf->af = af;
        gbuf->nValidN = snpinfo.N;
        gbuf->nValidAllele = snpinfo.AlCount;
        gbuf->mean = 2.0 * af;
        gbuf->sd = std;

        if(bMakeGeno){
            // get get genotype in lookup table
            
            if(isSexXY == 1){
                double weight;
                bool needWeight;
                setMaleWeight(weight, needWeight);
                if(needWeight){
                    for(int i = 0 ; i < keepMaleSampleCT; i++){
                        gbuf->geno[keepMaleExtractIndex[i]] *= weight;
                    }
                }
            }
        }
        if(bMakeMiss){
            gbuf->missing.resize(missPtrSize, 0);
        }

    }
}

void Geno::getGenoDouble_bed(uintptr_t *buf, int idx, GenoBufItem* gbuf){
    SNPInfo snpinfo;
    uintptr_t *cur_buf = buf + idx * bedRawGenoBuf1PtrSize;
    uint8_t isSexXY = isMarkersSexXYs[curBufferIndex];
    bool hasNoHET = true;
    if(isSexXY != 1){
        PgenReader::CountHardFreqMissExt(cur_buf, keepMaskInterPtr, rawSampleCT, keepSampleCT, &snpinfo, f_std);
    }else{
        string errmsg;
        hasNoHET = PgenReader::CountHardFreqMissExtX(cur_buf, keepMaskInterPtr, maleMaskInterPtr, rawSampleCT, keepSampleCT, keepMaleSampleCT, &snpinfo, errmsg, iDC==1, f_std);
        /*
            if(errmsg == "het"){
                LOGGER.e(0, "found heterozygote coding (=1) on ChrX (23) for male samples. GCTA treats ChrX (23) as non-PAR and codes male genotype on ChrX as either 0 or 2. \nPlease check the gender information in your genotype file, or use PLINK --split-x to separate the ChrX into non-PAR (chr23) and PAR (chr25) regions for further analysis.\nThe non-PAR region is assumed to have the dosage compensation issue, and you can specify it with the --dc flag, with 0: no compensation, and 1: dosage compensation (the default is 1 in GCTA).");
            }
        }
        */
    }
    uint32_t curExtractIndex = gbuf->extractedMarkerIndex;
    bool isEffRev = marker->isEffecRev(curExtractIndex);
    double af = isEffRev ? (1.0 - snpinfo.af) : snpinfo.af;
    if(bHasPreAF){
        af = AFA1[curExtractIndex];
        snpinfo.mean = 2 * af;
        snpinfo.std = 2 * af * (1.0 - af); 
    }
    double maf = std::min(af, 1.0 - af);
    if(maf >= min_maf && maf <= max_maf){
        if(snpinfo.nMissRate >= dFilterMiss){
            gbuf->valid = true;
            gbuf->af = af;
            gbuf->nValidN = snpinfo.N;
            gbuf->nValidAllele = snpinfo.AlCount;

            if(bGRM){
                gbuf->mean = 2.0 * af;
            }else{
                gbuf->mean = snpinfo.mean;
            }
            if(f_std){
                gbuf->sd = snpinfo.std;
            }else{
                gbuf->sd = gbuf->mean * (1.0 - af);
            }
            if(bMakeGeno){
                double mu = gbuf->mean;
                double sd = gbuf->sd;
                if(sd < 1.0e-50){
                    gbuf->valid = false;
                    return;
                }

                double center_value = 0.0;
                double rdev = 1.0;
                double a0, a1, a2, na;

                if(!bGRMDom){
                    if(bGenoCenter){
                        center_value = mu;
                    }
                    if(bGenoStd){
                        rdev = sqrt(1.0 / sd);
                    }
                    double aa0 = 0.0, aa2 = 2.0;
                    if(isEffRev){
                        double temp = aa0;
                        aa0 = aa2;
                        aa2 = aa0;
                    }

                    a0 = (aa0 - center_value) * rdev;
                    a1 = (1.0 - center_value) * rdev;
                    a2 = (aa2 - center_value) * rdev;
                    na = (mu - center_value) * rdev;
               }else{
                   double psq = 0.5 * mu * mu;
                   if(bGenoCenter)center_value = psq; // psq
                   if(bGenoStd){
                       rdev = 1.0 / sd;
                   }
                   double aa0 = 0.0, aa2 = 2.0 * mu - 2.0;
                   if(isEffRev){
                       double temp = aa0;
                       aa0 = aa2;
                       aa2 = temp;
                   }
                   a0 = (aa0 - center_value) * rdev;
                   a1 = (mu - center_value) * rdev;
                   a2 = (aa2 - center_value) * rdev;
                   na = (psq - center_value)*rdev;
                }

                const double lookup[32] __attribute__ ((aligned (16))) = GET_TABLE16(a0, a1, a2, na);
                gbuf->geno.resize(keepSampleCT);
                uintptr_t * pmiss = NULL;
                if(bMakeMiss){
                    gbuf->missing.resize(missPtrSize); 
                    pmiss = gbuf->missing.data();
                }
                PgenReader::ExtractDoubleExt(cur_buf, keepMaskPtr, rawSampleCT, keepSampleCT, lookup, gbuf->geno.data(), pmiss); 
                // adjust for chr X;
                if(isSexXY == 1){
                    /* don't set to missing
                    if(!hasNoHET){
                        for(int i = 0 ; i < keepMaleSampleCT; i++){
                            uint32_t curMaleIndex = keepMaleExtractIndex[i];
                            if(gbuf->geno[curMaleIndex] == a1){
                                gbuf->geno[curMaleIndex] = na;
                            }
                        }
                    }
                    */
                    double weight;
                    bool needWeight;
                    setMaleWeight(weight, needWeight);
                    if(needWeight){
                        if(bGRM){
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                gbuf->geno[keepMaleExtractIndex[i]] *= weight;
                            }
                        }else{
                            double correctWeight = (weight - 1) * rdev * center_value;
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                uint32_t curIndex = keepMaleExtractIndex[i];
                                gbuf->geno[curIndex] *= weight;
                                gbuf->geno[curIndex] += correctWeight;
                            }
                        }
                    }

                    
                }
            }
            return;
        }
    }
    gbuf->valid = false;
}

void Geno::readGeno_bgen(const vector<uint32_t> &extractIndex){
    const vector<uint32_t> raw_marker_index = marker->get_extract_index();
    vector<uint32_t> rawIndices(extractIndex.size());
    std::transform(extractIndex.begin(), extractIndex.end(), rawIndices.begin(), 
            [&raw_marker_index](size_t pos){return raw_marker_index[pos];});

    openGFiles();

    uintptr_t *g_buf = NULL;
    uint32_t numMarker = extractIndex.size();
    uint32_t finishedMarker = 0;
    uint32_t nextSize;
    int fileIndex = 0;
    bool chr_ends;
    uint8_t isSexXY;
    int curWriteBufIndex = 0;
    while(finishedMarker != numMarker && (nextSize = marker->getNextSize(rawIndices, finishedMarker, numMarkerBlock,fileIndex, chr_ends, isSexXY)) != 0){
        g_buf = asyncBuf64->start_write();
        FILE *bgenFile = gFiles[fileIndex];
        for(int i = 0; i < nextSize; i++){
            int processIndex = finishedMarker + i;
            int rawIndex = rawIndices[processIndex];

            uint64_t pos, size;
            marker->getStartPosSize(rawIndex, pos, size);
            fseek(bgenFile, pos, SEEK_SET);
            if(fread(g_buf, sizeof(char), size, bgenFile) != size){
                int lag_index = rawIndex - baseIndexLookup[fileIndex];
                LOGGER.e(0, "can't read " + to_string(lag_index) + "th SNP in [" + geno_files[fileIndex] + "].");
            }
            g_buf += bgenRawGenoBuf1PtrSize;
        }

        finishedMarker += nextSize;
        numMarkersReadBlocks[curWriteBufIndex] = nextSize;
        isMarkersSexXYs[curWriteBufIndex] = isSexXY;
        fileIndexBuf[curWriteBufIndex] = fileIndex;
        asyncBuf64->end_write();
        curWriteBufIndex = nextBufIndex(curWriteBufIndex);
    }
}

void Geno::setMaleWeight(double &weight, bool &needWeight){
    weight = 1.0;
    if(bGRM){ // GRM 
        weight = sqrt(0.5);
        if(iGRMdc == 1){
            weight *= sqrt(2.0);
        }else if(iGRMdc == 0){
            weight *= sqrt(0.5);
        }
    }else{
        if(iDC == 0){
            weight = 0.5;
        }
        /*
           if(!bGenoStd){
           if(iDC == 0){
           weight = 0.5;
           }
           }else{
           if(iDC == 0){
           weight = sqrt(0.5);
           }
           }
           */
    }
    if(std::abs(weight - 1.0) > 1e-6){
        needWeight = true;
    }else{
        needWeight = false;
    }
}

void calDosage_bgen(uint32_t prob1, uint32_t prob2, uint64_t &dosage, uint32_t &prob1d){
    prob1d = prob1 * 2;
    dosage = prob1d + prob2;
}

void calDosagePhase_bgen(uint32_t prob1, uint32_t prob2, uint64_t &dosage, uint32_t &prob1d){
    dosage = prob1 + prob2;
    prob1d = 2 * prob1 * prob2;
}


void Geno::getGenoDouble_bgen(uintptr_t *buf, int idx, GenoBufItem* gbuf){
    SNPInfo snpinfo;
    uintptr_t *cur_buf = buf + idx * bgenRawGenoBuf1PtrSize;
    // skip the header
    uint8_t *curbuf = (uint8_t*)cur_buf;
    int fileIndex = fileIndexBuf[curBufferIndex];

    int compressFormat = compressFormats[fileIndex];

    uint16_t L16;
    uint32_t L32;
    //skip Lid
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16) + L16;

    //skip rs
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16) + L16;

    //skip chr
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16) + L16;

    //skip pos
    curbuf += sizeof(uint32_t);

        // skip n allels
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16);
    for(int i = 0; i < L16; i++){
        memcpy(&L32, curbuf, sizeof(L32));
        curbuf += sizeof(L32) + L32;
    }

    uint32_t len_comp, len_decomp;
    memcpy(&len_comp, curbuf, sizeof(len_comp));
    curbuf += sizeof(len_comp);

    if(compressFormat == 0){
        len_decomp = len_comp;
    }else{
        len_comp -= 4;
        memcpy(&len_decomp, curbuf, sizeof(len_decomp));
        curbuf += sizeof(len_decomp);
    }

    string error_promp = to_string(gbuf->extractedMarkerIndex) + "th SNP of [" + geno_files[fileIndex] + "]."; 
    uint8_t *dec_data;
    if(compressFormat != 0){
        dec_data = new uint8_t[len_decomp + 8];
        uint32_t curCompSize = len_comp;
        if(compressFormat == 1){
            uint32_t Ldecomp = len_decomp;
            int z_result = uncompress((Bytef*)dec_data, (uLongf*)&Ldecomp, (Bytef*)curbuf, curCompSize);
            if(z_result != Z_OK || len_decomp != Ldecomp){
                LOGGER.e(0, "decompressing genotype data error in " + error_promp); 
            }
        }else if(compressFormat == 2){
            //zstd  
            uint64_t const rSize = ZSTD_getFrameContentSize((void*)curbuf, curCompSize);
            switch(rSize){
                case ZSTD_CONTENTSIZE_ERROR:
                    LOGGER.e(0, "not compressed by zstd in " + error_promp);
                    break;
                case ZSTD_CONTENTSIZE_UNKNOWN:
                    LOGGER.e(0, "original size unknown in " + error_promp);
                    break;
            }
            if(rSize != len_decomp){
                LOGGER.e(0, "size stated in the compressed file is different from " + error_promp);
            }
            size_t const dSize = ZSTD_decompress((void *)dec_data, len_decomp, (void*)curbuf, curCompSize); 

            if(ZSTD_isError(dSize)){
                LOGGER.e(0, "decompressing genotype error: " + string(ZSTD_getErrorName(dSize)) + " in " + error_promp);
            }
        }else{
            LOGGER.e(0, "unknown compress format in " + error_promp);
        }
    }else{
        dec_data = curbuf;
    }

    uint32_t n_sample = *(uint32_t *)dec_data;
    if(n_sample != rawCountSamples[fileIndex]){
        LOGGER.e(0, "inconsistent number of individuals in " + error_promp);
    }
    uint16_t num_alleles = *(uint16_t *)(dec_data + 4);
    if(num_alleles != 2){
        LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
    }

    uint8_t min_ploidy = *(uint8_t *)(dec_data + 6);//2
    uint8_t max_ploidy = *(uint8_t *)(dec_data + 7); //2
    uint8_t * sample_ploidy = (uint8_t *)(dec_data + 8);
    //check all ploidy are 2 or not
    if(min_ploidy != 2){
        LOGGER.e(0, "multiploidy detected in " + error_promp);
    }

    uint8_t *geno_prob = sample_ploidy + n_sample;
    uint8_t is_phased = *(geno_prob);
    uint8_t bits_prob = *(geno_prob+1);
    uint8_t* X_prob = geno_prob + 2;
    uint32_t len_prob = len_decomp - n_sample - 10;
    void (*calFunc)(uint32_t, uint32_t, uint64_t&, uint32_t &prob1d);
    if(is_phased){
        calFunc = &calDosagePhase_bgen;
    }else{
        calFunc = &calDosage_bgen;
    }

    uint8_t double_bits_prob = bits_prob * 2;
    vector<uint32_t> miss_index;

    uint8_t isSexXY = isMarkersSexXYs[curBufferIndex];

    uint64_t mask = (1U << bits_prob) - 1;
    uint64_t dosage_sum = 0, fij_sum = 0, dosage2_sum = 0;
    uint32_t validN = 0;
    uint32_t validAllele = 0;
    vector<uint32_t> dosages(keepSampleCT);
    uint32_t max_dos = mask * 2 + 1;
    bool has_miss = false;

    uint32_t curSampleCT = keepSampleCT;
    vector<uint32_t> *curSampleIndexPtr = &sampleKeepIndex;
    /* // use default female mode
    if(isSexXY == 1){
        curSampleCT = keepSexSampleCT;
        curSampleIndexPtr = &keepSexIndex;
    }
    */

    for(int j = 0; j < curSampleCT; j++){
        uint32_t sindex = (*curSampleIndexPtr)[j];
        uint8_t item_ploidy = sample_ploidy[sindex];
        if(item_ploidy > 128){
            miss_index.push_back(sindex);
            has_miss = true;
            dosages[j] = max_dos;
        }else if(item_ploidy == 2){
            uint32_t start_bits = sindex * double_bits_prob;
            uint64_t geno_temp;
            memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(geno_temp));
            geno_temp = geno_temp >> (start_bits % CHAR_BIT);
            uint32_t prob1 = geno_temp & mask;
            uint32_t prob2 = (geno_temp >> bits_prob) & mask;
            /*
            uint32_t prob1d = prob1 * 2;
            uint64_t dosage = prob1d + prob2;
            */
            uint32_t prob1d;
            uint64_t dosage;
            calFunc(prob1, prob2, dosage, prob1d);
            dosages[j] = dosage;
            dosage_sum += dosage;
            dosage2_sum += dosage * dosage;

            //uint64_t fij = dosage + prob1d;
            //fij_sum += fij;
            fij_sum += prob1d;
            validN++;
            validAllele += 2;
        }else{
            LOGGER.e(0, "multiploidy detected in " + error_promp);
        }
    }

    
    double dosage_sum_half = dosage_sum;
    double dosage2_sum_half = dosage2_sum;
    if(isSexXY == 1){
        for(int j = 0; j < keepMaleSampleCT; j++){
            uint32_t sindex = keepMaleIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(geno_temp));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;
                // check prob2 equal to 0 in male?
                /*
                uint32_t prob1d = prob1 * 2;
                uint64_t dosage = prob1d + prob2;
                */
                //dosages[j] = dosage;
                uint32_t prob1d;
                uint64_t dosage;
                calFunc(prob1, prob2, dosage, prob1d);
                uint32_t prob1_true = prob1d / 2;
                dosage_sum_half -= prob1_true;
                dosage2_sum_half -= (dosage * dosage - (uint64_t)prob1_true * prob1_true);
                //dosage_sum -= (prob1 + prob2);
                //dosage2_sum -= (dosage * dosage - (uint64_t)prob1 * prob1);
                //
                //fij_sum -= prob1_true;
                validAllele--;
            }
        }
    }


    if(compressFormat != 0){
        delete[] dec_data;
    }


    double maskd = (double)mask;
    double af = (double)dosage_sum_half / maskd / validAllele;
    double mean;
    bool bEffRev = this->marker->isEffecRev(gbuf->extractedMarkerIndex);
    if(bEffRev){
        af = 1.0 - af;
    }
    double std = 2.0 * af * (1.0 - af);
    double info = 0.0;
    double mask2 = mask * mask;
    if(std < 1e-50){
        info = 1.0;
    }else{
        double dos2_fij_sum = (double)(dosage_sum + fij_sum)/ maskd - (double)dosage2_sum / mask2; 
        info = 1.0 - dos2_fij_sum / (std * validN);
    }
    if(is_phased){
        info = 1;
    }

    if(bHasPreAF){
        af = AFA1[gbuf->extractedMarkerIndex];
        mean = 2.0 * af;
        std = 2.0 * af * (1.0 - af);
    }else{
        if(iDC == 1 && (!bGRM)){
            double dos_double = (double)dosage_sum / maskd;
            mean = dos_double / validN;
            std = ((double)dosage2_sum / mask2 - dos_double * mean)/(validN - 1);
        }else{
            double dos_double = (double)dosage_sum_half / maskd;
            mean = dos_double / validN;
            std = ((double)dosage2_sum_half / mask2 - dos_double * mean)/(validN - 1);
        }
 
    }

    double maf = std::min(af, 1.0 - af);
    if(maf >= min_maf && maf <= max_maf){
        double nMissRate = 1.0*validN / curSampleCT;
        //LOGGER << "dFilterMiss: " << dFilterMiss << "*" << nMissRate << ", info: " << dFilterInfo << "*" << info << ", " << std << std::endl; 
        if(nMissRate >= dFilterMiss && info >= dFilterInfo){
            gbuf->valid = true;
            gbuf->af = af;
            gbuf->nValidN = validN;
            gbuf->nValidAllele = validAllele;
            gbuf->info = is_phased ? (std::numeric_limits<double>::quiet_NaN()) : info;
            gbuf->mean = mean;
            gbuf->sd = std;
            if(bMakeGeno){
                double mu = gbuf->mean;
                if(std < 1.0e-50){
                    gbuf->valid = false;
                    return;
                }

                double* dos_lookup = new double[max_dos + 2];
                double center_value = 0.0;
                double rdev = 1.0;
                if(!bGRMDom){
                    if(bGenoCenter){
                        center_value = mu;
                    }
                    if(bGenoStd){
                        rdev = sqrt(1.0 / std);
                    }

                    for(uint32_t i = 0; i < max_dos; i++){
                        double tdos = (double)i / mask;
                        tdos = bEffRev ? (2.0 - tdos) : tdos;
                        dos_lookup[i] = (tdos - center_value) *rdev;
                    }
                    dos_lookup[max_dos] = ( mu - center_value) * rdev;
                }else{// dominance
                    uint32_t cut05 = ceil(0.5 * mask);
                    uint32_t cut15 = ceil(1.5 * mask);
                    double psq = 0.5 * mu * mu;
                    double dos00 = bEffRev ? (2.0 * mu - 2.0) : 0.0;
                    double dos01 = mu;
                    double dos10 = bEffRev ? 0.0 : (2.0 * mu - 2.0);
                    double dosna = psq;

                    if(bGenoCenter)center_value = psq; // psq
                    if(bGenoStd){
                        rdev = 1.0 / std;
                    }
                    double a0 = (dos00 - center_value) * rdev;
                    double a1 = (dos01 - center_value) * rdev;
                    double a2 = (dos10 - center_value) * rdev;
                    double na = (dosna - center_value) * rdev;
                    
                    for(uint32_t i = 0; i < cut05; i++){
                        dos_lookup[i] = a0;
                    }
                    for(uint32_t i = cut05; i < cut15; i++){
                        dos_lookup[i] = a1;
                    }
                    for(uint32_t i = cut15; i < max_dos; i++){
                        dos_lookup[i] = a2;
                    }

                    dos_lookup[max_dos] = (dosna - center_value) * rdev;
                }

                gbuf->geno.resize(curSampleCT);
                for(int j = 0; j < curSampleCT; j++){
                    gbuf->geno[j] = dos_lookup[dosages[j]];
                }
                delete[] dos_lookup;
                // adjust for chr X;
                if(isSexXY == 1){
                    double weight;
                    bool needWeight;
                    setMaleWeight(weight, needWeight);
                    if(needWeight){
                        if(bGRM){
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                gbuf->geno[keepMaleExtractIndex[i]] *= weight;
                            }
                        }else{
                            double correctWeight = (weight - 1) * rdev * center_value;
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                uint32_t curIndex = keepMaleExtractIndex[i];
                                gbuf->geno[curIndex] *= weight;
                                gbuf->geno[curIndex] += correctWeight;
                            }
 
                        }
                    }
                }
            }
            if(bMakeMiss){
                gbuf->missing.resize(missPtrSize, 0); 
                const int ptrsize = sizeof(uintptr_t) * CHAR_BIT;
                for(int j = 0; j < miss_index.size(); j++){
                    int cur_index = miss_index[j];
                    gbuf->missing[cur_index / ptrsize] |= (1UL << (cur_index % ptrsize));
                }
            }
            return;
        }
    }
    gbuf->valid = false;
 
}

void Geno::setGRMMode(bool grm, bool dominace){
    this->bGRM = grm;
    this->bGRMDom = dominace;
}

void Geno::endGenoDouble_bed(){
    delete asyncBuf64;
    delete[] keepMaskPtr;
    delete[] keepMaskInterPtr;

    delete[] sexMaskPtr;
    delete[] sexMaskInterPtr;
    delete[] maleMaskPtr;
    delete[] maleMaskInterPtr;
}

void Geno::endGenoDouble(){
    (this->*endGenoDoubleFuncs[genoFormat])();
}


void Geno::endGenoDouble_bgen(){
    delete asyncBuf64;
}

void Geno::endGenoDouble_pgen(){
    delete asyncBuf64;
}

bool Geno::getGenoHasInfo(){
    return hasInfo;
}

void Geno::loopDouble(const vector<uint32_t> &extractIndex, int numMarkerBuf, bool bMakeGeno, bool bGenoCenter, bool bGenoStd, bool bMakeMiss, vector<function<void (uintptr_t *buf, const vector<uint32_t> &exIndex)>> callbacks, bool showLog){
   
    preGenoDouble(numMarkerBuf, bMakeGeno, bGenoCenter, bGenoStd, bMakeMiss);
    thread read_thread([this, &extractIndex](){this->readGeno(extractIndex);});
    read_thread.detach();
    // main loop
    
    LOGGER.ts("LOOP_GENO_PRE");
    LOGGER.ts("LOOP_GENO_TOT");
    int nTMarker = extractIndex.size();
    uint32_t nFinishedMarker = 0;

    int pre_block = 0;
    curBufferIndex = 0;

    while(nFinishedMarker < nTMarker){
       uintptr_t *r_buf = NULL;
       bool isEOF = false;
       std::tie(r_buf, isEOF) = asyncBuf64->start_read();
       if(isEOF){
           LOGGER.e(0, "the reading process reached to the end of the genotype file but couldn't finish.");
       }
        
       int nMarker = numMarkersReadBlocks[curBufferIndex];
       uint32_t endIndex = nFinishedMarker + nMarker;
       endIndex = endIndex > nTMarker ? nTMarker : endIndex;
       vector<uint32_t> curExtractIndex(extractIndex.begin() + nFinishedMarker, 
              extractIndex.begin() + endIndex);

       for(auto callback : callbacks){
          callback(r_buf, curExtractIndex);
       }
       asyncBuf64->end_read();

       nFinishedMarker += nMarker;
       curBufferIndex = nextBufIndex(curBufferIndex);

        // show progress
       if(showLog){
           int cur_block = nFinishedMarker >> 14;
           if(cur_block > pre_block){
               pre_block = cur_block;
               float time_p = LOGGER.tp("LOOP_GENO_PRE");
               if(time_p > 300){
                   LOGGER.ts("LOOP_GENO_PRE");
                   float elapse_time = LOGGER.tp("LOOP_GENO_TOT");
                   float finished_percent = (float) nFinishedMarker / nTMarker;
                   float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                   std::ostringstream ss;
                   ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 

                   LOGGER.i(1, ss.str());
               }
           }
       }
    }
    //finished 
    if(showLog){
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << "100% finished in " << LOGGER.tp("LOOP_GENO_TOT") << " sec";
        LOGGER.i(1, ss.str());
        LOGGER << nFinishedMarker << " SNPs have been processed." << std::endl;
    }
    endGenoDouble();
}



void Geno::preProcess(GenoBuf *gbuf, int numMarkerBuf, vector<uint32_t> *rawMarkerIndex){
   sampleKeepIndex = pheno->get_index_keep();
   gbuf->n_sample = pheno->count_keep();
   gbuf->n_marker = numMarkerBuf;
   numMarkerBlock = gbuf->n_marker;
   if(rawMarkerIndex != NULL){
       asyncMode = true;
       rawMarkerIndexProceed = *rawMarkerIndex;
   }else{
       asyncMode = false;
   }
       
   openGFiles();
   (this->*preProcessFuncs[genoFormat])(gbuf);
   setGenoBufSize(gbuf, numMarkerBuf);
}



void Geno::endProcess(){
    (this->*endProcessFuncs[genoFormat])();
}

void Geno::getGenoArray(const vector<uint32_t>& raw_marker_index, GenoBuf *gbuf){
    (this->*getGenoArrayFuncs[genoFormat])(raw_marker_index, gbuf);
}

// this will fill the full buffer, exept chr end meet;
void Geno::getGenoArrayExtractFull(const vector<uint32_t> & fullExtractIndex, uint32_t start, uint32_t num, vector<uint32_t> &keptIndex, uint32_t &totalReadMarker, GenoBuf *gbuf){

    keptIndex.clear();
    keptIndex.reserve(num);
    uint32_t n_sample = gbuf->n_sample;
    uint32_t nByteGeno = n_sample * 8;

    bool chr_ends = false;
    totalReadMarker = 0;
    while(keptIndex.size() < num && (!chr_ends)){
        uint32_t startIndex = fullExtractIndex[start + totalReadMarker];
        bool isX;
        vector<uint32_t> indices = marker->getNextSizeIndex(startIndex, num - keptIndex.size(), chr_ends, isX);
        if(indices.size() == 0){
            break;
        }
        gbuf->bSex = isX;
        gbuf->indexStartMarker = keptIndex.size();
        getGenoArrayExtract(indices, gbuf);
        totalReadMarker += indices.size();
        int move_dis = 0;
        for(int i = gbuf->indexStartMarker; i < gbuf->indexStartMarker + indices.size(); i++){
            if(gbuf->usedIndex[i]){
                keptIndex.push_back(indices[i - gbuf->indexStartMarker]);
                if(move_dis != 0){
                    int moved_index = i - move_dis;
                    gbuf->af[moved_index] = gbuf->af[i];
                    gbuf->preAF[moved_index] = gbuf->preAF[i];
                    gbuf->nValidN[moved_index] = gbuf->nValidN[i];
                    gbuf->nValidAllele[moved_index] = gbuf->nValidAllele[i];
                    if(gbuf->hasInfo)gbuf->info[moved_index] = gbuf->info[i];
                    if(gbuf->saveGeno)memcpy(gbuf->geno.data() + moved_index * n_sample, 
                            gbuf->geno.data() + i * n_sample, nByteGeno);
                    if(gbuf->saveMiss)memcpy(gbuf->miss.data() + moved_index * gbuf->nBLMiss, 
                            gbuf->miss.data() + i * gbuf->nBLMiss, gbuf->nBLMiss);
                }
            }else{
                move_dis++;
            }
        }
    }
}

void Geno::getGenoArrayExtract(const vector<uint32_t>& extract_index, GenoBuf *gbuf){
    int cur_num_marker = extract_index.size();
    uint32_t saveStartMarker = gbuf->indexStartMarker;
    if((cur_num_marker + saveStartMarker) > gbuf->n_marker){
        LOGGER.e(0, "the size of the extracted SNPs is larger than the buffer size.");
    }
    vector<uint32_t> raw_marker_index(cur_num_marker);
    std::transform(extract_index.begin(), extract_index.end(), raw_marker_index.begin(), [this](uint32_t pos){
        return this->marker->getRawIndex(pos);
    });

    if(bHasPreAF){
        gbuf->bHasPreAF = true;
        for(int i = 0; i < cur_num_marker; i++){
            gbuf->preAF[saveStartMarker + i] = AFA1[extract_index[i]];
        }
    }else{
        gbuf->bHasPreAF = false;
    }

    (this->*getGenoArrayFuncs[genoFormat])(raw_marker_index, gbuf);
}



//function for BGEN
void Geno::preProcess_bgen(GenoBuf *gbuf){
    //check bgen files
    gbuf->hasInfo = true;
    compressFormats.clear();
    rawCountSamples.clear();
    rawCountSNPs.clear();
    for(int i = 0; i < gFiles.size(); i++){
        MarkerParam curParam = marker->getMarkerParams(i); 
        if(curParam.rawCountSample != pheno->count_raw()){
            LOGGER.e(0, "inconsistent sample sizes between the .bgen file [" + geno_files[i] + "] and the .sample file (specified by --sample).");
        }
        compressFormats.push_back(curParam.compressFormat);
        rawCountSamples.push_back(curParam.rawCountSample);
        rawCountSNPs.push_back(curParam.rawCountSNP);
    }
}


inline void bgen12ExtractVal(uint64_t val, uint8_t bit_prob, uint64_t mask, uint64_t &v1, uint64_t &v2){
    v1 = val & mask;
    v2 = (val >> bit_prob) & mask;
}

// call genotype
void dosageFunc(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss){
    gval = (double)(prob1 * 2 + prob2) / mask;
}


void dosageCallFunc(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss){
    uint32_t dos = prob1 * 2 + prob2;
    if(dos > A1U){
        gval = 2.0;
    }else if(dos < A1L){
        gval = 0.0;
    }else{
        gval = 1.0;
    }
}

void hardCallFunc(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss){
    miss = false;
    if(prob1 >= cutVal){
        gval = 2.0;
    }else if(prob2 >= cutVal){
        gval = 1.0;
    }else if(mask - prob1 - prob2 >= cutVal){
        gval = 0.0;
    }else{
        gval = 0.0;
        miss = true;
    }
};



void Geno::getGenoArray_bgen(const vector<uint32_t>& raw_marker_index, GenoBuf * gbuf){
    uint32_t indexStartMarker = gbuf->indexStartMarker;
    uint32_t n_res = gbuf->n_sample;
    if(gbuf->saveMiss){
        memset(((gbuf->miss).data() + indexStartMarker * gbuf->nBLMiss), 0, gbuf->nBLMiss * raw_marker_index.size() * 8);
    }
    // read bgen files
    //vector<uint8_t *> snp_data;
    int nMarker = raw_marker_index.size();
    vector<uint32_t> compSizes(nMarker), decSizes(nMarker);
    vector<int> compressMethods(nMarker), fileIDs(nMarker);
    vector<uint64_t> readStarts(nMarker);
    //vector<uint32_t> countSamples; 
    vector<uint32_t> bufStarts(nMarker + 1);
    vector<uint32_t> decDataStarts(nMarker + 1);
    uint32_t sumStarts = 0, sumDecDataStart = 0;
    bufStarts[0] = 0;
    decDataStarts[0] = 0;
    for(int i = 0; i < nMarker; i++){
        uint32_t raw_index = raw_marker_index[i];
        int fileID = marker->getMIndex(raw_index);
        fileIDs[i] = fileID;
        //countSamples.push_back(rawCountSamples[fileID]);

        int compressFormat = compressFormats[fileID];
        compressMethods[i] = compressFormat;

        FILE * curBgen = gFiles[fileID];
        uint64_t pos, size;
        marker->getStartPosSize(raw_index, pos, size);
        Marker::extractBgenMarkerInfo(curBgen, pos);
        uint32_t len_comp = read1Byte<uint32_t>(curBgen);
        uint32_t len_decomp;
        if(compressFormat == 0){
            len_decomp = len_comp;
        }else{
            len_comp -= 4;
            len_decomp = read1Byte<uint32_t>(curBgen);
        }
        compSizes[i] = len_comp; 
        sumStarts += len_comp;
        bufStarts[i + 1] = sumStarts;
        decSizes[i] = len_decomp;
        sumDecDataStart += len_decomp + 8;
        decDataStarts[i + 1] = sumDecDataStart;
        readStarts[i] = ftell(curBgen);
    }
    
    // read genotype
    uint8_t * snp_data = new uint8_t[bufStarts[nMarker]];
    for(int i = 0; i < nMarker; i++){
        int fileID = fileIDs[i];
        FILE * curBgen = gFiles[fileID];
        fseek(curBgen, readStarts[i], SEEK_SET);
        readBytes(curBgen, compSizes[i], snp_data + bufStarts[i]);
    }
        
    // processing each SNP.
    uint8_t *dec_data = new uint8_t[decDataStarts[nMarker]]; 
    //strange check!
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nMarker; i++){
        uint32_t cur_raw_index = raw_marker_index[i];
        uint32_t curFileID = fileIDs[i];
        string error_promp = to_string(cur_raw_index) + "th SNP of "
            + to_string(curFileID) + "th BGEN file.";

        uint32_t len_decomp = decSizes[i];
        int curCompressFormat = compressMethods[i];
        uint8_t * cur_snp_data = snp_data + bufStarts[i];
        uint8_t * cur_dec_data = dec_data + decDataStarts[i];
        if(curCompressFormat == 0){
            memcpy(cur_dec_data, cur_snp_data, compSizes[i]);
        }else{
            uint32_t curCompSize = compSizes[i];
            if(curCompressFormat == 1){
                int z_result = uncompress((Bytef*)cur_dec_data, (uLongf*)&len_decomp, (Bytef*)cur_snp_data, curCompSize);
                if(z_result != Z_OK || len_decomp != decSizes[i]){
                    LOGGER.e(0, "decompressing genotype data error in " + error_promp); 
                }
            }else if(curCompressFormat == 2){
                //zstd
                uint64_t const rSize = ZSTD_getFrameContentSize((void*)cur_snp_data, curCompSize);
                switch(rSize){
                case ZSTD_CONTENTSIZE_ERROR:
                    LOGGER.e(0, "not compressed by zstd in " + error_promp);
                    break;
                case ZSTD_CONTENTSIZE_UNKNOWN:
                    LOGGER.e(0, "original size unknown in " + error_promp);
                    break;
                }
                if(rSize != len_decomp){
                    LOGGER.e(0, "the size stated in compressed data is different from " + error_promp);
                }
                size_t const dSize = ZSTD_decompress((void *)cur_dec_data, len_decomp, (void*)cur_snp_data, curCompSize); 

                if(ZSTD_isError(dSize)){
                    LOGGER.e(0, "decompressing genotype error: " + string(ZSTD_getErrorName(dSize)) + " in " + error_promp);
                }
            }else{
                LOGGER.e(0, "unknown compress format in " + error_promp);
            }
        }
    }
    delete[] snp_data;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nMarker; i++){
        uint32_t cur_raw_index = raw_marker_index[i];
        uint32_t curFileID = fileIDs[i];
        string error_promp = to_string(cur_raw_index) + "th SNP of "
            + to_string(curFileID) + "th BGEN file.";

        //decompress data;
        uint32_t len_decomp = decSizes[i];
        //bool bFreeData = false;
        uint8_t * cur_dec_data = dec_data + decDataStarts[i];

        uint32_t n_sample = *(uint32_t *)cur_dec_data;
        if(n_sample != rawCountSamples[curFileID]){
            LOGGER.e(0, "inconsistent number of samples in " + error_promp );
        }
        uint16_t num_alleles = *(uint16_t *)(cur_dec_data + 4);
        if(num_alleles != 2){
            LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
        }

        uint8_t min_ploidy = *(uint8_t *)(cur_dec_data + 6);//2
        uint8_t max_ploidy = *(uint8_t *)(cur_dec_data + 7); //2
        uint8_t * sample_ploidy = (uint8_t *)(cur_dec_data + 8);
        //check all ploidy are 2 or not
        if(min_ploidy != 2){
            LOGGER.e(0, "multi-ploidy detected in " + error_promp);
        }

        uint8_t *geno_prob = sample_ploidy + n_sample;
        uint8_t is_phased = *(geno_prob);
        uint8_t bits_prob = *(geno_prob+1);
        uint8_t* X_prob = geno_prob + 2;
        uint32_t len_prob = len_decomp - n_sample - 10;
        if(is_phased){
            LOGGER.e(0, "the current version of GCTA does not support phased data in " + error_promp);
        }

        uint8_t double_bits_prob = bits_prob * 2;
        vector<uint32_t> miss_index;

        uint64_t mask = (1U << bits_prob) - 1;
        uint64_t dosage_sum = 0, fij_sum = 0, dosage2_sum = 0;
        uint32_t validN = 0;
        vector<uint32_t> dosages(n_res);
        uint32_t max_dos = mask * 2 + 1;
        bool has_miss = false;
        //LOGGER << "n_res: " << n_res << ", len_prob " << len_prob << ", bits " << (int)bits_prob << std::endl;
        for(int j = 0; j < n_res; j++){
            uint32_t sindex = sampleKeepIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy > 128){
                miss_index.push_back(sindex);
                has_miss = true;
                dosages[j] = max_dos;
            }else if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                //LOGGER << "j:" << sindex << ", bits: " << start_bits << ", byte: " << start_bits / CHAR_BIT << std::endl;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(uint64_t));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;
                uint32_t prob1d = prob1 * 2;
                uint64_t dosage = prob1d + prob2;
                dosages[j] = dosage;
                dosage_sum += dosage;
                dosage2_sum += dosage * dosage;

                uint64_t fij = dosage + prob1d;
                fij_sum += fij;
                validN++;
            }else{
               LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
            }
        }


        double af = (double)dosage_sum / mask / (2 * validN);
        bool bEffRev = this->marker->isEffecRevRaw(cur_raw_index);
        if(bEffRev){
            af = 1.0 - af;
        }
        double std = 2.0 * af * (1.0 - af);
        //double std = (1.0 * dosage2_sum - 1.0 * dosage_sum / validN * dosage_sum) /(validN - 1)/mask/mask;
        //LOGGER << 2.0 * af * (1.0 - af) << "\t" << std << std::endl;
        double info = 0.0;
        if(std < 1e-50){
            info = 1.0;
        }else{
            double dos2_fij_sum = (double)fij_sum / mask - (double)dosage2_sum/mask/mask; 
            info = 1.0 - dos2_fij_sum / (std * validN);
        }

        int gbufIndex = indexStartMarker + i;
        gbuf->af[gbufIndex] = af;
        gbuf->info[gbufIndex] = info;
        gbuf->nValidN[gbufIndex] = validN;
        gbuf->nValidAllele[gbufIndex] = validN * 2;
        //LOGGER << to_string(i) + "After AF calculation" << std::endl;

        //use pre Allele
        if(gbuf->bHasPreAF){
            af = gbuf->preAF[gbufIndex];
            std = 2.0 * af * (1.0 - af);
        }

        double maf = std::min(af, 1.0 - af);
        bool keep_snp = false;
        if(maf >= min_maf && maf <= max_maf){
            if(info >= dFilterInfo && ((1.0*validN / n_sample) >= dFilterMiss)){
                keep_snp = true;
            }
        }

        gbuf->usedIndex[gbufIndex] = keep_snp;
        //LOGGER << to_string(i) + "After filter" << std::endl;

        if(keep_snp){
            if(gbuf->saveGeno){
                double *base_geno = (gbuf->geno).data() + gbufIndex * n_res;
                
                double mu = 0.0;
                double rstd = 1.0;
                if(gbuf->center){
                    mu = 2.0 * af;
                }
                if(gbuf->std){
                    rstd = (std < 1e-50) ? 0 : sqrt(1.0/std);
                }
                double* dos_lookup = new double[max_dos + 2];
                for(uint32_t i = 0; i < max_dos; i++){
                    double tdos = (double)i / mask;
                    tdos = bEffRev ? (2.0 - tdos) : tdos;
                    dos_lookup[i] = (tdos - mu) *rstd;
                }
                dos_lookup[max_dos] = ( 2.0 * af - mu) * rstd;

                for(int j = 0; j < n_res; j++){
                    base_geno[j] = dos_lookup[dosages[j]];
                }
                delete[] dos_lookup;
            }
            //LOGGER << to_string(i) + "geno typed" << std::endl;
            if(has_miss){
                if(gbuf->saveMiss){
                    uint64_t* base_miss = (gbuf->miss).data() + gbufIndex * gbuf->nBLMiss; 
                    for(int j = 0; j < miss_index.size(); j++){
                        int cur_index = miss_index[j];
                        base_miss[cur_index / 64] |= (1UL << (cur_index % 64));
                    }
                }
            }
            //LOGGER << to_string(i) + "after miss all done" << std::endl;
        }
     }
     delete[] dec_data;

}
        /* obsoleted functions hard calling
        uint8_t double_bits_prob = bits_prob * 2;
        uint32_t A1U = floor(mask * 1.5);
        uint32_t A1L = ceil(mask * 0.5);
        uint32_t cutVal = ceil(mask * options_d["hard_call_thresh"]);

        void (*callFunc)(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss);
        if(options.find("dosage_call") != options.end()){
            callFunc = dosageCallFunc;
        }else if(options.find("dosage") != options.end()){
            callFunc = dosageFunc;
        }else{
            callFunc = hardCallFunc;
        }

        auto callFuncB = [callFunc, mask, cutVal, A1U, A1L](uint32_t prob1, uint32_t prob2, double &gval, bool&miss){
            callFunc(prob1, prob2, mask, cutVal, A1U, A1L, gval, miss);
        };
    
        //vector<uint32_t> prob1s(n_res), prob2s(n_res);
        vector<uint32_t> miss_index;
        double sum_geno = 0;
        double *base_geno = (gbuf->geno).data() + i * n_res;
        uint64_t *base_miss;
        uint32_t n_valid = 0;
        if(gbuf->saveMiss) base_miss = (gbuf->miss).data() + i * n_bmiss; 

        uint64_t dosage_sum = 0, fij_dosage2_sum = 0;
        bool *lookup_miss = new bool[mask + 1];


        for(int j = 0; j < n_res; j++){
            uint32_t sindex = sampleKeepIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy > 128){
                miss_index.push_back(sindex);
            }else if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(uint64_t));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;
                uint32_t prob1d = prob1 * 2;
                uint32_t dosage = prob1d + prob2;

                

                dosage_sum += dosage;

                uint32_t fij = dosage + prob1d;
                fij_dosage2_sum += (fij - dosage * dosage);
            }else{
               LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
            }
        }

        for(int j = 0; j < n_res; j++){
            uint32_t sindex = sampleKeepIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy > 128){
                miss_index.push_back(sindex);
            }else if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(uint64_t));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;

                //total_prob += ((prob1 << 1) + prob2);
                //prob1s[i] = prob1;
                //prob2s[i] = prob2;
                double gval;
                bool miss;
                callFuncB(prob1, prob2, gval, miss);
                base_geno[j] = gval;
                if(miss){
                    miss_index.push_back(j);
                    if(gbuf->saveMiss) base_miss[j / 64] |= (1UL << (j % 64)); 
                }else{
                    n_valid++;
                }
           }else{
               LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
           }
        }
        delete[] dec_data;
        Eigen::VectorXd Vgeno = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> > (base_geno, n_res);
        double af = Vgeno.sum() / (2*n_valid);
        (gbuf->validN)[i] = n_valid;
        (gbuf->af)[i] = af;
        */
        


void Geno::endProcess_bgen(){
    closeGFiles();
    if(asyncMode){
        if(asyncBufn) delete asyncBufn;
    }
}

// functions for BED
void Geno::preProcess_bed(GenoBuf *gbuf){
    gbuf->hasInfo = false;
    uint32_t rawCountSample = pheno->count_raw();
    numBytePerMarker = (rawCountSample + 3) / 4; 
    compressFormats.clear();
    rawCountSamples.clear();
    rawCountSNPs.clear();
    //checking the bed files
    bool has_error = false;
    uint64_t f_size;
    uint8_t buffer[3];
    string message;
 
    for(int i = 0; i < gFiles.size(); i++){
        MarkerParam curParam = marker->getMarkerParams(i); 
        compressFormats.push_back(curParam.compressFormat);
        rawCountSamples.push_back(rawCountSample);
        rawCountSNPs.push_back(curParam.rawCountSNP);
    }

    bedIndexLookup.clear();
    bedIndexLookup.push_back(0);
    int32_t sumIndex = 0;
    for(int i = 0; i < gFiles.size() - 1; i++){
        sumIndex += rawCountSNPs[i];
        bedIndexLookup.push_back(sumIndex);
    }

    //check all files
    for(int i = 0; i < gFiles.size(); i++){
        string bed_file = geno_files[i];
        uint32_t cur_size = rawCountSNPs[i];

        FILE *pFile = gFiles[i];
        if(pFile == NULL){
            has_error = true;
            message += "Can't open [" + bed_file + "] to read.\n";
            continue;
        }
        fseek(pFile, 0, SEEK_END);
        f_size = ftell(pFile);
        rewind(pFile);

        if((f_size - 3) != ((uint64_t)numBytePerMarker) * cur_size){
            has_error = true;
            message += "Invalid bed file [" + bed_file +
                "]. The sample and SNP number in bed file are different from bim and fam file.\n";
            continue;
        }

        size_t read_count = fread(buffer, 1, 3, pFile);

        if((read_count != 3) &&
                (*buffer != 0x6c) &&
                (*(buffer+1) != 0x1b) &&
                (*(buffer+2) != 0x01)){
            has_error = true;
            message += "Invalid bed file [" + bed_file +
                "], please convert it into new format (SNP major).\n";
        }
    }

    //delete[] buffer;
    if(has_error){
        LOGGER.e(0, message);
    }

    bedGenoBuf1Size = (sampleKeepIndex.size() + 31) /32;
    keepMask64 = new uint64_t[(rawCountSamples[0] + 63)/64](); 
    pheno->getMaskBit(keepMask64); 
    //if(gbuf->bSex){
    maleMask64 = new uint64_t[(sampleKeepIndex.size() + 63) / 64](); 
    pheno->getMaskBitMale(maleMask64); 
    countMale = pheno->count_male();

    if(asyncMode){
        int n_marker = gbuf->n_marker;
        asyncBufn = new AsyncBuffer<uint8_t>(numBytePerMarker * n_marker);
   }
    //}
    //no improve
    //bed_buf8 = new uint8_t[gbuf->n_marker * numBytePerMarker];
    //bed_buf64 = new uint64_t[gbuf->n_marker * bedGenoBuf1Size];
}

void Geno::read_bed2(const vector<uint32_t> &raw_marker_index, int numMarkerBlock){
    uint8_t *g_buf = NULL;
    //read genotype
    uint32_t numMarker = raw_marker_index.size();
    bool bNewWrite = true;
    int numMarkerRead = 0;
    for(int i = 0; i < numMarker; i++){
        if(bNewWrite){
            g_buf = asyncBufn->start_write();
            bNewWrite = false;
        }
        uint32_t curRawIndex = raw_marker_index[i];
        int curFileID = marker->getMIndex(curRawIndex);
        FILE *pFile = gFiles[curFileID];

        int64_t lag_index = curRawIndex - bedIndexLookup[curFileID];
        if(lag_index < 0)LOGGER.e(0, "strange index in " + to_string(curRawIndex) 
                + "th SNP of " + to_string(curFileID)  + "th BED file.");

        fseek(pFile, (uint64_t) lag_index * numBytePerMarker + 3, SEEK_SET);
        if(fread(g_buf, 1, numBytePerMarker, pFile) != numBytePerMarker){
            perror("Errors:");
            LOGGER << "Error index: " << i << ", raw index: " << curRawIndex << std::endl;
            LOGGER << "Error buffer:" << static_cast<void *>(g_buf) << std::endl;
            LOGGER.e(0, "error in reading [" + geno_files[curFileID] + "].\nThere might be some problems with your storage, or the file has been changed.");
        }
        g_buf += numBytePerMarker;
        numMarkerRead += 1;

        if(numMarkerRead == numMarkerBlock){
            asyncBufn->end_write();
            bNewWrite = true;
            numMarkerRead = 0;
        }
    }
    if(!bNewWrite){
        asyncBufn->end_write();
    }
}

void Geno::bed64ToDouble(uint64_t * geno_buf, GenoBuf * gbuf, const vector<uint32_t> &raw_marker_index){
    //freq
    int indexStartMarker = gbuf->indexStartMarker;
    uint32_t numMarker = raw_marker_index.size();
    if((!gbuf->bHasPreAF) || dFilterMiss != 0){
        if(!gbuf->bSex){
            freq64_bed(geno_buf, raw_marker_index, gbuf);
        }else{
            freq64_bedX(geno_buf, raw_marker_index, gbuf);
        }
    }

    vector<double> &AF = gbuf->bHasPreAF ? (gbuf->preAF) : (gbuf->af);
    uint32_t n_actual_sample = (gbuf->n_sub_sample == 0) ? gbuf->n_sample : gbuf->n_sub_sample;
    uint32_t num_item_n_sample = (n_actual_sample + 31) / 32;
    uint32_t last_sample = (n_actual_sample % 32 == 0) ? 32 : (n_actual_sample % 32);
    uint32_t last_8block = last_sample / 4;
    uint32_t last_2block = last_sample % 4; 

    static std::array<uint8_t, 65536> miss8_lookup16 = [](){
        uint8_t g1_lookup[4];
        g1_lookup[0] = 0;
        g1_lookup[1] = 1;
        g1_lookup[2] = 0;
        g1_lookup[3] = 0;
        uint16_t mask = 3;
        std::array<uint8_t, 65536> miss8_lookup16;
        for(uint32_t b = 0; b < 65536; b++){
            miss8_lookup16[b] = 0;
            for(int j = 0; j < 8; j++){
                miss8_lookup16[b] |= g1_lookup[(b >> (j * 2)) & mask] << j;
            }
        }
        return miss8_lookup16;
    }();


    //code each SNP
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numMarker; i++){
        int gbufIndex = indexStartMarker + i;
        double af = AF[gbufIndex];
        double maf = std::min(af, 1.0 - af);
        bool keep_snp = false;
        if(maf >= min_maf && maf <= max_maf){
            if( (0.5 * gbuf->nValidAllele[gbufIndex] / gbuf->n_sample) >= dFilterMiss){
                keep_snp = true;
            }
        }

        gbuf->usedIndex[gbufIndex] = keep_snp;
        if(keep_snp){
            uint64_t *cur_buf = geno_buf + i * bedGenoBuf1Size;
            bool bEffRev = this->marker->isEffecRevRaw(raw_marker_index[i]);
            if(gbuf->saveGeno){
                double mu = 2.0 * af;
                double center_value = 0.0;
                if(gbuf->center){
                    center_value = mu;
                }
                double rdev = 1.0;
                if(gbuf->std){
                    double dev = mu * (1.0 - af);
                    rdev = (dev < 1.0e-50) ? 0 : sqrt(1.0 / dev);
                }

                double g1_lookup[4];
                g1_lookup[0] = ((bEffRev ? 0.0 : 2.0) - center_value) * rdev;
                g1_lookup[1] = (mu - center_value) * rdev;
                g1_lookup[2] = (1.0 - center_value) * rdev;
                g1_lookup[3] = ((bEffRev ? 2.0 : 0.0)- center_value) * rdev;

                double g_lookup[256][4];
                for(uint16_t i = 0; i <= 255; i++){
                    for(uint16_t j = 0; j < 4; j++){
                        g_lookup[i][j] = g1_lookup[(i >> (2 * j)) & 3];
                    }
                }

                int sub_index = 0;
                uint32_t index = 0;
                double *w_buf = gbuf->geno.data() + (gbufIndex) * gbuf->n_sample;

                for(; index < num_item_n_sample - 1; index++){
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

            if(gbuf->saveMiss){
                uint16_t *g16buf = (uint16_t*)cur_buf;
                uint8_t *curMiss8 = (uint8_t*)(gbuf->miss.data() + gbuf->nBLMiss * gbufIndex);
                int miss8block = (n_actual_sample + 7) / 8;
                int nLast8block = n_actual_sample % 8;
                
                for(int j = 0; j < miss8block; j++){
                    *curMiss8 = miss8_lookup16[*g16buf];
                    curMiss8++;
                    g16buf++;
                }
                uint8_t final8Mask = (nLast8block == 0) ? ((1<<8) - 1) : ((1 << nLast8block) - 1);
                *(--curMiss8) &= final8Mask;
           }
        }
    }
 
}

void Geno::getGenoArray_bed(const vector<uint32_t>& raw_marker_index, GenoBuf* gbuf){
    int numMarker = raw_marker_index.size();
    uint64_t *geno_buf = new uint64_t[bedGenoBuf1Size * numMarker];
    if(asyncMode){
        uint8_t *r_buf = NULL;
        bool isEOF = false;
        std::tie(r_buf, isEOF) = asyncBufn->start_read();
        if(isEOF){
            LOGGER.e(0, "the reading process reached to the end of the BED file but couldn’t finish.");
        }
        move_geno(r_buf, keepMask64, rawCountSamples[0], gbuf->n_sample, numMarker, geno_buf); 
        asyncBufn->end_read();
    }else{
        int indexStartMarker = gbuf->indexStartMarker;
        if(gbuf->saveMiss){
            memset(((gbuf->miss).data() + indexStartMarker * gbuf->nBLMiss), 0, gbuf->nBLMiss * raw_marker_index.size() * 8);
        }

        if(numMarker > (gbuf->n_marker)){
            LOGGER.e(0, "requested marker number exceeded the buffer size.");
        }
        //uint8_t *g_buf = bed_buf8;
        uint8_t *g_buf = new uint8_t[numMarker * numBytePerMarker];
        //read genotype
        for(int i = 0; i < numMarker; i++){
            uint32_t curRawIndex = raw_marker_index[i];
            int curFileID = marker->getMIndex(curRawIndex);
            FILE *pFile = gFiles[curFileID];

            int64_t lag_index = curRawIndex - bedIndexLookup[curFileID];
            if(lag_index < 0)LOGGER.e(0, "strange index in " + to_string(curRawIndex) 
                    + "th SNP of " + to_string(curFileID)  + "th BED file.");

            fseek(pFile, (uint64_t) lag_index * numBytePerMarker + 3, SEEK_SET);
            if(fread(g_buf + i * numBytePerMarker, 1, numBytePerMarker, pFile) != numBytePerMarker){
                LOGGER.e(0, "error in reading [" + geno_files[curFileID] + "].\nThere might be some problems with your storage, or the file has been changed.");
            }
        }

        move_geno(g_buf, keepMask64, rawCountSamples[0], gbuf->n_sample, numMarker, geno_buf); 
        delete[] g_buf;
    }

    //avoid the removed samples
    //uint64_t *geno_buf = bed_buf64;

    /*
    FILE *debug = fopen("test_gbuf1.bin", "wb");
    uint64_t bsize = numMarker * numBytePerMarker;
    fwrite((void *)g_buf, 1, bsize, debug);
    fclose(debug);

    FILE *masks = fopen("test_mask1.bin", "wb");
    uint64_t msize = (rawCountSamples[0] + 63) / 64;
    fwrite((void *)keepMask64, 8, msize, masks);
    fclose(masks);
    
    LOGGER << "1.gbuf size: " << bsize << ", mask size: " << msize << ", raw ct: " << rawCountSamples[0] << ", keep ct: " << gbuf->n_sample << ", nMarker: " << numMarker << std::endl;
    */

     bed64ToDouble(geno_buf, gbuf, raw_marker_index);

     delete[] geno_buf;
}

    
void Geno::freq64_bed(uint64_t *buf, const vector<uint32_t> &markerIndex, GenoBuf *gbuf) {
    const static uint64_t MASK = 6148914691236517205UL; 
    int num_marker = markerIndex.size();

    #pragma omp parallel for schedule(dynamic) 
    for(int cur_marker_index = 0; cur_marker_index < num_marker; ++cur_marker_index){
        //uint32_t curA1A1, curA1A2, curA2A2;
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0;
        uint64_t *p_buf = buf + cur_marker_index * bedGenoBuf1Size;
        for(int index = 0; index < bedGenoBuf1Size ; index++){
            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            odd_ct += popcount(g_buf & MASK);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf & g_buf_h);
        }

        //curA1A1 = num_keep_sample + both_ct - even_ct - odd_ct;
        //curA1A2 = even_ct - both_ct;
        //curA2A2 = both_ct;

        //countA1A1[raw_index_marker] = curA1A1;
        //countA1A2[raw_index_marker] = curA1A2;
        //countA2A2[raw_index_marker] = curA2A2;
        uint32_t cur_total_markers = (gbuf->n_sample - (odd_ct - both_ct)) * 2;
        double cur_af = 1.0*(even_ct + both_ct) / cur_total_markers;
        if(!marker->isEffecRevRaw(markerIndex[cur_marker_index])){
            cur_af = 1.0 - cur_af;
        }
        int curGbufIndex = gbuf->indexStartMarker + cur_marker_index;
        gbuf->af[curGbufIndex] = cur_af;
        gbuf->nValidAllele[curGbufIndex] = cur_total_markers;
        gbuf->nValidN[curGbufIndex] = cur_total_markers / 2;
    }
}

void Geno::freq64_bedX(uint64_t *buf, const vector<uint32_t> &markerIndex, GenoBuf * gbuf){
    const static uint64_t MASK = 6148914691236517205UL; 

    int num_marker = markerIndex.size();
    uint32_t *gender_mask = (uint32_t *)maleMask64;
    uint32_t totalMakers = 2 * gbuf->n_sample - countMale;
    
    #pragma omp parallel for schedule(dynamic) 
    for(int cur_marker_index = 0; cur_marker_index < num_marker; ++cur_marker_index){
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0, odd_ct_m = 0, both_ct_m = 0;
        uint64_t *p_buf = buf + cur_marker_index * bedGenoBuf1Size;
        for(int index = 0; index < bedGenoBuf1Size; index++){
            uint64_t mask_gender = *(gender_mask + index);
            mask_gender = ~fill_inter_zero(mask_gender);

            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            uint64_t g_buf_l = g_buf & MASK;
            uint64_t g_buf_b = g_buf & g_buf_h;
            odd_ct += popcount(g_buf_l);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf_b);

            odd_ct_m += popcount(g_buf_l & mask_gender);
            both_ct_m += popcount(g_buf_b & mask_gender);
        }

        uint32_t cur_total_markers = totalMakers - odd_ct_m - odd_ct + both_ct_m + both_ct;

        double cur_af = 1.0 * (even_ct + both_ct_m) / cur_total_markers;
        if(!marker->isEffecRevRaw(markerIndex[cur_marker_index])){
            cur_af = 1.0 - cur_af;
        }

        int curGbufIndex = gbuf->indexStartMarker + cur_marker_index;
        gbuf->af[curGbufIndex] = cur_af;
        gbuf->nValidAllele[curGbufIndex] = cur_total_markers;
        gbuf->nValidN[curGbufIndex] = cur_total_markers / 2;
    }
}


void Geno::endProcess_bed(){
    closeGFiles();
    if(keepMask64)delete[] keepMask64;
    if(maleMask64)delete[] maleMask64;
    if(bed_buf8) delete[] bed_buf8;
    if(bed_buf64) delete[] bed_buf64;
    if(asyncMode){
        if(asyncBufn) delete asyncBufn;
    }
}

void Geno::openGFiles(){
    closeGFiles();
    gFiles.resize(geno_files.size());
    for(int i = 0; i < geno_files.size(); i++){
        string cur_filename = geno_files[i];
        gFiles[i] = fopen(cur_filename.c_str(), "rb");
        if(gFiles[i] == NULL) {
            LOGGER.e(0, "failed to open genotype [" + cur_filename + "], " + string(strerror(errno)));
        }
    }
}

void Geno::closeGFiles(){
    for(int i = 0; i < gFiles.size(); i++){
        if(gFiles[i])fclose(gFiles[i]);
    }
    gFiles.resize(0);
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

/*
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
*/

void Geno::sum_geno_x(uint64_t *buf, int num_marker) {
    static bool inited = false;
    static std::ofstream out;
    if(!inited){
        inited = true;
        out.open((options["out"] + ".sum").c_str());
        if (!out) { LOGGER.e(0, "failed to open the file [" + options["out"] + ".sum" + "] to write"); }
        out << "CHR\tSNP\tPOS\tA1\tA2\tAAm\tABm\tBBm\tMm\tAAf\tABf\tBBf\tMf" << std::endl;
    }

    const static uint64_t MASK = 6148914691236517205UL; 
    if(num_marker_freq >= marker->count_extract()) return;

    int cur_num_marker_read = num_marker;
    uint32_t *gender_mask = (uint32_t *)keep_male_mask;

    vector<string> out_contents;
    out_contents.resize(cur_num_marker_read);
    
    #pragma omp parallel for schedule(dynamic) 
    for(int cur_marker_index = 0; cur_marker_index < cur_num_marker_read; ++cur_marker_index){
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0, even_ct_m = 0, odd_ct_m = 0, both_ct_m = 0;
        uint64_t *p_buf = buf + cur_marker_index * num_item_1geno;

        for(int index = 0; index < num_item_1geno ; index++){
            uint64_t mask_gender = *(gender_mask + index);
            mask_gender = fill_inter_zero(mask_gender);
            mask_gender += (mask_gender << 1);

            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            uint64_t g_buf_l = g_buf & MASK;
            uint64_t g_buf_b = g_buf & g_buf_h;
            odd_ct += popcount(g_buf_l);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf_b);

            even_ct_m += popcount(g_buf_h & mask_gender);
            odd_ct_m += popcount(g_buf_l & mask_gender);
            both_ct_m += popcount(g_buf_b & mask_gender);
        }

        int all_BB = both_ct;
        int all_AB = even_ct - both_ct;
        int all_miss = odd_ct - both_ct;
        int all_AA = num_keep_sample - odd_ct - even_ct + both_ct;

        int m_BB = both_ct_m;
        int m_AB = even_ct_m - both_ct_m;
        int m_miss = odd_ct_m - both_ct_m;
        int m_AA = num_male_keep_sample - odd_ct_m - even_ct_m + both_ct_m;
        
        std::ostringstream os;
        int raw_index_marker = num_marker_freq + cur_marker_index;
        os << marker->get_marker(marker->getRawIndex(raw_index_marker)) << "\t";
        os << m_AA << "\t" << m_AB << "\t" << m_BB << "\t" << m_miss << "\t";
        os << all_AA - m_AA << "\t" << all_AB - m_AB << "\t" << all_BB - m_BB << "\t" << all_miss - m_miss;

        out_contents[cur_marker_index] = os.str();
 
    }

    std::copy(out_contents.begin(), out_contents.end(), std::ostream_iterator<string>(out, "\n"));
 
    num_marker_freq += num_marker;

}


void Geno::freq64_x(uint64_t *buf, int num_marker) {
    const static uint64_t MASK = 6148914691236517205UL; 
    if(num_marker_freq >= marker->count_extract()) return;

    int cur_num_marker_read = num_marker;
    uint32_t *gender_mask = (uint32_t *)keep_male_mask;
    
    #pragma omp parallel for schedule(dynamic) 
    for(int cur_marker_index = 0; cur_marker_index < cur_num_marker_read; ++cur_marker_index){
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0, odd_ct_m = 0, both_ct_m = 0;
        uint64_t *p_buf = buf + cur_marker_index * num_item_1geno;
        for(int index = 0; index < num_item_1geno ; index++){
            uint64_t mask_gender = *(gender_mask + index);
            mask_gender = ~fill_inter_zero(mask_gender);

            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            uint64_t g_buf_l = g_buf & MASK;
            uint64_t g_buf_b = g_buf & g_buf_h;
            odd_ct += popcount(g_buf_l);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf_b);

            odd_ct_m += popcount(g_buf_l & mask_gender);
            both_ct_m += popcount(g_buf_b & mask_gender);
        }

        int raw_index_marker = num_marker_freq + cur_marker_index;
        uint32_t cur_total_markers = total_markers - odd_ct_m - odd_ct + both_ct_m + both_ct;

        double cur_af = 1.0 * (even_ct + both_ct_m) / cur_total_markers;
        if(!marker->isEffecRev(raw_index_marker)){
            cur_af = 1.0 - cur_af;
        }
        AFA1[raw_index_marker] = cur_af;
        //RDev[raw_index_marker] = 2.0 * cur_af * (1.0 - cur_af);
        countMarkers[raw_index_marker] = cur_total_markers;
    }
    num_marker_freq += num_marker;

}

union Geno_prob{
    char byte[4];
    uint32_t value = 0;
};


void Geno::bgen2bed(const vector<uint32_t> &raw_marker_index){
    LOGGER << "Old bgen 2 bed" << std::endl;
    LOGGER.ts("LOOP_BGEN_BED");
    LOGGER.ts("LOOP_BGEN_TOT");
    vector<uint32_t>& index_keep = pheno->get_index_keep();

    /*
    std::ofstream out((options["out"] + "_sample_index.txt").c_str());
    for(auto & item : index_keep){
        out << item << std::endl;
    }
    out.close();
    */
    
    auto buf_size = (num_raw_sample + 31) / 32;
    size_t buf_size_byte = buf_size * 8;

    int num_marker = 1;
    int num_markers = raw_marker_index.size();
    LOGGER << "samples: " << num_raw_sample << ", keep_sample: " << index_keep.size() << std::endl;
    LOGGER << "Markers: " << num_markers << std::endl;

    FILE * h_bgen = fopen(options["bgen_file"].c_str(), "rb");
    /*
    std::ofstream infos("exmaple.txt");
    infos << "index\traw_index\tpos\tLen_comp\tLen_decomp" << std::endl;
    */
    #pragma omp parallel for schedule(static) ordered
    for(uint32_t index = 0; index < num_markers; index++){
        //LOGGER.i(0, to_string(index) + "NUM_thread: " + to_string(omp_get_max_threads()));
        auto raw_index = raw_marker_index[index];
        uint64_t *buf = new uint64_t[buf_size]();
        uint64_t byte_pos, byte_size;
        this->marker->getStartPosSize(raw_index, byte_pos, byte_size);
        uint32_t len_comp, len_decomp;
        char * snp_data;

        #pragma omp ordered
        {
            fseek(h_bgen, byte_pos, SEEK_SET);
            len_comp = read1Byte<uint32_t>(h_bgen) - 4;
            len_decomp = read1Byte<uint32_t>(h_bgen);
            snp_data = new char[len_comp];
            readBytes(h_bgen, len_comp, snp_data);
            //infos << index << "\t" << raw_index << "\t" << byte_pos << "\t" << len_comp << "\t" << len_decomp << std::endl;
        }
        uLongf dec_size = len_decomp;

        char * dec_data =  new char[len_decomp];
        int z_result = uncompress((Bytef*)dec_data, &dec_size, (Bytef*)snp_data, len_comp);
        delete[] snp_data;
        if(z_result == Z_MEM_ERROR || z_result == Z_BUF_ERROR || dec_size != len_decomp){
            LOGGER.e(0, "decompressing genotype data error in " + to_string(raw_index) + "th SNP."); 
        }

        uint32_t n_sample = *(uint32_t *)dec_data;
        if(n_sample != num_raw_sample){
            LOGGER.e(0, "inconsistent number of samples in " + to_string(raw_index) + "th SNP." );
        }
        uint16_t num_alleles = *(uint16_t *)(dec_data + 4);
        if(num_alleles != 2){
            LOGGER.e(0, "multi-allelic SNPs detected likely because the bgen file is malformed.");
        }

        uint8_t min_ploidy = *(uint8_t *)(dec_data + 6);//2
        uint8_t max_ploidy = *(uint8_t *)(dec_data + 7); //2
        uint8_t * sample_ploidy = (uint8_t *)(dec_data + 8);

        uint8_t *geno_prob = sample_ploidy + n_sample;
        uint8_t is_phased = *(geno_prob);
        uint8_t bits_prob = *(geno_prob+1);
        uint8_t* X_prob = geno_prob + 2;
        uint32_t len_prob = len_decomp - n_sample - 10;
        if(is_phased){
            LOGGER.e(0, "GCTA does not support phased data currently.");
        }

        int byte_per_prob = bits_prob / 8;
        int double_byte_per_prob = byte_per_prob * 2;
        if(bits_prob % 8 != 0){
            LOGGER.e(0, "GCTA does not support probability bits other than in byte units.");
        }

        if(len_prob != double_byte_per_prob * n_sample){
            LOGGER.e(0, "malformed data in " + to_string(raw_index) + "th SNP.");
        }
        /*
        infos << index << "_2\t" << raw_index << "\t" << byte_pos << "\t" << len_comp << "\t" << len_prob << std::endl;
        FILE *obgen = fopen((to_string(index) + ".bin").c_str(), "wb");
        fwrite(snp_data, sizeof(char), len_comp, obgen);
        fwrite(X_prob, sizeof(char), len_prob, obgen);
        fclose(obgen);
        */

        uint32_t base_value = (1 << bits_prob) - 1;

        uint8_t *buf_ptr = (uint8_t *)buf;
        if(options.find("dosage_call") == options.end()){
            uint32_t cut_value = ceil(base_value * options_d["hard_call_thresh"]);
            for(uint32_t i = 0; i < num_keep_sample; i++){
                uint32_t item_byte = i >> 2;
                uint32_t move_byte = (i & 3) << 1;

                uint32_t sindex = index_keep[i];
                uint8_t item_ploidy = sample_ploidy[sindex];

                uint8_t geno_value;
                if(item_ploidy > 128){
                    geno_value = 1;
                }else if(item_ploidy == 2){
                    auto base = sindex * double_byte_per_prob;
                    auto base1 = base + byte_per_prob;
                    Geno_prob prob_item;
                    Geno_prob prob_item1;
                    /*
                       memcpy(prob_item.byte, X_prob + base,  byte_per_prob); 
                       memcpy(prob_item1.byte, X_prob + base1, byte_per_prob); 
                       */
                    for(int i = 0 ; i != byte_per_prob; i++){
                        prob_item.byte[i] = X_prob[base + i];
                        prob_item1.byte[i] = X_prob[base1 + i];
                    }

                    uint32_t t1 = prob_item.value;
                    uint32_t t2 = prob_item1.value;
                    uint32_t t3 = base_value - t1 - t2;
                    if(t1 >= cut_value){
                        geno_value = 0;
                    }else if(t2 >= cut_value){
                        geno_value = 2;
                    }else if(t3 >= cut_value){
                        geno_value = 3;
                    }else{
                        geno_value = 1;
                    }
                }else{
                    LOGGER.e(0, "multi-allelic SNPs detected in the " + to_string(raw_index) + "th SNP.");
                }
                buf_ptr[item_byte] += geno_value << move_byte;
            }
        }else{
            uint32_t A1U = floor(base_value * 1.5);
            uint32_t A1L = ceil(base_value * 0.5);


            for(uint32_t i = 0; i < num_keep_sample; i++){
                uint32_t item_byte = i >> 2;
                uint32_t move_byte = (i & 3) << 1;

                uint32_t sindex = index_keep[i];
                uint8_t item_ploidy = sample_ploidy[sindex];

                uint8_t geno_value;
                if(item_ploidy > 128){
                    //missing
                    geno_value = 1;
                }else if(item_ploidy == 2){
                    auto base = sindex * double_byte_per_prob;
                    auto base1 = base + byte_per_prob;
                    Geno_prob prob_item;
                    Geno_prob prob_item1;
                    /*
                       memcpy(prob_item.byte, X_prob + base,  byte_per_prob); 
                       memcpy(prob_item1.byte, X_prob + base1, byte_per_prob); 
                       */
                    for(int i = 0 ; i != byte_per_prob; i++){
                        prob_item.byte[i] = X_prob[base + i];
                        prob_item1.byte[i] = X_prob[base1 + i];
                    }

                    uint32_t t1 = prob_item.value;
                    uint32_t t2 = prob_item1.value;
                    uint32_t dosageA = 2 * t1 + t2;
                    if(dosageA > A1U){
                        geno_value = 0;
                    }else if(dosageA < A1L){
                        geno_value = 3;
                    }else{
                        geno_value = 2;
                    }
                }else{
                    LOGGER.e(0, " multi-allelic SNPs detected in the " + to_string(raw_index) + "th SNP.");
                }
                buf_ptr[item_byte] += geno_value << move_byte;
            }
 
        }
        //LOGGER.i(0, "MIDDLE: " + to_string(index) + "NUM_thread: " + to_string(omp_get_max_threads()));

        #pragma omp ordered
        save_bed(buf, num_marker);
        delete[] buf;
        delete[] dec_data;
        //#pragma omp ordered
        //LOGGER.i(0, "Finished " + to_string(index) + "NUM_thread: " + to_string(omp_get_max_threads()));
        if(index % 10000 == 0){
            float time_p = LOGGER.tp("LOOP_BGEN_BED");
            if(time_p > 300){
                LOGGER.ts("LOOP_BGEN_BED");
                float elapse_time = LOGGER.tp("LOOP_BGEN_TOT");
                float finished_percent = (float) index / num_markers;
                float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                std::ostringstream ss;
                ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 
                
                LOGGER.i(1, ss.str());
            }
        }

    }
    //infos.close();
    closeOut();
    fclose(h_bgen);
}



void Geno::save_bed(uint64_t *buf, int num_marker){
    static string err_string = "can't write to [" + options["out"] + ".bed].";
    static bool inited = false;
    if(!inited){
        hOut = fopen((options["out"] + ".bed").c_str(), "wb");
        if(hOut == NULL){
            LOGGER.e(0, err_string);
        }
        uint8_t * buffer = new uint8_t[3];
        buffer[0] = 0x6c;
        buffer[1] = 0x1b;
        buffer[2] = 0x01;
        if(3 != fwrite(buffer, sizeof(uint8_t), 3, hOut)){
            LOGGER.e(0, err_string);
        }
        inited = true;
        delete[] buffer;
    }

    uint64_t base_buffer = 0;
    for(int i = 0; i < num_marker; i++){
        uint8_t * buffer = (uint8_t *) (buf + base_buffer);
        if(fwrite(buffer, sizeof(uint8_t), num_byte_keep_geno1,hOut) != num_byte_keep_geno1){
            LOGGER.e(0, err_string);
        }
        base_buffer += num_item_1geno;
    }

}

void Geno::closeOut(){
    fclose(hOut);
}

void Geno::resetFreq(){
    num_marker_freq = 0;
    bFreqFiltered = false;
}

void Geno::resetLoop(){
    resetFreq();
    num_finished_markers = 0;
}

void Geno::freq64(uint64_t *buf, int num_marker) {
    //pheno->mask_geno_keep(buf, num_marker);
    const static uint64_t MASK = 6148914691236517205UL; 
    if(bFreqFiltered) return;
    if(isX){
        freq64_x(buf, num_marker);
        return;
    }

    int cur_num_marker_read = num_marker;
    
    #pragma omp parallel for schedule(dynamic) 
    for(int cur_marker_index = 0; cur_marker_index < cur_num_marker_read; ++cur_marker_index){
        //uint32_t curA1A1, curA1A2, curA2A2;
        uint32_t even_ct = 0, odd_ct = 0, both_ct = 0;
        uint64_t *p_buf = buf + cur_marker_index * num_item_1geno;
        for(int index = 0; index < num_item_1geno ; index++){
            uint64_t g_buf = p_buf[index];
            uint64_t g_buf_h = MASK & (g_buf >> 1);
            odd_ct += popcount(g_buf & MASK);
            even_ct += popcount(g_buf_h);
            both_ct += popcount(g_buf & g_buf_h);
        }

        //curA1A1 = num_keep_sample + both_ct - even_ct - odd_ct;
        //curA1A2 = even_ct - both_ct;
        //curA2A2 = both_ct;

        int raw_index_marker = num_marker_freq + cur_marker_index;

        //countA1A1[raw_index_marker] = curA1A1;
        //countA1A2[raw_index_marker] = curA1A2;
        //countA2A2[raw_index_marker] = curA2A2;
        uint32_t cur_total_markers = (total_markers - 2 * (odd_ct - both_ct));
        double cur_af = 1.0*(even_ct + both_ct) / cur_total_markers;
        if(!marker->isEffecRev(raw_index_marker)){
            cur_af = 1.0 - cur_af;
        }
        AFA1[raw_index_marker] = cur_af;
        //RDev[raw_index_marker] = 2.0 * cur_af * (1.0 - cur_af);
        countMarkers[raw_index_marker] = cur_total_markers;
    }
    num_marker_freq += num_marker;
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

    static uint32_t num_byte_keep_geno = (num_keep_sample + 3) / 4;
    static uint32_t num_byte_per_marker = (num_raw_sample + 3) / 4;
    static uint32_t num_qword_per_marker = (num_byte_keep_geno + 7) / 8;

    static int remain_bit = (num_byte_per_marker % 8) * 8;
    static int move_bit = 64 - remain_bit;
    static uint64_t MASK = (UINT64_MAX << move_bit) >> move_bit;

    const int move_markers = 4;

    #pragma omp parallel for schedule(dynamic) 
    for(uint32_t index = 0; index < num_marker; index += move_markers){
        uint64_t *pbuf[move_markers], *gbuf[move_markers];
        int num_marker_move = 0;
        for(uint32_t i = 0; i < move_markers; i++){
            uint32_t startIndex = index + i;
            if(startIndex < num_marker){ 
                pbuf[i] = (uint64_t *) (buf + startIndex * num_byte_per_marker);
                gbuf[i] = geno_buf + startIndex * num_qword_per_marker; 
                num_marker_move++;
            }
        }
        copy_quaterarr_nonempty_subset(pbuf, keep_list, num_raw_sample, num_keep_sample, gbuf, num_marker_move);
    }
}

void Geno::loop_64block(const vector<uint32_t> &raw_marker_index, vector<function<void (uint64_t *buf, int num_marker)>> callbacks, bool showLog) {
    // show log also act as function that 
    if(showLog){
        LOGGER.i(0, "Reading PLINK BED file(s) in SNP-major format...");
        num_finished_markers = 0;
        if(!bFreqFiltered){
            num_marker_freq = 0;
        }
    }
    //LOGGER << "DEBUG num_finished_markers: " << num_finished_markers << std::endl;
    thread read_thread([this, &raw_marker_index](){this->read_bed(raw_marker_index);});
    read_thread.detach();

    uint8_t *r_buf = NULL;
    bool isEOF = false;
    int cur_num_marker_read;

    if(showLog){
        LOGGER.ts("LOOP_GENO_TOT");
        LOGGER.ts("LOOP_GENO_PRE");
    }
    int cur_num_blocks = (raw_marker_index.size() + Constants::NUM_MARKER_READ - 1) / Constants::NUM_MARKER_READ;

    uint64_t *geno_buf = new uint64_t[num_item_geno_buffer];
    for(int cur_block = 0; cur_block < cur_num_blocks; ++cur_block){
        std::tie(r_buf, isEOF) = asyncBuffer->start_read();
        //LOGGER.i(0, "time get buffer: " + to_string(LOGGER.tp("LOOP_GENO_PRE")));

        LOGGER.d(0, "Process block " + std::to_string(cur_block));
        if(isEOF && cur_block != (cur_num_blocks - 1)){
            LOGGER.e(0, "the reading process reached to the end of the BED file but couldn’t finish.");
        }
        //correct the marker read;
        if(cur_block == (cur_num_blocks - 1)){
            cur_num_marker_read = raw_marker_index.size() - Constants::NUM_MARKER_READ * cur_block;
        }else{
            cur_num_marker_read = Constants::NUM_MARKER_READ;
        }

        //FILE * f1 = fopen(("b" + to_string(cur_block) + ".bin").c_str(), "wb");
        //fwrite(r_buf, sizeof(uint8_t), num_byte_buffer, f1);
        //fclose(f1);
        /*
        if(cur_block == 0){
            FILE *debug = fopen("test_gbuf0.bin", "wb");
            uint64_t bsize = num_byte_buffer;
            fwrite((void *)r_buf, 1, bsize, debug);
            fclose(debug);

            FILE *masks = fopen("test_mask0.bin", "wb");
            uint64_t msize = (num_raw_sample + 63)/64;
            fwrite((void *)keep_mask, 8, msize, masks);
            fclose(masks);

            LOGGER << "0. gbuf size: " << bsize << ", mask size: " << msize << ", raw ct: " << num_raw_sample << ", keep ct: " << num_keep_sample << ", nMarker: " << cur_num_marker_read << std::endl;
        }
        */
 
        move_geno(r_buf, keep_mask, num_raw_sample, num_keep_sample, cur_num_marker_read, geno_buf);
        asyncBuffer->end_read();

        //FILE * f2 = fopen(("b" + to_string(cur_block) + ".mbin").c_str(), "wb");
        //fwrite(geno_buf, sizeof(uint64_t), num_item_geno_buffer, f2);
        //fclose(f2);


        for(auto callback : callbacks){
            callback(geno_buf, cur_num_marker_read);
        }


        num_finished_markers += cur_num_marker_read;
        if(showLog && cur_block % 100 == 0){
            float time_p = LOGGER.tp("LOOP_GENO_PRE");
            if(time_p > 300){
                LOGGER.ts("LOOP_GENO_PRE");
                float elapse_time = LOGGER.tp("LOOP_GENO_TOT");
                float finished_percent = (float) cur_block / cur_num_blocks;
                float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                std::ostringstream ss;
                ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 
                
                LOGGER.i(1, ss.str());
            }
        }
    }
    delete[] geno_buf;
    if(showLog){
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << "100% finished in " << LOGGER.tp("LOOP_GENO_TOT") << " sec";
        LOGGER.i(1, ss.str());
    }
}


void Geno::makeMarkerX(uint64_t *buf, int cur_marker, double *w_buf, bool center, bool std, uint32_t num_sample){
    static uint32_t n_actual_sample = (num_sample == 0) ? num_keep_sample : num_sample;
    uint32_t num_item_n_sample = (n_actual_sample + 31) / 32;
    static uint32_t last_sample = (n_actual_sample % 32 == 0) ? 32 : (n_actual_sample % 32);
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
        double dev = mu * (1.0 - af);
        rdev = (dev < 1.0e-50) ? 0 : sqrt(1.0 / dev);
        //use RDev
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
    for(; index < num_item_n_sample - 1; index++){
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

void Geno::setMAF(double val){
    if(val < 0){
        LOGGER.e(0, "MAF can't be negative: " + to_string(val));
    }
    this->min_maf = val * (1.0 - Constants::SMALL_EPSILON);
    deterFilterMAF();
}

double Geno::getMAF(){
    return this->min_maf;
}

double Geno::getFilterInfo(){
    return this->dFilterInfo;
}

double Geno::getFilterMiss(){
    return this->dFilterMiss;
}

void Geno::deterFilterMAF(){
    if(max_maf < min_maf){
        LOGGER.e(0, "the value specified for --max-maf can't be smaller than that for --min-maf");
    }

    if(std::abs(min_maf) <= 1e-10 && std::abs(max_maf - 0.5) <= 1e-10){
        this->bFilterMAF = false;
    }else{
        this->bFilterMAF = true;
    }
}

void Geno::setMaxMAF(double val){
    if(val <= 0 || val > 0.5){
        LOGGER.e(0, "the value specified for --max-maf can't be negative or larger than 0.5");
    }
    this->max_maf = val * (1.0 + Constants::SMALL_EPSILON);
    deterFilterMAF();
}

void Geno::setFilterInfo(double val){
    this->dFilterInfo = val;
}

void Geno::setFilterMiss(double val){
    this->dFilterMiss = val;
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
    addOneFileOption("geno_file", ".bed", "--bfile", options_in);
    addOneFileOption("bgen_file", "", "--bgen", options_in);
    addOneFileOption("pgen_file", ".pgen", "--pfile", options_in);
    addOneFileOption("pgen_file", ".pgen", "--bpfile", options_in);

    addMFileListsOption("m_file", ".bed", "--mbfile", options_in, options);
    addMFileListsOption("mbgen_file", ".bgen", "--mbgen", options_in, options);
    addMFileListsOption("mpgen_file", ".pgen", "--mbpfile", options_in, options);
    addMFileListsOption("mpgen_file", ".pgen", "--mpfile", options_in, options);

    options_d["min_maf"] = 0.0;
    options_d["max_maf"] = 0.5;
    if(options_in.find("--maf") != options_in.end()){
        auto option = options_in["--maf"];
        if(option.size() == 1){
            try{
                options_d["min_maf"] = std::stod(option[0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "invalid value for --maf");
            }
            if(options_d["min_maf"]<0.0 || options_d["max_maf"]>0.5){
                LOGGER.e(0, "value specified for--maf can't be smaller than 0 or larger than 0.5");
            }else if(options_d["min_maf"] == 0.0){
                options_in["--nofilter"] = {};
            }

        }else{
            LOGGER.e(0, "GCTA does not support multiple values for --maf currently");
        }
        options_in.erase("--maf");
    }

     if(options_in.find("--max-maf") != options_in.end()){
        auto option = options_in["--max-maf"];
        if(option.size() == 1){
            try{
                options_d["max_maf"] = std::stod(option[0]);
           }catch(std::invalid_argument&){
                LOGGER.e(0, "invalid value for --maf");
           }
           if(options_d["max_maf"] < 0.0 || options_d["max_maf"] > 0.5){
               LOGGER.e(0, "the value specified for --max-maf can't be smaller than 0 or larger than 0.5");
           }
        }else{
            LOGGER.e(0, " GCTA does not support multiple values for --maf currently ");
        }
        options_in.erase("--max-maf");
    }

    if(options_d["min_maf"] > options_d["max_maf"]){
        LOGGER.e(0, "value specified for --max-maf can't be smaller than that for --min-maf");
    }


    addOneValOption<double>("geno_rate", "--geno", options_in, options_d, 1.0, 0.0, 1.0);
    if(options_in.find("--geno") != options_in.end()){
        if(options_d["geno_rate"] == 0.0){
            options_in["--nofilter"] = {};
        }
    }

    addOneValOption<double>("info_score", "--info", options_in, options_d, 0.0, 0.0, 1.0);
    addOneValOption<double>("dos_dc", "--dc", options_in, options_d, -1.0, -1.0, 1.0);



    options_d["hard_call_thresh"] = 0.9;
    string flag = "--hard-call-thresh";
    if(options_in.find(flag) != options_in.end()){
        auto option = options_in[flag];
        if(options.size() == 1){
            try{
                options_d["hard_call_thresh"] = std::stod(option[0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "invalid value in " + flag);
            }
        }else{
            LOGGER.e(0, " GCTA does not support multiple values for " + flag + " currently.");
        }
        options_in.erase(flag);
    }

    flag = "--dosage-call";
    if(options_in.find(flag) != options_in.end()){
        options["dosage_call"] = "true";
        options_in.erase(flag);
    }

    if(options_in.find("--freq") != options_in.end()){
        processFunctions.push_back("freq");
        if(options_in["--freq"].size() != 0){
            LOGGER.w(0, "--freq should not be followed by any parameter, if you want to calculate the allele frequencies in founders only, "
                    "please specify by --founders option");
        }
        options_in.erase("--freq");

        options["out"] = options_in["--out"][0];

        return_value++;
    }

    if(options_in.find("--freqx") != options_in.end()){
        processFunctions.push_back("freqx");
        if(options_in["--freqx"].size() != 0){
            LOGGER.w(0, "--freq should not be followed by any other parameter, if you want to calculate the allele frequencies in founders only, "
                    "please specify by --founders option");
        }
        options_in.erase("--freqx");

        options["out"] = options_in["--out"][0];

        return_value++;
    }

    if(options_in.find("--make-bed") != options_in.end()){
        if(options.find("bgen_file") == options.end()){
            processFunctions.push_back("make_bed");
        }else{
            processFunctions.push_back("make_bed_bgen");
        } 

        options_in.erase("--make-bed");
        options["out"] = options_in["--out"][0];

        return_value++;
    }

    if(options_in.find("--recodet") != options_in.end()){
        processFunctions.push_back("recodet");
        options["recode_method"] = "nomiss";
        if(options_in["--recodet"].size() > 0){
            string cop = options_in["--recodet"][0];
            vector<string> ops = {"nomiss", "raw", "std"};
            if(std::find(ops.begin(), ops.end(), cop) != ops.end()){
                options["recode_method"] = cop;
            }else{
                LOGGER.e(0, "can't recognize recode method: " + options_in["--recodet"][0]);
            }
        }
        options_in.erase("--recodet");
        options["out"] = options_in["--out"][0];
        return_value++;
    }


    addOneFileOption("update_freq_file", "", "--update-freq", options_in);

    if(options_in.find("--filter-sex") != options_in.end()){
        options["sex"] = "yes";
    }

    if(options_in.find("--sum-geno-x") != options_in.end()){
        processFunctions.push_back("sum_geno_x");
        options["sex"] = "yes";
        std::map<string, vector<string>> t_option;
        t_option["--chrx"] = {};
        t_option["--filter-sex"] = {};
        Pheno::registerOption(t_option);
        Marker::registerOption(t_option);
        options["out"] = options_in["--out"][0];
        return_value++;
    }



    return return_value;
}

void Geno::setGenoBufSize(GenoBuf *gbuf, uint32_t n_marker){
    uint32_t n_sample = gbuf->n_sample;
    gbuf->n_marker = n_marker;

    gbuf->usedIndex.resize(n_marker);
    gbuf->af.resize(n_marker);
    gbuf->nValidN.resize(n_marker);
    gbuf->nValidAllele.resize(n_marker);
    gbuf->preAF.resize(n_marker);

    if(gbuf->hasInfo)gbuf->info.resize(n_marker);

    if(gbuf->saveGeno)gbuf->geno.resize(n_marker * n_sample);
    if(gbuf->saveMiss){
        gbuf->nBLMiss = (n_sample + 63)/64;
        gbuf->miss.resize(n_marker * gbuf->nBLMiss);
    }
}

void Geno::loopDoublePre(const vector<uint32_t> &extractMarkerIndex, GenoBuf *gbuf, vector<function<void (GenoBuf *gbuf, const vector<uint32_t> &markerIndex)>> callbacks){
    bool showLog = true;
    if(showLog){
        LOGGER.ts("LOOP_GENO_TOT");
        LOGGER.ts("LOOP_GENO_PRE");
    }

    int nMarker = gbuf->n_marker;
    thread read_thread([this, nMarker](){this->read_bed2(this->rawMarkerIndexProceed, nMarker);});
    read_thread.detach();

 
    int nTMarker = extractMarkerIndex.size();
    
    uint32_t nFinishedMarker = 0;

    int pre_block = 0;
    while(nFinishedMarker < nTMarker){
        uint32_t endIndex = nFinishedMarker + nMarker;
        endIndex = endIndex > nTMarker ? nTMarker : endIndex;
        vector<uint32_t> curExtractIndex(extractMarkerIndex.begin() + nFinishedMarker, 
                extractMarkerIndex.begin() + endIndex);
        getGenoArrayExtract(curExtractIndex, gbuf);

        for(auto callback : callbacks){
            callback(gbuf, curExtractIndex);
        }
        nFinishedMarker += nMarker;

        // show progress
        if(showLog){
            int cur_block = nFinishedMarker >> 14;
            if(cur_block > pre_block){
                pre_block = cur_block;
                float time_p = LOGGER.tp("LOOP_GENO_PRE");
                if(time_p > 300){
                    LOGGER.ts("LOOP_GENO_PRE");
                    float elapse_time = LOGGER.tp("LOOP_GENO_TOT");
                    float finished_percent = (float) nFinishedMarker / nTMarker;
                    float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                    std::ostringstream ss;
                    ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 

                    LOGGER.i(1, ss.str());
                }
            }
        }
    }
    if(showLog){
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << "100% finished in " << LOGGER.tp("LOOP_GENO_TOT") << " sec";
        LOGGER.i(1, ss.str());
        LOGGER << nFinishedMarker << " SNPs have been processed." << std::endl;
    }
}

void Geno::loopDoubleFull(const vector<uint32_t> &extractMarkerIndex, GenoBuf *gbuf, vector<function<void (GenoBuf *gbuf, const vector<uint32_t> &markerIndex)>> callbacks){
    bool showLog = true;
    if(showLog){
        LOGGER.ts("LOOP_GENO_TOT");
        LOGGER.ts("LOOP_GENO_PRE");
    }
 
    int nMarker = gbuf->n_marker;
    int nTMarker = extractMarkerIndex.size();
    
    uint32_t nFinishedMarker = 0;
    uint32_t num_valid_marker_processed = 0;
    vector<uint32_t> keptIndex;
    int pre_block = 0;
    while(nFinishedMarker < nTMarker){
        uint32_t endIndex = nFinishedMarker + nMarker;
        endIndex = endIndex > nTMarker ? nTMarker : endIndex;
        uint32_t totalReadMarker;
        //vector<uint32_t> curExtractIndex(extractMarkerIndex.begin() + nFinishedMarker, 
         //       extractMarkerIndex.begin() + endIndex);
        getGenoArrayExtractFull(extractMarkerIndex, nFinishedMarker, nMarker, keptIndex, totalReadMarker, gbuf);
        if(totalReadMarker == 0){
            break;
        }

        for(auto callback : callbacks){
            callback(gbuf, keptIndex);
        }
        nFinishedMarker += totalReadMarker;
        num_valid_marker_processed += keptIndex.size();

        // show progress
        if(showLog){
            int cur_block = nFinishedMarker >> 14;
            if(cur_block > pre_block){
                pre_block = cur_block;
                float time_p = LOGGER.tp("LOOP_GENO_PRE");
                if(time_p > 300){
                    LOGGER.ts("LOOP_GENO_PRE");
                    float elapse_time = LOGGER.tp("LOOP_GENO_TOT");
                    float finished_percent = (float) nFinishedMarker / nTMarker;
                    float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                    std::ostringstream ss;
                    ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 

                    LOGGER.i(1, ss.str());
                }
            }
        }
    }
    if(showLog){
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << "100% finished in " << LOGGER.tp("LOOP_GENO_TOT") << " sec";
        LOGGER.i(1, ss.str());
        LOGGER << num_valid_marker_processed << " SNPs have been processed." << std::endl;
    }
}

void Geno::processFreq(){
    string name_out = options["out"] + ".frq";
    int buf_size = 23068672;
    osBuf.resize(buf_size);
    osOut.rdbuf()->pubsetbuf(&osBuf[0], buf_size);
 
    osOut.open(name_out.c_str());
    if (!osOut) { LOGGER.e(0, "cannot open the file [" + name_out + "] to write."); }
    osOut << "CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS";
    if(hasInfo){
        osOut << "\tINFO";
    }
    osOut << "\n";

    LOGGER << "Computing allele frequencies and saving them to [" << name_out << "]..." << std::endl;

    int nMarker = 128;
    vector<uint32_t> extractIndex(marker->count_extract());
    std::iota(extractIndex.begin(), extractIndex.end(), 0);
    
    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    callBacks.push_back(bind(&Geno::freq_func, this, _1, _2));

    numMarkerOutput = 0;
    loopDouble(extractIndex, nMarker, false, false, false, false, callBacks);

    osOut.flush();
    osOut.close();
    LOGGER << "Saved " << numMarkerOutput << " SNPs." << std::endl;
}

void Geno::processRecodet(){
    string name_out = options["out"] + ".xmat";
    int buf_size = 23068672;
    osBuf.resize(buf_size);
    osOut.rdbuf()->pubsetbuf(&osBuf[0], buf_size);
 
    osOut.open(name_out.c_str());
    if (!osOut) { LOGGER.e(0, "cannot open the file [" + name_out + "] to write."); }
    osOut << "CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS";
    if(hasInfo){
        osOut << "\tINFO";
    }

    LOGGER << "Recoding genotypes and saving them to [" << name_out << "]..." << std::endl;
    uint32_t n_sample = pheno->count_keep();
    vector<string> phenoID = pheno->get_id(0, n_sample - 1, "|");

    for(auto & phenItem : phenoID){
        osOut << "\t" << phenItem;
    }
    osOut << "\n";

    int nMarker = 128;
    bool center, std, saveMiss;
    if(options["recode_method"] == "std"){
        center = true;
        std = true;
        saveMiss = false;
    }else if(options["recode_method"] == "nomiss"){
        center = false;
        std = false;
        saveMiss = false;
    }else if(options["recode_method"] == "raw"){
        center = false;
        std = false;
        saveMiss = true;
    }
    bRecodeSaveMiss = saveMiss;

    vector<uint32_t> extractIndex(marker->count_extract());
    std::iota(extractIndex.begin(), extractIndex.end(), 0);
    
    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    callBacks.push_back(bind(&Geno::recode_func, this, _1, _2));

    numMarkerOutput = 0;
    loopDouble(extractIndex, nMarker, true, center, std, saveMiss, callBacks);

    osOut.flush();
    osOut.close();
    LOGGER << "Saved " << numMarkerOutput << " SNPs." << std::endl;

}

void Geno::freq_func(uintptr_t* genobuf, const vector<uint32_t> &markerIndex){
    int num_marker = markerIndex.size();
    vector<uint8_t> isValids(num_marker);
    vector<double> af(num_marker);
    vector<uint32_t> nValidAllele(num_marker);
    vector<double> info(num_marker);
    /*
    vector<string> outs(num_marker);
    int bufsize = keepSampleCT * 10;
    for(int i = 0; i < num_marker; i++){
        outs[i].reserve(bufsize);
    }
    */
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;

        getGenoDouble(genobuf, i, &item);
        isValids[i] = item.valid;
        if(item.valid){
            af[i] = item.af;
            nValidAllele[i] = item.nValidAllele;
            info[i] = item.info;
        }
    }
    //output
    for(int i = 0; i != num_marker; i++){
        if(isValids[i]){
            numMarkerOutput++;
            osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t"  << af[i]
                <<"\t" << nValidAllele[i];
            if(hasInfo)osOut << "\t" << info[i];
            osOut << "\n";
        }
    }

}

void Geno::recode_func(uintptr_t* genobuf, const vector<uint32_t> &markerIndex){
    int num_marker = markerIndex.size();
    //vector<uint8_t> isValids(num_marker);
    /*
    vector<string> outs(num_marker);
    int bufsize = keepSampleCT * 10;
    for(int i = 0; i < num_marker; i++){
        outs[i].reserve(bufsize);
    }
    */
    #pragma omp parallel for ordered schedule(static,1)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;

        getGenoDouble(genobuf, i, &item);

        #pragma omp ordered
        {
            if(item.valid) {
                numMarkerOutput++;
                osOut << marker->getMarkerStrExtract(cur_marker) << "\t"
                    << item.af << "\t" << item.nValidAllele;
                if(hasInfo) osOut << "\t" << item.info;
                for(int j = 0; j < keepSampleCT; j++){
                    osOut << "\t";
                    if(bRecodeSaveMiss && (item.missing[j/64] & (1UL << (j %64)))){
                        osOut << "NA";
                    }else{
                        osOut << item.geno[j];
                    }
                }
                osOut << "\n";
            }
        }
    }
}

void Geno::processMain() {
    //vector<function<void (uint8_t *, int)>> callBacks;
    vector<function<void (uint64_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "freq"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.processFreq();
        }
            /*
            if(!geno.filterMAF() ){
                LOGGER.i(0, "Computing allele frequencies...");
                callBacks.push_back(bind(&Geno::freq64, &geno, _1, _2));
                geno.loop_64block(marker.get_extract_index(), callBacks);
            }
            geno.out_freq(options["out"]);
            */

            /* // Full manner;
            string name_frq = options["out"] + ".frq";
           


            int nMarker = 128;
            int nTMarker = marker.count_extract();
            int nBMarker = (nTMarker + nMarker - 1) / nMarker;
            GenoBuf gbuf;
            gbuf.center = false;
            gbuf.std = false;
            gbuf.saveMiss = false;
            gbuf.saveGeno = false;
            geno.preProcess(&gbuf, nMarker);

            std::ofstream o_freq(name_frq.c_str());
            if (!o_freq) { LOGGER.e(0, "cannot open the file [" + name_frq + "] to write"); }
            o_freq << "CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS";
            if(gbuf.hasInfo){
                o_freq << "\tINFO\n";
            }else{
                o_freq << "\n";
            }


            uint32_t nFinishedMarker = 0;
            vector<uint32_t> fullExtractIndex(nTMarker);
            std::iota(fullExtractIndex.begin(), fullExtractIndex.end(), 0);
            vector<uint32_t> keptIndex;
            int nValidMarker = 0;
            while(nFinishedMarker < nTMarker){
                uint32_t totalReadMarker;
                geno.getGenoArrayExtractFull(fullExtractIndex, nFinishedMarker, nMarker, keptIndex, totalReadMarker, &gbuf);
                //LOGGER << totalReadMarker << " SNPs read." << std::endl;

                for(int i = 0; i != keptIndex.size(); i++){
                    o_freq << marker.getMarkerStrExtract(keptIndex[i]) << "\t"  << to_string(gbuf.af[i])
                        <<"\t" << gbuf.nValidAllele[i];
                    if(gbuf.hasInfo)o_freq << "\t" << gbuf.info[i];
                    o_freq << "\n";
                }
                if(totalReadMarker == 0){
                    break;
                }
                nFinishedMarker += totalReadMarker;
                nValidMarker += keptIndex.size();
            }
             //final 
            geno.endProcess();
            o_freq.close();
            LOGGER.i(0, "Allele frequencies of " + to_string(nValidMarker) + " SNPs have been saved in the file [" + name_frq + "]");

            callBacks.clear();
        }
        */


        /*  extract mode example
            uint32_t nFinishedMarker = 0;
            while(nFinishedMarker < markerIndices.size()){
                bool isX;
                bool chr_ends;
                vector<uint32_t> indices = marker.getNextSizeIndex(nFinishedMarker, nMarker, chr_ends, isX);
                gbuf.bSex = isX;
                if(indices.size() == 0){
                    break;
                }
                geno.getGenoArrayExtract(indices, &gbuf);
                //output marker
                for(int i = 0; i != indices.size(); i++){
                    o_freq << marker.getMarkerStrExtract(indices[i]) << "\t"  << to_string(gbuf.af[i])
                        <<"\t" << gbuf.nValidAllele[i];
                    if(gbuf.hasInfo)o_freq << "\t" << gbuf.info[i];
                    o_freq << "\n";
                }
                nFinishedMarker += indices.size();
            }
            */
/* 
            for(int bindex = 0; bindex < nBMarker; bindex++){
                int baseIndex = bindex * nMarker; 
                int startIndex = baseIndex;
                int endIndex = (bindex + 1) * nMarker;
                endIndex = endIndex > nTMarker ? nTMarker : endIndex;
                vector<uint32_t> curMarkerIndices(markerIndices.begin() + startIndex,
                    markerIndices.begin() + endIndex);
                uint32_t n_marker = curMarkerIndices.size();
                geno.getGenoArray(curMarkerIndices, &gbuf);

                //output marker
                for(int i = 0; i != n_marker; i++){
                    o_freq << marker.get_marker(curMarkerIndices[i]) << "\t"  << to_string(gbuf.af[i])
                        <<"\t" << gbuf.nValidAllele[i];

                    if(gbuf.hasInfo)o_freq << "\t" << gbuf.info[i];

                    o_freq << "\n";

                }
            }
            */
        if(process_function == "recodet"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.processRecodet();
        }

/*
        if(process_function == "freqx"){
            Geno::setSexMode();

            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            if(!geno.filterMAF()){
                LOGGER.i(0, "Computing allele frequencies...");
                callBacks.push_back(bind(&Geno::freq64, &geno, _1, _2));
                geno.loop_64block(marker.get_extract_index(), callBacks);
            }
            geno.out_freq(options["out"]);
            callBacks.clear();
        }
        */

        if(process_function == "make_bed"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.filterMAF();
            string filename = options["out"];
            pheno.save_pheno(filename + ".fam");
            marker.save_marker(filename + ".bim");
            LOGGER.i(0, "Saving genotype to PLINK binary PED format [" + filename + ".bed]...");
            callBacks.push_back(bind(&Geno::save_bed, &geno, _1, _2));
            geno.loop_64block(marker.get_extract_index(), callBacks);
            geno.closeOut();
            LOGGER.i(0, "Genotype has been saved.");
        }

        if(process_function == "make_bed_bgen"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            string filename = options["out"];
            pheno.save_pheno(filename + ".fam");
            marker.save_marker(filename + ".bim");
            LOGGER.i(0, "Converting bgen to PLINK binary PED format [" + filename + ".bed]...");
            geno.bgen2bed(marker.get_extract_index());
            LOGGER.i(0, "Genotype has been saved.");
        }


        if(process_function == "sum_geno_x"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.filterMAF();
            LOGGER.i(0, "Summing up genotypes based on sex"); 
            callBacks.push_back(bind(&Geno::sum_geno_x, &geno, _1, _2));
            geno.loop_64block(marker.get_extract_index(), callBacks);
            LOGGER.i(0, "Summary has been saved.");
        }


    }

}
