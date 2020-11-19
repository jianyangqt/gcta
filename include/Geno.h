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

#ifndef GCTA2_GENO_H
#define GCTA2_GENO_H
#include <string>
#include <map>
#include <vector>
#include "Pheno.h"
#include "Marker.h"
#include "Logger.h"
#include "AsyncBuffer.hpp"
#include <functional>
#include "tables.h"
#include <unordered_map>

using std::function;
using std::string;
using std::map;
using std::vector;
using std::placeholders::_1;
using std::placeholders::_2;
using std::bind;

typedef struct GenoBuffer{
    bool center; //center or not //IN
    bool std;  //std or not  //IN
    bool saveMiss; // save missing or not //IN
    bool saveGeno; // save geno or no //IN

    uint32_t n_sample; // number of filtered keep samples //preIN
    uint32_t n_marker; // number of largest buffer size; //preIN
    uint32_t nBLMiss; // missing data 1 marker size (number of uint64_t); //preIN

    uint32_t n_sub_sample = 0; // subset samples, 0 mean all; //IN
    uint32_t indexStartMarker = 0; // which index to start.  //IN

    vector<uint8_t> usedIndex; //true the index pass the filter.
    vector<double> af; // allele frequencies
    vector<float> info; // info scores
    vector<uint32_t> nValidN; // non-missing samples in a marker
    vector<uint32_t> nValidAllele; // non-missing SNPs in a marker;
    vector<double> geno;  // genotype double per marker
    vector<uint64_t> miss; // missing bits per marker, 1 means miss.

    bool bSex = false; // sex mode (chrX) //IN
    bool hasInfo = false; // determine has info score or not //preIN
    bool bHasPreAF = false; //use Pre AF? //preIN
    vector<double> preAF; //pre AF //preIN
} GenoBuf;

typedef struct GenoBufItem{
    // in
    uint32_t extractedMarkerIndex;   // for allele lookup

    // out
    bool valid;
  //  int refAllele; //0 as 1, 1 as GCTA
    uint8_t isSexXY;   // 0: no, 1: X, 2: Y
    vector<double> geno;
    vector<uintptr_t> missing;
    double af;
    double mean;
    double sd; //std^2
    double info;
    uint32_t nValidN;
    uint32_t nValidAllele;
} GenoBufItem;


class Geno {
public:
    Geno(Pheno* pheno, Marker* marker);
    ~Geno();
    void loop_64block(const vector<uint32_t>& raw_marker_index, vector<function<void (uint64_t *buf, int num_marker)>> callbacks = vector<function<void (uint64_t *buf, int num_marker)>>(), bool showLog = true);
    void freq(uint8_t *buf, int num_marker);
    void freq2(uint8_t *buf, int num_marker);
    void freq64(uint64_t *buf, int num_marker);
    void save_bed(uint64_t *buf, int num_marker);
    void sum_geno(uint64_t *buf, int num_marker);
    void sum_geno_x(uint64_t *buf, int num_marker);
    void closeOut();
    void init_keep();
    void freq64_x(uint64_t *buf, int num_marker);
    void resetFreq();
    void resetLoop();
    void setMAF(double val);
    void setMaxMAF(double val);
    void setFilterInfo(double val);
    void setFilterMiss(double val);
    double getFilterInfo();
    double getMAF();
    double getFilterMiss();
            
    bool check_bed();
    void out_freq(string filename);
    void makeMarkerX(uint64_t *buf, int cur_marker, double *w_buf, bool center, bool std, uint32_t num_sample=0);
    static void move_geno(uint8_t *buf, uint64_t *keep_list, uint32_t num_raw_sample, 
            uint32_t num_keep_sample, uint32_t num_marker, uint64_t *geno_buf);
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    bool filterMAF();
    static void setSexMode();
    uint32_t getTotalMarker();

// seq read start;
    void preProcess(GenoBuf *gbuf, int numMarkerBuf, vector<uint32_t> *rawMarkerIndexProceed = NULL);
    void getGenoArray(const vector<uint32_t>& raw_marker_index, GenoBuf *gbuf);
    void getGenoArrayExtract(const vector<uint32_t>& extract_index, GenoBuf *gbuf);
    void getGenoArrayExtractFull(const vector<uint32_t> & fullExtractIndex, uint32_t start, uint32_t num, vector<uint32_t> &keptIndex, uint32_t &totalReadMarker, GenoBuf *gbuf);
    void endProcess();

    void loopDoubleFull(const vector<uint32_t> &extractMarkerIndex, GenoBuf *gbuf, vector<function<void (GenoBuf *gbuf, const vector<uint32_t> &markerIndex)>> callbacks = vector<function<void (GenoBuf *gbuf, const vector<uint32_t> &markerIndex)>>());
    void loopDoublePre(const vector<uint32_t> &extractMarkerIndex, GenoBuf *gbuf, vector<function<void (GenoBuf *gbuf, const vector<uint32_t> &markerIndex)>> callbacks = vector<function<void (GenoBuf *gbuf, const vector<uint32_t> &markerIndex)>>());

    // new loop subset manner
    void preGenoDouble(int numMarkerBuf, bool bMakeGeno, bool bGenoCenter, bool bGenoStd, bool bMakeMiss);
    void getGenoDouble(uintptr_t *buf, int bufIndex, GenoBufItem* gbuf);
    void endGenoDouble();

    void loopDouble(const vector<uint32_t> &extractIndex, int numMarkerBuf, bool bMakeGeno, bool bGenoCenter, bool bGenoStd, bool bMakeMiss, vector<function<void (uintptr_t *buf, const vector<uint32_t> &exIndex)>> callbacks = vector<function<void (uintptr_t *buf, const vector<uint32_t> &exIndex)>>());

    bool getGenoHasInfo();

    void setGRMMode(bool grm, bool dominace);
    void setGenoItemSize(uint32_t &genoSize, uint32_t &missSize);
 
private:
    Pheno* pheno;
    Marker* marker;
    FILE* hOut = NULL;
    int64_t num_byte_per_marker;
    int last_byte_NA_sample;
    int64_t num_byte_buffer;
    int num_blocks;
    int num_finished_markers = 0;
    int num_marker_freq = 0;
    bool bFreqFiltered = false;
    AsyncBuffer<uint8_t>* asyncBuffer = NULL;

    GBitCountTable g_table;
    vector<double> AFA1;
    
    //vector<uint32_t> countA1A1;
    //vector<uint32_t> countA1A2;
    //vector<uint32_t> countA2A2;
    vector<uint32_t> countMarkers;
    //vector<double> RDev;
    uint32_t total_markers;

    
    static map<string, string> options;
    static map<string, double> options_d;
    static vector<string> processFunctions;
    static void addOneFileOption(string key_store, string append_string, string key_name,
                                 map<string, vector<string>> options_in);

    void read_bed(const vector<uint32_t> &raw_marker_index);

    void init_AF(string alleleFileName);
    void init_AsyncBuffer();
    uint64_t num_item_1geno;
    uint64_t num_raw_sample;
    uint64_t num_keep_sample;
    uint64_t num_byte_keep_geno1;
    uint64_t num_male_keep_sample;
    uint64_t num_item_geno_buffer;
    uint64_t *keep_mask = NULL;
    uint64_t *keep_male_mask = NULL;
    bool isX;
    void bgen2bed(const vector<uint32_t> &raw_marker_index);

    friend class LD;

// seq read start;
    string genoFormat = "";
    vector<string> geno_files;
    vector<FILE *> gFiles;
    vector<int> compressFormats;
    vector<uint32_t> rawCountSamples;
    vector<uint32_t> rawCountSNPs;
    vector<uint32_t> sampleKeepIndex;

    bool bHasPreAF = false;
    int numMarkerBlock;
    vector<uint32_t> rawMarkerIndexProceed;

    double min_maf = 0.0;
    double max_maf = 0.5;
    double dFilterInfo = 0;
    bool bFilterMAF = false;
    double dFilterMiss = 0; // 1 - missingrate
    void deterFilterMAF();

    // define map to calls
    typedef void (Geno::*PreProcessFunc)(GenoBuf *gbuf);
    typedef std::unordered_map<string, PreProcessFunc> PreProcessFuncs;

    typedef void (Geno::*GetGenoArrayFunc)(const vector<uint32_t>& raw_marker_index, GenoBuf * gbuf);
    typedef std::unordered_map<string, GetGenoArrayFunc> GetGenoArrayFuncs;

    typedef void (Geno::*EndProcessFunc)(void);
    typedef std::unordered_map<string, EndProcessFunc> EndProcessFuncs;


    AsyncBuffer<uint8_t>* asyncBufn = NULL;
    bool asyncMode = false;

    PreProcessFuncs preProcessFuncs;
    GetGenoArrayFuncs getGenoArrayFuncs;
    EndProcessFuncs endProcessFuncs;

    void openGFiles();
    void closeGFiles();

    void setGenoBufSize(GenoBuf *gbuf, uint32_t n_marker);

    //bgen format
    void preProcess_bgen(GenoBuf *gbuf);
    void getGenoArray_bgen(const vector<uint32_t>& raw_marker_index, GenoBuf * gbuf);
    void endProcess_bgen();

    //BED format
    void preProcess_bed(GenoBuf *gbuf);
    void getGenoArray_bed(const vector<uint32_t>& raw_marker_index, GenoBuf * gbuf);
    void endProcess_bed();
    int64_t numBytePerMarker; // (raw_sample + 3) / 4
    vector<int32_t> bedIndexLookup; // in new define
    uint64_t bedGenoBuf1Size; // how many 64bit geno of keep sample save in 64bit// bedGenoBufSize = (sampleKeepIndex.size() + 31) /32;
    uint64_t *keepMask64 = NULL;
    uint64_t *maleMask64 = NULL; 
    uint32_t countMale = 0; 
    
    void freq64_bed(uint64_t *buf, const vector<uint32_t> &markerIndex, GenoBuf * gbuf);
    void freq64_bedX(uint64_t *buf, const vector<uint32_t> &markerIndex, GenoBuf * gbuf);
    uint8_t *bed_buf8 = NULL;
    uint64_t *bed_buf64 = NULL;

    //pre read bed
    void read_bed2(const vector<uint32_t> &raw_marker_index, int numMakersBlock);
    void bed64ToDouble(uint64_t * g64buf, GenoBuf * gbuf, const vector<uint32_t> &raw_marker_index);

    //TODO:  X chromosome not work
    
    /**** for subset reading ****
     * ****/

    typedef void (Geno::*PreGenoDoubleFunc)();
    typedef std::unordered_map<string, PreGenoDoubleFunc> PreGenoDoubleFuncs;

    typedef void (Geno::*GetGenoDoubleFunc)(uintptr_t *buf, int idx, GenoBufItem* gbuf);
    typedef std::unordered_map<string, GetGenoDoubleFunc> GetGenoDoubleFuncs;

    typedef void (Geno::*EndGenoDoubleFunc)(void);
    typedef std::unordered_map<string, EndGenoDoubleFunc> EndGenoDoubleFuncs;

    typedef void (Geno::*ReadGenoFunc)(const vector<uint32_t> &extractIndex);
    typedef std::unordered_map<string, ReadGenoFunc> ReadGenoFuncs;

    PreGenoDoubleFuncs preGenoDoubleFuncs;
    GetGenoDoubleFuncs getGenoDoubleFuncs;
    EndGenoDoubleFuncs endGenoDoubleFuncs;
    ReadGenoFuncs readGenoFuncs;
    
    void readGeno(const vector<uint32_t> &extractIndex);

    bool hasInfo = false;
    AsyncBuffer<uintptr_t>* asyncBuf64 = NULL;

    bool bMakeGeno;
    bool bGenoCenter;
    bool bGenoStd;
    bool bMakeMiss;
    bool bGRM = false;
    bool bGRMDom = false;
    int iGRMdc = -1; // 0 no male dosage comp; 1 full comp; //default value shall be -1, equal variance
    int iDC = 1;
    bool f_std = false;
    void setMaleWeight(double &weight, bool &needWeight); // set the male weight by bGRM, dc specity

    int8_t alleModel = 1; // 1: add; 2: Dom; 3: Reces; 4: Het; //currently unused affect a0 a1 a2 na;

    int curBufferIndex;
    vector<int> numMarkersReadBlocks;
    vector<uint8_t> isMarkersSexXYs;
    vector<int> fileIndexBuf;
    vector<int32_t> baseIndexLookup;

    //BED format;
    void preGenoDouble_bed();
    void getGenoDouble_bed(uintptr_t *buf, int idx, GenoBufItem* gbuf);
    void endGenoDouble_bed();
    void readGeno_bed(const vector<uint32_t> &extractIndex);
    //BGEN format;
    void preGenoDouble_bgen();
    void getGenoDouble_bgen(uintptr_t *buf, int idx, GenoBufItem* gbuf);
    void endGenoDouble_bgen();
    void readGeno_bgen(const vector<uint32_t> &extractIndex);
    //PGEN format;
    void preGenoDouble_pgen();
    void getGenoDouble_pgen(uintptr_t *buf, int idx, GenoBufItem* gbuf);
    void endGenoDouble_pgen();
    void readGeno_pgen(const vector<uint32_t> &extractIndex);
 
    //BED
    int bedRawGenoBuf1PtrSize; // how many 64bit geno of raw sample save 
    int maskPtrSize;
    int missPtrSize;
    uint32_t rawSampleCT;
    uint32_t keepSampleCT;

    uintptr_t *keepMaskPtr = NULL;
    uintptr_t *keepMaskInterPtr = NULL; 
    
    uint32_t keepSexSampleCT;
    uint32_t keepMaleSampleCT;

    vector<uint32_t> keepSexIndex;  // in the raw index of fam
    uintptr_t *sexMaskPtr = NULL;
    uintptr_t *sexMaskInterPtr = NULL;

    vector<uint32_t> keepMaleIndex; // in raw index of fam
    vector<uint32_t> keepMaleExtractIndex;
    uintptr_t *maleMaskPtr = NULL;
    uintptr_t *maleMaskInterPtr = NULL; 
    //uintptr_t *maleMaskExtractPtr = NULL;

    //BGEN
    int bgenRawGenoBuf1PtrSize;

    //PGEN
    int pgenGenoBuf1PtrSize;
    int pgenGenoPtrSize;
    int pgenDosagePresentPtrSize;
    int pgenDosageMainPtrSize;

    std::ofstream osOut;
    FILE * bOut = NULL;
    vector<char> osBuf;
    uint32_t numMarkerOutput = 0;

    // main funcs
    void processRecodet();
    bool bRecodeSaveMiss = false;
    void recode_func(uintptr_t* genobuf, const vector<uint32_t> &markerIndex);

    void processFreq();
    void freq_func(uintptr_t * genobuf, const vector<uint32_t> &markerIndex);

 };


#endif //GCTA2_GENO_H
