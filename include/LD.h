
#ifndef GCTA_LD_H
#define GCTA_LD_H
#include "Geno.h"
#include <string>
#include <vector>
#include <map>
#include <memory>

using std::map;
using std::string;
using std::vector;
using std::unique_ptr;

struct LDHeader{
    char magic[3]; // ZLD
    uint8_t version; //0: orignal version
    uint8_t valueType; //0: r; 1: r2;
    uint8_t matrixType; // 0: triangular;  1: full
    uint32_t numMarker; 
    uint8_t compressType; // 0: zfp; 1: SZ
    uint8_t order;  // 0: colMajor; 1: rowMajor;
    uint32_t window;  //number of bp
    uint64_t fileSize; // size of file in bype;
    // reserve 256btye for further extension
    uint64_t markerInfoStart;  // the start of marker = 256; 
    //no marker info if LDinfostart==markerInfoStart
    uint64_t LDInfoStart; // the start of LD infomation
    // For check the data corrupted or not
    uint32_t headerCRC; // the CRC of above bytes
    uint32_t resourceCRC;  // the CRC from the 256byte to end;
};

struct LDBlockInfo{
    uint32_t index;
};
    

struct LDInfoStart{
    uint32_t startIndex;  // marker index in the marker info or .gld.esi
    uint32_t numMarker; 
    // marker Index 50 seems resonalbe 30K * 8 * 1024/1024  * 50 = 114MB
    uint32_t* numRelAfter; // how many marker relation saved; size = 1
    uint32_t* numRelBefore; // how many marker before;
    uint64_t startByte; // the byte of start
    uint64_t compressSize; // compressed size;
    uint64_t decompressSize; // decompressed size;
};

class LD{
public:
    LD(Geno *geno);
    ~LD();
    void readGeno(uint64_t *buf, int num_marker);
    void readLD();
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    
private:
    Geno *geno;
    static uint32_t num_indi;
    static map<string, string> options;
    static map<string, int> options_i;
    static vector<string> processFunctions;

    static bool chr_ends;
    static unique_ptr<double[]> geno_buffer[2];
    static int cur_buffer;
    static uint64_t cur_buffer_offset[2];
    void calcLD();
    uint32_t ld_window;
    bool is_r2;
    uint32_t cur_process_marker_index;

    FILE *h_ld;

};
    
    






#endif //GCTA_LD_H:w


