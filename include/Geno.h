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

using std::function;
using std::string;
using std::map;
using std::vector;
using std::placeholders::_1;
using std::placeholders::_2;
using std::bind;

class Geno {
public:
    Geno(Pheno* pheno, Marker* marker);
    ~Geno();
    void loop_64block(const vector<uint32_t>& raw_marker_index, vector<function<void (uint64_t *buf, int num_marker)>> callbacks = vector<function<void (uint64_t *buf, int num_marker)>>());
    void freq(uint8_t *buf, int num_marker);
    void freq2(uint8_t *buf, int num_marker);
    void freq64(uint64_t *buf, int num_marker);
    bool check_bed();
    void out_freq(string filename);
    void makeMarkerX(uint64_t *buf, int cur_marker, double *w_buf, bool center, bool std);
    static void move_geno(uint8_t *buf, uint64_t *keep_list, uint32_t num_raw_sample, 
            uint32_t num_keep_sample, uint32_t num_marker, uint64_t *geno_buf);
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();

private:
    Pheno* pheno;
    Marker* marker;
    vector<string> bed_files;
    int64_t num_byte_per_marker;
    int last_byte_NA_sample;
    int64_t num_byte_buffer;
    int num_blocks;
    int num_finished_markers = 0;
    int num_marker_freq = 0;
    AsyncBuffer<uint8_t>* asyncBuffer = NULL;

    GBitCountTable g_table;
    vector<double> AFA1;
    vector<uint32_t> countA1A1;
    vector<uint32_t> countA1A2;
    vector<uint32_t> countA2A2;
    vector<double> RDev;

    static map<string, string> options;
    static map<string, double> options_d;
    static vector<string> processFunctions;
    static void addOneFileOption(string key_store, string append_string, string key_name,
                                 map<string, vector<string>> options_in);

    void read_bed(const vector<uint32_t> &raw_marker_index);

    void init_AF(string alleleFileName);
    void init_AsyncBuffer();
    void filter_MAF();
    uint64_t num_item_1geno;
    uint64_t num_raw_sample;
    uint64_t num_keep_sample;
    uint64_t num_item_geno_buffer;
    uint64_t *keep_mask;

    friend class GRM;
    friend class FastFAM;
 };


#endif //GCTA2_GENO_H
