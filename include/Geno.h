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
    void loop_block(vector<function<void (uint8_t *buf, int num_marker)>> callbacks
                    = vector<function<void (uint8_t *buf, int num_marker)>>());
    void freq(uint8_t *buf, int num_marker);
    bool check_bed();
    void out_freq(string filename);
    static bool registerOption(map<string, vector<string>>& options_in);
    static void processMain();

private:
    Pheno* pheno;
    Marker* marker;
    string bed_file;
    int64_t num_byte_per_marker;
    int last_byte_NA_sample;
    int64_t num_byte_buffer;
    int num_blocks;
    AsyncBuffer<uint8_t>* asyncBuffer;

    GBitCountTable g_table;
    vector<double> AFA1;
    vector<uint32_t> countA1A1;
    vector<uint32_t> countA1A2;
    vector<uint32_t> countA2A2;

    static map<string, string> options;
    static map<string, double> options_d;
    static vector<string> processFunctions;
    static void addOneFileOption(string key_store, string append_string, string key_name,
                                 map<string, vector<string>> options_in);

    void read_bed();

    void init_AF();
    void init_AsyncBuffer();
    void filter_MAF();

    friend class GRM;
 };


#endif //GCTA2_GENO_H
