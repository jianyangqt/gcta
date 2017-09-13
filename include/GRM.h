/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: generate, read and process GRM.

   Depends on the class of genotype

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GCTA2_GRM_H
#define GCTA2_GRM_H
#include "Geno.h"
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include "constants.hpp"
#include "mem.hpp"

using std::string;
using std::vector;
using std::pair;

class GRM {
public:
    GRM(Geno *geno);
    GRM();
    ~GRM() {
        posix_mem_free(grm);
        posix_mem_free(N);
        delete[] lookup_GRM_table;
        delete[] sub_miss;
        posix_mem_free(geno_buf);
        posix_mem_free(mask_buf);
        posix_mem_free(cmask_buf);
    };

    void calculate_GRM(uint8_t *buf, int num_marker);
    void grm_thread(int grm_index_from, int grm_index_to);
    void N_thread(int grm_index_from, int grm_index_to);
    void deduce_GRM();
    vector<uint32_t> divide_parts(uint32_t from, uint32_t to, uint32_t num_parts);
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    void loop_block(vector<function<void (double *buf, int num_block)>> callbacks
                    = vector<function<void (double *buf, int num_block)>>());
    void cut_rel(float thresh, bool no_grm = false);
    void prune_fam(float thresh, bool isSparse = true);

private:
    Geno *geno;
    vector<uint32_t> index_keep;
    uint32_t part;
    uint32_t num_parts;
    uint32_t num_individual;
    uint64_t num_grm;
    pair<uint32_t, uint32_t> part_keep_indices; //the index in pheno keep index;
    vector<pair<int, int>> index_grm_pairs;

    double *grm;
    uint32_t *N;
    uint32_t *sub_miss;

    uint16_t *geno_buf;
    uint16_t *mask_buf;
    uint64_t *cmask_buf;
    int cur_num_block = 0;
    int finished_marker = 0;
    const static int num_marker_block = 5;
    const static int num_lookup_table = 8 * 8 * 8 * 8 * 8;

    const static int num_cmask_block = 60;


    double GRM_table[Constants::NUM_MARKER_READ][8];
    //double lookup_GRM_table[Constants::NUM_MARKER_READ / num_marker_block][num_lookup_table];
    double (*lookup_GRM_table)[num_lookup_table];

    vector<uint16_t> elements = {0, 2, 3, 4, 5, 6, 7};
    uint64_t num_byte_geno;
    uint64_t num_byte_cmask;
    int reg_bit_width;
    int num_grm_handle;
    int num_count_handle;
    int num_block_handle;
    int num_marker_process_block;

    void output_id();

    string o_name;

    static map<string, string> options;
    static map<string, double> options_d;
    static map<string, bool> options_b;
    static vector<string> processFunctions;

    string grm_file;
    AsyncBuffer<double>* asyncBuffer;
    void init_AsyncBuffer();
    int num_byte_buffer;
    uint64_t num_subjects;
    vector<string> grm_ids;
    const int num_byte_GRM_read = 100 * 1024 * 1024;
    vector<uint64_t> byte_part_grms;

    void outBinFile(FILE *sFile, FILE *dFile);

    //Just for testing
#ifndef NDEBUG
    FILE * o_geno0;
    FILE * o_mask0;
#endif

};


#endif //GCTA2_GRM_H
