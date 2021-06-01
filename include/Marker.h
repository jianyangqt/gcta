/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: SNPs information in plink format.

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GCTA2_MARKER_H
#define GCTA2_MARKER_H
#include <string>
#include <vector>
#include <map>
#include <utility>
using std::vector;
using std::string;
using std::to_string;
using std::pair;
using std::map;

struct MarkerInfo{
    string chr;
    string name;
    uint32_t pd;
    vector<string> alleles;
    bool A_rev;
};

struct MarkerParam{
    uint32_t rawCountSNP;
    uint32_t rawCountSample; //if unknown, 0;
    int compressFormat; // -1, don't know;
    uint64_t posGenoDataStart;
};

class Marker {

public:
    Marker();
    uint32_t count_raw(int part = -1);
    uint32_t count_extract();
    bool isInExtract(uint32_t index);
    uint32_t getRawIndex(uint32_t extractedIndex);
    vector<uint32_t>& get_extract_index(); // return the raw index for autosome;
    vector<uint32_t> get_extract_index_autosome(); // return the extract index for autosome; //not raw index
    vector<uint32_t> get_extract_index_X();
    int getMIndex(uint32_t raw_index);
    //uint64_t getStartPosSize(uint32_t raw_index);
    void getStartPosSize(uint32_t raw_index, uint64_t &pos, uint64_t &size);
    bool isEffecRev(uint32_t extractedIndex);
    bool isEffecRevRaw(uint32_t rawIndex);
    string get_marker(int rawindex, bool bflip=false);
    string getMarkerStrExtract(int extractindex, bool bflip=false);
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    static MarkerInfo extractBgenMarkerInfo(FILE *h_bgen, uint64_t &pos);
    static MarkerParam getBgenMarkerParam(FILE *h_bgen, string &outputs);
    void extract_marker(vector<string> markers, bool isExtract);
    void reset_exclude();
    void keep_raw_index(const vector<uint32_t>& keep_index);
    void keep_extracted_index(const vector<uint32_t>& keep_index);
    void matchSNPListFile(string filename, int num_min_fields, const vector<int>& field_return, vector<string> &fields, vector<bool>& a_rev, bool update_a_rev = false);
    void save_marker(string filename);
    vector<uint32_t> getNextWindowIndex(uint32_t cur_marker_index, uint32_t window, bool& chr_ends, bool& isX, bool retRaw = true);
    vector<uint32_t> getNextSizeIndex(uint32_t cur_marker_index, uint32_t num, bool& chr_ends, bool& isX, bool retRaw = false);
    uint32_t getNextSize(const vector<uint32_t> &rawRef, uint32_t curExtractIndex, uint32_t num, int &fileIndex, bool &chr_ends, uint8_t &isSexXY);
    uint32_t getNextWindowSize(uint32_t cur_marker_index, uint32_t window);

    MarkerParam getMarkerParams(int part_num);
    uint8_t mapCHR(string chr_str, bool &success);

    uint64_t getMaxGenoMarkerUptrSize();
    vector<pair<string, vector<uint32_t>>> read_gene(string gfile);

private:
    vector<uint8_t> chr;
    vector<string> name;
    vector<float> gd;
    vector<uint32_t> pd;
    vector<string> a1;
    vector<string> a2;
    vector<bool> A_rev; //effect allele;
    vector<uint64_t> byte_start;
    vector<uint64_t> byte_size;
    vector<uint32_t> raw_limits;

    //bgen
    uint64_t maxGeno1ByteSize = 0;

    vector<uint32_t> index_extract;
    vector<uint32_t> index_exclude;
    uint32_t num_marker;
    uint32_t num_extract;
    uint32_t num_exclude;

    void read_bim(string bim_file);
    void read_mbim(string bim_file);
    void read_bgen(string bgen_file);
    void read_mbgen(string mbgen_file);

    void read_pvar(string pvar_file);
    void read_mpvar(string mpvar_file);

    static map<string, string> options;
    static map<string, int> options_i;
    static void addOneFileOption(string key_store, string append_string, string key_name,
                                 map<string, vector<string>> options_in);
    vector<string> read_snplist(string snplist_file);
    void read_bgen_index(string bgen_file);
    map<string, uint8_t> chr_maps;
    vector<MarkerParam> markerParams;
};


#endif //GCTA2_MARKER_H
