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
using std::vector;
using std::string;
using std::map;

class Marker {

public:
    Marker();
    uint32_t count_raw();
    uint32_t count_extract();
    bool isInExtract(uint32_t index);
    uint32_t getExtractIndex(uint32_t extractedIndex);
    string get_marker(int extract_index);
    static bool registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    void extract_marker(vector<string> markers, bool isExtract);
    void remove_extracted_index(vector<int> remove_index);

private:
    vector<uint8_t> chr;
    vector<string> name;
    vector<double> gd;
    vector<uint32_t> pd;
    vector<string> a1;
    vector<string> a2;
    vector<string> A; //effect allele;

    vector<uint32_t> index_extract;
    vector<uint32_t> index_exclude;
    uint32_t num_marker;
    uint32_t num_extract;
    uint32_t num_exclude;

    void read_bim(string bim_file);

    static map<string, string> options;
    static map<string, int> options_i;
    static void addOneFileOption(string key_store, string append_string, string key_name,
                                 map<string, vector<string>> options_in);
    vector<string> read_snplist(string snplist_file);
    map<string, uint8_t> chr_maps;
};


#endif //GCTA2_MARKER_H
