/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: holds covar information

   Developed by Zhili Zheng<zhilizheng@outlook.com>, 2017

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef gcta2_covar_h
#define gcta2_covar_h
#include <vector>
#include <string>
#include <map>
#include "Pheno.h"

using std::map;
using std::vector;
using std::string;

class Covar{
public:
    Covar();
    bool getCommonSampleIndex(const vector<string> &sampleIDs, vector<uint32_t> &keep_index, vector<uint32_t> &covar_index);
    bool getCovarXRaw(const vector<string> &sampleIDs, vector<double> &X, vector<uint32_t> &keep_index);
    bool getCovarX(const vector<string> &sampleIDs, vector<double> &X, vector<uint32_t> &keep_index);
    bool setCovarMapping(bool is_rcovar, vector<vector<string>> *map_order = NULL);
    const vector<string>& getSampleID() const;

    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    static void read_covar(string filename, vector<string>& sub_list, vector<vector<double>>* covar = NULL, vector<map<string, int>>* labels = NULL, vector<int>* keep_row_p = NULL);

private:
    static map<string, string> options;
    vector<string> sample_id;
    vector<map<string,int>> labels_covar;
    vector<map<string,int>> labels_rcovar;
    vector<vector<double>> labels_covar_mapping;
    vector<vector<double>> labels_rcovar_mapping;
    vector<vector<double>> covar;
    vector<vector<double>> qcovar;
    vector<vector<double>> rcovar;
};

#endif //gcta2_covar_h

