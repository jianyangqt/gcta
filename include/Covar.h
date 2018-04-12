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

    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();

private:
    static map<string, string> options;
};

#endif //gcta2_covar_h

