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

#include "Covar.h"
#include "Logger.h"

map<string, string> Pheno::options;

Covar::Covar(){

}

int Covar::registerOption(map<string, vector<string>>& options_in){
    string curFlag = "--qcovar";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["qcovar"] = options_in[curFlag][0];
        }else{
            LOGGER.e(0, curFlag + "can't deal with more than 1 --qcovar file, you shall concat them together.");
        }
    }

    curFlag = "--covar";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["covar"] = options_in[curFlag][0];
        }else{
            LOGGER.e(0, curFlag + "can't deal with more than 1 --covar file, you shall concat them together.");
        }
    }

    return 0;
}

void Covar::processMain(){
    LOGGER.e(0, "No main function in covariate yet.");
}
