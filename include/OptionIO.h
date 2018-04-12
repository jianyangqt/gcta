/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   IO or option helper functions;

   Developed by Zhili Zheng<zhilizheng@outlook.com>, 2017

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef gcta2_optionio_h
#define gcta2_optionio_h
#include <vector>
#include <string>
#include <map>
using std::map;
using std::vector;
using std::string;

void addOneFileOption(string key_store, string append_string, string key_name, map<string, vector<string>> options_in, map<string,string>& options);

bool checkFileReadable(string filename);
#endif  //gcta2_optionio_h
