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
#include <cstdio>
#include "Logger.h"
using std::map;
using std::vector;
using std::string;
using std::to_string;

void addOneFileOption(string key_store, string append_string, string key_name, map<string, vector<string>> options_in, map<string,string>& options);
void addMFileListsOption(string key_store, string append_string, string key_name, map<string, vector<string>> options_in, map<string,string>& options);

// get the header of the text file; true: get success, false: only string start with #
bool getListHeaderNew(std::ifstream *input, bool &ioerr, int &ncol, vector<string> &head, int &nHeader);
bool readTxtList(string fileName, int minFields, vector<string> &head, int &nHeader, map<int, vector<string>> &lists);

bool checkFileReadable(string filename);

uint64_t getFileSize(FILE * file);

template <typename T>
void readBytes(FILE * file, int num_item, T* buffer){
    if(num_item != fread(buffer, sizeof(T), num_item, file)){
        LOGGER.e(0, "on reading bgen file on position " + to_string(ftell(file)));
    }
}

template <typename T>
T read1Byte(FILE * file){
    T buffer;
    if(1 != fread(&buffer, sizeof(T), 1, file)){
        LOGGER.e(0, "on reading bgen file on position " + to_string(ftell(file)));
    }
    return buffer;
}

template<typename T>
void addOneValOption(string key_store, string key_name, map<string, vector<string>> options_in, map<string,T>& options, T defVal, T lowVal, T upVal){
    string curFlag = key_name;
    T val = defVal;
    if(options_in.find(curFlag) != options_in.end()){
        auto option = options_in[curFlag];
        if(option.size() == 1){
            try{
                val = std::stod(option[0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "illegal value in " + curFlag + ".");
            }
            if(val < lowVal || val > upVal){
                LOGGER.e(0, curFlag + " can't be smaller than " + to_string(lowVal) + " or larger than " + to_string(upVal) + ".");
            }
        }else{
            LOGGER.e(0, "multiple value in " + curFlag + " are not supported currently");
        }
    }
    options[key_store] = val;
}



#endif  //gcta2_optionio_h
