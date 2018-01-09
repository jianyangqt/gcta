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

#include <sstream>
#include "Marker.h"
#include "Logger.h"
#include "constants.hpp"
#include <iterator>
#include <algorithm>
#include <numeric>

using std::to_string;

map<string, string> Marker::options;
map<string, int> Marker::options_i;

Marker::Marker() {
    for(uint8_t index = 0; index <= options_i["last_chr"]; index++){
        chr_maps[to_string(index)] = index;
    }
    int last_chr_autosome = options_i["last_chr_autosome"];
    chr_maps["X"] = last_chr_autosome + 1;
    chr_maps["x"] = last_chr_autosome + 1;
    chr_maps["Y"] = last_chr_autosome + 2;
    chr_maps["y"] = last_chr_autosome + 2;
    chr_maps["XY"] = last_chr_autosome + 3;
    chr_maps["xy"] = last_chr_autosome + 3;
    chr_maps["MT"] = last_chr_autosome + 4;
    chr_maps["mt"] = last_chr_autosome + 4;


    if(options.find("marker_file") != options.end()){
        read_bim(options["marker_file"]);
    }else{
        LOGGER.e(0, "No marker file present");
    }

    if(options.find("extract_file") != options.end()){
        vector<string> extractlist = read_snplist(options["extract_file"]);
        extract_marker(extractlist, true);
    }

    if(options.find("exclude_file") != options.end()){
        vector<string> excludelist = read_snplist(options["exclude_file"]);
        extract_marker(excludelist, false);
    }

    if(options.find("update_ref_allele_file") != options.end()){


    }
    

}

void Marker::read_bim(string bim_file) {
    LOGGER.i(0, "Reading PLINK BIM file from [" + bim_file + "]...");
    std::ifstream bim(bim_file.c_str());
    if(!bim){
        LOGGER.e(0, "can not open the file [" + bim_file + "] to read");
    }

    int line_number = 0;
    int last_length = 0;
    string line;
    uint8_t chr_item;
    while(std::getline(bim, line)){
        line_number++;
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() < Constants::NUM_BIM_COL) {
            LOGGER.e(0, "the bim file [" + bim_file + "], line " + to_string(line_number)
                   + " has elements less than " + to_string(Constants::NUM_BIM_COL));
        }
        if(line_number > 1 && line_elements.size() != last_length){
            LOGGER.w(0, "the bim file [" + bim_file + "], line " + to_string(line_number) +
                    " have different elements, take care");
        }
        // replace strings

        try{
            chr_item = chr_maps.at(line_elements[0]);
        }catch(std::out_of_range&){
            LOGGER.e(0, "Line " + to_string(line_number) + " of [" + bim_file +
                    "] contains illegal chr number, please check");
        }
        chr.push_back(chr_item);
        name.push_back(line_elements[1]);
        try{
            gd.push_back(std::stod(line_elements[2]));
            pd.push_back(std::stoi(line_elements[3]));
        }catch(std::invalid_argument&){
            LOGGER.e(0, "Line " + to_string(line_number) + " of [" + bim_file +
                        "] contains illegal distance value, please check");
        }
        a1.push_back(line_elements[4]);
        a2.push_back(line_elements[5]);
        A.push_back(line_elements[4]);
        last_length = line_elements.size();
        if(chr_item >= options_i["start_chr"] && chr_item <= options_i["end_chr"]){
            index_extract.push_back(line_number - 1);
        }
    }
    num_marker = name.size();
    num_extract = index_extract.size();
    LOGGER.i(0, to_string(num_marker) + " SNPs to be included from BIM file.");
    if(num_marker != num_extract){
        LOGGER.i(0, to_string(num_extract) + " SNPs remained after filtered by valid chromosome");
    }
    bim.close();
    if(num_marker == 0){
        LOGGER.e(0, "0 SNPs remained");
    }
}

uint32_t Marker::count_raw() {
    return num_marker;
}

uint32_t Marker::count_extract(){
    return index_extract.size();
}

void Marker::extract_marker(vector<string> markers, bool isExtract) {
    std::sort(markers.begin(), markers.end());
    vector<uint32_t> index_in;
    uint32_t counter = 0;
    for(auto &item_mark : name){
        if(std::binary_search(markers.begin(), markers.end(), item_mark)){
            index_in.push_back(counter);
        }
        counter++;
    }

    if(isExtract){
        auto common = [&index_in](const int key) ->bool{
            return std::find(index_in.begin(), index_in.end(), key) == index_in.end();
        };
        index_extract.erase(std::remove_if(index_extract.begin(), index_extract.end(), common), index_extract.end());
    }else{
        auto diff = [&index_in](const int key) ->bool{
            return std::find(index_in.begin(), index_in.end(), key) != index_in.end();
        };
        index_extract.erase(std::remove_if(index_extract.begin(), index_extract.end(), diff), index_extract.end());
    }

    vector<uint32_t> whole_index(num_marker);
    std::iota(whole_index.begin(), whole_index.end(), 0);
    auto diff = [this](const int key) ->bool{
        return std::find(this->index_extract.begin(), this->index_extract.end(), key) != this->index_extract.end();
    };
    whole_index.erase(std::remove_if(whole_index.begin(), whole_index.end(), diff), whole_index.end());

    index_exclude = whole_index;

    num_extract = index_extract.size();
    num_exclude = index_exclude.size();

    LOGGER.i(0, string("After ") + (isExtract? "extracting" : "excluding") +  " SNP, " + to_string(num_exclude) + " SNPs removed, " + to_string(num_extract) + " SNPs remained.");

}

void Marker::remove_extracted_index(vector<int> remove_index) {
    auto sorted_remove_index = remove_index;
    std::sort(sorted_remove_index.begin(), sorted_remove_index.end());

    for(auto i = sorted_remove_index.rbegin(); i != sorted_remove_index.rend(); ++i){
        index_exclude.push_back(index_extract[*i]);
        index_extract.erase(index_extract.begin() + (*i));
    }
    std::sort(index_exclude.begin(), index_exclude.end());
    std::sort(index_extract.begin(), index_extract.end());
    num_exclude = index_exclude.size();
    num_extract = index_extract.size();
}

string Marker::get_marker(int extract_index){
    return std::to_string(chr[extract_index]) + "\t" + name[extract_index] + "\t" + 
            std::to_string(pd[extract_index]) + "\t" + a1[extract_index] + "\t" + a2[extract_index];
}

bool Marker::isInExtract(uint32_t index) {
    if(index_extract.size() > index_exclude.size()){
        return !std::binary_search(index_exclude.begin(), index_exclude.end(), index);
    }else{
        return std::binary_search(index_extract.begin(), index_extract.end(), index);
    }
}

uint32_t Marker::getExtractIndex(uint32_t extractedIndex){
    return index_extract[extractedIndex];
}

//TODO support multiple SNP list, currently only take the first SNP list
// If multiple column provided, it will miss the correct SNP name.
vector<string> Marker::read_snplist(string snplist_file) {
    std::ifstream if_snplist(snplist_file.c_str());
    vector<string> snplist;

    string line;
    int line_number = 0;
    int last_length = 0;
    while(std::getline(if_snplist, line)){
        line_number++;
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() < 1){
            LOGGER.e(0, "the SNP list file [" + snplist_file + "], line " + to_string(line_number) +
                        " has elements less than 1");
        }

        if(line_number > 1 && line_elements.size() != last_length){
            LOGGER.w(0, "the SNP list file [" + snplist_file + "], line " + to_string(line_number) +
                        " has different elements");
        }
        snplist.push_back(line_elements[0]);
        last_length = line_elements.size();
    }
    if_snplist.close();
    LOGGER.i(0, "Get " + to_string(line_number) + " SNPs from [" + snplist_file + "]");
    return snplist;
}

void Marker::addOneFileOption(string key_store, string append_string, string key_name,
                                    map<string, vector<string>> options_in) {
    if(options_in.find(key_name) != options_in.end()){
        if(options_in[key_name].size() == 1){
            options[key_store] = options_in[key_name][0] + append_string;
        }else if(options_in[key_name].size() > 1){
            options[key_store] = options_in[key_name][0] + append_string;
            LOGGER.w(0, "Marker: multiple " + key_name + ", use the first one only" );
        }else{
            LOGGER.e(0, "no " + key_name + " parameter found");
        }
        std::ifstream f(options[key_store].c_str());
        if(!f.good()){
            LOGGER.e(0, key_name + " " + options[key_store] + " not found");
        }
        f.close();
    }
}

int Marker::registerOption(map<string, vector<string>>& options_in){
    addOneFileOption("marker_file", ".bim", "--bfile", options_in);
    addOneFileOption("marker_file", "", "--bim", options_in);
    addOneFileOption("extract_file", "", "--extract", options_in);
    addOneFileOption("exclude_file", "", "--exclude", options_in);
    addOneFileOption("update_ref_allele_file", "", "--update-ref-allele", options_in);

    if(options_in.find("--autosome-num") != options_in.end()){
        if(options_in["--autosome-num"].size() == 1){
            try{
                options_i["last_chr_autosome"] = std::stoi(options_in["--autosome-num"][0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "Invalid autosome number: " + options_in["--autosome-num"][0]);
            }
        }else{
            LOGGER.e(0, "Multiple --autosome-num is not supported");
        }
    }else if(options_i.find("last_chr_autosome") == options_i.end()){
        options_i["last_chr_autosome"] = 22;
    }
    // include the x y xy MT
    options_i["last_chr"] = options_i["last_chr_autosome"] + 4;
    if(options_i.find("start_chr") == options_i.end()){
        options_i["start_chr"] = 0;
        options_i["end_chr"] = options_i["last_chr"];
    }

    bool filterChrFlag = false;
    static bool specifiedChrFlag = false;

    if(options_in.find("--autosome") != options_in.end()){
        if(specifiedChrFlag){
            if(options_i["start_chr"] < 1 || options_i["end_chr"] > options_i["last_chr_autosome"]){
                LOGGER.e(0, "Chromosome range has been fixed by --chr flag, however it was not in autosome range");
            }
        }else{
            options_i["start_chr"] = 1;
            options_i["end_chr"] = options_i["last_chr_autosome"];
            filterChrFlag = true;
        }
    }

    if(options_in.find("--autosome-x-y") != options_in.end()){
        if(filterChrFlag){
            LOGGER.w(0, "One of the chromosome filter has been applied, it has been overrided by --autosome-x-y");
        }
        if(specifiedChrFlag){
            if(options_i["start_chr"] < 1 || options_i["end_chr"] > options_i["last_chr_autosome"] + 2){
                LOGGER.e(0, "Chromosome range has been fixed by --chr flag, however it was not in autosome range");
            }
        }else{
            options_i["start_chr"] = 1;
            options_i["end_chr"] = options_i["last_chr_autosome"] + 2;
            filterChrFlag = true;
        }
    }

    if(options_in.find("--chr") != options_in.end()){
        specifiedChrFlag = true;
        if(filterChrFlag){
            LOGGER.w(0, "One of the chromosome filter has been applied, it has been overrided by --chr");
        }
        if(options_in["--chr"].size() == 1){
            try{
                options_i["start_chr"] = std::stoi(options_in["--chr"][0]);
                options_i["end_chr"] = options_i["start_chr"];
            }catch(std::invalid_argument&){
                LOGGER.e(0, "--chr contains no numeric value");
            }

        }else if(options_in["--chr"].size() == 2){
            try{
                options_i["start_chr"] = std::stoi(options_in["--chr"][0]);
                options_i["end_chr"] = std::stoi(options_in["--chr"][1]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "--chr contains no numeric value");
            }
        }else{
            LOGGER.e(0, "multiple --chr is not supported currently");
        }

        if(options_i["start_chr"] < 0 || options_i["end_chr"] > options_i["last_chr"]){
            LOGGER.e(0, "--chr is out of chromosome range");
        }

    }

    return 0;
}


void Marker::processMain(){
    LOGGER.e(0, "Marker has no main process this time");
}
