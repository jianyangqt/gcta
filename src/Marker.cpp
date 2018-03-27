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
#include <boost/algorithm/string.hpp>
#include "utils.hpp"

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
    reset_exclude();

    if(options.find("update_ref_allele_file") != options.end()){
        LOGGER.i(0, "Reading reference alleles of SNPs from [" + options["update_ref_allele_file"] + "]...");
        vector<int> field_return;
        vector<string> fields;
        vector<bool> a_rev;
        matchSNPListFile(options["update_ref_allele_file"], 2, field_return, fields, a_rev, true);

        LOGGER.i(0, "Reference alleles are updated."); 
    }

}

void Marker::matchSNPListFile(string filename, int num_min_fields, const vector<int>& field_return, vector<string> &fields, vector<bool>& a_rev, bool update_a_rev){
    vector<string> temp_fields;
    vector<bool> temp_a_rev;
    std::ifstream allele_file(filename.c_str());
    if(allele_file.fail()){
        LOGGER.e(0, "can't open [" + filename + "] to read]");
    }

    int max_fields = num_min_fields;
    int require_fields = (field_return.size() > 0) ? (field_return[field_return.size() - 1] + 1): 0;
    require_fields = (max_fields > require_fields) ? max_fields : require_fields;

    vector<string> marker_name, ref_allele;
    string line;
    vector<uint32_t> bad_lines;
    int num_line = 0;
    while(std::getline(allele_file, line)){
        num_line++;
        vector<string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        if(line_elements.size() < require_fields){
            bad_lines.push_back(num_line);
        }
        if(line_elements.size() >= 1){
            marker_name.push_back(line_elements[0]);
        }

        if(line_elements.size() >= 2){
            ref_allele.push_back(line_elements[1]);
        }
        for(auto index : field_return){
            temp_fields.push_back(line_elements[index]);
        }

    }

    // match snp name
    vector<uint32_t> marker_index, ref_index;
    vector_commonIndex_sorted1(name, marker_name, marker_index, ref_index);

    // match alleles
    vector<uint32_t> index_remained;
    index_remained.reserve(marker_index.size());
    if(ref_allele.size()){
        temp_a_rev.reserve(marker_index.size());
        for(int i = 0; i < marker_index.size(); i++){
            uint32_t cur_marker_index = marker_index[i];
            string cur_ref = ref_allele[ref_index[i]];
            std::transform(cur_ref.begin(), cur_ref.end(), cur_ref.begin(), toupper);
            if(a1[cur_marker_index] == cur_ref){
                temp_a_rev.push_back(false);
                index_remained.push_back(cur_marker_index);
            }else if(a2[cur_marker_index] == cur_ref){
                temp_a_rev.push_back(true);
                index_remained.push_back(cur_marker_index);
            }
        }
    }else{
        index_remained = marker_index;
    }

    if(update_a_rev){
        for(int i = 0; i < index_remained.size(); i++){
            A_rev[index_remained[i]] = temp_a_rev[i];
        }
    }

    // match original kept list
    vector<uint32_t> extract_index, export_index;
    vector_commonIndex_sorted1(index_extract, index_remained, extract_index, export_index);

    vector<uint32_t> index_common;
    index_common.resize(extract_index.size());
    std::transform(extract_index.begin(), extract_index.end(), index_common.begin(), [this](size_t pos){return (this->index_extract)[pos];});

    int rm_snps = index_extract.size() - index_common.size();
    if(rm_snps){
        LOGGER.w(0, to_string(rm_snps) + " SNPs are removed due to mismatching SNP name or alleles.");
    }

    keep_raw_index(index_common);

    if(ref_allele.size()){
        a_rev.resize(export_index.size());
        std::transform(export_index.begin(), export_index.end(), a_rev.begin(), [&temp_a_rev](size_t pos){return temp_a_rev[pos];});
    }

    int field_num = field_return.size();
    if(field_num){
        fields.reserve(export_index.size() * field_num);
        for(auto i : export_index){
            for(int j = i; j < i + field_num; j++){
                fields.push_back(temp_fields[j]);
            }
        }
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
        //std::istringstream line_buf(line);
        //std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements;
        //vector<string> line_elements(begin, end);
        boost::split(line_elements, line, boost::is_any_of("\t "));
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
        std::transform(line_elements[4].begin(), line_elements[4].end(), line_elements[4].begin(), toupper);
        std::transform(line_elements[5].begin(), line_elements[5].end(), line_elements[5].begin(), toupper);
        a1.push_back(line_elements[4]);
        a2.push_back(line_elements[5]);
        A_rev.push_back(false);
        last_length = line_elements.size();
        if(chr_item >= options_i["start_chr"] && chr_item <= options_i["end_chr"]){
            index_extract.push_back(line_number - 1);
        }
    }
    num_marker = name.size();
    num_extract = index_extract.size();
    LOGGER.i(0, to_string(num_marker) + " SNPs to be included from BIM file.");
    if(num_marker != num_extract){
        LOGGER.i(0, to_string(num_extract) + " SNPs to be included from valid chromosome number");
    }
    bim.close();
    if(num_marker == 0){
        LOGGER.e(0, "0 SNP remain.");
    }
}

uint32_t Marker::count_raw() {
    return num_marker;
}

uint32_t Marker::count_extract(){
    return index_extract.size();
}

void Marker::extract_marker(vector<string> markers, bool isExtract) {
    vector<uint32_t> ori_index, marker_index;
    vector_commonIndex_sorted1(name, markers, ori_index, marker_index);
    vector<uint32_t> remain_index;
    if(isExtract){
        std::set_intersection(index_extract.begin(), index_extract.end(),
                              ori_index.begin(), ori_index.end(),
                              std::back_inserter(remain_index));
    }else{
        std::set_difference(index_extract.begin(), index_extract.end(),
                            ori_index.begin(), ori_index.end(),
                            std::back_inserter(remain_index));
    }

    index_extract = remain_index;
    num_extract = index_extract.size();

    if(num_extract == 0){
        LOGGER.e(0, "0 SNP remain.");
    }

    LOGGER.i(0, string("After ") + (isExtract? "extracting" : "excluding") +  " SNP, " +  to_string(num_extract) + " SNPs remain.");

}

void Marker::reset_exclude(){
    vector<uint32_t> whole_index(num_marker);
    std::iota(whole_index.begin(), whole_index.end(), 0);

    index_exclude.resize(whole_index.size() - index_extract.size());
    std::set_difference(whole_index.begin(), whole_index.end(), index_extract.begin(), 
            index_extract.end(), index_exclude.begin());
    num_exclude = index_exclude.size();
}

// didn't check whether in extracted list or not;
void Marker::keep_raw_index(const vector<uint32_t>& keep_index) {
    index_extract.resize(keep_index.size());
    index_extract = keep_index;

    num_extract = index_extract.size();
    if(num_extract == 0){
        LOGGER.e(0, "0 SNP remain.");
    }
    reset_exclude();
}

void Marker::keep_extracted_index(const vector<uint32_t>& keep_index) {
    vector<uint32_t> raw_index;
    raw_index.reserve(keep_index.size());
    for(auto index : keep_index){
        raw_index.push_back(index_extract[index]);
    }
    std::sort(raw_index.begin(), raw_index.end());
    keep_raw_index(raw_index);
}

string Marker::get_marker(int extract_index){ // raw index
    string return_string = std::to_string(chr[extract_index]) + "\t" + name[extract_index] + "\t" + 
        std::to_string(pd[extract_index]) + "\t";
    if(A_rev[extract_index]){
        return return_string + a2[extract_index] + "\t" + a1[extract_index];
    }else{
        return return_string + a1[extract_index] + "\t" + a2[extract_index];
    }
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

bool Marker::isEffecRev(uint32_t extractedIndex){
    return A_rev[index_extract[extractedIndex]];
}

//TODO support multiple column SNP list, currently only take the first SNP list
// If multiple column provided, it will miss the correct SNP name.
vector<string> Marker::read_snplist(string snplist_file) {
    std::ifstream if_snplist(snplist_file.c_str());
    vector<string> snplist;

    string line;
    int line_number = 0;
    int last_length = 0;
    while(std::getline(if_snplist, line)){
        line_number++;
        //std::istringstream line_buf(line);
        //std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
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
