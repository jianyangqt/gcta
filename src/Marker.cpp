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
#include <boost/algorithm/string/join.hpp>
#include "utils.hpp"
#include "OptionIO.h"
#include <memory>

using std::to_string;
using std::unique_ptr;

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
 
    bool has_marker = false;

    if(options.find("marker_file") != options.end()){
        read_bim(options["marker_file"]);
        has_marker = true;
        raw_limits.push_back(num_marker);
    }
    
    if(options.find("m_file") != options.end()){
        has_marker = true;
        read_mbim(options["m_file"]);
    }

    if(options.find("bgen_file") != options.end()){
        has_marker = true;
        read_bgen(options["bgen_file"]);
        raw_limits.push_back(num_marker);
    }

    if(!has_marker){
        LOGGER.e(0, "no marker exist");
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
        uint16_t num_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        num_elements = line_elements.size();
        boost::replace_all(line_elements[num_elements - 1], "\r", "");
        if(num_elements < require_fields){
            bad_lines.push_back(num_line);
        }
        if(num_elements >= 1){
            marker_name.push_back(line_elements[0]);
        }

        if(num_elements >= 2){
            ref_allele.push_back(line_elements[1]);
        }
        for(auto index : field_return){
            temp_fields.push_back(line_elements[index]);
        }

    }
    allele_file.close();

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
    // bug here?
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

void Marker::read_mbim(string mfile){
    vector<string> mfiles;
    boost::split(mfiles, mfile, boost::is_any_of("\t "));
    raw_limits.clear();
    for(auto & mfile : mfiles){
        read_bim(mfile + ".bim");
        raw_limits.push_back(num_marker);
    }
}

vector<uint32_t>& Marker::get_extract_index(){
    return this->index_extract;
}

//cur_marker_index:  is the index point to position of index_extract
vector<uint32_t> Marker::getNextWindowIndex(uint32_t cur_marker_index, uint32_t window, bool& chr_ends){
    vector<uint32_t> indices;
    if(cur_marker_index >= index_extract.size()){
        chr_ends = true;
        return indices;
    }
    uint32_t cur_index = index_extract[cur_marker_index];

    uint32_t cur_pd = pd[cur_index];
    uint8_t cur_chr = chr[cur_index];
    uint32_t final_pd = window + cur_pd;

    chr_ends = false;

    for(uint32_t marker_index = cur_marker_index; marker_index < index_extract.size(); marker_index++){
        uint32_t temp_index = index_extract[marker_index];
        if(chr[temp_index] != cur_chr){
            chr_ends = true;
            break;
        }
        if(pd[temp_index] <= final_pd){
            indices.push_back(temp_index);
        }else{
            break;
        }
    }
    return indices;
}

uint32_t Marker::getNextWindowSize(uint32_t cur_marker_index, uint32_t window){
    if(cur_marker_index >= index_extract.size()){
        return 0;
    }
    uint32_t cur_index = index_extract[cur_marker_index];

    uint32_t cur_pd = pd[cur_index];
    uint8_t cur_chr = chr[cur_index];
    uint32_t final_pd = window + cur_pd;

    uint32_t count = 0;

    for(uint32_t marker_index = cur_marker_index; marker_index < index_extract.size(); marker_index++){
        uint32_t temp_index = index_extract[marker_index];
        if(chr[temp_index] != cur_chr){
            break;
        }
        if(pd[temp_index] <= final_pd){
            count++;
        }else{
            break;
        }
    }
    return count;
}


int Marker::getMIndex(uint32_t raw_index){
    for(int i = 0; i < raw_limits.size(); i++){
        if(raw_index < raw_limits[i]){
            return i;
        }
    }
    LOGGER.e(0, "too large SNP index " + to_string(raw_index));
    return 0;
}

void Marker::save_marker(string filename){
    LOGGER.i(0, "Saving SNP information to [" + filename + "]...");
    std::ofstream out(filename.c_str());
    for(auto & i : index_extract){
        out << (int)chr[i] << "\t" << name[i] << "\t" << gd[i] 
            << "\t" << pd[i] << "\t" << a1[i] << "\t"
            << a2[i] << std::endl;
    }
    LOGGER.i(0, to_string(index_extract.size()) + " SNPs saved.");

}

void Marker::read_bim(string bim_file) {
    LOGGER.i(0, "Reading PLINK BIM file from [" + bim_file + "]...");
    std::ifstream bim(bim_file.c_str());
    if(!bim){
        LOGGER.e(0, "can not open the file [" + bim_file + "] to read");
    }

    int line_number = name.size();
    int last_length = 0;
    string line;
    int start_line_number = name.size() + 1;
    uint8_t chr_item;
    while(std::getline(bim, line)){
        line_number++;
        //std::istringstream line_buf(line);
        //std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements;
        //vector<string> line_elements(begin, end);
        boost::split(line_elements, line, boost::is_any_of("\t "));
        boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
        if(line_elements.size() < Constants::NUM_BIM_COL) {
            LOGGER.e(0, "the bim file [" + bim_file + "], line " + to_string(line_number)
                   + " has elements less than " + to_string(Constants::NUM_BIM_COL));
        }
        if(line_number > start_line_number && line_elements.size() != last_length){
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
            gd.push_back(std::stof(line_elements[2]));
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
        //TODO refractor the code.
        byte_start.push_back(1);
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
    if(num_extract == 0){
        LOGGER.e(0, "0 SNP remain.");
    }
}

uint64_t Marker::getStartPos(uint32_t raw_index){
    return byte_start[raw_index];
}


void Marker::read_bgen(string bgen_file){
    LOGGER << "Extracting biallelic SNPs from bgen [" << bgen_file << "]..." << std::endl;
    FILE* h_bgen = fopen(bgen_file.c_str(), "rb");
    if(!h_bgen){
        LOGGER.e(0, "can't open bgen file to read");
    }
    uint32_t start_byte = read1Byte<uint32_t>(h_bgen);
    uint32_t start_data_block = start_byte + 4;

    uint32_t len_header = read1Byte<uint32_t>(h_bgen);

    uint32_t n_variants = read1Byte<uint32_t>(h_bgen);
    LOGGER << n_variants << " SNPs, ";

    auto n_sample = read1Byte<uint32_t>(h_bgen);
    LOGGER << n_sample << " samples in the bgen file." << std::endl;

    char magic[5];
    readBytes<char>(h_bgen, 4, magic);

    int32_t skip_byte = (int32_t)len_header - 20;
    if(skip_byte > 0){
        fseek(h_bgen, skip_byte, SEEK_CUR); 
    }else if(skip_byte < 0){
        LOGGER.e(0, "strange header length, might be an invalid bgen.");
    }

    uint32_t flags = read1Byte<uint32_t>(h_bgen);
    int compress_block = flags & 3; // first 2 bit
    int layout = (flags >> 2) & 15; //2-5 bit
    int has_sample = (flags >> 31) & 1;

    LOGGER << "bgen version ";
    switch(layout){
        case 1:
            LOGGER << "1.1, ";
            break;
        case 2:
            LOGGER << "1.2, ";
            break;
        default:
            LOGGER << "unkown, ";
    }
    switch(compress_block){
        case 0:
            LOGGER << "no compress";
            break;
        case 1:
            LOGGER << "compressed by zlib";
            break;
        case 2:
            LOGGER << "compressed by zstd";
            break;
        default:
            LOGGER << "unknown format";
    }
    LOGGER << "." << std::endl;
               
    if(compress_block != 1 || layout != 2){
        LOGGER.e(0, "bgen < 1.2, compress not by zlib is not supported currently");
    }

    LOGGER << "Looking for bialleric allels..." << std::endl;
    fseek(h_bgen, start_data_block, SEEK_SET);
    uint32_t count_chr_error = 0, count_multi_alleles = 0;

    for(int index = 0; index < n_variants; index++){
        auto Lid = read1Byte<uint16_t>(h_bgen);
        fseek(h_bgen, Lid, SEEK_CUR);

        auto len_rs = read1Byte<uint16_t>(h_bgen);
		unique_ptr<char[]> rsid(new char[len_rs + 1]());
        readBytes(h_bgen, len_rs, rsid.get());

        auto len_chr = read1Byte<uint16_t>(h_bgen);
		unique_ptr<char[]> snp_chr(new char[len_chr + 1]());
        readBytes(h_bgen, len_chr, snp_chr.get());
        uint8_t chr_item = 0;
        bool keep_snp = true;
        try{
            chr_item = chr_maps.at(snp_chr.get());
        }catch(std::out_of_range&){
            count_chr_error++;
            keep_snp = false;
        }

        if(chr_item < options_i["start_chr"] || chr_item > options_i["end_chr"]){
            count_chr_error++;
            keep_snp = false;
        }

        auto snp_pos = read1Byte<uint32_t>(h_bgen);

        auto n_alleles = read1Byte<uint16_t>(h_bgen);
        if(n_alleles != 2){
            count_multi_alleles++;
            keep_snp = false;
        }

        auto len_a1 = read1Byte<uint32_t>(h_bgen);
		unique_ptr<char[]> snp_a1(new char[len_a1 + 1]());
        readBytes(h_bgen, len_a1, snp_a1.get());

        auto len_a2 = read1Byte<uint32_t>(h_bgen);
		unique_ptr<char[]> snp_a2(new char[len_a2 + 1]());
        readBytes(h_bgen, len_a2, snp_a2.get());

        uint64_t snp_start = ftell(h_bgen);

        for(int index_allele = 2; index_allele < n_alleles; index_allele++){
            auto len_allele = read1Byte<uint32_t>(h_bgen);
            fseek(h_bgen, len_allele, SEEK_CUR);
        }

        auto len_comp = read1Byte<uint32_t>(h_bgen); 
        fseek(h_bgen, len_comp, SEEK_CUR);
        
        if(keep_snp){
            chr.push_back(chr_item);
            name.push_back(rsid.get());
            gd.push_back(0);
            pd.push_back(snp_pos);
            std::transform(snp_a1.get(), snp_a1.get() + len_a1, snp_a1.get(), toupper);
            std::transform(snp_a2.get(), snp_a2.get() + len_a2, snp_a2.get(), toupper);
            a1.push_back(snp_a1.get());
            a2.push_back(snp_a2.get());
            A_rev.push_back(false);
            byte_start.push_back(snp_start);
        }
    }
    num_marker = name.size();
    index_extract.resize(name.size());
    std::iota(index_extract.begin(), index_extract.end(), 0);
    num_extract = index_extract.size();
    LOGGER.i(0, to_string(num_marker) + " SNPs to be included from bgen file.");
    LOGGER << count_chr_error << " SNPs excluded due to filter on chromosome." << std::endl;
    LOGGER << count_multi_alleles << " SNPs excluded due to multiple alleles." << std::endl; 
    fclose(h_bgen);
    if(num_extract == 0){
        LOGGER.e(0, "0 SNP remain for further analysis.");
    }
}

uint32_t Marker::count_raw(int part) {
    if(part == -1){
        return num_marker;
    }else{
        return raw_limits[part];
    }
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
        boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
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
    addOneFileOption("bgen_file", "", "--bgen", options_in);

    if(options_in.find("m_file") != options_in.end()){
        for(auto & item : options_in["m_file"]){
            std::ifstream file_item((item + ".bim").c_str());
            if(file_item.fail()){
                LOGGER.e(0, "can't read BIM file in [" + item + "].");
            }
            file_item.close();
        }
        options["m_file"] = boost::algorithm::join(options_in["m_file"], "\t");
        boost::replace_all(options["m_file"], "\r", "");
    }

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
    static map<string, string> options;
    static map<string, int> options_i;
    
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
            LOGGER.w(0, "One of the CHR filter has been applied, it has been overrided by --chr");
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
        filterChrFlag = true;
    }

    if(options_in.find("--chrx") != options_in.end()){
        specifiedChrFlag = true;
        if(filterChrFlag){
            LOGGER.w(0, "One of the CHR filter has been applied, it has been overrided by --chrx");
        }
        options_i["start_chr"] = options_i["last_chr_autosome"] + 1;
        options_i["end_chr"] = options_i["start_chr"];
        filterChrFlag = true;
    }

    return 0;
}


void Marker::processMain(){
    LOGGER.e(0, "Marker has no main process this time");
}
