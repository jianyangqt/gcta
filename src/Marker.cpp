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
#include <utility>
#include <sqlite3.h>

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

    if(options.find("pvar_file") != options.end()){
        has_marker = true;
        read_pvar(options["pvar_file"]);
        raw_limits.push_back(num_marker);
    }

    if(options.find("m_pvar") != options.end()){
        has_marker = true;
        read_mpvar(options["m_pvar"]);
    }

    if(options.find("bgen_file") != options.end()){
        has_marker = true;
        read_bgen_index(options["bgen_file"]);
        raw_limits.push_back(num_marker);
        index_extract.resize(num_marker);
        std::iota(index_extract.begin(), index_extract.end(), 0); 

    }

    if(options.find("mbgen_file") != options.end()){
        has_marker = true;
        read_mbgen(options["mbgen_file"]);
    }

    if(!has_marker){
        LOGGER.e(0, "no marker exist");
    }


    if(options.find("extract_file") != options.end()){
        vector<string> extractlist = read_snplist(options["extract_file"]);
        LOGGER << "Get " << extractlist.size() << " SNPs from list [" << options["extract_file"] << "]." << std::endl;
        extract_marker(extractlist, true);
    }

    if(options.find("exclude_file") != options.end()){
        vector<string> excludelist = read_snplist(options["exclude_file"]);
        LOGGER << "Get " << excludelist.size() << " SNPs from list [" << options["exclude_file"] << "]." << std::endl;
        extract_marker(excludelist, false);
    }
    reset_exclude();

    if(options.find("update_ref_allele_file") != options.end()){
        LOGGER.i(0, "Reading reference alleles of SNPs from [" + options["update_ref_allele_file"] + "]...");
        vector<int> field_return;
        vector<string> fields;
        vector<bool> a_rev;
        matchSNPListFile(options["update_ref_allele_file"], 2, field_return, fields, a_rev, true);

        LOGGER.i(0, to_string(index_extract.size()) + " reference alleles are updated."); 
    }

    if(index_extract.size() == 0){
        LOGGER.e(0, "0 SNP remained!");
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
    vector<uint32_t> ref_index_remained;
    ref_index_remained.reserve(marker_index.size());
    if(ref_allele.size()){
        temp_a_rev.reserve(marker_index.size());
        for(int i = 0; i < marker_index.size(); i++){
            uint32_t cur_marker_index = marker_index[i];
            uint32_t cur_ref_index = ref_index[i];
            string cur_ref = ref_allele[cur_ref_index];
            std::transform(cur_ref.begin(), cur_ref.end(), cur_ref.begin(), toupper);
            if(a1[cur_marker_index] == cur_ref){
                temp_a_rev.push_back(A_rev[cur_marker_index]);
                index_remained.push_back(cur_marker_index);
                ref_index_remained.push_back(cur_ref_index);
            }else if(a2[cur_marker_index] == cur_ref){
                temp_a_rev.push_back(!A_rev[cur_marker_index]);
                index_remained.push_back(cur_marker_index);
                ref_index_remained.push_back(cur_ref_index);
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
            uint32_t t_index = ref_index_remained[i];
            uint32_t t_index_from = t_index * field_num;
            for(int j = t_index_from; j < t_index_from + field_num; j++){
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
        read_bim(mfile);
        raw_limits.push_back(num_marker);
    }
}

void Marker::read_mpvar(string mfile){
    vector<string> mfiles;
    boost::split(mfiles, mfile, boost::is_any_of("\t "));
    raw_limits.clear();
    for(auto & mfile : mfiles){
        read_pvar(mfile);
        raw_limits.push_back(num_marker);
    }
}


void Marker::read_mbgen(string mbgen_file){
    vector<string> mfiles;
    boost::split(mfiles, mbgen_file, boost::is_any_of("\t"));
    raw_limits.clear();
    for(auto & mfile : mfiles){
        read_bgen_index(mfile);
        raw_limits.push_back(num_marker);
    }

    index_extract.resize(num_marker);
    std::iota(index_extract.begin(), index_extract.end(), 0); 
}


//cur_marker_index:  is the index point to position of index_extract
vector<uint32_t> Marker::getNextWindowIndex(uint32_t cur_marker_index, uint32_t window, bool& chr_ends, bool& isX, bool retRaw){
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
    bool success;
    isX = (cur_chr == mapCHR("X", success));

    for(uint32_t marker_index = cur_marker_index; marker_index < index_extract.size(); marker_index++){
        uint32_t temp_index = index_extract[marker_index];
        if(chr[temp_index] != cur_chr){
            chr_ends = true;
            break;
        }
        if(pd[temp_index] <= final_pd){
            indices.push_back(retRaw ? temp_index : marker_index);
        }else{
            break;
        }
    }
    return indices;
}

uint32_t Marker::getNextSize(const vector<uint32_t> &rawRef, uint32_t curExtractIndex, uint32_t num, int &fileIndex, bool &chr_ends, uint8_t &isSexXY){
    static bool success;
    static uint8_t Xchr = mapCHR("X", success);
    static uint8_t Ychr = mapCHR("Y", success);
    if(curExtractIndex >= rawRef.size()){
        chr_ends = true;
        return 0;
    }
    uint32_t cur_index = rawRef[curExtractIndex];
    uint8_t cur_chr = chr[cur_index];


    chr_ends = false;
    
    if(cur_chr == Xchr){
        isSexXY = 1;
    }else if(cur_chr == Ychr){
        isSexXY = 2;
    }else{
        isSexXY = 0;
    }

    //isSexXY = (cur_chr == mapCHR("X", success)); // 1 X
    //isSexXY = (cur_chr == mapCHR("Y", success)? 2 : 0); // 2 Y;

    fileIndex = getMIndex(cur_index);
    uint32_t rawIndexLimits = raw_limits[fileIndex];

    uint32_t retSize = 0; 
    for(uint32_t marker_index = curExtractIndex; marker_index < curExtractIndex + num; marker_index++){
        if(marker_index >= rawRef.size()){
            chr_ends = true;
            break;
        }
        uint32_t temp_index = rawRef[marker_index];
        if(chr[temp_index] != cur_chr){
            chr_ends = true;
            break;
        }
        if(temp_index >= rawIndexLimits){
            break;
        }
        retSize++;
    }
    return retSize;
}



vector<uint32_t> Marker::getNextSizeIndex(uint32_t cur_marker_index, uint32_t num, bool& chr_ends, bool& isX, bool retRaw){
    vector<uint32_t> indices;
    indices.reserve(num);
    if(cur_marker_index >= index_extract.size()){
        chr_ends = true;
        return indices;
    }
    uint32_t cur_index = index_extract[cur_marker_index];
    uint8_t cur_chr = chr[cur_index];

    chr_ends = false;
    bool success;
    isX = (cur_chr == mapCHR("X", success));

    for(uint32_t marker_index = cur_marker_index; marker_index < cur_marker_index + num; marker_index++){
        if(marker_index >= index_extract.size()){
            chr_ends = true;
            break;
        }
        uint32_t temp_index = index_extract[marker_index];
        if(chr[temp_index] != cur_chr){
            chr_ends = true;
            break;
        }
        indices.push_back(retRaw ? temp_index : marker_index);
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

void Marker::read_pvar(string pvar_file){
    LOGGER.i(0, "Reading PLINK2 PVAR file from [" + pvar_file + "]...");
    vector<string> head;
    map<int, vector<string>> lists;
    int nHeader = 0;
    if(readTxtList(pvar_file, 5, head, nHeader, lists)){
        int ncol = lists.size();
        int iChr, iPOS, iID, iA1, iA2;
        if(head.empty()){
            if(ncol==6){
                iChr = 0;
                iID = 1;
                iPOS = 3;
                iA1 = 4;
                iA2 = 5;
            }else if(ncol == 5){
                iChr = 0;
                iID = 1;
                iPOS = 2;
                iA1 = 3;
                iA2 = 4;
            }else{
                LOGGER.e(0, "can't read variant list file with " + to_string(ncol) + " columns, not a valid BIM file.");
            }
        }else{
            if(head[0]=="#CHROM"){
                iChr = 0;
                int validHeadCT = 1;
                bool found;
                iPOS = findElementVector(head, string("POS"), found);
                if(found)validHeadCT++;

                iID = findElementVector(head, string("ID"), found);
                if(found) validHeadCT++;

                iA2 = findElementVector(head, string("REF"), found);
                if(found) validHeadCT++;

                iA1 = findElementVector(head, string("ALT"), found);
                if(found) validHeadCT++;

                if(validHeadCT != 5){
                    LOGGER.e(0, "can't find enough essential columns in PVAR file.");
                }
            }else{
                LOGGER.e(0, "invaild PVAR file, should start with #CHROM.");
            }
        }


        uint32_t nrows = lists[0].size();
        uint32_t oriSize = chr.size();
        uint32_t newSize = oriSize + nrows;
        chr.resize(newSize);
        name.resize(newSize);
        pd.resize(newSize);
        a1.resize(newSize);
        a2.resize(newSize);
        A_rev.resize(newSize, false);
        byte_start.resize(newSize, 1);
        vector<uint8_t> validSNP(nrows);
        #pragma omp parallel for
        for(int i = 0; i < nrows; i++){
            uint32_t curRow = i + oriSize;
            bool success;
            chr[curRow] = mapCHR(lists[iChr][i], success);
            if(success){
                validSNP[i] = 1;
            }else{
                validSNP[i] = 0;
            }
            name[curRow] = lists[iID][i];
            try{
                //gd here if need;
                pd[curRow] = std::stoi(lists[iPOS][i]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "Line " + to_string(nHeader + i + 1) + " contains illegal distance value.");
            }
            auto &curA1 = lists[iA1][i];
            std::transform(curA1.begin(), curA1.end(), curA1.begin(), toupper);
            a1[curRow] = curA1;

            auto &curA2 = lists[iA2][i];
            std::transform(curA2.begin(), curA2.end(), curA2.begin(), toupper);
            a2[curRow] = curA2;
        }

        index_extract.reserve(newSize);
        uint32_t nValidSNP = 0;
        for(int i = 0; i < nrows; i++){
            uint32_t curRow = i + oriSize;
            if(validSNP[i]){
                index_extract.push_back(curRow);
                nValidSNP++;
            }
        }
        index_extract.shrink_to_fit();

        num_marker = name.size();
        num_extract = index_extract.size();

        LOGGER.i(0, to_string(num_marker) + " SNPs to be included from PVAR file(s).");
        if(nrows != nValidSNP){
            LOGGER.i(0, to_string(num_extract) + " SNPs to be included from valid chromosome number");
        }

        MarkerParam markerParam;
        markerParam.rawCountSNP = nrows;
        markerParam.rawCountSample = 0; // dummy 
        markerParam.compressFormat = 0;
        markerParam.posGenoDataStart = 3; // just dummy, pgen don't start with 3
        markerParams.push_back(markerParam);
    }else{
        LOGGER.e(0, "invalid PVAR file.");
    }
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
    LOGGER.i(0, to_string(num_marker) + " SNPs to be included from BIM file(s).");
    if(num_marker != num_extract){
        LOGGER.i(0, to_string(num_extract) + " SNPs to be included from valid chromosome number");
    }
    bim.close();
    MarkerParam markerParam;
    markerParam.rawCountSNP = num_marker - start_line_number + 1;
    markerParam.rawCountSample = 0;
    markerParam.compressFormat = 0;
    markerParam.posGenoDataStart = 3;
    markerParams.push_back(markerParam);
}

vector<pair<string, vector<uint32_t>>> Marker::read_gene(string gfile){
    std::ifstream genein(gfile.c_str());
    if(!genein){
        LOGGER.e(0, "can not open the file [" + gfile + "] to read");
    }

    string line;
    std::getline(genein, line);
    vector<string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
    int nElements;
    if(line_elements[0] == "gene" && line_elements[1] == "chr" && 
            line_elements[2] == "start" && line_elements[3] == "end" && line_elements[4] == "strand"){
        nElements = line_elements.size();
    }else{
        LOGGER.e(0, "invalid gene list file");
    }

    int line_number  = 1;
    uint8_t chr_item;
    uint8_t curChrItem = 0;
    bool bFind = false;
    uint32_t chrStartIndex = 0;
    uint32_t chrEndIndex = 0;
    vector<uint32_t> curIndices;

    vector<pair<string, vector<uint32_t>>> gene;
    while(std::getline(genein, line)){
        line_number++;
        vector<string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
        if(line_elements.size() != nElements) {
            LOGGER.e(0, "  line " + to_string(line_number)
                   + " has elements not equal to " + to_string(nElements));
        }
        // replace strings
        try{
            chr_item = chr_maps.at(line_elements[1]);
        }catch(std::out_of_range&){
            LOGGER.e(0, "  Line " + to_string(line_number) + " of [" + gfile +
                    "] contains invalid chr, please check");
        }
        
        string curGeneName = line_elements[0];
        int curStart, curEnd; 
        try{
            curStart = std::stoi(line_elements[2]);
            curEnd = std::stoi(line_elements[3]);
        }catch(std::invalid_argument&){
           LOGGER.e(0, "  Line " + to_string(line_number) + " of [" + gfile +
                    "] contains invalid position value, please check");
        }

        if(chr_item != curChrItem){
            bool isSetStart = false;
            bool isSetEnd = false;
            for(uint32_t marker_index = 0; marker_index < index_extract.size(); marker_index++){
                uint32_t temp_index = index_extract[marker_index];
                //LOGGER << temp_index << "\t" << (int)chr_item << "\t" << (int)chr[temp_index] << std::endl;
                if(chr[temp_index] == chr_item){
                    if((!isSetStart)){
                        chrStartIndex = marker_index;
                        isSetStart = true;
                    }
                }else{
                    if(isSetStart){
                        chrEndIndex = marker_index - 1;
                        isSetEnd = true;
                        break;
                    }
                }
            }
            if(isSetStart){
                if(!isSetEnd) chrEndIndex = index_extract.size() - 1;
            }else{
                //LOGGER << "gene " << curGeneName << " can not find, chr: " <<  int(chr_item) << std::endl;
                //curChrItem = chr_item;
                continue;
                //LOGGER.e(0, " Line " + to_string(line_number) + ", chr information can't be found in the genotype");
            }
            curChrItem = chr_item;
        }

        vector<uint32_t> indicies;
        for(uint32_t i = chrStartIndex; i <= chrEndIndex; i++){
            uint32_t temp_index = index_extract[i];
            uint32_t cur_pd = pd[temp_index];
            if(cur_pd >= curStart && cur_pd <= curEnd){
                indicies.push_back(i);
            }
        }
        if(indicies.size() == 0) continue;
        pair<string, vector<uint32_t>> pair_gene = std::make_pair(curGeneName+"\t" + to_string(curStart) + "\t" + to_string(curEnd), indicies);

        gene.emplace_back(pair_gene);
    }
    return gene;
}


void Marker::getStartPosSize(uint32_t raw_index, uint64_t &pos, uint64_t &size){
    pos = byte_start[raw_index];
    size = byte_size[raw_index];
}

MarkerParam Marker::getBgenMarkerParam(FILE *h_bgen, string &outputs){
    // check file formats
    uint32_t start_byte = read1Byte<uint32_t>(h_bgen);
    uint32_t start_data_block = start_byte + 4;

    uint32_t len_header = read1Byte<uint32_t>(h_bgen);

    uint32_t n_variants_total = read1Byte<uint32_t>(h_bgen);

    auto n_sample = read1Byte<uint32_t>(h_bgen);
    outputs = outputs + to_string(n_sample) + " samples, " + to_string(n_variants_total) + " SNPs.\n";

    char magic[5];
    readBytes<char>(h_bgen, 4, magic);
    if(strncmp(magic, "bgen", 4) != 0){
        LOGGER.e(0, "bad magic number in the bgen file.");
    }

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

    outputs += "BGEN version ";
    switch(layout){
        case 1:
            outputs += "1.1, ";
            break;
        case 2:
            outputs += "1.2, ";
            break;
        default:
            outputs += "unkown, ";
    }
    switch(compress_block){
        case 0:
            outputs += "no compress";
            break;
        case 1:
            outputs += "compressed by zlib";
            break;
        case 2:
            outputs += "compressed by zstd";
            break;
        default:
            outputs += "unknown format";
    }
    outputs += ".";

    if(layout != 2){
        LOGGER.e(0, "GCTA only support bgen version 1.2, 1.3. Use QCTOOLv2 to convert to new version.");
    }

    if(has_sample){
        LOGGER.w(0, "GCTA reads sample information from '--sample' inputs, but ignore the built-in sample data.");
    }

    MarkerParam param;
    param.rawCountSNP = n_variants_total;
    param.rawCountSample = n_sample;
    param.compressFormat = compress_block;
    param.posGenoDataStart = start_data_block;
    return param;
}


#define RCSTR(TEMP) std::string(reinterpret_cast<const char*>(TEMP))

void Marker::read_bgen_index(string bgen_file){
    sqlite3 *db;
    int rc;
    string index_fname = bgen_file + ".bgi";
    string query_file = "file:" + index_fname + "?nolock=1";
    LOGGER.i(0, "Loading bgen index from [" + index_fname + "]...");
    rc = sqlite3_open_v2(query_file.c_str(), &db, SQLITE_OPEN_READONLY | SQLITE_OPEN_URI, NULL);

    //string prompt_index = "'gcta64 --bgen test.bgen --bgen-index --out test.bgen.bgi' or 'bgenix -g test.bgen -index'";
    string prompt_index = "'bgenix -g test.bgen -index'";
    if(rc){
        LOGGER.e(0, "can't open index file: " + string(sqlite3_errmsg(db)) + 
                "\nIndex the bgen file by " + prompt_index + ".");
    }

    sqlite3_stmt *stmt;

    //check the bgen file
    const char * sql_meta = "SELECT * FROM Metadata";
    rc = sqlite3_prepare_v2(db, sql_meta, -1, &stmt, NULL);
    if(rc != SQLITE_OK){
        LOGGER.e(0, "bad index file: " + string(sqlite3_errmsg(db)) +
                "\nTry to regenerate the index by " + prompt_index + ".");
    }
    string filename;
    uint64_t file_size;
    const char * first1000bytes;
    
    rc = sqlite3_step(stmt);
    if(rc == SQLITE_ROW){
        filename = RCSTR(sqlite3_column_text(stmt, 0));
        file_size = sqlite3_column_int64(stmt, 1);
        first1000bytes = reinterpret_cast< const char* >(sqlite3_column_blob(stmt, 3));
    }

    //check 1000 bytes and file size;
    int num_firstNbytes = file_size < 1000 ? file_size : 1000;

    FILE * h_bgen = fopen(bgen_file.c_str(), "rb");
    if(h_bgen == NULL){
        LOGGER.e(0, "can't read bgen file [" + bgen_file + "], " + string(strerror(errno)));
    }

    char * geno_first1000bytes = new char[num_firstNbytes];
    if(num_firstNbytes != fread(geno_first1000bytes, sizeof(char), num_firstNbytes, h_bgen)){
        LOGGER.e(0, "can't read bgen file [" + bgen_file + "], " + string(strerror(errno)));
    }

    if(strncmp(first1000bytes, geno_first1000bytes, num_firstNbytes) != 0){
        LOGGER.e(0, "bad index file, first 1000 bytes aren't consistent."
                "\nTry to regenerate the index by " + prompt_index + ".");
    }

    fseek(h_bgen, 0L, SEEK_END);
    if(file_size != ftell(h_bgen)){
        LOGGER.e(0, "bad index file, file size isn't consistent."
                "\nTry to regenerate the index by " + prompt_index + ".");
    }
    rewind(h_bgen);
    rc = sqlite3_reset(stmt);

    // check variants;
    const char * sql_allvar = "SELECT count(*) FROM Variant";
    const char * sql_var1 = "SELECT * FROM Variant";
    rc = sqlite3_prepare_v2(db, sql_allvar, -1, &stmt, NULL);
    if(rc != SQLITE_OK){
        LOGGER.e(0, "bad index file: " + string(sqlite3_errmsg(db)) +
                "\nTry to regenerate the index by " + prompt_index + ".");
    }
    rc = sqlite3_step(stmt);

    int n_variants_total_index;
    if(rc == SQLITE_ROW){
        n_variants_total_index = sqlite3_column_int(stmt, 0);
    }
    rc = sqlite3_reset(stmt);

    string outputs;
    MarkerParam markerParam = getBgenMarkerParam(h_bgen, outputs);
    LOGGER << outputs << std::endl;
    if(markerParam.rawCountSNP != n_variants_total_index){
        LOGGER.e(0, "bad index file, the indexed SNPs are different from bgen file."
                "\nTry to regenerate the index by " + prompt_index + ".");
    }
    markerParams.push_back(markerParam);

    // load the index from bgi
    const char *sql_count = "SELECT count(*) FROM Variant WHERE number_of_alleles=2";
    const char *sql = "SELECT chromosome,position,rsid,allele1,allele2,file_start_position,size_in_bytes FROM Variant WHERE number_of_alleles=2";

    // count the alleles
    rc = sqlite3_prepare_v2(db, sql_count, -1, &stmt, NULL);
    if(rc != SQLITE_OK){
        LOGGER.e(0, "bad index file: " + string(sqlite3_errmsg(db)) +
                "\nTry to regenerate the index by " + prompt_index + ".");
    }
    int n_variants = 0;
    while((rc = sqlite3_step(stmt)) == SQLITE_ROW){
        n_variants = sqlite3_column_int(stmt, 0);
    }
    if(rc != SQLITE_DONE){
       LOGGER.e(0, "bad index file: " + string(sqlite3_errmsg(db)) +
                "\nTry to regenerate the index by " + prompt_index + ".");
    }
    rc = sqlite3_reset(stmt);
    if(n_variants == 0){
        LOGGER.w(0, "None variant with bialleric SNPs.");
        return;
    }
    int cur_total_num_variants = chr.size() + n_variants;
    chr.reserve(cur_total_num_variants);
    name.reserve(cur_total_num_variants);
    gd.reserve(cur_total_num_variants);
    pd.reserve(cur_total_num_variants);
    a1.reserve(cur_total_num_variants);
    a2.reserve(cur_total_num_variants);
    A_rev.reserve(cur_total_num_variants);
    byte_start.reserve(cur_total_num_variants);
    byte_size.reserve(cur_total_num_variants);
    LOGGER.i(0, to_string(n_variants) + " SNPs to be included from bgen index file.");

    // retrieve each variants
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if(rc != SQLITE_OK){
        LOGGER.e(0, "bad index file: " + string(sqlite3_errmsg(db)) +
                "\nTry to regenerate the index by " + prompt_index + ".");
    }

    int count_chr_error = 0;
    while((rc = sqlite3_step(stmt)) == SQLITE_ROW){
        string snp_chr = RCSTR(sqlite3_column_text(stmt, 0));
        uint8_t chr_item = 0;
        bool keep_snp = true;
        try{
            chr_item = std::stoi(snp_chr);
        }catch(std::invalid_argument&){
            try{
                chr_item = chr_maps.at(snp_chr);
            }catch(std::out_of_range&){
                count_chr_error++;
                keep_snp = false;
            }
        }

        if(chr_item < options_i["start_chr"] || chr_item > options_i["end_chr"]){
            count_chr_error++;
            keep_snp = false;
        }
        if(keep_snp){
            chr.push_back(chr_item);
            name.push_back(RCSTR(sqlite3_column_text(stmt, 2)));
            gd.push_back(0);
            pd.push_back(sqlite3_column_int(stmt, 1));
            string snp_a1 = RCSTR(sqlite3_column_text(stmt, 3));
            string snp_a2 = RCSTR(sqlite3_column_text(stmt, 4));
            std::transform(snp_a1.begin(), snp_a1.end(), snp_a1.begin(), toupper);
            std::transform(snp_a2.begin(), snp_a2.end(), snp_a2.begin(), toupper);
            a1.push_back(snp_a1);
            a2.push_back(snp_a2);

            A_rev.push_back(false);
            byte_start.push_back(sqlite3_column_int64(stmt, 5));
            uint64_t cur_size = sqlite3_column_int64(stmt, 6);
            byte_size.push_back(cur_size);
            if(maxGeno1ByteSize < cur_size){
                maxGeno1ByteSize = cur_size;
            }
        }
    }

    if (rc != SQLITE_DONE) {
        LOGGER.e(0, string(sqlite3_errmsg(db)) + "\nBad index file, regenerate by " + prompt_index + ".");
    }

    sqlite3_finalize(stmt);
    sqlite3_close(db);
    num_marker = cur_total_num_variants - count_chr_error;
    int num_var_added = n_variants - count_chr_error;

    // check the first and last SNPs
    if(num_var_added > 0){
        int index1 = chr.size() - num_var_added;
        int index2 = index1 + num_var_added - 1;
        uint64_t pos1 = byte_start[index1];
        uint64_t pos2 = byte_start[index2];
        MarkerInfo m1 = extractBgenMarkerInfo(h_bgen, pos1);
        MarkerInfo m2 = extractBgenMarkerInfo(h_bgen, pos2);
        vector<string> allele1 = {a1[index1], a2[index1]};
        vector<string> allele2 = {a1[index2], a2[index2]};
        bool success;
        if(mapCHR(m1.chr, success) != chr[index1] || m1.name != name[index1] || m1.pd != pd[index1] || m1.alleles != allele1){
            LOGGER.e(0, "The first variant in bgen file aren't consistent with index file."
                "\nTry to regenerate the index by " + prompt_index + ".");
        }

        if(mapCHR(m2.chr, success) != chr[index2] || m2.name != name[index2] || m2.pd != pd[index2] || m2.alleles != allele2){
            LOGGER.e(0, "The variants in bgen file aren't consistent with index file."
                "\nTry to regenerate the index by " + prompt_index + ".");
        }
    }

    fclose(h_bgen);

    if(count_chr_error > 0){
        LOGGER << count_chr_error << " SNPs excluded due to filter on chromosome. " << std::endl;
    }

    LOGGER << "Total SNPs included: " << num_var_added  << "/" <<  num_marker << "." << std::endl;
 
}

uint64_t Marker::getMaxGenoMarkerUptrSize(){
    int suptr = sizeof(uintptr_t);
    return (maxGeno1ByteSize + suptr - 1) / suptr; 
}

MarkerInfo Marker::extractBgenMarkerInfo(FILE *h_bgen, uint64_t &pos){
    fseek(h_bgen, pos, SEEK_SET);
    auto Lid = read1Byte<uint16_t>(h_bgen);
    fseek(h_bgen, Lid, SEEK_CUR);

    auto len_rs = read1Byte<uint16_t>(h_bgen);
    unique_ptr<char[]> rsid(new char[len_rs + 1]());
    readBytes(h_bgen, len_rs, rsid.get());

    auto len_chr = read1Byte<uint16_t>(h_bgen);
    unique_ptr<char[]> snp_chr(new char[len_chr + 1]());
    
    readBytes(h_bgen, len_chr, snp_chr.get());

    string chr(snp_chr.get());

    auto snp_pos = read1Byte<uint32_t>(h_bgen);

    auto n_alleles = read1Byte<uint16_t>(h_bgen);
    vector<string> alleles;
    for(int i = 0; i < n_alleles; i++){
        auto len_a1 = read1Byte<uint32_t>(h_bgen);
        unique_ptr<char[]> snp_a1(new char[len_a1 + 1]());
        readBytes(h_bgen, len_a1, snp_a1.get());
        std::transform(snp_a1.get(), snp_a1.get() + len_a1, snp_a1.get(), toupper);
        alleles.push_back(snp_a1.get());
    }

    uint64_t snp_start = ftell(h_bgen);

    //auto len_comp = read1Byte<uint32_t>(h_bgen); 
    //fseek(h_bgen, len_comp, SEEK_CUR);

    MarkerInfo markerInfo;
    markerInfo.chr = chr;
    markerInfo.name = string(rsid.get());
    markerInfo.pd = snp_pos;
    markerInfo.alleles = alleles;
    markerInfo.A_rev = false;

    pos = snp_start; // set the position easy for seek

    return markerInfo;
}
    

MarkerParam Marker::getMarkerParams(int part_num){
    if(part_num < markerParams.size()){
        return markerParams[part_num];
    }else{
        LOGGER.e(0, "Get Marker Params out of range.");
        return markerParams[0]; // dummy return
    }
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

    if(layout != 2){
        LOGGER.e(0, "GCTA only support bgen version 1.2, 1.3. Use QCTOOL to convert to new version.");
    }
               
    if(compress_block != 1){
        LOGGER.e(0, "Compress not by zlib is not supported currently");
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
            chr_item = std::stoi(string(snp_chr.get()));
        }catch(std::invalid_argument&){
            try{
                chr_item = chr_maps.at(snp_chr.get());
            }catch(std::out_of_range&){
                count_chr_error++;
                keep_snp = false;
            }
        }
        //if(errno == ERANGE){
        //if(errno != 0){

        if(chr_item < options_i["start_chr"] || chr_item > options_i["end_chr"]){
            count_chr_error++;
            keep_snp = false;
        }

        //LOGGER << options_i["start_chr"] << " " << options_i["end_chr"] << " " << (int)chr_item << " count_chr_error:" << count_chr_error << std::endl;
        //LOGGER.e(0, "test")


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

uint8_t Marker::mapCHR(string chr_str, bool &success){
    uint8_t chr_item = 0;
    bool keep_snp = true;
    try{
        chr_item = std::stoi(chr_str);
    }catch(std::invalid_argument&){
        try{
            chr_item = chr_maps.at(chr_str);
        }catch(std::out_of_range&){
            keep_snp = false;
        }
    }

    if(chr_item < options_i["start_chr"] || chr_item > options_i["end_chr"]){
        keep_snp = false;
    }
    success = keep_snp;
    return chr_item;
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

vector<uint32_t>& Marker::get_extract_index(){ // raw index
    return this->index_extract;
}

vector<uint32_t> Marker::get_extract_index_autosome(){ // returns the index in keep list not raw
    uint8_t last_auto_chr = options_i["last_chr_autosome"] + 1;
    uint8_t xy_chr = options_i["last_chr_autosome"] + 3;
    vector<uint32_t> auto_index;
    for(int i = 0; i < index_extract.size(); i++){
        uint32_t curIndex = index_extract[i];
        if(chr[curIndex] < last_auto_chr || chr[curIndex] == xy_chr){
            auto_index.push_back(i);
        }
    }
    return auto_index;
}

vector<uint32_t> Marker::get_extract_index_X(){
   uint8_t chrx = options_i["last_chr_autosome"] + 1;
   vector<uint32_t> auto_index;
   for(int i = 0; i < index_extract.size(); i++){
       uint32_t curIndex = index_extract[i];
       if(chr[curIndex] == chrx){
           auto_index.push_back(i);
       }
   }
   return auto_index;
}


void Marker::extract_marker(vector<string> markers, bool isExtract) {
   vector<uint32_t> ori_index, marker_index;
   size_t numOriMarker = markers.size();
   removeDuplicateSort(markers);
   int numDup = numOriMarker - markers.size();
   if(numDup != 0) LOGGER.w(0, to_string(numDup) + " duplicated SNPs were ignored in the list." );

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

string Marker::get_marker(int rawindex, bool bflip){ // raw index
    string return_string = std::to_string(chr[rawindex]) + "\t" + name[rawindex] + "\t" + 
        std::to_string(pd[rawindex]) + "\t";
    if(A_rev[rawindex] ^ bflip){
        return return_string + a2[rawindex] + "\t" + a1[rawindex];
    }else{
        return return_string + a1[rawindex] + "\t" + a2[rawindex];
    }
}

string Marker::getMarkerStrExtract(int extractindex, bool bflip){ // extract index
    return get_marker(getRawIndex(extractindex), bflip);
}
 

bool Marker::isInExtract(uint32_t index) {
    if(index_extract.size() > index_exclude.size()){
        return !std::binary_search(index_exclude.begin(), index_exclude.end(), index);
    }else{
        return std::binary_search(index_extract.begin(), index_extract.end(), index);
    }
}

uint32_t Marker::getRawIndex(uint32_t extractedIndex){
    return index_extract[extractedIndex];
}

bool Marker::isEffecRev(uint32_t extractedIndex){
    return A_rev[index_extract[extractedIndex]];
}

bool Marker::isEffecRevRaw(uint32_t rawIndex){
    return A_rev[rawIndex];
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

    addOneFileOption("pvar_file", ".pvar", "--pfile", options_in);
    addOneFileOption("marker_file", ".bim", "--bpfile", options_in);

    addMFileListsOption("mbgen_file", ".bgen", "--mbgen", options_in, options);
    addMFileListsOption("m_pvar", ".pvar", "--mpfile", options_in, options);
    addMFileListsOption("m_file", ".bim", "--mbpfile", options_in, options);
    addMFileListsOption("m_file", ".bim", "--mbfile", options_in, options);
        
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
    }
    if(options_i.find("last_chr_autosome") == options_i.end()){
        options_i["last_chr_autosome"] = 22;
    }
    // include the x y xy MT
    options_i["last_chr"] = options_i["last_chr_autosome"] + 4;
    if(options_i.find("start_chr") == options_i.end()){
        options_i["start_chr"] = 0;
        options_i["end_chr"] = options_i["last_chr"];
    }

    bool filterChrFlag = false;
    //static map<string, string> options;
    //static map<string, int> options_i;
    
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
