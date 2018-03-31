/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: holds phenotype information in plink format

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#include "Pheno.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "constants.hpp"
#include "Logger.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include "utils.hpp"
#include <cstring>
#include <boost/algorithm/string.hpp>
#include <set>

using std::to_string;

map<string, string> Pheno::options;

const uintptr_t k1LU = (uintptr_t)1;

Pheno::Pheno() {
    if(options.find("m_file") != options.end()){
   }
            
        
    if(options.find("pheno_file") != options.end()){
        this->read_fam(options["pheno_file"]);
    }else{
        LOGGER.e(0, "no phenotype file presents");
    }

    if(options.find("keep_file") != options.end()){
        vector<string> keep_subjects = read_sublist(options["keep_file"]);
        set_keep(keep_subjects, mark, index_keep,  true);
    }

    if(options.find("remove_file") != options.end()){
        vector<string> remove_subjects = read_sublist(options["remove_file"]);
        set_keep(remove_subjects, mark, index_keep, false);
    }

    if(options.find("qpheno_file") != options.end()){
        vector<vector<double>> phenos;
        LOGGER.i(0, "Reading phenotype data from [" + options["qpheno_file"] + "]...");
        vector<string> pheno_subjects = read_sublist(options["qpheno_file"], &phenos);
        int cur_pheno = 1;
        if(options.find("mpheno") != options.end()){
            try{
                cur_pheno = std::stoi(options["mpheno"]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "--mpheno isn't a numberic value");
            }
        }

        if(cur_pheno <= 0 || cur_pheno > phenos.size()){
            LOGGER.e(0, "selected pheno column can't be less than 0 or larger than --pheno columns");
        }

        cur_pheno -= 1;

        update_pheno(pheno_subjects, phenos[cur_pheno]);
        LOGGER.i(0, to_string(index_keep.size()) + " overlapped individuals with non-missing data to be included from the phenotype file.");
        
    }

    if(options.find("sex_file") != options.end()){
        vector<vector<double>> phenos;
        LOGGER.i(0, "Reading gender information from [" + options["sex_file"] + "]...");
        vector<string> subjects = read_sublist(options["sex_file"], &phenos);
        vector<double> sex_info = phenos[0];
        update_sex(subjects, sex_info);
        LOGGER.i(0, to_string(index_keep.size()) + " individuals with valid sex information to be included from the phenotype file.");
    }

    if(options.find("filter_sex") != options.end()){
        filter_sex();
    }

    reinit();
}

void Pheno::filter_keep_index(vector<uint32_t>& k_index){
    if(k_index.size() == index_keep.size()){return;}
    vector<uint32_t> index_keep2(k_index.size(), 0);
    std::transform(k_index.begin(), k_index.end(), index_keep2.begin(), 
            [this](size_t pos){return this->index_keep[pos];});
    index_keep = index_keep2;

    reinit();
}

void Pheno::filter_sex(){
    vector<uint32_t> new_index;
    vector<uint32_t> male_index;
    new_index.reserve(index_keep.size());
    uint32_t valid_index = 0;
    for(int i = 0; i < index_keep.size(); i++){
        uint32_t index = index_keep[i];
        int8_t sex_item = sex[index];
        switch(sex_item){
            case 1:
                male_index.push_back(valid_index);
                new_index.push_back(index);
                valid_index++;
                break;
            case 2:
                new_index.push_back(index);
                valid_index++;
                break;
        }
    }

    new_index.shrink_to_fit();
    index_keep = new_index;
    index_keep_male = male_index;
    LOGGER.i(0, to_string(new_index.size()) + " individuals have gender information.");
}


void Pheno::reinit(){
    reinit_rm(index_keep, index_rm, num_ind);
    num_keep = index_keep.size();
    num_rm = index_rm.size();
    if(num_keep == 0){
        LOGGER.e(0, "0 individual remain for further analysis.");
    }
    //init_mask_block();
    
}


// TODO filter the non-number strings other than nan
vector<string> Pheno::read_sublist(string sublist_file, vector<vector<double>> *phenos, vector<int> *keep_row_p) {
    vector<string> subject_list;
    std::ifstream sublist(sublist_file.c_str());
    if(!sublist.good()){
        LOGGER.e(0, "cann't read [" + sublist_file + "]");
    }
    vector<int> keep_row;
    string err_file = "the subject list file [" + sublist_file + "]";

    string line;
    int line_number = 0;
    int last_length = 0;
    bool init_pheno = true;
    int large_elements = 0;
    while(std::getline(sublist, line)){
        line_number++;
        vector<string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        int num_elements = line_elements.size();
        if(num_elements < 2){
            LOGGER.e(0, err_file + ", line " + to_string(line_number) +
                        " has elements less than 2");
        }
        if(phenos && init_pheno){
            if(keep_row_p){
                keep_row = *keep_row_p;
                phenos->resize(keep_row.size());
            }else{
                phenos->resize(num_elements - 2);
                keep_row.resize(phenos->size());
                std::iota(keep_row.begin(), keep_row.end(), 0);
            }
            large_elements = keep_row[keep_row.size() - 1] + 1 + 2;
            if(large_elements > num_elements){
                LOGGER.e(0, err_file + " has not enough column to read");
            }
            init_pheno = false;
        }

        if(line_number > 1 && num_elements != last_length){
            string errmsg = err_file + ", line " + to_string(line_number) +
                        " has different elements";
            LOGGER.w(0, errmsg);
        }

        subject_list.push_back(line_elements[0] + "\t" + line_elements[1]);
        if(phenos){
            if(large_elements > num_elements){
                LOGGER.e(0, err_file + ", line " + to_string(line_number) +
                        " has not enough elements");
            }
            for(int index = 0; index != keep_row.size(); index++){
                const char *temp_str = line_elements[index+2].c_str();
                char* pEnd;
                double temp_double = strtod(temp_str, &pEnd);
                if(strlen(temp_str) != pEnd - temp_str){ 
                    temp_double = strtod("nan", NULL);
                }
                (*phenos)[index].push_back(temp_double);
                
            }
        }

        last_length = num_elements;
    }
    sublist.close();
    return subject_list;
}

void Pheno::read_fam(string fam_file) {
    LOGGER.i(0, "Reading PLINK FAM file from [" + fam_file + "]...");
    std::ifstream fam(fam_file.c_str());
    if(!fam){
        LOGGER.e(0, "can not open the file [" + fam_file + "] to read");
    }

    int line_number = 0;
    int last_length = 0;
    string line;
    while(std::getline(fam, line)){
        line_number++;
        vector<string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t ")); 
        if(line_elements.size() < Constants::NUM_FAM_COL) {
            LOGGER.e(0, "the fam file [" + fam_file + "], line " + to_string(line_number)
                   + " has elements less than " + to_string(Constants::NUM_FAM_COL));
        }
        if(line_number > 1 && line_elements.size() != last_length){
            LOGGER.w(0, "the fam file [" + fam_file + "], line " + to_string(line_number) + " have different elements");
        }
        fid.push_back(line_elements[0]);
        pid.push_back(line_elements[1]);
        mark.push_back(line_elements[0] + "\t" + line_elements[1]);
        fa_id.push_back(line_elements[2]);
        mo_id.push_back(line_elements[3]);
        sex.push_back(std::stoi(line_elements[4]));
        pheno.push_back(strtod("nan", NULL)); //5
        last_length = line_elements.size();

        index_keep.push_back(line_number - 1);
        //LOGGER.p(0, to_string(line_number) + " subjects in fam processed");
    }
    num_ind = fid.size();
    num_bytes = (num_ind + 3) / 4;
    num_keep = index_keep.size();
    LOGGER.i(0, to_string(num_ind) + " individuals to be included from FAM file.");
    fam.close();
}

vector<string> Pheno::get_id(int from_index, int to_index){
    vector<string> out_id;
    int raw_index;
    out_id.reserve(num_keep);
    for(int index = from_index; index <= to_index; index++){
        raw_index = index_keep[index];
        //out_id.push_back(fid[raw_index] + "\t" + pid[raw_index]);
        out_id.push_back(mark[raw_index]);
    }
    return out_id;
}

// ori_index, original raw count in fam file
PhenoMask Pheno::get_indi_mask(uint32_t ori_index){
    uint32_t pos_byte = ori_index % 4;
    switch(pos_byte){
        case 0:
            return PhenoMask(3, 0);
        case 1:
            return PhenoMask(12, 2);
        case 2:
            return PhenoMask(48, 4);
        default:
            return PhenoMask(192, 6);
    }
}

vector<uint32_t>& Pheno::get_index_keep() {
    return this->index_keep;
}

void Pheno::get_pheno(vector<string>& ids, vector<double>& pheno){
    ids.clear();
    ids.reserve(index_keep.size());
    pheno.clear();
    pheno.reserve(index_keep.size());
    for(auto& index : index_keep){
        if(std::isfinite(this->pheno[index])){
            ids.push_back(mark[index]);
            pheno.push_back(this->pheno[index]);
        }
    }
}


uint8_t Pheno::extract_genobit(uint8_t *const buf, int index_in_keep) {
    //return 3;
    uint32_t raw_index = index_keep[index_in_keep];
    //PhenoMask mask = get_indi_mask(raw_index);
    uint8_t mask, shift;
    switch(raw_index % 4){
        case 0:
            mask = 3; shift = 0; break;
        case 1:
            mask = 12; shift = 2; break;
        case 2:
            mask = 48; shift = 4; break;
        default:
            mask = 192; shift = 6;
    }
    uint8_t extract_byte = *(buf + raw_index / 4);
    return (extract_byte & mask) >> shift;
}

uint32_t Pheno::count_raw(){
    return num_ind;
}

uint32_t Pheno::count_keep(){
    return num_keep;
}

uint32_t Pheno::count_male(){
    return index_keep_male.size();
}

int8_t Pheno::get_sex(uint32_t index){
    return sex[index_keep[index]];
}

// remove have larger priority than keep, once the SNP has been removed, it
// will never be kept again
void Pheno::set_keep(vector<string>& indi_marks, vector<string>& marks, vector<uint32_t>& keeps, bool isKeep) {
    vector<uint32_t> keep_index, indi_index;
    vector_commonIndex_sorted1(marks, indi_marks, keep_index, indi_index);

    vector<uint32_t> remain_index;
    if(isKeep){
        std::set_intersection(keeps.begin(), keeps.end(),
                               keep_index.begin(), keep_index.end(),
                               std::back_inserter(remain_index));
    }else{
        std::set_difference(keeps.begin(), keeps.end(),
                            keep_index.begin(), keep_index.end(),
                            std::back_inserter(remain_index));
    }
    keeps = remain_index;

    LOGGER.i(0, string("After ") + (isKeep?"keeping":"removing") +  " individuals, " + to_string(keeps.size()) + " subjects remain.");
/* 
    std::sort(indi_marks.begin(), indi_marks.end());
    vector<uint32_t> pIN;
    uint32_t counter = 0;
    for(auto& item_mark : marks){
        if(std::binary_search(indi_marks.begin(), indi_marks.end(), item_mark)){
            pIN.push_back(counter);
        }
        counter++;
    }

    if(isKeep){
        auto common = [&pIN](int key) ->bool{
            return std::find(pIN.begin(), pIN.end(), key) == pIN.end();
        };
        keeps.erase(std::remove_if(keeps.begin(), keeps.end(), common), keeps.end());
    }else{
        auto diff = [&pIN](int key) ->bool{
            return std::find(pIN.begin(), pIN.end(), key) != pIN.end();
        };
        keeps.erase(std::remove_if(keeps.begin(), keeps.end(), diff), keeps.end());
    }
    */

}

void Pheno::update_sex(vector<string>& indi_marks, vector<double>& phenos){
    vector<uint32_t> pheno_index, update_index;
    vector_commonIndex_sorted1(mark, indi_marks, pheno_index, update_index);
    
    vector<uint32_t> common_index, pheno_index2;
    vector_commonIndex(index_keep, pheno_index, common_index, pheno_index2); 

    vector<uint32_t> indicies;
    indicies.reserve(common_index.size());
    for(int i = 0; i < common_index.size(); i++){
        uint32_t raw_index = index_keep[common_index[i]];
        int temp_update_pheno = std::round(phenos[update_index[pheno_index2[i]]]);
        if(temp_update_pheno > 0 && temp_update_pheno < 3){
            indicies.push_back(raw_index);
            sex[raw_index] = temp_update_pheno;
        }
    }
    index_keep = indicies;
}


void Pheno::update_pheno(vector<string>& indi_marks, vector<double>& phenos){
    vector<uint32_t> pheno_index, update_index;
    vector_commonIndex_sorted1(mark, indi_marks, pheno_index, update_index);
    
    vector<uint32_t> common_index, pheno_index2;
    vector_commonIndex(index_keep, pheno_index, common_index, pheno_index2); 

    vector<uint32_t> indicies;
    indicies.reserve(common_index.size());
    for(int i = 0; i < common_index.size(); i++){
        uint32_t raw_index = index_keep[common_index[i]];
        double temp_update_pheno = phenos[update_index[pheno_index2[i]]];
        if(!std::isnan(temp_update_pheno)){
            indicies.push_back(raw_index);
            pheno[raw_index] = temp_update_pheno;
        }
    }
    index_keep = indicies;

    /*
    for(int i = 0; i < pheno_index.size(); i++){
        double temp_update_pheno = phenos[update_index[i]];
        int raw_index = pheno_index[i];
        if(std::binary_search(index_keep.begin(), index_keep.end(), raw_index) &&
                (!std::isnan(temp_update_pheno))){
            pIN.push_back(raw_index);
            pheno[raw_index] = temp_update_pheno;
        }
    }

    index_keep = pIN;

    for(auto& index : index_keep){
        auto lower = std::lower_bound(indi_marks.begin(), indi_marks.end(), mark[index]);
        int lower_index = lower - indi_marks.begin();
        if(lower != indi_marks.end() && (*lower) == mark[index] && ! std::isnan(phenos[lower_index])){
            pIN.push_back(index);
            pheno[index] = phenos[lower_index];
        }
    }

    index_keep = pIN;
    */
}



void Pheno::reinit_rm(vector<uint32_t> &keeps, vector<uint32_t> &rms, int total_sample_number) {
    vector<uint32_t> whole_index(total_sample_number);
    std::iota(whole_index.begin(), whole_index.end(), 0);
    rms.clear();
    if(whole_index.size() == keeps.size()){
        return;
    }else{
        std::set_difference(whole_index.begin(), whole_index.end(), keeps.begin(), keeps.end(),
                std::inserter(rms, rms.begin()));
    }
}


void Pheno::init_bmask_block(){
    int max_block = (index_keep[index_keep.size()-1] + 31) / 32;
    block_num = max_block;
    int max_index = max_block * 32 - 1;

    int start_rm_index;
    if(index_rm.size() == 0){
        start_rm_index = index_keep[index_keep.size()-1] + 1;
    }else{
        start_rm_index = index_rm[index_rm.size() - 1] + 1;
    }

    for(int rm_index = start_rm_index; rm_index <= max_index; rm_index++){
        index_rm.push_back(rm_index);
    }

    int cur_block = 0;
    auto last_it_keep = index_keep.begin();
    auto last_it_rm = index_rm.begin();
    while(cur_block < max_block){
        uint32_t start_val = cur_block * 32;
        uint32_t end_val = start_val + 31;
        auto it_keep = std::upper_bound(last_it_keep, index_keep.end(), end_val); 
        bool keep_found = true;
        uint32_t cur_index = *(it_keep - 1);
        if(cur_index >= start_val && cur_index <= end_val){
            keep_block_index.push_back(cur_block);
            uint64_t mask_item = 0xFFFFFFFFFFFFFFFF;
            uint64_t mask_add_item = 0;
            auto begin_it_rm = std::upper_bound(last_it_rm, index_rm.end(), start_val);
            auto end_it_rm = std::upper_bound(begin_it_rm, index_rm.end(), end_val);
            if((*begin_it_rm) >= start_val && (*end_it_rm) <= end_val){
                for(auto it = begin_it_rm; it <= end_it_rm; it++){
                    uint32_t cur_offset = (*it) - start_val;
                    int cur_byte = cur_offset / 4;
                    uint8_t *mask_pos = (uint8_t *) ((&mask_item) + cur_byte);
                    uint8_t *mask_add_pos = (uint8_t *)((&mask_add_item) + cur_byte);
                    int cur_bitpos = cur_offset % 4;
                    uint8_t mask_piece;
                    uint8_t mask_add_piece;
                    switch(cur_bitpos){
                        case 0:
                            mask_piece = 0b11111100;
                            mask_add_piece = 0b00000001;
                            break;
                        case 1:
                            mask_piece = 0b11110011;
                            mask_add_piece = 0b00000100;
                            break;
                        case 2:
                            mask_piece = 0b11001111;
                            mask_add_piece = 0b00010000;
                            break;
                        case 3:
                            mask_piece = 0b00111111;
                            mask_add_piece = 0b01000000;
                    }
                    (*mask_pos) &= mask_piece;
                    (*mask_add_pos) += mask_add_piece;
                }
            }
            mask_items.push_back(mask_item);
            mask_add_items.push_back(mask_add_item);
            last_it_rm = end_it_rm;
        }
        last_it_keep = it_keep;
        cur_block++;
    }
}


void Pheno::init_mask_block(){
    mask_block.clear();
    mask_add_block.clear();
    for(const auto &index : index_rm){
        int byte_pos = index / 4;
        Mask_t::iterator cur_item = mask_block.lower_bound(byte_pos);
        Mask_t::iterator cur_add_item = mask_add_block.lower_bound(byte_pos);
        uint8_t mask_item = 0;
        uint8_t mask_add_item = 0;
        switch(index % 4){
            case 0:
                mask_item = 0b11111100;
                mask_add_item = 0b00000001;
                break;
            case 1:
                mask_item = 0b11110011;
                mask_add_item = 0b00000100;
                break;
            case 2:
                mask_item = 0b11001111;
                mask_add_item = 0b00010000;
                break;
            case 3:
                mask_item = 0b00111111;
                mask_add_item = 0b01000000;
        }

        if(cur_item != mask_block.end() && !(mask_block.key_comp()(byte_pos, cur_item->first))){
            cur_item->second = cur_item->second & mask_item;
            cur_add_item->second = cur_add_item->second + mask_add_item;
        }else{
            mask_block.insert(cur_item, Mask_t::value_type(byte_pos, mask_item));
            mask_add_block.insert(cur_add_item, Mask_t::value_type(byte_pos, mask_add_item));
        }

    }
}

void Pheno::mask_geno_keep(uint8_t *const geno_1block, int num_blocks) {
    if(mask_block.size() == 0){
        return;
    }
    #pragma omp parallel for schedule(dynamic) 
    for(int cur_block = 0; cur_block < num_blocks; cur_block++){
        uint8_t *cur_pos =  geno_1block + cur_block * num_bytes;
        for(const auto &item : mask_block){
            uint8_t *mask_pos = cur_pos + item.first;
            *mask_pos = ((*mask_pos) & item.second) + mask_add_block[item.first];
        }
    }
}

void Pheno::getMaskBit(uint64_t *maskp){
    for(auto keep_item : index_keep){
        uint32_t cur_qword = keep_item / 64;
        uint32_t cur_offset = keep_item % 64;
        maskp[cur_qword] |= k1LU << cur_offset;
    }
}

void Pheno::getMaskBitMale(uint64_t *maskp){
    for(auto keep_item : index_keep_male){
        uint32_t cur_qword = keep_item / 64;
        uint32_t cur_offset = keep_item % 64;
        maskp[cur_qword] |= k1LU << cur_offset;
    }
}



void Pheno::addOneFileOption(string key_store, string append_string, string key_name,
                                    map<string, vector<string>> options_in, map<string,string>& options) {
    if(options_in.find(key_name) != options_in.end()){
        if(options_in[key_name].size() == 1){
            options[key_store] = options_in[key_name][0] + append_string;
        }else if(options_in[key_name].size() > 1){
            options[key_store] = options_in[key_name][0] + append_string;
            LOGGER.w(0, "There are multiple " + key_name + ". Only the first one will be used in the analysis" );
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

int Pheno::registerOption(map<string, vector<string>>& options_in){
    addOneFileOption("pheno_file", ".fam", "--bfile", options_in, options);
    addOneFileOption("pheno_file", "", "--fam", options_in, options);
    options_in.erase("--fam");

    addOneFileOption("m_file", "", "--mbfile", options_in, options);
    if(options.find("m_file") != options.end()){
       vector<string> m_files;
        std::ifstream m_file(options["m_file"]);
        string line;
        vector<string> initial_ids;
        bool inited = false;
        while(std::getline(m_file,line)){
            vector<string> line_elements;
            boost::split(line_elements, line, boost::is_any_of("\t "));
            for(auto elements : line_elements){
                if(!elements.empty()){
                    vector<string> ids = read_sublist(elements+".fam");
                    m_files.push_back(elements);
                    if(!inited){
                        initial_ids = ids;
                        inited = true;
                        continue;
                    }
                    if(ids == initial_ids){
                        continue;
                    }else{
                        LOGGER.e(0, "fam is different in [" + elements + "].");
                    }
                }
            }
        }
        m_file.close();
        if(m_files.size() > 0){
            options["pheno_file"] = m_files[0] + ".fam";
            std::set<string> unique_mfiles;
            auto ends = std::remove_if(m_files.begin(), m_files.end(),
                    [&unique_mfiles](const string& item){
                        if(unique_mfiles.find(item) != unique_mfiles.end()){
                            return true;
                        }
                        unique_mfiles.insert(item);
                        return false;
                    });
            m_files.erase(ends, m_files.end());
            LOGGER.i(0, to_string(m_files.size()) + " genotype(s) in [" + options["m_file"] + "].");
            options_in.erase("--mbfile");
            options_in["m_file"] = m_files;
        }else{
            LOGGER.e(0, "no item in [" + options["m_file"] + "].");
        }
    }
  

    addOneFileOption("keep_file", "", "--keep", options_in, options);
    //options_in.erase("--keep");
    addOneFileOption("remove_file", "", "--remove", options_in,options);
    //options_in.erase("--remove"); // also may use in the GRM
    
    addOneFileOption("sex_file", "", "--update-sex", options_in, options);

    if(options_in.find("--pheno") != options_in.end()){
        addOneFileOption("qpheno_file", "", "--pheno", options_in, options);
        options_in.erase("--pheno");
    }

    if(options_in.find("--mpheno") != options_in.end()){
        if(options.find("qpheno_file") == options.end()){
            LOGGER.e(0, "--mpheno has to combine with --pheno");
        }
        options["mpheno"] = options_in["--mpheno"][0];
        options_in.erase("--mpheno");
    }

    if(options_in.find("--filter-sex") != options_in.end()){
        options["filter_sex"] = "yes";
    }


    // no main
    return 0;
}

void Pheno::processMain(){
    LOGGER.e(0, "Phenotype has no main process this time");
}
