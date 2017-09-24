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

using std::to_string;

map<string, string> Pheno::options;

Pheno::Pheno() {
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
        
    }

    reinit_rm(index_keep, index_rm, num_ind);
    num_keep = index_keep.size();
    num_rm = index_rm.size();

    init_mask_block();
}

// TODO filter the non-number strings other than nan
vector<string> Pheno::read_sublist(string sublist_file, vector<vector<double>> *phenos) {
    vector<string> subject_list;
    std::ifstream sublist(sublist_file.c_str());

    string line;
    int line_number = 0;
    int last_length = 0;
    bool init_pheno = true;
    while(std::getline(sublist, line)){
        line_number++;
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);
        int num_elements = line_elements.size();
        if(num_elements < 2){
            LOGGER.e(0, "the subject list file [" + sublist_file + "], line " + to_string(line_number) +
                        " has elements less than 2");
        }
        if(phenos && init_pheno){
            phenos->resize(num_elements - 2);
            init_pheno = false;
        }

        if(line_number > 1 && num_elements != last_length){
            string errmsg = "the subject list file [" + sublist_file + "], line " + to_string(line_number) +
                        " has different elements";
            if(phenos){
                LOGGER.e(0, errmsg);
            }else{
                LOGGER.w(0, errmsg);
            }
        }

        subject_list.push_back(line_elements[0] + "\t" + line_elements[1]);
        if(phenos){
            for(int index = 2; index != num_elements; index++){
                (*phenos)[index - 2].push_back(strtod(line_elements[index].c_str(), NULL));
            }
        }

        last_length = num_elements;
    }
    sublist.close();
    LOGGER.i(0, "Get " + to_string(line_number) + " subjects from [" + sublist_file + "]");
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
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);
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
    LOGGER.i(0, to_string(num_ind) + " individuals to be included from [" + fam_file + "].");
    fam.close();
}

vector<string> Pheno::get_id(int from_index, int to_index){
    vector<string> out_id;
    int raw_index;
    out_id.reserve(num_keep);
    for(int index = from_index; index <= to_index; index++){
        raw_index = index_keep[index];
        out_id.push_back(fid[raw_index] + "\t" + pid[raw_index]);
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
        ids.push_back(mark[index]);
        pheno.push_back(this->pheno[index]);
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

// remove have larger priority than keep, keep in mind, once the SNP has been removed, it
// will never be kept again
void Pheno::set_keep(vector<string>& indi_marks, vector<string>& marks, vector<uint32_t>& keeps, bool isKeep) {
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

    LOGGER.i(0, string("After ") + (isKeep?"keeping":"removing") +  " subjects, " + to_string(keeps.size()) + " subjects remained.");
}

void Pheno::update_pheno(vector<string>& indi_marks, vector<double>& phenos){
    std::sort(indi_marks.begin(), indi_marks.end());
    vector<uint32_t> pIN;
    for(auto& index : index_keep){
        auto lower = std::lower_bound(indi_marks.begin(), indi_marks.end(), mark[index]);
        int lower_index = lower - indi_marks.begin();
        if(lower != indi_marks.end() && (*lower) == mark[index] && ! std::isnan(phenos[lower_index])){
            pIN.push_back(index);
            pheno[index] = phenos[lower_index];
        }
    }

    index_keep = pIN;

    LOGGER.i(0, "After updating phenotypes, " + to_string(index_keep.size()) + " subjects remained.");
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
    uint8_t *cur_pos = NULL, *mask_pos = NULL;
    int start_byte;
    for(int cur_block = 0; cur_block != num_blocks; cur_block++){
        cur_pos =  geno_1block + cur_block * num_bytes;
        for(const auto &item : mask_block){
            mask_pos = cur_pos + item.first;
            *mask_pos = ((*mask_pos) & item.second) + mask_add_block[item.first];
        }
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
    addOneFileOption("keep_file", "", "--keep", options_in, options);
    //options_in.erase("--keep");
    addOneFileOption("remove_file", "", "--remove", options_in,options);
    //options_in.erase("--remove"); // also may use in the GRM

    if(options_in.find("--update-sex") != options_in.end()){
       // LOGGER.w(0, "--update-sex didn't work this time");
    }

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

    // no main
    return 0;
}

void Pheno::processMain(){
    LOGGER.e(0, "Phenotype has no main process this time");
}
