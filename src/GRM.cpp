/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: generate, read and process GRM.

   Depends on the class of genotype

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#include "GRM.h"
#include "Logger.h"
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <numeric>
#include <unordered_set>
#include "utils.hpp"
#include "ThreadPool.h"
#include "AsyncBuffer.hpp"
#include "utils.hpp"

using std::to_string;

map<string, string> GRM::options;
map<string, double> GRM::options_d;
vector<string> GRM::processFunctions;


uint64_t numByteGRMFile(FILE * file) {
    uint64_t f_size;
    if(file == NULL){
        LOGGER.e(0, "Read GRM file fail");
    }
    fseek(file, 0, SEEK_END);
    f_size = ftell(file);
    rewind(file);
    return f_size;
}

GRM::GRM(){
    if(options.find("grm_file") == options.end()){
        LOGGER.e(0, "You must specify --grm to invoke this function");
    }
    this->grm_file = options["grm_file"];
    grm_ids = Pheno::read_sublist(grm_file + ".grm.id");
    num_subjects = grm_ids.size();

    uint64_t num_grm = num_subjects * (num_subjects + 1) / 2;
    FILE *file = fopen((grm_file + ".grm.bin").c_str(), "rb");
    if(num_grm * 4 != numByteGRMFile(file)){
        LOGGER.e(0, "The GRM id and GRM binary is not matching [" + grm_file + "]");
    }
    fclose(file);

    uint64_t num_grm_byte = num_grm * 4;
    uint64_t num_parts = (num_grm_byte + num_byte_GRM_read - 1) / num_byte_GRM_read;
    vector<uint32_t> parts = divide_parts(0, num_subjects - 1, num_parts);
    index_grm_pairs.reserve(parts.size());
    byte_part_grms.reserve(parts.size());

    index_grm_pairs.push_back(std::make_pair(0, parts[0]));
    byte_part_grms.push_back((0 + parts[0] + 2) * (parts[0] + 1) / 2);
    for(int index = 1; index != parts.size(); index++){
        index_grm_pairs.push_back(std::make_pair(parts[index - 1] + 1, parts[index]));
        byte_part_grms.push_back(((uint64_t)parts[index - 1] + parts[index] + 3) * (parts[index] - parts[index - 1]) / 2);
    }

    num_byte_buffer = *std::max_element(byte_part_grms.begin(), byte_part_grms.end());

    index_keep.resize(num_subjects);
    std::iota(index_keep.begin(), index_keep.end(), 0);

    if(options.find("keep_file") != options.end()){
        vector<string> keep_lists = Pheno::read_sublist(options["keep_file"]);
        Pheno::set_keep(keep_lists, grm_ids, index_keep, true);
    }

    if(options.find("remove_file") != options.end()){
        vector<string> rm_lists = Pheno::read_sublist(options["remove_file"]);
        Pheno::set_keep(rm_lists, grm_ids, index_keep, false);
    }

}

void GRM::prune_FAM(float thresh){
    LOGGER.i(0, "Pruning the FAM with a cut off of " + to_string(thresh) + "...");
    LOGGER.i(0, "Total number of parts to proceed: " + to_string(index_grm_pairs.size()));
    FILE *grmFile = fopen((grm_file + ".grm.bin").c_str(), "rb");

    std::ofstream o_id((options["out"] + ".grm.id").c_str());
    std::ofstream o_fam((options["out"] + ".grm.fam").c_str());
    if(!o_id) LOGGER.e(0, "can't write to [" + options["out"] + ".grm.id]");
    if(!o_fam) LOGGER.e(0, "can't write to [" + options["out"] + ".grm.fam]");

    //Save the kept IDs, which may change by --keep and --remove
    vector<string> keep_ID;
    keep_ID.reserve(index_keep.size());
    for(auto & index : index_keep){
        keep_ID.emplace_back(grm_ids[index]);
    }
    LOGGER.i(2, "Saving " + to_string(keep_ID.size()) + " individual IDs");
    std::copy(keep_ID.begin(), keep_ID.end(), std::ostream_iterator<string>(o_id, "\n"));
    o_id.close();
    keep_ID.clear();
    keep_ID.shrink_to_fit();


    // Save pair1 par2 GRM
    std::unordered_set<int> keeps_ori(index_keep.begin(), index_keep.end());
    float *grm_buf = new float[num_byte_buffer];
    float cur_grm, *cur_grm_pos0;
    vector<float> rm_grm;
    vector<int> rm_grm_ID1, rm_grm_ID2;
    for(int part_index = 0; part_index != index_grm_pairs.size(); part_index++){
        if(fread(grm_buf, sizeof(float), byte_part_grms[part_index], grmFile) != byte_part_grms[part_index]){
            LOGGER.e(0, "Read GRM failed between line " + to_string(index_grm_pairs[part_index].first + 1) + " and "
                        + to_string(index_grm_pairs[part_index].second + 1));
        };
        if(part_index % 50 == 0){
            LOGGER.i(2, "Processing part " + to_string(part_index + 1));
        }
        cur_grm_pos0 = grm_buf;
        for(int id1 = index_grm_pairs[part_index].first; id1 != index_grm_pairs[part_index].second + 1; id1++){
            if (keeps_ori.find(id1) == keeps_ori.end()) {
                cur_grm_pos0 += id1 + 1;
                continue;
            }
            for(int id2 : index_keep){
                if(id2 > id1) break;
                cur_grm = *(cur_grm_pos0 + id2);
                if(cur_grm > thresh){
                    rm_grm_ID1.push_back(id1);
                    rm_grm_ID2.push_back(id2);
                    rm_grm.push_back(cur_grm);
                }
            }
            cur_grm_pos0 += id1 + 1;
        }
    }
    delete [] grm_buf;

    LOGGER.i(0, "Saving FAM (" + to_string(rm_grm.size()) + " pairs) to [" + options["out"] + ".grm.fam]");
    auto sorted_index = sort_indexes(rm_grm_ID2, rm_grm_ID1);
    for(auto index : sorted_index){
        o_fam << rm_grm_ID1[index] << "\t" << rm_grm_ID2[index] << "\t" << rm_grm[index] << std::endl;
    }
    o_fam.close();
    LOGGER.i(0, "Success:", "FAM pruning finished");

}



void GRM::cut_rel(float thresh, bool no_grm){
    LOGGER.i(0, "Pruning the GRM with a cutoff of " + to_string(thresh) + "...");
    LOGGER.i(0, "Total number of parts to proceed: " + to_string(index_grm_pairs.size()));
    FILE *grmFile = fopen((grm_file + ".grm.bin").c_str(), "rb");
    // put this first to avoid unwritable disk
    std::ofstream o_keep((options["out"] + ".grm.id").c_str());
    if((!o_keep)){
        LOGGER.e(0, "can't write to [" + options["out"] + ".grm.id]");
    }

    bool detail_flag = false;
    std::ofstream out_grm, o_remove;
    if(options.find("cutoff_detail") != options.end()){
        detail_flag = true;
        out_grm.open((options["out"] + ".rel.pair").c_str());
        o_remove.open((options["out"] + ".rel.id").c_str());
        if((!out_grm) || (!o_remove)){
            LOGGER.e(0, "can't write to [" + options["out"] + ".rel.pair, rel.id]");
        }
    }
 

    std::unordered_set<int> keeps_ori(index_keep.begin(), index_keep.end());

    float *grm_buf = new float[num_byte_buffer];
    float cur_grm, *cur_grm_pos0;
    vector<float> rm_grm;
    vector<int> rm_grm_ID1, rm_grm_ID2;
    for(int part_index = 0; part_index != index_grm_pairs.size(); part_index++){
        if(fread(grm_buf, sizeof(float), byte_part_grms[part_index], grmFile) != byte_part_grms[part_index]){
            LOGGER.e(0, "Read GRM failed between line " + to_string(index_grm_pairs[part_index].first + 1) + " and "
                        + to_string(index_grm_pairs[part_index].second + 1));
        };
        if(part_index % 50 == 0){
            LOGGER.i(2, "Processing part " + to_string(part_index + 1));
        }
        cur_grm_pos0 = grm_buf;
        for(int id1 = index_grm_pairs[part_index].first; id1 != index_grm_pairs[part_index].second + 1; id1++){
            if (keeps_ori.find(id1) == keeps_ori.end()) {
                cur_grm_pos0 += id1 + 1;
                continue;
            }
            for(int id2 : index_keep){
                if(id2 >= id1) break;
                cur_grm = *(cur_grm_pos0 + id2);
                if(cur_grm > thresh){
                    rm_grm_ID1.push_back(id1);
                    rm_grm_ID2.push_back(id2);
                    rm_grm.push_back(cur_grm);
                }
            }
            cur_grm_pos0 += id1 + 1;
        }
    }
    delete [] grm_buf;

    if(detail_flag){
        LOGGER.i(0, "Saving related pairs");
        for(int index = 0; index != rm_grm.size(); index++){
            out_grm << grm_ids[rm_grm_ID1[index]] << "\t" <<  grm_ids[rm_grm_ID2[index]] << "\t" << rm_grm[index] << std::endl;
        }
        out_grm.close();
        LOGGER.i(0, "Related sample pairs have been saved to " + options["out"] + ".rel.pair");
    }

    int i_buf;
    // copy from GCTA 1.26
    vector<int> rm_uni_ID(rm_grm_ID1);
    rm_uni_ID.insert(rm_uni_ID.end(), rm_grm_ID2.begin(), rm_grm_ID2.end());
    stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    std::map<int, int> rm_uni_ID_count;
    for (int i = 0; i < rm_uni_ID.size(); i++) {
        i_buf = count(rm_grm_ID1.begin(), rm_grm_ID1.end(), rm_uni_ID[i]) + count(rm_grm_ID2.begin(), rm_grm_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(std::pair<int, int>(rm_uni_ID[i], i_buf));
    }

    // swapping
    std::map<int, int>::iterator iter1, iter2;
    for (int i = 0; i < rm_grm_ID1.size(); i++) {
        iter1 = rm_uni_ID_count.find(rm_grm_ID1[i]);
        iter2 = rm_uni_ID_count.find(rm_grm_ID2[i]);
        if (iter1->second < iter2->second) {
            i_buf = rm_grm_ID1[i];
            rm_grm_ID1[i] = rm_grm_ID2[i];
            rm_grm_ID2[i] = i_buf;
        }
    }

    stable_sort(rm_grm_ID1.begin(), rm_grm_ID1.end());
    rm_grm_ID1.erase(unique(rm_grm_ID1.begin(), rm_grm_ID1.end()), rm_grm_ID1.end());
    vector<string> removed_ID;
    for (auto &index : rm_grm_ID1) removed_ID.push_back(grm_ids[index]);

    auto diff = [&rm_grm_ID1](int value) ->bool{
        return std::binary_search(rm_grm_ID1.begin(), rm_grm_ID1.end(), value);
    };
    index_keep.erase(std::remove_if(index_keep.begin(), index_keep.end(), diff), index_keep.end());

    vector<string> keep_ID;
    keep_ID.reserve(index_keep.size());
    for(auto & index : index_keep){
        keep_ID.emplace_back(grm_ids[index]);
    }

    LOGGER.i(0, "After pruning the GRM, there are " + to_string(removed_ID.size()) + " individuals (" + to_string(keep_ID.size()) + " individuals removed).");
    std::copy(keep_ID.begin(), keep_ID.end(), std::ostream_iterator<string>(o_keep, "\n"));
    o_keep.close();
    LOGGER.i(2, "Pruned unrelated IDs have been saved to " + options["out"] + ".grm.id");

    if(detail_flag){
        std::copy(removed_ID.begin(), removed_ID.end(), std::ostream_iterator<string>(o_remove, "\n"));
        o_remove.close();
        LOGGER.i(2, "IDs removed have been saved to " + options["out"] + ".rel.id");
    }

    if(no_grm) {
        fclose(grmFile);
        return;
    }

    LOGGER.i(0, "Pruning GRM values, total parts " + std::to_string(index_grm_pairs.size()));
    FILE *grm_out_file = fopen((options["out"] + ".grm.bin").c_str(), "wb");
    if(!grm_out_file){
        LOGGER.e(0, "can't open [" + options["out"] + ".grm.bin] for write");
    }
    fseek(grmFile, 0, SEEK_SET);
    outBinFile(grmFile, grm_out_file);
    fclose(grmFile);
    fclose(grm_out_file);
    LOGGER.i(2, "GRM values has been saved to [" + options["out"] + ".grm.bin]");

    LOGGER.i(0, "Pruning number of SNPs to calculate GRM, total parts " + std::to_string(index_grm_pairs.size()));
    FILE *NFile = fopen((grm_file + ".grm.N.bin").c_str(), "rb");
    if(!NFile){
        LOGGER.w(2, "There is no [" + grm_file + ".grm.N.bin]");
        return;
    }
    FILE *N_out_file = fopen((options["out"] + ".grm.N.bin").c_str(), "wb");
    if(!N_out_file){
        LOGGER.w(2, "can't open [" + options["out"] + ".grm.N.bin] to write, ignore this step");
        fclose(NFile);
        return;
    }
    outBinFile(NFile, N_out_file);
    fclose(NFile);
    fclose(N_out_file);
    LOGGER.i(2, "number of SNPs has been saved to [" + options["out"] + ".grm.N.bin]");
}

void GRM::outBinFile(FILE *sFile, FILE *dFile) {
    std::unordered_set<int> keeps(index_keep.begin(), index_keep.end());
    float *grm_buf = new float[num_byte_buffer];
    float *grm_out_buffer = new float[num_byte_buffer];
    float *cur_grm_pos, *cur_out_pos;
    for(int part_index = 0; part_index != index_grm_pairs.size(); part_index++){
        if(fread(grm_buf, sizeof(float), byte_part_grms[part_index], sFile) != byte_part_grms[part_index]){
            LOGGER.e(0, "Read failed between line " + to_string(index_grm_pairs[part_index].first + 1) + " and "
                        + to_string(index_grm_pairs[part_index].second + 1));
        };
        if(part_index % 50 == 0){
            LOGGER.i(2, "Processing part " + to_string(part_index + 1));
        }
        cur_grm_pos = grm_buf;
        cur_out_pos = grm_out_buffer;
        for(int id1 = index_grm_pairs[part_index].first; id1 != index_grm_pairs[part_index].second + 1; id1++) {
            if (keeps.find(id1) == keeps.end()) {
                cur_grm_pos += id1 + 1;
                continue;
            }
            for (int id2 : index_keep) {
                if(id2 > id1) break;
                *(cur_out_pos++) = *(cur_grm_pos + id2);
            }
            cur_grm_pos += id1 + 1;
        }
        uint64_t size_write = cur_out_pos - grm_out_buffer;
        if(fwrite(grm_out_buffer, sizeof(float), size_write, dFile) != size_write){
            LOGGER.e(0, "Write output file failed, please check the disk condition or permission");
        };
    }
    delete [] grm_buf;
    delete [] grm_out_buffer;
}

void GRM::init_AsyncBuffer(){
    if(asyncBuffer){
        delete asyncBuffer;
    }
    asyncBuffer = new AsyncBuffer<double>(num_byte_buffer);
}


GRM::GRM(Geno *geno) {
    //clock_t begin = t_begin();
    this->geno = geno;
    this->index_keep = geno->pheno->get_index_keep();
    this->p_index_keep = this->index_keep.data();
    this->part = std::stoi(options["cur_part"]);
    this->num_parts = std::stoi(options["num_parts"]);

    this->num_byte_geno = sizeof(uint16_t) * Constants::NUM_MARKER_READ / num_marker_block * index_keep.size();

    this->num_byte_cmask = sizeof(uint64_t) * Constants::NUM_MARKER_READ * index_keep.size() / num_cmask_block;

    lookup_GRM_table = new double[Constants::NUM_MARKER_READ / num_marker_block][num_lookup_table];

    int ret1 = posix_memalign((void **) &geno_buf, 32, num_byte_geno);
    int ret2 = posix_memalign((void **) &mask_buf, 32, num_byte_geno);
    int ret3 = posix_memalign((void **) &cmask_buf, 32, num_byte_cmask); 
    if(ret1 != 0 | ret2 != 0 | ret3 != 0){
        LOGGER.e(0, "can't allocate enough memory to genotype buffer.");
    }

    reg_bit_width = 64;

    this->num_grm_handle = reg_bit_width / 64;
    this->num_count_handle = reg_bit_width / 32;
    this->num_block_handle = reg_bit_width / 16;
    this->num_marker_process_block = this->num_marker_block * this->num_block_handle;

    vector<uint32_t> parts = divide_parts(0, this->geno->pheno->count_keep() - 1, num_parts);

    if (parts.size() != num_parts) {
        LOGGER.w(0, "can not divide into " + to_string(num_parts) + ". Use " + to_string(parts.size()) + " instead.");
        this->num_parts = parts.size();
    }

    if (part > this->num_parts) {
        LOGGER.e(0, "can not calculate" + to_string(part) + " larger than " + to_string(this->num_parts));
    } else {
        if (part == 1) {
            part_keep_indices = std::make_pair(0, parts[0]);
        } else {
            part_keep_indices = std::make_pair(parts[part - 2] + 1, parts[part - 1]);
        }
    }

    num_individual = part_keep_indices.second - part_keep_indices.first + 1;
    num_grm = ((uint64_t) part_keep_indices.first + part_keep_indices.second + 2) * num_individual / 2;

    uint64_t fill_grm = (num_grm + num_count_handle - 1) / num_count_handle * num_count_handle;
    int ret_grm = posix_memalign((void **)&grm, 32, fill_grm * sizeof(double));
    if(ret_grm){
        LOGGER.e(0, "can't allocate enough memory to (parted) GRM: " + to_string(fill_grm*sizeof(double) / 1024.0/1024/1024) + "GB");
    }
    memset(grm, 0, fill_grm * sizeof(double));

    int ret_N = posix_memalign((void **)&N, 32, fill_grm * sizeof(uint32_t));
    if(ret_N){
        LOGGER.e(0, "can't allocate enough memory to (parted) N: " + to_string(fill_grm*sizeof(uint32_t) / 1024.0/1024/1024) + "GB");
    }
    memset(N, 0, fill_grm * sizeof(uint32_t));

    sub_miss = new uint32_t[index_keep.size()]();

    //calculate each index in pair thread;
    int num_thread = THREADS.getThreadCount() + 1;
    index_grm_pairs.reserve(num_thread);
    vector<uint32_t> thread_parts = divide_parts(part_keep_indices.first, part_keep_indices.second, num_thread);
    if(num_thread != thread_parts.size()){
        LOGGER.w(0, "can not run in " + to_string(num_thread) + " threads. Use " + to_string(thread_parts.size()) + " instead");
    }

    index_grm_pairs.push_back(std::make_pair(part_keep_indices.first, thread_parts[0]));
    for(int index = 1; index != thread_parts.size(); index++){
        index_grm_pairs.push_back(std::make_pair(thread_parts[index - 1] + 1, thread_parts[index]));
    }

    //t_print(begin, "  INIT finished");

    LOGGER.i(0, "Start: ", "part " + to_string(part) + ", Subjects ID: " + to_string(part_keep_indices.first) + "-" + to_string(part_keep_indices.second) + ". Total parts " + to_string(num_parts) + ". " + to_string(thread_parts.size()) + " threads");
    LOGGER.i(2, to_string(num_individual) + " samples, " + to_string(geno->marker->count_extract()) + " markers, " + to_string(num_grm) + " GRM");

    o_name = options["out"];

    output_id();

#ifndef NDEBUG
    o_geno0 = fopen("./test.bin", "wb");
    o_mask0 = fopen("./test_mask.bin", "wb");
#endif
}


void GRM::output_id() {
    vector<string> out_id = geno->pheno->get_id(part_keep_indices.first, part_keep_indices.second);

    string o_grm_id = o_name + ".grm.id";
    std::ofstream grm_id(o_grm_id.c_str());

    if (!grm_id) { LOGGER.e(0, "can not open the file [" + o_grm_id + "] to write"); }
    std::copy(out_id.begin(), out_id.end(), std::ostream_iterator<string>(grm_id, "\n"));
    grm_id.close();
    LOGGER.i(0, "IDs for the GRM file has been saved in the file [" + o_grm_id + "]");

}

void GRM::calculate_GRM(uint8_t *buf, int num_marker) {
    //LOGGER.i(0, "GRM: ", "logged " + to_string(num_marker));
    clock_t begin;// = t_begin();

    //prepare the 7 freq arrays in current range;
    for (int index_marker = 0; index_marker != num_marker; index_marker++) {
        double af = geno->AFA1[finished_marker + index_marker];
        double mu = af * 2.0;
        double rdev = 1.0 / (mu * (1.0 - af));
        GRM_table[index_marker][0] = (2.0 - mu) * (2.0 - mu) * rdev;  // 00 00
        GRM_table[index_marker][2] = (2.0 - mu) * (1.0 - mu) * rdev;  // 00 10
        GRM_table[index_marker][3] = (2.0 - mu) * (0.0 - mu) * rdev;  // 00 11
        GRM_table[index_marker][4] = (1.0 - mu) * (1.0 - mu) * rdev;  // 10 10
        GRM_table[index_marker][5] = (0.0 - mu) * (1.0 - mu) * rdev;  // 11 10
        GRM_table[index_marker][6] = (0.0 - mu) * (0.0 - mu) * rdev;  // 11 11
        GRM_table[index_marker][7] = 0.0;                             // 01 01
    }

    //init the remained arrays to 0;
    for(int index_marker = num_marker; index_marker != Constants::NUM_MARKER_READ; index_marker++){
        for(int i = 0; i != 8; i++){
            GRM_table[index_marker][i] = 0.0;
        }
    }
    //t_print(begin, "  FREQ init");

    //begin = t_begin();

    int num_block = (num_marker + num_marker_block - 1) / num_marker_block;

    int length = elements.size();

    //init GRM lookup table
    for(int cur_block = 0; cur_block != num_block; cur_block++) {
        int cur_block_base = cur_block * num_marker_block;
        for (int e1 = 0; e1 != length; e1++) {
            for (int e2 = 0; e2 != length; e2++) {
                for (int e3 = 0; e3 != length; e3++) {
                    for (int e4 = 0; e4 != length; e4++) {
                        for (int e5 = 0; e5 != length; e5++) {
                            lookup_GRM_table[cur_block][elements[e5]
                                             | (elements[e4] << 3)
                                             | (elements[e3] << 6)
                                             | (elements[e2] << 9)
                                             | (elements[e1] << 12)] = GRM_table[cur_block_base][elements[e5]] +
                                                                       GRM_table[cur_block_base + 1][elements[e4]] +
                                                                       GRM_table[cur_block_base + 2][elements[e3]] +
                                                                       GRM_table[cur_block_base + 3][elements[e2]] +
                                                                       GRM_table[cur_block_base + 4][elements[e1]];
                        }
                    }
                }
            }
        }

    }

    //t_print(begin, "  FREQ LOOKUPtable");

   // begin = t_begin();

    memset(this->geno_buf, 0, num_byte_geno);
    memset(this->mask_buf, 0, num_byte_geno);
    memset(this->cmask_buf, 0, num_byte_cmask);

    int num_process_block = (num_marker + num_marker_process_block - 1) / num_marker_process_block;
    this->cur_num_block = num_process_block;
    int num_marker_in_process_block;

    uint8_t *cur_buf;
    uint16_t *cur_geno;
    uint16_t temp_geno;
    uint16_t *cur_mask;
    uint32_t raw_index;
    uint16_t * p_geno;
    uint16_t * p_mask;
    int cur_block, cur_block_pos, indi_block;
    int base_marker_index, base_buf_pos;
    for(int index_process_block = 0; index_process_block != num_process_block; index_process_block++){
        base_marker_index = index_process_block * num_marker_process_block;
        cur_buf = buf + base_marker_index * geno->num_byte_per_marker;
        base_buf_pos = num_block_handle * index_process_block * index_keep.size();
        cur_geno = this->geno_buf + base_buf_pos;
        cur_mask = this->mask_buf + base_buf_pos;
        num_marker_in_process_block = index_process_block == num_process_block - 1 ?
                                      (num_marker - base_marker_index): num_marker_process_block;
        // recode genotype
        for(int index_indi = 0; index_indi != part_keep_indices.second + 1; index_indi += 1) {
            raw_index = index_keep[index_indi];
            indi_block = index_indi * num_block_handle;

            p_geno = cur_geno + indi_block;
            p_mask = cur_mask + indi_block;

            for (int marker_index = 0; marker_index != num_marker_in_process_block; marker_index++) {
                cur_block = marker_index / num_marker_block;
                cur_block_pos = (marker_index % num_marker_block) * 3;
                temp_geno = (*(cur_buf + marker_index * geno->num_byte_per_marker + raw_index / 4)
                        >> ((raw_index % 4) * 2)) & 3LU;
                *(p_geno + cur_block) |= (temp_geno << cur_block_pos);
                if (temp_geno == 1) {
                    *(p_mask + cur_block) |= (7LU << cur_block_pos);
                    sub_miss[index_indi] += 1;
                    //cmask
                    cmask_buf[index_indi + ((base_marker_index + marker_index) / num_cmask_block)*index_keep.size()] |=
                            (1LLU << ((base_marker_index + marker_index) % num_cmask_block));
                }
            }

            for (int marker_index = num_marker_in_process_block;
                 marker_index != num_marker_process_block; marker_index++) {
                cur_block = marker_index / num_marker_block;
                cur_block_pos = (marker_index % num_marker_block) * 3;
                *(p_mask + cur_block) |= (7LU << cur_block_pos);

            }
        }
    }
 //   t_print(begin, "  init finished");

    //submit thread
    //begin = t_begin();
    int pair_size = index_grm_pairs.size();
    for(int index = 0; index != pair_size - 1; index++){
        THREADS.AddJob(std::bind(&GRM::grm_thread, this, index_grm_pairs[index].first, index_grm_pairs[index].second));
    }

    for(int index = 0; index != pair_size - 1; index++){
        THREADS.AddJob(std::bind(&GRM::N_thread, this, index_grm_pairs[index].first, index_grm_pairs[index].second));
    }

    auto last_pair = index_grm_pairs[pair_size - 1];
    grm_thread(last_pair.first, last_pair.second);
    N_thread(last_pair.first, last_pair.second);

    THREADS.WaitAll();

    finished_marker += num_marker;
    //std::stringstream out_message;
    //out_message << std::fixed << std::setprecision(2) << finished_marker * 100.0 / geno->marker->count_extract();
    //LOGGER.i(0, out_message.str() + "% has been finished");

    if(finished_marker % 30000 == 0){
        LOGGER.i(4, to_string(finished_marker) + " SNPs finished");
    }
}

void GRM::deduce_GRM(){
    LOGGER.i(0, "Success", "GRM calculations finished");
    //Just for test
#ifndef NDEBUG
    fclose(o_geno0);
    fclose(o_mask0);
#endif

    //clock_t begin = t_begin();
    string grm_name = o_name + ".grm.bin";
    string N_name = o_name + ".grm.N.bin";
    FILE *grm_out = fopen(grm_name.c_str(), "wb");
    FILE *N_out = fopen(N_name.c_str(), "wb");
    if((!grm_out) || (!N_out)){
        LOGGER.e(0, "can't open " + o_name + ".grm.bin or .grm.N.bin to write");
    }

    uint32_t num_sample = index_keep.size();
    float *w_grm = new float[num_sample];
    float *w_N = new float[num_sample];

    double *po_grm = grm;
    uint32_t *po_N = N;
    uint32_t sub_miss1;
    uint32_t sub_N;
    for(int pair1 = part_keep_indices.first; pair1 != part_keep_indices.second + 1; pair1++){
        sub_miss1 = finished_marker - sub_miss[pair1];
        for(int pair2 = 0; pair2 != pair1 + 1; pair2++){
            sub_N = *po_N + sub_miss1 - sub_miss[pair2];
            w_N[pair2] = (float)sub_N;

            if(sub_N){
                w_grm[pair2] = (float)((*po_grm)/sub_N);
            }else{
                w_grm[pair2] = 0.0;
            }
            po_N++;
            po_grm++;
        }
        fwrite(w_grm, sizeof(float), pair1 + 1, grm_out);
        fwrite(w_N, sizeof(float), pair1 + 1, N_out);

    }
    fclose(grm_out);
    fclose(N_out);
    delete[] w_grm;
    delete[] w_N;
    //t_print(begin, "  GRM deduce finished");
    LOGGER.i(0, "GRM has been saved in the file [" + o_name + ".grm.bin]");
    LOGGER.i(0, "Number of SNPs in each pair of individuals has been saved in the file [" + o_name + ".grm.N.bin]");
}

void GRM::N_thread(int grm_index_from, int grm_index_to){
    //LOGGER.i(0, "Nthread: ", to_string(grm_index_from) + " " + to_string(grm_index_to));
    uint32_t *po_N;
    uint64_t cmask1, cmask2, cmask;
    uint64_t *p_cmask1, *p_cmask2;
    uint64_t *cur_cmask;
    uint32_t *po_N_start = N + (grm_index_from + 1 + part_keep_indices.first) * (grm_index_from - part_keep_indices.first) / 2;
    for(int cur_block = 0; cur_block != Constants::NUM_MARKER_READ / num_cmask_block; cur_block++){
        cur_cmask = cmask_buf + cur_block * index_keep.size();
        po_N = po_N_start;
        p_cmask1 = cur_cmask + grm_index_from;
        for(int index_pair1 = grm_index_from; index_pair1 != grm_index_to + 1; index_pair1++){
            p_cmask2 = cur_cmask;
            cmask1 = *p_cmask1++;
            if(cmask1){
                for(int index_pair2 = 0; index_pair2 != index_pair1 + 1; index_pair2++){
                    cmask2 = *p_cmask2++;
                    if(cmask2){
                        cmask = cmask1 & cmask2;
                        if(cmask){
                            *po_N += popcount(cmask);
                        }
                    }
                    po_N++;
                }
            }else{
                po_N += index_pair1 + 1;
            }
        }
    }
}

void GRM::grm_thread(int grm_index_from, int grm_index_to) {
    //LOGGER.i(0, "GRMthread: ", to_string(grm_index_from) + " " + to_string(grm_index_to));

    double *po_grm;
    uint32_t *po_N;
    uint64_t geno, mask;
    uint64_t *p_geno1, *p_geno2, geno1;
    uint64_t *p_mask1, *p_mask2, mask1;
    uint64_t *cur_geno, *cur_mask;
    uint16_t geno_p1, geno_p2, geno_p3, geno_p4;
    double *table0, *table1, *table2, *table3;
    int cur_geno_block, start_geno_block;
    int grm_pos_offset = (grm_index_from + 1 + part_keep_indices.first) * (grm_index_from - part_keep_indices.first) / 2;
    uint32_t *po_N_start = N + grm_pos_offset;
    double *po_grm_start = grm + grm_pos_offset;
    for (int cur_block = 0; cur_block != cur_num_block; cur_block++) {
        cur_geno_block = cur_block * num_block_handle;
        start_geno_block = cur_geno_block * index_keep.size();

        table0 = lookup_GRM_table[cur_geno_block];
        table1 = lookup_GRM_table[cur_geno_block + 1];
        table2 = lookup_GRM_table[cur_geno_block + 2];
        table3 = lookup_GRM_table[cur_geno_block + 3];

        cur_geno = (uint64_t *) (geno_buf + start_geno_block);
        cur_mask = (uint64_t *) (mask_buf + start_geno_block);

        po_grm = po_grm_start;
        po_N = po_N_start;

        p_geno1 = cur_geno + grm_index_from;
        p_mask1 = cur_mask + grm_index_from;

        for (int index_pair1 = grm_index_from; index_pair1 != grm_index_to + 1; index_pair1++) {
            geno1 = *(p_geno1);
            mask1 = *(p_mask1);
            //Just for test
            #ifndef NDEBUG
            if(index_pair1 == 0){
                fwrite(p_geno1, sizeof(uint64_t), 1, o_geno0);
                fwrite(p_mask1, sizeof(uint64_t), 1, o_mask0);
            }
            #endif
            p_geno2 = cur_geno;
            p_mask2 = cur_mask;
            if(mask1) {
                for (int index_pair2 = 0; index_pair2 != index_pair1 + 1; index_pair2++) {
                    geno = geno1 + *(p_geno2);
                    mask = mask1 | *(p_mask2);

                    geno |= mask;

                    geno_p1 = (uint16_t) geno;
                    geno_p2 = (uint16_t) (geno >> 16);
                    geno_p3 = (uint16_t) (geno >> 32);
                    geno_p4 = (uint16_t) (geno >> 48);

                    *po_grm += table0[geno_p1] +
                               table1[geno_p2] +
                               table2[geno_p3] +
                               table3[geno_p4];


                    p_geno2++;
                    p_mask2++;
                    po_grm++;
                    po_N++;
                }
            }else{
                for (int index_pair2 = 0; index_pair2 != index_pair1 + 1; index_pair2++) {
                    geno = geno1 + *(p_geno2);

                    geno |= *(p_mask2);

                    geno_p1 = (uint16_t) geno;
                    geno_p2 = (uint16_t) (geno >> 16);
                    geno_p3 = (uint16_t) (geno >> 32);
                    geno_p4 = (uint16_t) (geno >> 48);

                    *po_grm += table0[geno_p1] +
                               table1[geno_p2] +
                               table2[geno_p3] +
                               table3[geno_p4];

                    p_geno2++;
                    p_mask2++;
                    po_grm++;
                    po_N++;
                }
            }

            p_geno1++;
            p_mask1++;
        }

    }
}


vector<uint32_t> GRM::divide_parts(uint32_t from, uint32_t to, uint32_t num_parts){
    vector<uint64_t> num_indi_grms;
    num_indi_grms.reserve(to - from + 1);
    for(uint64_t index = from; index != to + 1; index++){
        num_indi_grms.push_back((index + 2 + from) * (index + 1 -  from) / 2);
    }

    uint64_t indi_parts = num_indi_grms[to - from] / num_parts;
    vector<uint32_t> parts;
    parts.reserve(num_parts);
    vector<uint64_t>::iterator upper;
    for(uint32_t index = 1; index <= num_parts; index++){
        upper = std::upper_bound(num_indi_grms.begin(), num_indi_grms.end(), indi_parts * index);
        if(upper != num_indi_grms.begin()){
            upper--;
            if(*upper < indi_parts * index){
                upper++;
            }
        }
        parts.push_back(upper - num_indi_grms.begin() + from );
    }
    parts.erase(unique(parts.begin(), parts.end()), parts.end());
    return parts;
}

int GRM::registerOption(map<string, vector<string>>& options_in) {
    int return_value = 0;
    options["out"] = options_in["out"][0];

    if(options_in.find("--grm") != options_in.end()){
        if(options_in["--grm"].size() == 1){
            options["grm_file"] = options_in["--grm"][0];
        }else{
            LOGGER.e(0, "can't handle multiple GRM files");
        }
        options_in.erase("--grm");
        if(options["grm_file"] == options["out"]){
            LOGGER.e(0, "can't set the name of input same with output name");
        }
    }

    Pheno::addOneFileOption("keep_file", "", "--keep", options_in, options);
    Pheno::addOneFileOption("remove_file", "", "--remove", options_in, options);

    int num_parts = 1;
    int cur_part = 1;

    string part_grm = "--make-grm-part";
    if(options_in.find(part_grm) != options_in.end()){
        if(options_in[part_grm].size() == 2){
            try{
                num_parts = std::stoi(options_in[part_grm][0]);
                cur_part = std::stoi(options_in[part_grm][1]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, part_grm + " can only deal with integer value");
            }
            if(num_parts <= 0 || cur_part <=0){
                LOGGER.e(0, part_grm + "arguments should >= 1");
            }
            if(num_parts < cur_part){
                LOGGER.e(0, part_grm + "1st parameter (number of parts) can't less than 2nd parameter");
            }

            std::string s_parts = std::to_string(num_parts);
            std::string c_parts = std::to_string(cur_part);
            options["out"] = options["out"] + ".part_" + s_parts + "_" + std::string(s_parts.length() - c_parts.length(), '0') + c_parts;
            options_in["out"][0] = options["out"];
            processFunctions.push_back("make_grm");
            options_in.erase(part_grm);
            std::map<string, vector<string>> t_option;
            t_option["--autosome"] = {};
            Marker::registerOption(t_option);
 
            return_value++;
        }else{
            LOGGER.e(0, part_grm + " takes two arguments, total parts and part to calculate currently");
        }
    }

    options["num_parts"] = std::to_string(num_parts);
    options["cur_part"] = std::to_string(cur_part);


    if(options_in.find("--make-grm") != options_in.end()){
        if(options.find("grm_file") == options.end()){
            processFunctions.push_back("make_grm");
            options_in.erase("--make-grm");
            std::map<string, vector<string>> t_option;
            t_option["--autosome"] = {};
            Marker::registerOption(t_option);
            return_value++;
        }else{
            if(options_in.find("--grm-cutoff") != options_in.end()){
                if(options_in["--grm-cutoff"].size() == 1){
                    options_d["grm_cutoff"] = std::stod(options_in["--grm-cutoff"][0]);
                }else{
                    LOGGER.e(0, "--grm-cutoff can't deal with more than one value currently");
                }
                if(options.find("grm_file") == options.end()){
                    LOGGER.e(0, "can't find --grm flag that is essential to --grm-cutoff");
                }

                processFunctions.push_back("grm_cutoff");
                options_in.erase("--grm-cutoff");

                return_value++;
            }

            if(options_in.find("--grm-no-relative") != options_in.end()){
                if(options_in["--grm-no-relative"].size() == 1){
                    options_d["grm_cutoff"] = std::stod(options_in["--grm-no-relative"][0]);
                }else{
                    LOGGER.e(0, "--grm-no-relative can't deal with more than one value currently");
                }
                if(options.find("grm_file") == options.end()){
                    LOGGER.e(0, "can't find --grm flag that is essential to --grm-no-relative");
                }

                processFunctions.push_back("grm_cutoff");
                options_in.erase("--grm-no-relative");
                options["no_grm"] = "true";

                return_value++;
            }


            if(options_in.find("--cutoff-detail") != options_in.end()){
                options["cutoff_detail"] = "true";
            }
        }

    }

    if(options_in.find("--make-fam") != options_in.end()){
        if(options.find("grm_file") == options.end()){
            LOGGER.e(0, "can't find --grm flag that is essential to --make-fam");
        }
        if(options_in["--make-fam"].size() == 1){
            options_d["grm_cutoff"] = std::stod(options_in["--make-fam"][0]);
            processFunctions.push_back("make_fam");
        }else{
            LOGGER.e(0, "--make-fam can't deal with more than one value currently");
        }
        return_value++;
    }

    return return_value;

}

void GRM::processMain() {
    vector<function<void (uint8_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "make_grm"){
            LOGGER.i(0, "Note: GRM is calculated based on autosome or the chromosome specified by --chr");
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            GRM grm(&geno);
            callBacks.push_back(bind(&Geno::freq, &geno, _1, _2));
            callBacks.push_back(bind(&GRM::calculate_GRM, &grm, _1, _2));
            geno.loop_block(callBacks);
            grm.deduce_GRM();
            return;
        }

        if(process_function == "grm_cutoff"){
            GRM grm;
            bool no_grm = false;
            if(options.find("no_grm") != options.end()) {
                no_grm = true;
            }
            grm.cut_rel(options_d["grm_cutoff"], no_grm);
           // grm.loop_block(callBacks);

        }

        if(process_function == "make_fam"){
            GRM grm;
            grm.prune_FAM(options_d["grm_cutoff"]);
        }
    }

}
