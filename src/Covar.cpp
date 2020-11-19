/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: covar information

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
#include "OptionIO.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <numeric>
#include "utils.hpp"
#include "StatLib.h"

map<string, string> Covar::options;
using std::to_string;

bool common_3vector(vector<string> **idListP, vector<vector<double>> **valueP,
        vector<string> *target_id, vector<vector<double>> **targetP){

    // can't be (*idlist)[0], this get the first element 
    vector<string> &s1 = *(idListP[0]);
    vector<string> &s2 = *(idListP[1]);
    vector<string> &s3 = *(idListP[2]);
    if(s1.size() == 0 || s2.size() == 0){
        return false;
    }

    vector<int> s1_index, s2_index;
    vector_commonIndex_sorted1(s1, s2, s1_index, s2_index);
    int remain_size = s1_index.size();
    if(remain_size == 0){
        return false;
    }
    vector<string> sample_id;
    sample_id.reserve(remain_size);
    std::transform(s1_index.begin(), s1_index.end(), 
            std::back_inserter(sample_id), [&s1](int pos){return s1[pos];}); 

    vector<int> s3_index;
    if(s3.size() != 0){
        vector<int> sa_index;
        vector_commonIndex_sorted1(sample_id, s3, sa_index, s3_index);
        remain_size = sa_index.size();
        if(remain_size == 0){
            return false;
        }
        vector<int> temp_s1, temp_s2;
        temp_s1.reserve(remain_size);
        temp_s2.reserve(remain_size);
        std::transform(sa_index.begin(), sa_index.end(), 
                std::back_inserter(temp_s1), [&s1_index](int pos){return s1_index[pos];});
        s1_index.resize(remain_size);
        s1_index = temp_s1;

        std::transform(sa_index.begin(), sa_index.end(), 
                std::back_inserter(temp_s2), [&s2_index](int pos){return s2_index[pos];});
        s2_index.resize(remain_size);
        s2_index = temp_s2;
    }

    target_id->resize(remain_size);
    std::transform(s1_index.begin(), s1_index.end(), target_id->begin(), 
            [&s1](int pos){return s1[pos];});

    auto &value1 = *valueP[0];
    auto &value2 = *valueP[1];
    auto &value3 = *valueP[2];
    auto &targetvalue1 = *targetP[0];
    auto &targetvalue2 = *targetP[1];
    auto &targetvalue3 = *targetP[2];

    vector<vector<double>> *values[] = {&value1, &value2, &value3};
    vector<int> *indicies[] = {&s1_index, &s2_index, &s3_index};
    vector<vector<double>> *targetvalues[] = {&targetvalue1, &targetvalue2, &targetvalue3};

    for(int i = 0; i < 3; i++){
        auto &temp_value = *(values[i]);
        auto &temp_index = *(indicies[i]);
        auto &temp_target = *(targetvalues[i]);
        vector<double> t_value(temp_index.size());
        for(int j = 0; j < temp_value.size(); j++){
            auto &pvalue = temp_value[j];
            std::transform(temp_index.begin(), temp_index.end(), 
                    t_value.begin(), [&pvalue](int pos){return pvalue[pos];}); 
            temp_target[j] = t_value;
        }
    }

    return true;
}



Covar::Covar(){
    vector<string> samples_qcovar, samples_covar, samples_rcovar;
    vector<vector<double>> qcovar, covar, rcovar;

    if(options.find("qcovar") != options.end()){
        string filename = options["qcovar"];
        LOGGER.i(0, "Reading quantitative covariates from [" + filename + "].");
        read_covar(filename, samples_qcovar, &qcovar, NULL, NULL);  
        if(hasVectorDuplicate(samples_qcovar)){
            LOGGER.e(0, "covariates can't have duplicate FID+IID.");
        }
        LOGGER.i(0, to_string(qcovar.size()) + " covariates of " + to_string(samples_qcovar.size()) + " samples to be included.");
    }

    if(options.find("covar") != options.end()){
        string filename = options["covar"];
        LOGGER.i(0, "Reading discrete covariates from [" + filename + "].");
        read_covar(filename, samples_covar, &covar, &labels_covar, NULL);  
        if(hasVectorDuplicate(samples_covar)){
            LOGGER.e(0, "covariates can't have duplicate FID+IID.");
        }
        LOGGER.i(0, to_string(covar.size()) + " covariates of " + to_string(samples_covar.size()) + " samples to be included.");
 
    }

    if(options.find("rcovar") != options.end()){
        string filename = options["rcovar"];
        LOGGER.i(0, "Reading ranked covariates from [" + filename + "].");
        read_covar(filename, samples_rcovar, &rcovar, &labels_rcovar, NULL);  
        if(hasVectorDuplicate(samples_rcovar)){
            LOGGER.e(0, "covariates can't have duplicate FID+IID.");
        }
        LOGGER.i(0, to_string(rcovar.size()) +  " covariates of " + to_string(samples_rcovar.size()) + " samples to be included.");
    }

    this->rcovar.resize(rcovar.size());
    this->covar.resize(covar.size());
    this->qcovar.resize(qcovar.size());

    int flags_samples = (samples_qcovar.size() != 0) + (samples_covar.size() != 0) * 10 + 
        (samples_rcovar.size() != 0) * 100;


    vector<string> *idList[3];
    vector<vector<double>> *values[3];
    vector<vector<double>> *targets[3];
    
    switch(flags_samples){
        case 0:
            return;
        case 1:
            sample_id = samples_qcovar;
            this->qcovar = qcovar;
            return;
        case 10:
            sample_id = samples_covar;
            this->covar = covar;
            setCovarMapping(false);
            return;
        case 100:
            sample_id = samples_rcovar;
            this->rcovar = rcovar;
            setCovarMapping(true);
            return;
        case 11:
            idList[0] = &samples_qcovar;
            idList[1] = &samples_covar;
            idList[2] = &samples_rcovar;
            values[0] = &qcovar;
            values[1] = &covar;
            values[2] = &rcovar;
            targets[0] = &(this->qcovar);
            targets[1] = &(this->covar);
            targets[2] = &(this->rcovar);
            break;
        case 101:
            idList[0] = &samples_qcovar;
            idList[2] = &samples_covar;
            idList[1] = &samples_rcovar;
            values[0] = &qcovar;
            values[2] = &covar;
            values[1] = &rcovar;
            targets[0] = &(this->qcovar);
            targets[2] = &(this->covar);
            targets[1] = &(this->rcovar);
            break;
        case 110:
            idList[2] = &samples_qcovar;
            idList[1] = &samples_covar;
            idList[0] = &samples_rcovar;
            values[2] = &qcovar;
            values[1] = &covar;
            values[0] = &rcovar;
            targets[2] = &(this->qcovar);
            targets[1] = &(this->covar);
            targets[0] = &(this->rcovar);
            break;
        case 111:
            idList[0] = &samples_qcovar;
            idList[1] = &samples_covar;
            idList[2] = &samples_rcovar;
            values[0] = &qcovar;
            values[1] = &covar;
            values[2] = &rcovar;
            targets[0] = &(this->qcovar);
            targets[1] = &(this->covar);
            targets[2] = &(this->rcovar);
            break;
        default:
            LOGGER.e(0, "something impossile happened.");
    }
    if(common_3vector(idList, values, &sample_id, targets)){
        LOGGER.i(0, to_string(qcovar.size()) + " qcovar, " + to_string(covar.size()) + " covar and " + to_string(rcovar.size()) + " rcovar to be included.");
        setCovarMapping(false);
        setCovarMapping(true);
    }else{
        LOGGER.e(0, "0 covartiate to be included.");
    }
    LOGGER.i(0, to_string(sample_id.size()) + " common samples in covariates to be included.");
}

const vector<string>& Covar::getSampleID() const{
    return this->sample_id;
}

bool Covar::setCovarMapping(bool is_rcovar, vector<vector<string>> *map_order){
    string err_string;
    if(is_rcovar){
        err_string = "can't specify inconsistent or single ranked covariate labels.";
    }else{
        err_string = "can't specify inconsistent or single covariate labels.";
    }
    if(sample_id.size() == 0 || covar.size() == 0){
        return false;
    }
    if(map_order){
        if(map_order->size() != covar.size()){
            LOGGER.e(0, err_string);
            return false;
        }
    }

    int expand_col_covar = 0;
    vector<int> start_col_X;
    start_col_X.push_back(0);

    auto *plabels_covar = &this->labels_covar;
    auto *pcovar = &this->covar;
    auto *plabels_covar_mapping = &this->labels_covar_mapping;
    if(is_rcovar){
        plabels_covar = &this->labels_rcovar;
        pcovar = &this->rcovar;
        plabels_covar_mapping = &this->labels_rcovar_mapping;
    }

    auto &labels_covar = *plabels_covar;
    auto &covar = *pcovar;
    auto &labels_covar_mapping = *plabels_covar_mapping;

    for(auto & label : labels_covar){
        expand_col_covar += label["LABEL_MAX_VALUE"];
        start_col_X.push_back(expand_col_covar);
    }

    labels_covar_mapping.resize(expand_col_covar);

    for(int i = 0; i < covar.size(); i++){
        auto &map_item = labels_covar[i];
        vector<string> elements;
        for(auto const& t_map : map_item){
            const string &key = t_map.first;
            if(key != "LABEL_MAX_VALUE"){
                elements.push_back(key);
            }
        }
        std::sort(elements.begin(), elements.end());
        if(elements.size() == 1){
            LOGGER.e(0, "column " + to_string(i) + ", " + err_string);
            return false;
        }
        if(map_order){
            vector<string> &pelements = (*map_order)[i];
            vector<string> mapper;
            for(auto const &element : pelements){
                if(std::binary_search(elements.begin(), elements.end(), element)){
                    mapper.push_back(element);
                }
            }
            if(elements.size() != mapper.size()){
                LOGGER.e(0, "column " + to_string(i) + ", " + err_string);
                return false;
            }
            elements = mapper;
        }

        if(!is_rcovar){
            for(int j = 1; j < elements.size(); j++){
                vector<double> item(elements.size(), 0);
                item[map_item[elements[j]]] = 1.0;
                labels_covar_mapping[start_col_X[i] + j - 1] = item;
            }
        }else{
            int n = elements.size();
            double *Z = new double[n*n];
            if(StatLib::rankContrast(n, Z)){
                for(int j = 1; j < n; j++){
                    vector<double> item(n);
                    int base_index_z = j * n;
                    std::transform(elements.begin(), elements.end(), item.begin(),
                            [&map_item, Z, &base_index_z](string t){return Z[base_index_z + map_item[t]];});

                    labels_covar_mapping[start_col_X[i] + j - 1] = item;
                }
            }else{
                delete[] Z;
                return false;
            }
            delete[] Z;
        }
    }


    return true;
}

bool Covar::getCommonSampleIndex(const vector<string> &sampleIDs, vector<uint32_t> &keep_index, vector<uint32_t> &covar_index){
    if(sampleIDs.size() == 0 || (!hasCovar())){
        return false;
    }
    vector_commonIndex_sorted1(sampleIDs, sample_id, keep_index, covar_index);
    if(keep_index.size() == 0){
        return false;
    }
    return true;
}

bool Covar::hasCovar(){
    int total_col_covar = qcovar.size() + covar.size() + rcovar.size();
    if(sample_id.size() == 0 || total_col_covar == 0){
        return false;
    }else{
        return true;
    }
}


bool Covar::getCovarX(const vector<string> &sampleIDs, vector<double> &X, vector<uint32_t> &keep_index){
    vector<uint32_t> covar_index;
    if(!getCommonSampleIndex(sampleIDs, keep_index, covar_index)){
        return false;
    }

    int common_sample_size = covar_index.size();

    vector<int> start_col_X;
    start_col_X.push_back(0);
    int expand_col_covar = 0;
    
    for(int i = 0; i < qcovar.size(); i++){
        expand_col_covar++;
        start_col_X.push_back(expand_col_covar);
    }

    for(auto & label : labels_covar){
        expand_col_covar += label["LABEL_MAX_VALUE"];
        start_col_X.push_back(expand_col_covar);
    }

    for(auto & label : labels_rcovar){
        expand_col_covar += label["LABEL_MAX_VALUE"];
        start_col_X.push_back(expand_col_covar);
    }

    X.resize(expand_col_covar * common_sample_size);

    for(int i = 0; i < qcovar.size(); i++){
        auto &covarp = this->qcovar[i];
        int base_pos = i * common_sample_size;
        std::transform(covar_index.begin(), covar_index.end(), X.begin() + base_pos,
                [&covarp](int pos){return covarp[pos];});
    }

    int base_covar = qcovar.size();
    int map_index = 0;
    for(int i = 0; i < covar.size(); i++){
        auto &covarp = this->covar[i];
        auto &labelp = this->labels_covar[i];
        int base_pos = (start_col_X[base_covar + i]) * common_sample_size;
        for(int j = 0; j < labelp["LABEL_MAX_VALUE"]; j++){
            auto &cur_table = labels_covar_mapping[map_index];
            std::transform(covar_index.begin(), covar_index.end(), X.begin() + base_pos + j * common_sample_size, [&cur_table, &covarp](int pos){ return cur_table[(int)covarp[pos]];});
            map_index++;
        }
    }

    int base_rcovar = qcovar.size() + covar.size();
    int map_rindex = 0;
    for(int i = 0; i < rcovar.size(); i++){
        auto &covarp = this->rcovar[i];
        auto &labelp = this->labels_rcovar[i];
        int base_pos = (start_col_X[base_rcovar + i]) * common_sample_size;
        for(int j = 0; j < labelp["LABEL_MAX_VALUE"]; j++){
            auto &cur_table = labels_rcovar_mapping[map_rindex];
            std::transform(covar_index.begin(), covar_index.end(), X.begin() + base_pos + j * common_sample_size, [&cur_table, &covarp](int pos){ return cur_table[(int)covarp[pos]];});
            //std::transform(covarp.begin(), covarp.end(), X.begin() + base_pos + j * common_sample_size, [&cur_table](double covar){ return cur_table[(int)covar];});
            map_rindex++;
        }
    }

    return true;
}

bool Covar::getCovarXRaw(const vector<string> &sampleIDs, vector<double> &X, vector<uint32_t> &keep_index){
    vector<uint32_t> covar_index;
    if(!getCommonSampleIndex(sampleIDs, keep_index, covar_index)){
        return false;
    }

    uint64_t common_sample_size = covar_index.size();
    uint64_t total_col_covar = qcovar.size() + covar.size() + rcovar.size();
    X.resize(total_col_covar * common_sample_size);

    vector<vector<double>> *covars[] = {&qcovar, &covar, &rcovar};
    vector<uint64_t> base_X_positions = {0, qcovar.size(), qcovar.size() + covar.size()};
    vector<uint64_t> look_X_positions = {qcovar.size(), qcovar.size() + covar.size(), total_col_covar};


    #pragma omp parallel for
    for(int i = 0; i < total_col_covar; i++){
        //get which covar
        auto const it = std::lower_bound(look_X_positions.begin(), look_X_positions.end(), i + 1);
        int j = std::distance(look_X_positions.begin(), it);

        auto &pcovar = *covars[j];
        int base_pos = i * common_sample_size;
        auto &t_covar = pcovar[i - base_X_positions[j]];
        std::transform(covar_index.begin(), covar_index.end(), X.begin() + base_pos, 
                [&t_covar](int pos){return t_covar[pos];});
    }

    return true;
/*  // can't be paralleled
    for(int j = 0; j < sizeof(covars) / sizeof(covars[0]); j++){
        auto &pcovar = *covars[j];
        int base_position = base_X_positions[j] * common_sample_size;
        #pragma omp parallel for
        for(int i = 0; i < pcovar.size(); i++){ 
            auto &t_covar = pcovar[i];
            int temp_begin_pos = base_position + i * common_sample_size;  
            std::transform(covar_index.begin(), covar_index.end(), X.begin() + temp_begin_pos, 
                    [&t_covar](int pos){return t_covar[pos];});
        }
    }
*/
}

void Covar::read_covar(string filename, vector<string>& sub_list, vector<vector<double>>* covar, vector<map<string, int>>* labels, vector<int>* keep_row_p){
    string err_string = "[" + filename + "].";
    std::ifstream hcovar(filename.c_str());
    if(!hcovar.good()){
        LOGGER.e(0, "can't read " + err_string);
    }
    string line;
    std::getline(hcovar, line);

    vector<string> line_elements;

    boost::split(line_elements, line, boost::is_any_of("\t "));
    int ncol = line_elements.size();
    boost::replace_all(line_elements[ncol - 1], "\r", "");
    int nkeep = 0;
    int last_keep = 0;
    if(keep_row_p){
        nkeep = keep_row_p->size();
        last_keep = (*keep_row_p)[nkeep - 1];
    }
    int least_col = 2;
    bool has_covar = false;
    if(covar){
        least_col = 3;
        has_covar = true;
    }

    last_keep += 2;
    if(ncol < least_col){
        LOGGER.e(0, "less than " + to_string(least_col) + " columns in " + err_string);
    }
    if(last_keep > ncol){
        LOGGER.e(0, "can't read " + to_string(last_keep) + "th column from " + err_string);
    }

    int line_number = 0;
    string first_element = line_elements[0];
    boost::to_upper(first_element);
    if(first_element[0] != '#' && first_element != "FID"){
        hcovar.clear();
        hcovar.seekg(0);
        line_number = 0;
    }else{
        line_number = 1;
    }

    vector<int> keep_col;
    if(nkeep == 0){
        keep_col.resize(ncol - 2);
        std::iota(keep_col.begin(), keep_col.end(), 2);
        nkeep = keep_col.size();
    }else{
        keep_col.resize(nkeep);
        std::transform(keep_row_p->begin(), keep_row_p->end(), keep_col.begin(), [](int value){return ++value;});
    }

    if(has_covar){
        covar->resize(nkeep);
    }else{
        nkeep = 0;
    }
    bool is_factor = false;
    if(labels){
        is_factor = true;
        labels->resize(nkeep);
        for(int i = 0; i< nkeep; i++){
            (*labels)[i]["LABEL_MAX_VALUE"] = -1;
        }
    }

    int n_item = atoi(options["covar_maxlevel"].c_str());

    while(std::getline(hcovar, line)){
        line_number++;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
        if(line_elements.size() < last_keep){
            LOGGER.e(0, "can't read " + to_string(last_keep) + "th column of line " + to_string(line_number) + " from " + err_string);
        }else if(line_elements.size() != ncol){
            LOGGER.w(0, "inconsistent column number in line " + to_string(line_number) + " from " + err_string);
        }

        bool is_skip = false;
        vector<double> covar_temp(nkeep);
        if(has_covar){

           // get every column in the keep
            for(int col_index = 0; col_index < nkeep; col_index++){
                string temp_string = line_elements[keep_col[col_index]];
                boost::to_upper(temp_string);
                double temp_item;
                if(temp_string == "NA" || temp_string == "NAN" || temp_string == "." || temp_string == "-9"){
                    is_skip = true;
                    break;
                }

                if(is_factor){
                    auto label_p = &(*labels)[col_index];
                    if(label_p->find(temp_string) != label_p->end()){
                        temp_item = (*label_p)[temp_string];
                    }else{ // new factor
                        (*label_p)["LABEL_MAX_VALUE"] += 1;
                        temp_item = (*label_p)["LABEL_MAX_VALUE"];
                        if(temp_item > n_item){
                            LOGGER.e(0, "too many levels in covariate #" + to_string(col_index + 1) + ". You may fit it as a quantitative covariate using --qcovar.");
                        }
                        (*label_p)[temp_string] = temp_item;
                    }
                }else{
                    try{
                        temp_item = boost::lexical_cast<double>(temp_string);
                    }catch(const boost::bad_lexical_cast &e){
                        LOGGER.e(0, "line " + to_string(line_number) + " contains non-numeric values in " + err_string);
                    }
                }
                covar_temp[col_index] = temp_item;
            }
        }

        // covar is complete or not
        if(!is_skip){
            sub_list.push_back(line_elements[0] + "\t" + line_elements[1]);
            for(int i = 0; i < nkeep; i++){
                (*covar)[i].push_back(covar_temp[i]);
            }
        }
    }
}

int Covar::registerOption(map<string, vector<string>>& options_in){
    addOneFileOption("qcovar", "", "--qcovar", options_in, options);
    addOneFileOption("covar", "", "--covar", options_in, options);
    addOneFileOption("rcovar", "", "--rcovar", options_in, options);

    string op = "--covar-maxlevel";
    options["covar_maxlevel"] = "33";
    if(options_in.find(op) != options_in.end()){
        auto options_tmp = options_in[op];
        if(options_tmp.size() >= 1){
            options["covar_maxlevel"] = options_tmp[0];
        }
    }

    if(options_in.find("--test-covar") != options_in.end()){
        string filename = options_in["--test-covar"][0];
        vector<string> samples;
        Covar::read_covar(filename, samples);
        LOGGER << "template: " << samples.size() << std::endl;

        Covar covar;
        vector<double> X;
        vector<uint32_t> sample_index;
        covar.getCovarXRaw(samples, X, sample_index);
        LOGGER<< "samples " << sample_index.size() << std::endl;
        LOGGER << "X size: " << X.size() << std::endl;

        for(int j = 0; j < X.size() / sample_index.size(); j++){
            for(int i = 0; i < sample_index.size(); i++){
                LOGGER << X[i + j * sample_index.size()] << " ";
            }
            LOGGER << std::endl;
        }

        covar.getCovarX(samples, X, sample_index);
        LOGGER<< "samples " << sample_index.size() << std::endl;
        LOGGER << "X size: " << X.size() << std::endl;

        for(int j = 0; j < X.size() / sample_index.size(); j++){
            for(int i = 0; i < sample_index.size(); i++){
                LOGGER << X[i + j * sample_index.size()] << " ";
            }
            LOGGER << std::endl;
        }

    }

    return 0;
}

void Covar::processMain(){
    LOGGER.e(0, "No main function in covariate yet.");
}


//template void Covar::read_covar(string filename, vector<string>& sub_list, vector<vector<int8_t>>* covar, vector<map<string, int8_t>>* labels, vector<int>* keep_row_p);
//template void Covar::read_covar(string filename, vector<string>& sub_list, vector<vector<double>>* covar, vector<map<string, double>>* labels, vector<int>* keep_row_p);
