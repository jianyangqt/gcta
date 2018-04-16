/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: holds covar information

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
    vector_commonIndex(s1, s2, s1_index, s2_index);
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
        vector_commonIndex(sample_id, s3, sa_index, s3_index);
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
        read_covar(filename, samples_qcovar, qcovar, NULL, NULL);  
        if(hasVectorDuplicate(samples_qcovar)){
            LOGGER.e(0, "covariates can't have duplicate FID+IID.");
        }
        LOGGER.i(0, to_string(qcovar.size()) + " covariates of " + to_string(samples_qcovar.size()) + " samples to be included.");
    }

    if(options.find("covar") != options.end()){
        string filename = options["covar"];
        LOGGER.i(0, "Reading discrete covariates from [" + filename + "].");
        read_covar(filename, samples_covar, covar, &labels_covar, NULL);  
        if(hasVectorDuplicate(samples_covar)){
            LOGGER.e(0, "covariates can't have duplicate FID+IID.");
        }
        LOGGER.i(0, to_string(covar.size()) + " covariates of " + to_string(samples_covar.size()) + " samples to be included.");
 
    }

    if(options.find("rcovar") != options.end()){
        string filename = options["rcovar"];
        LOGGER.i(0, "Reading ranked covariates from [" + filename + "].");
        read_covar(filename, samples_rcovar, rcovar, &labels_rcovar, NULL);  
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
            return;
        case 100:
            sample_id = samples_rcovar;
            this->rcovar = rcovar;
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
    }else{
        LOGGER.e(0, "0 covartiate to be included.");
    }
    LOGGER.i(0, to_string(sample_id.size()) + " common samples in covariates to be included.");
}

bool Covar::getCovarX(const vector<string> &sampleIDs, vector<double> &X, vector<int> &keep_index){
    uint64_t total_col_covar = qcovar.size() + covar.size() + rcovar.size();
    if(sampleIDs.size() == 0 || sample_id.size() == 0 || total_col_covar == 0){
        return false;
    }
    vector<int> covar_index;
    vector_commonIndex_sorted1(sampleIDs, sample_id, keep_index, covar_index);
    if(keep_index.size() == 0){
        return false;
    }
    uint64_t common_sample_size = covar_index.size();
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

void Covar::read_covar(string filename, vector<string>& sub_list, vector<vector<double>>& covar, vector<map<string, double>>* labels, vector<int>* keep_row_p){
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
    int nkeep = 0;
    int last_keep = 0;
    if(keep_row_p){
        nkeep = keep_row_p->size();
        last_keep = (*keep_row_p)[nkeep - 1];
    }
    last_keep += 2;
    if(ncol < 3){
        LOGGER.e(0, "less than 3 columns in " + err_string);
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

    covar.resize(nkeep);
    bool is_factor = false;
    if(labels){
        is_factor = true;
        labels->resize(nkeep);
        for(int i = 0; i< nkeep; i++){
            (*labels)[i]["LABEL_MAX_VALUE"] = -1.0;
        }
    }
    
    while(std::getline(hcovar, line)){
        line_number++;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        if(line_elements.size() < last_keep){
            LOGGER.e(0, "can't read " + to_string(last_keep) + "th column of line " + to_string(line_number) + " from " + err_string);
        }else if(line_elements.size() != ncol){
            LOGGER.w(0, "inconsistent column number in line " + to_string(line_number) + " from " + err_string);
        }

        bool is_skip = false;

        vector<double> covar_temp(nkeep);

        // get every column in the keep
        for(int col_index = 0; col_index < nkeep; col_index++){
            string temp_string = line_elements[keep_col[col_index]];
            boost::to_upper(temp_string);
            double temp_item;
            if(temp_string == "NA" || temp_string == "NAN"){
                is_skip = true;
                break;
            }

            if(is_factor){
                auto label_p = &(*labels)[col_index];
                if(label_p->find(temp_string) != label_p->end()){
                    temp_item = (*label_p)[temp_string];
                }else{ // new factor
                    (*label_p)["LABEL_MAX_VALUE"] += 1.0;
                    temp_item = (*label_p)["LABEL_MAX_VALUE"];
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

        // covar is complete or not
        if(!is_skip){
            sub_list.push_back(line_elements[0] + "\t" + line_elements[1]);
            for(int i = 0; i < nkeep; i++){
                covar[i].push_back(covar_temp[i]);
            }
        }
    }
}

int Covar::registerOption(map<string, vector<string>>& options_in){
    addOneFileOption("qcovar", "", "--qcovar", options_in, options);
    addOneFileOption("covar", "", "--covar", options_in, options);
    addOneFileOption("rcovar", "", "--rcovar", options_in, options);
    if(options_in.find("--test-covar") != options_in.end()){
        string filename = options_in["--test-covar"][0];
        vector<string> samples;
        vector<vector<double>> s_covar;
        Covar::read_covar(filename, samples, s_covar);

        Covar covar;
        vector<double> X;
        vector<int> sample_index;
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
