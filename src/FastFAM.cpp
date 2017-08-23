/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   FastFAM regression

   Depends on the class of genotype

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/
#include "FastFAM.h"
#include "main/StatFunc.h"
#include <cmath>
#include <algorithm>
#include <Eigen/SparseCholesky>
#include <sstream>
#include <iterator>
#include "utils.hpp"
#include "Logger.h"
#include "ThreadPool.h"


map<string, string> FastFAM::options;
map<string, double> FastFAM::options_d;
vector<string> FastFAM::processFunctions;

FastFAM::FastFAM(Geno *geno){
    this->geno = geno;
    num_indi = geno->pheno->count_keep();
    num_marker = geno->marker->count_extract();

    beta = new double[num_marker];
    se = new double[num_marker];
    p = new double[num_marker];

    double VG;
    double VR;
    bool flag_est_GE = true;
    if(options.find("G") != options.end()){
        VG = std::stod(options["G"]);
        VR = std::stod(options["E"]);
        flag_est_GE = false;
    }

    vector<string> ids;
    geno->pheno->get_pheno(ids, phenos);
    if(ids.size() != num_indi){
        LOGGER.e(0, "Phenotype is not equal, this shall be a flag bug");
    }
    phenoVec = Map<VectorXd> (phenos.data(), num_indi);

    SpMat fam;

    vector<string> grm_id;
    string ffam_file = options["grmsparse_file"];
    readFAM(ffam_file, fam, ids);

    vector<double> Aij;
    vector<double> Zij;
 
    if(flag_est_GE){
        LOGGER.i(0, "Estimate VG with HE regression");
        for(int k = 0; k < fam.outerSize(); ++k){
            for(SpMat::InnerIterator it(fam, k); it; ++it){
                if(it.row() < it.col()){
                    Aij.push_back(it.value());
                    Zij.push_back(phenos[it.row()] * phenos[it.col()]);
                }
            }
        }

        VG = HEreg(Zij, Aij);
        VR = 1.0 - VG;
        LOGGER.i(2, "VG=" + to_string(VG) + ", VE=" + to_string(VR));
    }

    inverseFAM(fam, VG, VR);
}

double FastFAM::HEreg(vector<double> &Zij, vector<double> &Aij){
    Map<VectorXd> ZVec(Zij.data(), Zij.size());
    Map<VectorXd> AVec(Aij.data(), Aij.size());

    double Zmean = ZVec.mean();
    double Amean = AVec.mean();
    ZVec -= (VectorXd::Ones(ZVec.size())) * Zmean;
    AVec -= (VectorXd::Ones(AVec.size())) * Amean;

    double A2v = (AVec.transpose() * AVec)(0, 0);
    if(A2v < 1e-6){
        LOGGER.e(0, "can't solve the regression");
    }
    double AZ = (AVec.transpose() * ZVec)(0, 0);
    return (1.0 / A2v) * AZ;
}
    

void FastFAM::readFAM(string filename, SpMat& fam, vector<string> &ids){
    uint32_t num_indi = ids.size();
    std::vector<string> sublist = Pheno::read_sublist(filename + ".grm.id");
    if(sublist != ids){
        LOGGER.e(0, "Sparse GRM has some difference to the kept samples in genotype, you can prune the GRM again with --keep");
    }

    std::ifstream pair_list((filename + ".grm.sp").c_str());
    if(!pair_list){
        LOGGER.e(0, "can't read [" + filename + ".grm.sp]");
    }

    string line;
    int line_number = 0;
    int last_length = 0;

    vector<uint32_t> id1;
    vector<uint32_t> id2;
    vector<double> grm;

    vector<uint32_t> num_elements(num_indi, 0);

    uint32_t tmp_id1 = 0, tmp_id2 = 0;
    double tmp_grm = 0.0;

    while(std::getline(pair_list, line)){
        line_number++;
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);

        tmp_id1 = (std::stoi(line_elements[0]));
        tmp_id2 = (std::stoi(line_elements[1]));
        tmp_grm = std::stod(line_elements[2]);
        id1.push_back(tmp_id1);
        id2.push_back(tmp_id2);
        num_elements[tmp_id2] += 1;
        grm.push_back(tmp_grm);
        if(tmp_id1 != tmp_id2){
            id1.push_back(tmp_id2);
            id2.push_back(tmp_id1);
            num_elements[tmp_id1] += 1;
            grm.push_back(tmp_grm);
        }
    }
    pair_list.close();

    auto sorted_index = sort_indexes(id2, id1);

    fam.resize(num_indi, num_indi);
    fam.reserve(num_elements);

    for(auto index : sorted_index){
        fam.insertBackUncompressed(id1[index], id2[index]) = grm[index];
    }
    fam.finalize();
    fam.makeCompressed();

}

void FastFAM::inverseFAM(SpMat& fam, double VG, double VR){
    SpMat eye(fam.rows(), fam.cols());
    eye.setIdentity();

    // V
    fam *= VG;
    fam += eye * VR;

    Eigen::SimplicialLDLT<SpMat> solver;
    solver.compute(fam);

    if(solver.info() != Eigen::Success){
        LOGGER.e(0, "can't inverse the FAM");
    }

    LOGGER.i(0, "Inversing the FAM, this may take long time");
    V_inverse = solver.solve(eye);
}


void FastFAM::calculate_fam(uint8_t *buf, int num_marker){
    // Memory fam_size * 2 * 4 + (N * 8 * 2 ) * thread_num + M * 3 * 8  B
    int num_thread = THREADS.getThreadCount() + 1; 

    int num_marker_part = (num_marker + num_thread - 1) / num_thread;
    for(int index = 0; index != num_thread - 1; index++){
        THREADS.AddJob(std::bind(&FastFAM::reg_thread, this, buf, index * num_marker_part, (index + 1) * num_marker_part));
    }

    reg_thread(buf, (num_thread - 1) * num_marker_part, num_marker);

    THREADS.WaitAll();

    num_finished_marker += num_marker;
    if(num_finished_marker % 30000 == 0){
        LOGGER.i(2, to_string(num_finished_marker) + " markers finished"); 
    }
}

void FastFAM::reg_thread(uint8_t *buf, int from_marker, int to_marker){
    double *w_buf = new double[num_indi];
    Map< VectorXd > xMat(w_buf, num_indi);
    MatrixXd XMat_V;
    for(int cur_marker = from_marker; cur_marker < to_marker; cur_marker++){
        geno->makeMarkerX(buf, cur_marker, w_buf, true, false);
        // Xt * V-1
        MatrixXd xMat_V = xMat.transpose() * V_inverse;
        // 
        double xMat_V_x = 1.0 / (xMat_V * xMat)(0, 0);
        double xMat_V_p = (xMat_V * phenoVec)(0, 0);
        
        double temp_beta =  xMat_V_x * xMat_V_p;
        double temp_se = sqrt(xMat_V_x);
        double temp_z = temp_beta / temp_se;

        uint32_t cur_raw_marker = num_finished_marker + cur_marker;

        beta[cur_raw_marker] = temp_beta; //* geno->RDev[cur_raw_marker]; 
        se[cur_raw_marker] = temp_se;
        {
            std::lock_guard<std::mutex> lock(chisq_lock);
            p[cur_raw_marker] = StatFunc::pchisq(temp_z * temp_z, 1); 
        } 
    }
    delete[] w_buf;
}


void FastFAM::output(string filename){
    //TODO get the real effect
    std::ofstream out(filename.c_str());
    vector<string> header{"CHR", "SNP", "POS", "A1", "A2", "AF1", "beta", "se", "p"};
    std::copy(header.begin(), header.end(), std::ostream_iterator<string>(out, "\t"));
    out << std::endl;
    for(int index = 0; index != num_marker; index++){
        out << geno->marker->get_marker(geno->marker->getExtractIndex(index)) << "\t" <<
            geno->AFA1[index] << "\t" << beta[index] << "\t" << se[index] << "\t" << p[index] << std::endl;
    }
    out.close();
    LOGGER.i(0, "Success:", "saved result to [" + filename +"]");
}

int FastFAM::registerOption(map<string, vector<string>>& options_in){
    int returnValue = 0;
    options["out"] = options_in["out"][0] + ".fastFAM.assoc";

    string curFlag = "--fastFAM";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        returnValue++;
        options_in.erase(curFlag);
    }

    curFlag = "--grm-sparse";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["grmsparse_file"] = options_in[curFlag][0];
        }else{
            LOGGER.e(0, curFlag + "can't deal with 0 or > 1 files");
        }
        options_in.erase(curFlag);
    }

    curFlag = "--ge";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 2){
            options["G"] = options_in[curFlag][0];
            options["E"] = options_in[curFlag][1];
        }else{
            LOGGER.e(0, curFlag + " can't handle other than 2 numbers");
        }
        options_in.erase(curFlag);
    }

    return returnValue;
}

void FastFAM::processMain(){
    vector<function<void (uint8_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "fast_fam"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            FastFAM ffam(&geno);

            callBacks.push_back(bind(&Geno::freq, &geno, _1, _2));
            callBacks.push_back(bind(&FastFAM::calculate_fam, &ffam, _1, _2));
            geno.loop_block(callBacks);

            ffam.output(options["out"]);
        }
    }
}


