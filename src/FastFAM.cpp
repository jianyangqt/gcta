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
#include "StatLib.h"
#include <cmath>
#include <algorithm>
#include <Eigen/SparseCholesky>
#include <Eigen/PardisoSupport>
#include <Eigen/IterativeLinearSolvers>
#include <sstream>
#include <iterator>
#include "utils.hpp"
#include "Logger.h"
#include "ThreadPool.h"
#include "omp.h"
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iomanip>

#include <iostream>

using std::to_string;
using Eigen::Matrix;

map<string, string> FastFAM::options;
map<string, double> FastFAM::options_d;
vector<string> FastFAM::processFunctions;

FastFAM::FastFAM(Geno *geno){
    //Eigen::setNbThreads(THREADS.getThreadCount() + 1);
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

    // read covar
    vector<uint32_t> remain_index, remain_index_covar;
    bool has_qcovar = false;
    Covar covar;
    if(covar.getCommonSampleIndex(ids, remain_index, remain_index_covar){
        has_covar = true;
        LOGGER.i(0, to_string(remain_index.size()) + " overlapped individuals with non-missing data to be included from the covariate file(s).");
    }else{
        remain_index.resize(ids.size());
        std::iota(remain_index.begin(), remain_index.end(), 0);
    }

    vector<string> remain_ids(remain_index.size());
    std::transform(remain_index.begin(), remain_index.end(), remain_ids.begin(), [&ids](size_t pos){return ids[pos];});

    // read fam
    string ffam_file = "";
    bool fam_flag = true;
    if(options.find("grmsparse_file") != options.end()){
        ffam_file = options["grmsparse_file"];
    }else{
        fam_flag = false;
    }

    // index in ids when merge to spare fam
    vector<uint32_t> remain_index_fam;
    SpMat fam;

    if(fam_flag){
        readFAM(ffam_file, fam, remain_ids, remain_index_fam);
    }else{
        remain_index_fam.resize(remain_ids.size());
        std::iota(remain_index_fam.begin(), remain_index_fam.end(), 0);
    }

    int n_remain_index_fam = remain_index_fam.size();

    //reorder phenotype, covar
    vector<double> remain_phenos(n_remain_index_fam);
    vector<string> remain_ids_fam(n_remain_index_fam);
    for(int i = 0; i != n_remain_index_fam; i++){
        remain_phenos[i] = phenos[remain_index[remain_index_fam[i]]];
        remain_ids_fam[i] = ids[remain_index[remain_index_fam[i]]];
    }
    LOGGER.i(0, "After matching all the files, " + to_string(remain_phenos.size()) + " individuals to be included in the analysis.");

    vector<double> remain_covar;
    vector<int> remain_inds_index;
    if(has_covar){
        covar.getCovarX(remain_ids_fam, remain_covar, remain_inds_index);
        remain_covar.resize(remain_covar.size() + n_remain_index_fam);
        std::fill(remain_covar.end() - n_remain_index_fam, remain_covar.end(), 1.0);
    }

    // standerdize the phenotype, and condition the covar
    phenoVec = Map<VectorXd> (remain_phenos.data(), remain_phenos.size());
    // condition the covar
    if(has_qcovar){
        MatrixXd concovar = Map<Matrix<double, Dynamic, Dynamic, Eigen::ColMajor>>(remain_covar.data(), remain_phenos.size(), 
                remain_covar.size() / remain_phenos.size());
        conditionCovarReg(phenoVec, concovar);
    }

    // Center
    double phenoVec_mean = phenoVec.mean();
    phenoVec -= VectorXd::Ones(phenoVec.size()) * phenoVec_mean;

    if(fam_flag){
        double Vpheno = phenoVec.array().square().sum() / (phenoVec.size() - 1);
        //phenoVec /= pheno_sd;

        LOGGER.i(0, "DEBUG: conditioned Pheno (first 5)");

        for(int i = 0; i < 5; i++){
            LOGGER.i(0, to_string(phenoVec[i]));
        }

        if(options.find("inv_file") == options.end()){
            vector<double> Aij;
            vector<double> Zij;

            if(flag_est_GE){
                LOGGER.i(0, "Estimating the genetic variance (Vg) by HE regression...");
                if(options["rel_only"] == "yes"){
                    LOGGER.i(0, "Use related pairs only.");
                    for(int k = 0; k < fam.outerSize(); ++k){
                        for(SpMat::InnerIterator it(fam, k); it; ++it){
                            if(it.row() < it.col()){
                                Aij.push_back(it.value());
                                Zij.push_back(phenoVec[it.row()] * phenoVec[it.col()]);
                            }
                        }
                    }

                    VG = HEreg(Zij, Aij);
                }else{
                    VG = HEreg(fam, phenoVec);
                }
                VR = Vpheno - VG;
                LOGGER.i(2, "Ve = " + to_string(VR));
                LOGGER.i(2, "Heritablity = " + to_string(VG/Vpheno));
            }

            inverseFAM(fam, VG, VR);
            if(options.find("save_inv") != options.end()){
                LOGGER.i(0, "Saving inverse of V for further analysis, use --load-inv for further analysis");
                std::ofstream inv_id((options["out"]+".grm.id").c_str());

                std::ofstream inv_out((options["out"]+".grm.inv").c_str());
                inv_out << std::setprecision( std::numeric_limits<double>::digits10+2); 
                for(int k = 0; k < V_inverse.outerSize(); ++k){
                    for(SpMat::InnerIterator it(V_inverse, k); it; ++it){
                        inv_out << it.row() << "\t" << it.col() << "\t" << it.value() << std::endl;
                    }
                }
                inv_out.close();
                LOGGER.i(0, "The inverse has been saved to [" + options["out"] + ".grm.inv]");
            }
        }else{
            V_inverse.resize(phenoVec.size(), phenoVec.size());
            string in_name = options["inv_file"] + ".grm.inv";
            LOGGER.i(0, "Loading inverse of V from " + in_name);
            LOGGER.ts("LOAD_INV");
            std::ifstream in_file(in_name.c_str());
            if(!in_file){
                LOGGER.e(0, "can't open the file");
            }
            string line;
            while(getline(in_file, line)){
                vector<string> line_elements;
                boost::split(line_elements, line, boost::is_any_of("\t "));
                if(line_elements.size() != 3){
                    LOGGER.e(0, "the inversed file seems to be incorrect");
                }
                V_inverse.insertBackUncompressed(stoi(line_elements[0]), stoi(line_elements[1])) 
                    = stod(line_elements[2]);
            }
            V_inverse.finalize();
            V_inverse.makeCompressed();
            LOGGER.i(0, "Inverse of V loaded in " + to_string(LOGGER.tp("LOAD_INV")) + " seconds");
        }
    }
}

void FastFAM::conditionCovarReg(VectorXd &pheno, MatrixXd &covar){
    MatrixXd t_covar = covar.transpose();
    VectorXd beta = (t_covar * covar).ldlt().solve(t_covar * pheno);
    LOGGER.i(0, "DEBUG: condition betas:");
    LOGGER << beta << std::endl;
    //VectorXd beta = covar.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(pheno);
    pheno -= covar * beta;
    //double pheno_mean = pheno.mean();
    //pheno -= (VectorXd::Ones(pheno.size())) * pheno_mean;
}

double FastFAM::HEreg(const Ref(const SpMat) fam, const Ref(const VectorXd) pheno){
    int num_covar = 1;
    int num_components = 1;
    int col_X = num_covar + num_components;
    MatrixXd XtX = MatrixXd::Zero(col_X, col_X);
    VectorXd XtY = VectorXd::Zero(col_X);
    double SSy = 0;

    uint64_t size = fam.cols();
    XtX(0, 0) = size;

    for(int i = 1; i < size; i++){
        double temp_pheno = pheno[i];
        auto fam_block = fam.block(0, i, i, 1);
        auto pheno_block = pheno.block(0, i);
        SSy += (pheno_block * temp_pheno).sum();


    }

    MatrixXd XtXi = XtX.inverse();
    VectorXd betas = XtXi * XtY;

    double sse = (SSy - betas.dot(XtY)) / (size - col_X);

    VectorXd SDs = sse * XtXi.diagonal();

    double hsq = betas(betas.size() - 1);
    double SD = SDs(SDs.size() - 1);
    double Zsq = hsq^2 / SD;
    double p = StatLib::pchisqd1(Zsq);

    LOGGER.i(2, "Vg = " + to_string(hsq) + ", se = " + to_string(sqrt(SD)) +  ", P = " + to_string(p));

    if(p > 0.05){
        LOGGER.e(0, "The estimate of Vg is not statistically significant. "
                "This is likely because the number of relatives is not large enough. "
                "We do not recommend to run fastFAM in this case.");
    }
    return hsq;
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
    double hsq = (1.0 / A2v) * AZ;

    VectorXd RZVec = ZVec - AVec * hsq;

    double delta = RZVec.array().square().sum() / (RZVec.size() - 2);
    double se = sqrt(delta / A2v);

    double z = hsq / se;

    double p = StatLib::pchisqd1(z * z);

    LOGGER.i(2, "Vg = " + to_string(hsq) + ", se = " + to_string(se) +  ", P = " + to_string(p));

    if(p > 0.05){
        LOGGER.e(0, "The estimate of Vg is not statistically significant. "
                "This is likely because the number of relatives is not large enough. "
                "We do not recommend to run fastFAM in this case.");
    }

    return hsq;
}
    

void FastFAM::readFAM(string filename, SpMat& fam, const vector<string> &ids, vector<uint32_t> &remain_index){
    LOGGER.i(0, "Reading the sparse GRM file from [" + filename + "]...");
    uint32_t num_indi = ids.size();
    vector<string> sublist = Pheno::read_sublist(filename + ".grm.id");
    vector<uint32_t> fam_index;
    vector_commonIndex(sublist, ids, fam_index, remain_index);
    LOGGER.i(0, "DEBUG: " + to_string(fam_index.size()) + " subjects remained");

    //Fix index order to outside, that fix the phenotype, covar order
    vector<size_t> index_list_order = sort_indexes(remain_index);
    vector<uint32_t> ordered_fam_index(remain_index.size(), 0);
    vector<uint32_t> ordered_remain_index(remain_index.size(), 0);
    std::transform(index_list_order.begin(), index_list_order.end(), ordered_fam_index.begin(), [&fam_index](size_t pos){
            return fam_index[pos];});
    std::transform(index_list_order.begin(), index_list_order.end(), ordered_remain_index.begin(), [&remain_index](size_t pos){
            return remain_index[pos];});
    remain_index = ordered_remain_index;

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

    map<uint32_t, uint32_t> map_index;
    for(uint32_t index = 0; index != ordered_fam_index.size(); index++){
        map_index[ordered_fam_index[index]] = index;
    }


    while(std::getline(pair_list, line)){
        line_number++;
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);

        uint32_t tmp_id1 = (std::stoi(line_elements[0]));
        uint32_t tmp_id2 = (std::stoi(line_elements[1]));
        if(map_index.find(tmp_id1) != map_index.end() &&
                map_index.find(tmp_id2) != map_index.end()){
            tmp_id1 = map_index[tmp_id1];
            tmp_id2 = map_index[tmp_id2];

            double tmp_grm = std::stod(line_elements[2]);
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
    }
    pair_list.close();

    auto sorted_index = sort_indexes(id2, id1);

    fam.resize(ordered_fam_index.size(), ordered_fam_index.size());
    fam.reserve(num_elements);

    for(auto index : sorted_index){
        fam.insertBackUncompressed(id1[index], id2[index]) = grm[index];
    }
    fam.finalize();
    fam.makeCompressed();

}

void FastFAM::inverseFAM(SpMat& fam, double VG, double VR){
    LOGGER.i(0, "Inverting the variance-covarinace matrix (This may take a long time).");
    LOGGER.i(0, string("Inverse method: ") + options["inv_method"]);
    LOGGER.i(0, "DEUBG: Inverse Threads " + to_string(Eigen::nbThreads()));
    LOGGER.ts("INVERSE_FAM");
    SpMat eye(fam.rows(), fam.cols());
    LOGGER.i(0, "FAM " + to_string(fam.rows()) + " * " + to_string(fam.cols()));
    eye.setIdentity();

    // V
    fam *= VG;
    fam += eye * VR;

    // inverse
    if(options["inv_method"] == "ldlt"){
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.compute(fam);
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }
        V_inverse = solver.solve(eye);
    }else if(options["inv_method"] == "cg"){
        Eigen::ConjugateGradient<SpMat> solver;
        solver.compute(fam);
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }
        V_inverse = solver.solve(eye);
    }else if(options["inv_method"] == "llt"){
        Eigen::SimplicialLLT<SpMat> solver;
        solver.compute(fam);
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }
        V_inverse = solver.solve(eye);
    }else if(options["inv_method"] == "pardiso1"){
        /*
        Eigen::PardisoLLT<SpMat> solver;
        solver.compute(fam);
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }
        V_inverse = solver.solve(eye);
        */
    }else if(options["inv_method"] == "tcg"){
        Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> solver;
        solver.compute(fam);
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }
        V_inverse = solver.solve(eye);
    }else if(options["inv_method"] == "lscg"){
        Eigen::LeastSquaresConjugateGradient<SpMat> solver;
        solver.compute(fam);
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }
        V_inverse = solver.solve(eye);
    }else{
        LOGGER.e(0, "Unknown inverse methods");
    }


    //solver.setTolerance(1e-3);
    ///solver.setMaxIterations(10);;

   //LOGGER.i(0, "# iteations: " + to_string(solver.iterations()));
    //LOGGER.i(0, "# error: " + to_string(solver.error()));

    LOGGER.i(0, "Inverted in " + to_string(LOGGER.tp("INVERSE_FAM")) + " seconds");
}

void FastFAM::calculate_gwa(uint64_t *buf, int num_marker){
    #pragma omp parallel for schedule(dynamic)
    for(int cur_marker = 0; cur_marker < num_marker; cur_marker++){
        double *w_buf = new double[num_indi];
        Map< VectorXd > xMat(w_buf, num_indi);
        MatrixXd XMat_V;

        geno->makeMarkerX(buf, cur_marker, w_buf, true, false);

        MatrixXd tMat_V = xMat.transpose();

        double xMat_V_x = 1.0 / (tMat_V * xMat)(0, 0);
        double xMat_V_p = (tMat_V * phenoVec)(0, 0);

        double temp_beta =  xMat_V_x * xMat_V_p;
        double temp_se = sqrt(xMat_V_x);
        double temp_z = temp_beta / temp_se;

        uint32_t cur_raw_marker = num_finished_marker + cur_marker;

        beta[cur_raw_marker] = temp_beta; //* geno->RDev[cur_raw_marker]; 
        se[cur_raw_marker] = temp_se;
        p[cur_raw_marker] = StatLib::pchisqd1(temp_z * temp_z); 
        delete[] w_buf;
    }

    num_finished_marker += num_marker;
    if(num_finished_marker % 30000 == 0){
        LOGGER.i(2, to_string(num_finished_marker) + " markers finished"); 
    }

}


void FastFAM::calculate_fam(uint64_t *buf, int num_marker){
    // Memory fam_size * 2 * 4 + (N * 8 * 2 ) * thread_num + M * 3 * 8  B
    //int num_thread = THREADS.getThreadCount() + 1; 
    #pragma omp parallel for schedule(dynamic)
    for(int cur_marker = 0; cur_marker < num_marker; cur_marker++){
        double *w_buf = new double[num_indi];
        Map< VectorXd > xMat(w_buf, num_indi);
        MatrixXd XMat_V;

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
        p[cur_raw_marker] = StatLib::pchisqd1(temp_z * temp_z); 
        delete[] w_buf;
    }
/*
    int num_thread = omp_get_max_threads();
    int num_marker_part = (num_marker + num_thread - 1) / num_thread;
    #pragma omp parallel for
    for(int index = 0; index <= num_thread; index++){
        if(index != num_thread){
            reg_thread(buf, index * num_marker_part, (index + 1) * num_marker);
        }else{
            reg_thread(buf, (num_thread - 1) * num_marker_part, num_marker);
        }
        //THREADS.AddJob(std::bind(&FastFAM::reg_thread, this, buf, index * num_marker_part, (index + 1) * num_marker_part));
    }

    //THREADS.WaitAll();
*/
    num_finished_marker += num_marker;
    if(num_finished_marker % 30000 == 0){
        LOGGER.i(2, to_string(num_finished_marker) + " markers finished"); 
    }
}
/*
void FastFAM::reg_thread(uint8_t *buf, int from_marker, int to_marker){
    //Eigen::setNbThreads(1);
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
            p[cur_raw_marker] = StatLib::pchisqd1(temp_z * temp_z); 
        } 
    }
    delete[] w_buf;
}
*/


void FastFAM::output(string filename){
    std::ofstream out(filename.c_str());
    vector<string> header{"CHR", "SNP", "POS", "A1", "A2", "AF1", "beta", "se", "p"};
    //std::copy(header.begin(), header.end(), std::ostream_iterator<string>(out, "\t"));
    string header_string = boost::algorithm::join(header, "\t");
    out << header_string << std::endl;
    for(int index = 0; index != num_marker; index++){
        out << geno->marker->get_marker(geno->marker->getExtractIndex(index)) << "\t" <<
            geno->AFA1[index] << "\t" << beta[index] << "\t" << se[index] << "\t" << p[index] << std::endl;
    }
    out.close();
    LOGGER.i(0, "The association results have been saved to [" + filename +"].");
}

int FastFAM::registerOption(map<string, vector<string>>& options_in){
    int returnValue = 0;
    //DEBUG: change to .fastFAM
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

    /*
    curFlag = "--qcovar";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["concovar"] = options_in[curFlag][0];
        }else{
            LOGGER.e(0, curFlag + "can't deal with covar other than 1");
        }
    }
    */

    options["inv_method"] = "ldlt";
    vector<string> flags = {"--cg", "--ldlt", "--llt", "--pardiso", "--tcg", "--lscg"};
    for(auto curFlag : flags){
        if(options_in.find(curFlag) != options_in.end()){
            boost::erase_all(curFlag, "--");
            options["inv_method"] = curFlag;
            options_in.erase(curFlag);
        }
    }

    curFlag = "--save-inv";
    if(options_in.find(curFlag) != options_in.end()){
        options["save_inv"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--load-inv";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["inv_file"] = options_in[curFlag][0];
        }else{
            LOGGER.e(0, "can't load multiple --load-inv files");
        }
    }

    curFlag = "--rel-only";
    if(options_in.find(curFlag) != options_in.end()){
        options["rel_only"] = "yes";
        options_in.erase(curFlag);
    }else{
        options["rel_only"] = "no";
    }

    return returnValue;
}

void FastFAM::processMain(){
    vector<function<void (uint64_t *, int)>> callBacks;
    //THREADS.JoinAll();
    for(auto &process_function : processFunctions){
        if(process_function == "fast_fam"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            FastFAM ffam(&geno);

            if(options.find("save_inv") != options.end()){
                LOGGER.i(0, "Using --load-inv to load the inversed file for fastFAM");
                return;
            }
            LOGGER.i(0, "Running fastFAM...");
            //Eigen::setNbThreads(1);
            callBacks.push_back(bind(&Geno::freq64, &geno, _1, _2));
            if(options.find("grmsparse_file") != options.end()){
                callBacks.push_back(bind(&FastFAM::calculate_fam, &ffam, _1, _2));
            }else{
                callBacks.push_back(bind(&FastFAM::calculate_gwa, &ffam, _1, _2));
            }
            geno.loop_64block(marker.get_extract_index(), callBacks);

            ffam.output(options["out"]);
        }
    }
}


