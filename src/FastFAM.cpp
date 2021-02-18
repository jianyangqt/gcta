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
#include "OptionIO.h"
#include "StatLib.h"
#include <cmath>
#include <algorithm>
#include <Eigen/SparseCholesky>
#include <Eigen/PardisoSupport>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <sstream>
#include <iterator>
#include "utils.hpp"
#include "Logger.h"
#include "ThreadPool.h"
#include "omp.h"
#include "Matrix.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include "Covar.h"
#include <cstdio>
#include <random>
#include <chrono>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <boost/math/distributions/beta.hpp>

struct InvItem{
    int32_t row;
    int32_t col;
    double val;
};

struct MDLHeader{
    char magic[5];
    uint32_t num_indi;
    bool fam_flag;
    bool isGrammar;
    uint64_t phenoVec_start;
    bool covarFlag;
    uint64_t covarVec_start;
    uint32_t covarVec_cols;
    uint64_t V_inverse_start;
    uint32_t V_inverse_items;
    double c_inf;
};
    


using std::to_string;
using Eigen::Matrix;
using Eigen::ArrayXd;

map<string, string> FastFAM::options;
map<string, double> FastFAM::options_d;
vector<string> FastFAM::processFunctions;


struct K1Res{
    double root;
    int nIter;
    bool bConverge;
};

void printVector(const vector<double> &v, string label){
    LOGGER << label << ": ";
    for(int i = 0; i < v.size(); i++){
        LOGGER << v[i] << " ";
    }
    LOGGER << std::endl;
}

class SPA{
private:
    double q, qinv;
    double gPos, gNeg;
    static ArrayXd mu;
    static int nSample;
    ArrayXd mu1muG;
    ArrayXd geno;
    bool bFast;
    ArrayXd muNZ;
    ArrayXd mu1muNZG;
    ArrayXd gNZ;
    double NAmu;
    double NAsigma;

    double K1adj(double t1, double q){
        //ArrayXd temp1 = mu1 * (-geno * t1).exp() + mu;
        //ArrayXd temp2 = mu * geno;
        if(bFast){
            double temp2dtemp1 = ((muNZ * gNZ)/((1.0 - muNZ) * (-gNZ * t1).exp() + muNZ)).sum();
            return temp2dtemp1 + NAmu + NAsigma * t1 - q;
        }else{
            return ((mu * geno)/((1.0 - mu) * (-geno * t1).exp() + mu)).sum() - q;
        }
    }

    double K2(double t1){
        if(bFast){
            ArrayXd genot1exp = (-gNZ * t1).exp();
            ArrayXd div21 = (mu1muNZG * genot1exp) / ((1.0 - muNZ) * genot1exp + muNZ).square();
            return ((!div21.isNaN()).select(div21, 0)).sum() + NAsigma;
        }else{
            ArrayXd genot1exp = (-geno * t1).exp();
            ArrayXd div21 = (mu1muG * genot1exp) / ((1.0 - mu) * genot1exp + mu).square();
            return ((!div21.isNaN()).select(div21, 0)).sum();
        }
    }

    double Korg(double t1){
        return (1.0 - mu + mu * (geno * t1).exp()).log().sum();
    }

    struct KS{ double k1; double k2;};

    inline void KorgK2(double t1, KS *ks){
        if(bFast){
            ArrayXd genot1exp = (-gNZ * t1).exp();
            ArrayXd div21 = (mu1muNZG * genot1exp) / ((1.0 - muNZ) * genot1exp + muNZ).square();
            ks->k2 = ((!div21.isNaN()).select(div21, 0)).sum() + NAsigma;
            ks->k1 = (1.0 - muNZ + muNZ / genot1exp).log().sum() + NAmu * t1 + 0.5 * NAsigma * t1 * t1;
        }else{
            ArrayXd genot1exp = (-geno * t1).exp();
            ArrayXd div21 = (mu1muG * genot1exp) / ((1.0 - mu) * genot1exp + mu).square();
            ks->k2 = ((!div21.isNaN()).select(div21, 0)).sum();
            ks->k1 = (1.0 - mu + mu / genot1exp).log().sum();
        }
    }

    double getSaddleProb(double zeta, double q){
        KS ks;
        KorgK2(zeta, &ks);
        //double k1 = Korg(zeta);
        //double k2 = K2(zeta);
        double k1 = ks.k1;
        double k2 = ks.k2;
        double pval = 0;
        if(std::isfinite(k1) && std::isfinite(k2)){
            double temp1 = zeta * q - k1;
            double w = sgn(zeta) * std::sqrt(2.0) * std::sqrt(temp1);
            double wi = 1.0 / w;
            double v = zeta * std::sqrt(k2);
            double z_test = w + wi * std::log(v * wi); 
            if(z_test > 0){
                pval = StatLib::pnorm(z_test, false);
            }else{
                pval = StatLib::pnorm(z_test, true);
            }
        }
        return(pval);
    }

    K1Res getRootK1(double init, double q, double thresh=0.0001220703, int maxIter=1000){
        K1Res res;
        if(q >= gPos || q <= gNeg){
            res.root = std::numeric_limits<double>::infinity();
            res.nIter = 0;
            res.bConverge = true;
            return(res);
        }else{
            double t1 = init;
            double K1eval = K1adj(t1, q);
            double prevJump = std::numeric_limits<double>::infinity();
            int nIter = 1;
            while(true){
                double K2eval = K2(t1);
                double tnew = t1 - K1eval / K2eval;
                if(!std::isfinite(tnew)){
                    res.bConverge = false;
                    break;
                }
                if(std::abs(tnew - t1) <= thresh){
                    res.bConverge = true;
                    break;
                }
                if(nIter == maxIter){
                    res.bConverge = false;
                    break;
                }
                double newK1 = K1adj(tnew, q);
                if(sgn(K1eval) != sgn(newK1)){
                    double absTnewT1 = std::abs(tnew - t1);
                    if(absTnewT1 > prevJump - thresh){
                        tnew = t1 + sgn(newK1 - K1eval) * prevJump * 0.5;
                        newK1 = K1adj(tnew, q);
                        prevJump = prevJump * 0.5;
                    }else{
                        prevJump = absTnewT1;
                    }
                }
                nIter++;
                t1 = tnew;
                K1eval = newK1;
            }
            res.root = t1;
            res.nIter = nIter;
            /*
               if(!res.bConverge){
               printVector(K1evals, "K1eval");
               printVector(prevJumps, "prevJump");
               printVector(K2evals, "K2eval");
               printVector(newK1s, "newK1");
               printVector(tnews, "tnew");
            //output the tempory
            }
            */
            return(res);
        }
    }

public:
    static void setMu(VectorXd mu){
        SPA::mu = mu.array();
        //SPA::mu1 = 1.0 - SPA::mu;
        //SPA::mu1mu = SPA::mu1 * SPA::mu;
        SPA::nSample = mu.size();
    }

    SPA(double q, double qinv, Ref<VectorXd> rgen, const vector<uint32_t> &index){
        this->q = q;
        this->qinv = qinv;
        //new (&geno) Map<VectorXd>(xp, nSample);
        this->geno = rgen.array();
        //LOGGER << "geno address: " << static_cast<void *>(geno.data()) << std::endl;
        //LOGGER << "orig address: " << static_cast<void *>(rgen.data()) << std::endl;
        //LOGGER.e(0, "debug");
        //double gPos = (geno.array()>0).select(geno, 0).sum();
        //double gNeg = (geno.array()<0).select(geno, 0).sum();
        gPos = 0;
        gNeg = 0;
        for(int i = 0; i < nSample; i++){
            double tgeno = geno[i];
            if(tgeno > 0){
                gPos+=tgeno;
            }else{
                gNeg+=tgeno;
            }
        }

        mu1muG = mu * (1.0 - mu) * geno.square();

        int nNZ = index.size();
        bFast = false;
        if((double)nNZ / nSample < 0.5){
            bFast = true;
            muNZ.resize(nNZ);
            mu1muNZG.resize(nNZ);
            gNZ.resize(nNZ);
            int curIndex = 0;
            for(int i=0; i < nNZ; i++){
                int tempIndex = index[i];
                muNZ[curIndex] = mu[tempIndex];
                mu1muNZG[curIndex] = mu1muG[tempIndex];
                gNZ[curIndex] = geno[tempIndex];
                curIndex++;
            }

            NAmu = (qinv + q) * 0.5 - (gNZ * muNZ).sum();
            NAsigma = mu1muG.sum() - mu1muNZG.sum();

        }
    }

    void saddleProb(SPARes *res){
        K1Res k1Res = getRootK1(0, q);
        K1Res k2Res = getRootK1(0, qinv);
        if(k1Res.bConverge && k2Res.bConverge){
            double p1 = getSaddleProb(k1Res.root, q);
            double p2 = getSaddleProb(k2Res.root, qinv);
            res->p_adj = std::abs(p1) + std::abs(p2);
            res->bConverge = true;
        }else{
            res->p_adj = std::numeric_limits<double>::quiet_NaN();
            res->bConverge = false;
        }

        if(res->p_adj != 0 && res->p / res->p_adj > 1000){
            res->p_adj = res->p;
            res->bConverge = false;
        }
    }

};
ArrayXd SPA::mu;
int SPA::nSample = 0;

void FastFAM::loadBinModel(){
        LOGGER << "Loading saved GLM model prefixed with [" << options["model_file"] << "]..." << std::endl;
        LOGGER << "Note: phenotype, covariates, sparse GRM and association test methods are provided by the model, thus these flags will be ignored." << std::endl;
        string model_file = options["model_file"] + ".mdl";

        string id_file = model_file + ".id";
        vector<string> modelIDs = Pheno::read_sublist(id_file);
        num_indi = modelIDs.size();

        uint32_t sample_keep_geno = pheno->count_keep();
        vector<string> sampleIDs = pheno->get_id(0, sample_keep_geno - 1, "\t");
        vector<uint32_t> id1, id2;
        vector_commonIndex_sorted1(sampleIDs, modelIDs, id1, id2);
        if(id2.size() != num_indi){
            LOGGER.e(0, "Some sample ID in model are not existed in genotype!");
        }

        pheno->filter_keep_index(id1);
        LOGGER << "  " << num_indi << " valid individuals to be included." << std::endl;

        fam_flag = true;
        covarFlag = true;
        if(fam_flag){
            options["grmsparse_file"] = "demo";
        }

        string bin_file = model_file + ".bin2";
        FILE* ibin = fopen(bin_file.c_str(), "rb");
        if(!ibin){
            LOGGER.e(0, "can't open file [" + bin_file + "] to read.");
        }
        char magic[5];
        if(fread(magic, sizeof(char), 5, ibin) != 5){
            LOGGER.e(0, "can't read magic number.");
        }
        if(strcmp(magic, "fGLM") != 0){
            LOGGER.e(0, "wrong header in [" + bin_file + "], this file can only generate from GCTA model.");
        }

        uint32_t num_indi2;
        if(fread(&num_indi2, sizeof(uint32_t), 1, ibin) != 1){
            LOGGER.e(0, "can't read sample size.");
        }
        if(num_indi != num_indi2){
            LOGGER.e(0, "wrong number of samples in [" + bin_file + "], this file can only generate from GCTA model.");
        }

        uint32_t covar_col;
        if(fread(&covar_col, sizeof(uint32_t), 1, ibin) != 1){
            LOGGER.e(0, "can't read number of covariates in [" + bin_file + "], this file can only generate from GCTA model.");
        }

        num_covar = covar_col;

        if(fread(&bPreciseCovar, sizeof(bool), 1, ibin) != 1){
            LOGGER.e(0, "can't read covariate flag");
        }

        if(fread(&taoVal, sizeof(double), 1, ibin) != 1){
            LOGGER.e(0, "can't tao value in [" + bin_file + "], this file can only generate from GCTA model.");
        }

        if(fread(&c_inf, sizeof(double), 1, ibin) != 1){
            LOGGER.e(0, "can't c_inf in [" + bin_file + "], this file can only generate from GCTA model.");
        }
        mu.resize(num_indi);
        phenoVec.resize(num_indi);
        covar.resize(num_indi, covar_col);
        H.resize(covar_col, num_indi);

        double *tempVec = new double[num_indi];
        if(fread(tempVec, sizeof(double), num_indi, ibin) != num_indi){
            LOGGER.e(0, "failed to read mu.");
        }
        //switch order;
        for(int i = 0; i < num_indi; i++){
            mu[i] = tempVec[id2[i]];
        }

        if(fread(tempVec, sizeof(double), num_indi, ibin) != num_indi){
            LOGGER.e(0, "failed to read phenotype.");
        }
        //switch order;
        for(int i = 0; i < num_indi; i++){
            phenoVec[i] = tempVec[id2[i]];
        }
        delete[] tempVec;
        uint64_t num_covar_elements = (uint64_t)num_indi * covar_col;
        double *tempCovar = new double[num_covar_elements];
        //covariate
        if(fread(tempCovar, sizeof(double), num_covar_elements, ibin) != num_covar_elements){
            LOGGER.e(0, "failed to read covariates from model.");
        }
        for(uint64_t i = 0; i < covar_col; i++){
            uint64_t cur_base = i * num_indi;
            for(uint64_t j = 0; j < num_indi; j++){
                covar(j, i) = tempCovar[cur_base + id2[j]];
            }
        }
        // H
        if(fread(tempCovar, sizeof(double), num_covar_elements, ibin) != num_covar_elements){
            LOGGER.e(0, "failed to read H from model.");
        }
        for(uint64_t i = 0; i < num_indi; i++){
            uint64_t cur_base = id2[i] * covar_col;
            for(uint64_t j = 0; j < covar_col; j++){
                H(j, i) = tempCovar[cur_base + j];
            }
        }
        delete[] tempCovar;

        initVar();

        fclose(ibin);
        LOGGER << "  loaded successfully." << std::endl;
}

void FastFAM::loadModel(){
        bBinary = false;
        LOGGER << "Loading saved model prefixed with [" << options["model_file"] << "]..." << std::endl;
        LOGGER << "Note: phenotype, covariates, sparse GRM and association test methods are provided by the model, thus these flags will be ignored." << std::endl;
        string model_file = options["model_file"] + ".mdl";

        string id_file = model_file + ".id";
        vector<string> modelIDs = Pheno::read_sublist(id_file);
        num_indi = modelIDs.size();

        uint32_t sample_keep_geno = pheno->count_keep();
        vector<string> sampleIDs = pheno->get_id(0, sample_keep_geno - 1, "\t");
        vector<uint32_t> id1, id2;
        vector_commonIndex_sorted1(sampleIDs, modelIDs, id1, id2);
        if(id2.size() != num_indi){
            LOGGER.e(0, "Some sample ID in model are not existed in genotype!");
        }

        pheno->filter_keep_index(id1);
        LOGGER << "  " << num_indi << " valid individuals to be included." << std::endl;

        MDLHeader head;
        string bin_file = model_file + ".bin";
        FILE* ibin = fopen(bin_file.c_str(), "rb");
        if(!ibin){
            LOGGER.e(0, "can't open file [" + bin_file + "] to read.");
        }
        if(fread(&head, sizeof(MDLHeader), 1, ibin) != 1){
            LOGGER.e(0, "can't read header in [" + bin_file + "].");
        }
        if(strcmp(head.magic, "fGWA") != 0){
            LOGGER.e(0, "wrong header in [" + bin_file + "], this file can only generate from GCTA model.");
        }
        if(head.num_indi != modelIDs.size()){
            LOGGER.e(0, "the sample in model is different from id file!");
        }
        fam_flag = head.fam_flag;
        if(fam_flag){
            options["grmsparse_file"] = "demo";
        }

        bGrammar = head.isGrammar;
        if(fam_flag && bGrammar){
            options["grammar"] = "yes";
        }

        covarFlag = head.covarFlag;
        c_inf = head.c_inf;

        // read the phenotype
        LOGGER << "  loading phenotypes..." << std::endl;
        double *phenoP;
        double *tempP = new double[num_indi];
        if((fam_flag && (!bGrammar)) || (!fam_flag)){
            phenoVec.resize(num_indi);
            phenoP = phenoVec.data();       
        }else{
            Vi_y_cinf.resize(num_indi);
            phenoP = Vi_y_cinf.data();
        }
        fseek(ibin, head.phenoVec_start, SEEK_SET);
        if(fread(tempP, sizeof(double), num_indi, ibin) != num_indi){
            LOGGER.e(0, "failed to read phenotype from model.");
        }
        // switch order;
        for(int i = 0; i < num_indi; i++){
            phenoP[i] = tempP[id2[i]];
        }
        delete[] tempP;

        // read covariates
        LOGGER << "  loading covariates..." << std::endl;
        uint64_t num_covar_elements = (uint64_t)num_indi * head.covarVec_cols;
        fseek(ibin, head.covarVec_start, SEEK_SET);
        this->covar.resize(num_indi, head.covarVec_cols);
        double *tempCovar = new double[num_covar_elements];
        if(fread(tempCovar, sizeof(double), num_covar_elements, ibin) != num_covar_elements){
            LOGGER.e(0, "failed to read covariates from model.");
        }
        for(uint64_t i = 0; i < head.covarVec_cols; i++){
            uint64_t cur_base = i * num_indi;
            for(uint64_t j = 0; j < num_indi; j++){
                covar(j, i) = tempCovar[cur_base + id2[j]];
            }
        }

        if(covarFlag){
            LOGGER << "  loading " << head.covarVec_cols << " covariates..." << std::endl;
            //TODO  H seems incorrect;
            H.resize(num_indi, head.covarVec_cols);
            if(fread(tempCovar, sizeof(double), num_covar_elements, ibin) != num_covar_elements){
                LOGGER.e(0, "failed to read covariates from model.");
            }
            for(uint64_t i = 0; i < head.covarVec_cols; i++){
                uint64_t cur_base = i * num_indi;
                for(uint64_t j = 0; j < num_indi; j++){
                    H(j, i) = tempCovar[cur_base + id2[j]];
                }
            }
        }
        delete[] tempCovar;

        if(fam_flag && (!bGrammar)){
            LOGGER << "  loading V matrix..." << std::endl;
            vector<size_t> lookup = sort_indexes(id2);

            vector<uint32_t> vi, vj;
            vector<double> val;
            vi.resize(head.V_inverse_items);
            vj.resize(head.V_inverse_items);
            val.resize(head.V_inverse_items);

            V_inverse.resize(num_indi, num_indi);

            uint32_t cur_item = 0;
            InvItem item;
            while(cur_item < head.V_inverse_items){
                if(fread(&item, sizeof(item), 1, ibin) == 1){
                    //V_inverse.insert(item.row, item.col) = item.val;
                    vi[cur_item] = lookup[item.row];
                    vj[cur_item] = lookup[item.col];
                    val[cur_item] = item.val; 

                    cur_item++;
                }else{
                    LOGGER.e(0, "failed to read the V inverse from model.");
                }
            }
            //reorder
            vector<size_t> ordIndex = sort_indexes(vj, vi);
            for(int i = 0; i < head.V_inverse_items; i++){
                uint32_t idx = ordIndex[i]; 
                V_inverse.insert(vi[idx], vj[idx]) = val[idx];
            }
        }

        V_inverse.finalize();
        V_inverse.makeCompressed();
        fclose(ibin);
        LOGGER << "  loaded successfully." << std::endl;
}

FastFAM::FastFAM(){
    //Eigen::setNbThreads(THREADS.getThreadCount() + 1);
    //Eigen::setNbThreads(1);
    
    pheno = new Pheno();
    marker = new Marker();
    this->geno = new Geno(pheno, marker);
    hasInfo = geno->getGenoHasInfo();

    num_indi = pheno->count_keep();

    if(options.find("binary") != options.end()){
        bBinary = true;
    }
 
    if(options["model_file"] != ""){
        if(bBinary){
            loadBinModel();
        }else{
            loadModel();
        }
        return;
    }

    if(options.find("grammar") !=options.end()){
        bGrammar = true;
    }


    if(options_d["seed"] == 0){
        seed = pheno->getSeed();
    }else{
        seed = options_d["seed"];
        LOGGER << "Use random seed: " << seed << std::endl;
    }

    double VG;
    double VR;
    bool flag_est_GE = true;
    if(options.find("G") != options.end()){
        VG = std::stod(options["G"]);
        VR = std::stod(options["E"]);
        flag_est_GE = false;
    }

    vector<string> ids;
    pheno->get_pheno(ids, phenos);
    if(ids.size() != num_indi){
        LOGGER.e(0, "Did you forget to specify --pheno?");
    }

    // read covar
    vector<uint32_t> remain_index, remain_index_covar;
    bool has_covar = false;
    Covar covar;
    if(covar.getCommonSampleIndex(ids, remain_index, remain_index_covar)){
        has_covar = true;
        LOGGER.i(0, to_string(remain_index.size()) + " overlapping individuals with non-missing data to be included from the covariate file(s).");
    }else{
        // have covariates, no overlapped
        if(covar.hasCovar()){
            LOGGER.e(0, "No overlappping individual with non-missing data to be included from the covariate file(s).");
        }else{
            // no covariates
            remain_index.resize(ids.size());
            std::iota(remain_index.begin(), remain_index.end(), 0);
        }
    }

    vector<string> remain_ids(remain_index.size());
    std::transform(remain_index.begin(), remain_index.end(), remain_ids.begin(), [&ids](size_t pos){return ids[pos];});

    // read fam
    string ffam_file = "";
    fam_flag = true;
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
    vector<uint32_t> total_remain_index(n_remain_index_fam);
    for(int i = 0; i != n_remain_index_fam; i++){
        uint32_t temp_index = remain_index[remain_index_fam[i]];
        remain_phenos[i] = phenos[temp_index];
        remain_ids_fam[i] = ids[temp_index];
        total_remain_index[i] = temp_index;
    }
    //fix pheno keep
    pheno->filter_keep_index(total_remain_index);
    //geno->init_keep();
    num_indi = pheno->count_keep();
    LOGGER.i(0, "After matching all the files, " + to_string(remain_phenos.size()) + " individuals to be included in the analysis.");


    // standerdize the phenotype, and condition the covar
    phenoVec = Map<VectorXd> (remain_phenos.data(), remain_phenos.size());
    rawPhenoVec = phenoVec;

    // condition the covar
    if(has_covar){
        vector<double> remain_covar;
        vector<uint32_t> remain_inds_index;
        covar.getCovarX(remain_ids_fam, remain_covar, remain_inds_index);

        // get rid of duplicates
        int num_remain = remain_phenos.size();
        int num_remain_col = remain_covar.size() / num_remain;
        vector<int> remain_cols;
        for(int i = 0; i < num_remain_col; i++){
            int baseIndex = i * num_remain;
            double temp = remain_covar[baseIndex]; 
            int num_ident = 0;
            for(int j = 0; j < num_remain; j++){
                if(std::abs(temp - remain_covar[j + baseIndex]) < 1e-8){
                    num_ident++;
                }
            }
            if(num_ident != num_remain){
                remain_cols.push_back(i);
            }
        }

        int num_remain_col2 = remain_cols.size();
        if(num_remain_col2 != num_remain_col){
            vector<double> remain_covar2;
            remain_covar2.reserve(num_remain_col2 * num_remain);
            LOGGER.w(0, "" + to_string(num_remain_col - num_remain_col2) + " covariates which contained identical values were removed from further analyais.");
            for(int i = 0; i < num_remain_col2; i++){
                int index2 = remain_cols[i] * num_remain;
                for(int j = 0; j < num_remain; j++){
                    remain_covar2.push_back(remain_covar[index2 + j]);
                }
            }
            remain_covar = remain_covar2;
        }

        remain_covar.resize(remain_covar.size() + num_remain);
        std::fill(remain_covar.end() - num_remain, remain_covar.end(), 1.0);

        MatrixXd concovar = Map<Matrix<double, Dynamic, Dynamic, Eigen::ColMajor>>(remain_covar.data(), num_remain, remain_covar.size() / num_remain);

        /*
        std::ofstream covar_w(options["out"] + "_aln_covar.txt"), pheno_w(options["out"] + "_aln_phen.txt"), pheno_w2(options["out"] + "_adj_phen.txt");
        pheno_w << phenoVec << std::endl;
        */
        /*
        std::ofstream covar_w(options["out"] + "_aln_covar.txt");
        for(int i = 0; i < phenoVec.size(); i++){
            covar_w << phenoVec(i);
            for(int j = 0; j < concovar.cols(); j++){
                covar_w << "\t" << concovar(i, j);
            }
            covar_w << std::endl;
        }
        covar_w.close();
        */
        /*
        FILE * out = fopen((options["out"] + ".covar.bin").c_str(), "wb");
        fwrite(concovar.data(), sizeof(double), remain_covar.size(), out);
        fclose(out);

        FILE * pout = fopen((options["out"] + ".phenob").c_str(), "wb");
        fwrite(phenoVec.data(), sizeof(double), phenoVec.size(), pout);
        fclose(pout);
        */
        covarFlag = true;
        makeIH(concovar);
        if(!bBinary)conditionCovarReg(phenoVec);
        if(options.find("save_pheno") != options.end()){
            std::ofstream pheno_w((options["out"] + ".cphen").c_str());
            if(!pheno_w) LOGGER.e(0, "failed to write " + options["out"]+".cphen");
            for(int k = 0; k < remain_ids_fam.size(); k++){
                pheno_w << remain_ids_fam[k] << "\t" << phenoVec[k] << std::endl;
            }
            pheno_w.close();
        }
        /*
        FILE * p1out = fopen((options["out"] + ".phenoa1").c_str(), "wb");
        fwrite(phenoVec.data(), sizeof(double), phenoVec.size(), p1out);
        fclose(p1out);
        */


        /*
        pheno_w2 << phenoVec << std::endl;
        covar_w.close();
        pheno_w.close();
        pheno_w2.close();
        */
    }else{
        this->covar.resize(phenoVec.size(), 1);
        this->covar.setZero();
        this->covar = this->covar.array() + 1;
    }

   // joint covar
    if(options.find("adj_covar") != options.end() && covarFlag){
        bPreciseCovar = true;
        LOGGER.i(0, "Fitting covariates jointly in association.");
        covarFlag = true;
    }else if(bBinary){
        covarFlag = true;
    }else{
        covarFlag = false;
        // not usefull now but useful to REML
        this->covar.resize(phenoVec.size(), 1);
        this->covar.setZero();
        this->covar = this->covar.array() + 1;
    }

    if(bBinary){
        initBinary(fam);
        return;
        //goto saveRes;
    }

    // Center
    double phenoVec_mean = phenoVec.mean();
    phenoVec -= VectorXd::Ones(phenoVec.size()) * phenoVec_mean;

    double raw_mean = rawPhenoVec.mean();
    rawPhenoVec = rawPhenoVec.array() - raw_mean;

    double Vpheno = phenoVec.array().square().sum() / (phenoVec.size() - 1);
    if(Vpheno < 1e-5){
        LOGGER.e(0, "Vp < 1e-5, which is not quite right. Please check: 1. Phenotype scale issue? 2. The covariates could explain all Vp; 3. Very rare prevalence of disease.");
    }


    //LOGGER << "DEBUG: samples size: " << geno->pheno->count_keep() << std::endl;
    /*
    FILE * p1out = fopen((options["out"] + ".phenoa2").c_str(), "wb");
    fwrite(phenoVec.data(), sizeof(double), phenoVec.size(), p1out);
    fclose(p1out);
    */

    std::map<string, string> mtdString;
    mtdString["REML"] = "fastGWA-REML (grid search)";
    mtdString["HE"] = "Haseman-Elston regression";

    if(fam_flag){

        if(options.find("inv_file") == options.end()){
            // no saved inv file
            vector<double> Aij;
            vector<double> Zij;

            if(flag_est_GE){
                LOGGER.i(0, "Estimating the genetic variance (Vg) by " + mtdString[options["VgEstMethod"]] + "...");
                LOGGER.ts("HE");
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

                    VG = HEreg(Zij, Aij, fam_flag);
                }else{
                    string vgEstMethod = options["VgEstMethod"];
                    if(vgEstMethod == "HE"){
                        VG = HEreg(fam, phenoVec, fam_flag);
                    }else if(vgEstMethod == "MCREML"){
                        VG = MCREML(fam, phenoVec, fam_flag);
                    }else if (vgEstMethod == "REML"){
                        VG = spREML(fam, phenoVec, fam_flag);
                    }else{
                        LOGGER.e(0, "Unknown method to estimate the Vg");
                    }
                }
                if(options.find("force_gwa") != options.end()){
                    if(!fam_flag){
                        LOGGER.w(0, "  Forced to run fastGWA MLM");
                    }
                    fam_flag = true;
                }

                if(fam_flag){
                    VR = Vpheno - VG;
                    if(VG < 0 && options.find("force_gwa") == options.end()){
                        LOGGER.w(0, "Constrain Vg to 0.");
                        fam_flag = false;
                    }else if(VG > Vpheno && fam_flag){
                        if(options.find("no_constrain") != options.end()){
                            LOGGER.w(0, "Vg larger than Vp");
                        }else{
                            VG = 0.99 * Vpheno;
                            LOGGER.w(0, "Constrain Vg to 0.99 * Vp: " + to_string(VG) + ".");
                        }
                    }
                }
                if(options["VgEstMethod"]=="REML"){
                    LOGGER << "fastGWA-REML runtime: ";
                }else if(options["VgEstMethod"]=="HE"){
                    LOGGER << "HE regression runtime: ";
                }
                LOGGER << LOGGER.tp("HE") << " sec." << std::endl;

            }
            if(!fam_flag){
                LOGGER.w(0, "the estimate of Vg is not statistically significant. "
                        "This is likely because the number of closely related individuals in the sample is not large enough. "
                        "\nIn this case, the program will use linear regression for association test.");
 
                fam.resize(0, 0);
                V_inverse.resize(0, 0);
                if(options.find("save_inv") != options.end()){
                    std::ofstream inv_id((options["out"]+".grm.id").c_str());
                    if(!inv_id) LOGGER.e(0, "failed to write " + options["out"]+".grm.id");
                    inv_id << "--fastGWA" << std::endl;
                    inv_id.close();
                }
                goto saveRes;
            }

            if(!bGrammar){
                // ordinary inverse way.
                inverseFAM(fam, VG, VR);
                if(options.find("save_inv") != options.end()){
                    LOGGER.i(0, "Saving inverse of V for further analysis, use --load-inv for further analysis");
                    std::ofstream inv_id((options["out"]+".grm.id").c_str());
                    if(!inv_id) LOGGER.e(0, "failed to write " + options["out"]+".grm.id");
                    for(int k = 0; k < remain_ids_fam.size(); k++){
                        inv_id << remain_ids_fam[k] << std::endl;
                    }

                    FILE * inv_out = fopen((options["out"]+".grm.inv").c_str(), "wb");
                    InvItem item;
                    for(int k = 0; k < V_inverse.outerSize(); ++k){
                        for(SpMat::InnerIterator it(V_inverse, k); it; ++it){
                            item.row = it.row();
                            item.col = it.col();
                            item.val = it.value();
                            if(fwrite(&item, sizeof(item), 1, inv_out) != 1){
                                LOGGER.e(0, "can't write to [" + options["out"] + ".grm.inv]");
                            }
                        }
                    }
                    fclose(inv_out);

                    LOGGER.i(0, "The inverse has been saved to [" + options["out"] + ".grm.inv]");
                }
            }else{
                //grammar approx.
                //LOGGER.i(0, "Performing GRAMMAR-gamma approximation...");
                grammar(fam, VG, VR);
           }
        }else{
            // use saved inv 
            string id_file = options["inv_file"] + ".grm.id";
            std::ifstream h_id(id_file.c_str());
            if(!h_id){
                LOGGER.e(0, "can't read file [" + id_file + "].");
            }
            string line;
            uint32_t cur_index = 0;
            getline(h_id, line);
            if(line == "--fastGWA") {
                fam_flag = false;
                LOGGER.w(0, "The estimate of Vg is not statistically significant. "
                        "This is likely because the number of relatives is not large enough. "
                        "\nPerforming simple regression via removing --grm-sparse instead...");
 
                fam.resize(0, 0);
                V_inverse.resize(0, 0);
                return;
            }
            // go to the original position
            h_id.clear(); 
            h_id.seekg(0);
            while(getline(h_id, line)){
                if(line != remain_ids_fam[cur_index]){
                    LOGGER.e(0, "samples are not same from line " + to_string(cur_index + 1) + " in [" + id_file + "].");
                }
                cur_index++;
            }
            h_id.close();
            if(cur_index == remain_ids_fam.size()){
                LOGGER.i(0, to_string(cur_index) + " samples are checked identical in inverse V [" + id_file + "].");
            }else{
                LOGGER.e(0, "Empty file or lines not consistent in inverse V [" + id_file + "].");
            }

            V_inverse.resize(phenoVec.size(), phenoVec.size());
            string in_name = options["inv_file"] + ".grm.inv";
            LOGGER.i(0, "Loading inverse of V from " + in_name + "...");
            LOGGER.ts("LOAD_INV");
            FILE * in_file = fopen(in_name.c_str(), "rb");
            if(!in_file){
                LOGGER.e(0, "can't open the file.");
            }
            fseek(in_file, 0L, SEEK_END);
            size_t file_size = ftell(in_file);
            rewind(in_file);

            InvItem item;
            size_t cur_pos = 0;
            while(cur_pos < file_size){
                if(fread(&item, sizeof(item), 1, in_file) == 1){
                    V_inverse.insert(item.row, item.col) = item.val;
                    cur_pos += sizeof(item);
                }else{
                    LOGGER.e(0, "can't read file in pos: " + to_string(cur_pos));
                }
            }
            fclose(in_file);

            V_inverse.finalize();
            V_inverse.makeCompressed();
            LOGGER.i(0, "Inverse of V loaded in " + to_string(LOGGER.tp("LOAD_INV")) + " seconds.");
        }
    }

saveRes:
    if(options.find("model_only") != options.end()){
        LOGGER << "Saving fastGWA model information..." << std::endl;
        string id_file = options["out"] + ".mdl.id";
        std::ofstream inv_id(id_file.c_str());
        if(!inv_id) LOGGER.e(0, "failed to write " + options["out"]+".grm.id");
        vector<string> ids = pheno->get_id(0, num_indi - 1, "\t");
        for(auto & t_id : ids){
            inv_id << t_id << "\n";
        }
        inv_id.close();
        LOGGER << "Sample information has been saved to [" << id_file << "]." << std::endl;

        MDLHeader header;
        header.magic[0] = 'f';
        header.magic[1] = 'G';
        header.magic[2] = 'W';
        header.magic[3] = 'A';
        header.magic[4] = '\0';
        header.num_indi = num_indi;
        header.fam_flag = fam_flag;
        header.isGrammar = bGrammar;
        header.phenoVec_start = 512;
        header.covarFlag = covarFlag;
        header.covarVec_start = header.phenoVec_start + header.num_indi * 8;
        header.covarVec_cols = this->covar.cols();
        uint64_t covar_size = header.covarVec_cols * header.num_indi * 8;
        header.V_inverse_start = header.covarVec_start + covar_size + (header.covarFlag ? covar_size : 0);
        if(header.fam_flag & (!header.isGrammar)) header.V_inverse_items = V_inverse.nonZeros();
        header.c_inf = c_inf;

        string bin_file = options["out"] + ".mdl.bin";
        FILE* obin = fopen(bin_file.c_str(), "wb");
        if(!obin){
            LOGGER.e(0, "can't open " + bin_file + " to write.");
        }

        if(fwrite(&header, sizeof(header), 1, obin) != 1){
            LOGGER.e(0, "can't write header to " + bin_file + ".");
        }

        fseek(obin, header.phenoVec_start, SEEK_SET);
        if(fam_flag && header.isGrammar){
            if(fwrite(Vi_y_cinf.data(), sizeof(double), header.num_indi, obin) != header.num_indi){
                LOGGER.e(0, "can't write phenotype to " + bin_file + ".");
            }
        }else{
            if(fwrite(this->phenoVec.data(), sizeof(double), header.num_indi, obin) != header.num_indi){
                LOGGER.e(0, "can't write phenotype to " + bin_file + ".");
            }
        }

        uint64_t num_covar_elements = (uint64_t)header.num_indi * header.covarVec_cols;
        fseek(obin, header.covarVec_start, SEEK_SET);
        if(fwrite(this->covar.data(), sizeof(double), num_covar_elements, obin)
                != num_covar_elements){
            LOGGER.e(0, "can't write covariates to " + bin_file + ".");
        }

        if(header.covarFlag){
            if(fwrite(this->H.data(), sizeof(double), num_covar_elements, obin)
                    != num_covar_elements){
                LOGGER.e(0, "can't write covariates to " + bin_file + ".");
            }
        }

        fseek(obin, header.V_inverse_start, SEEK_SET);
        if(header.fam_flag && (!header.isGrammar)){
            InvItem item;
            for(int k = 0; k < V_inverse.outerSize(); ++k){
                for(SpMat::InnerIterator it(V_inverse, k); it; ++it){
                    item.row = it.row();
                    item.col = it.col();
                    item.val = it.value();
                    if(fwrite(&item, sizeof(item), 1, obin) != 1){
                        LOGGER.e(0, "can't write V inverse to " + bin_file + ".");
                    }
                }
            }
        }
        fclose(obin);
        LOGGER << "Model has been saved to [" << bin_file << "]." << std::endl;
    }

 }
/*
void FastFAM::initMarkerVars(){
    num_marker = marker->count_extract();
    if(beta)delete[] beta;
    if(se) delete[] se;
    if(p) delete[] p;
    beta = new float[num_marker];
    se = new float[num_marker];
    p = new double[num_marker];
}
*/

FastFAM::~FastFAM(){
    delete pheno;
    delete marker;
    delete geno;
    //if(beta)delete[] beta;
    //if(se) delete[] se;
    //if(p) delete[] p;
    //if(countMarkers) delete[] countMarkers;
    //if(af) delete[] af;
}

void FastFAM::makeIH(MatrixXd &X){
    MatrixXd ctc = X.transpose() * X;
    //MatrixXd cty = covar.transpose() * pheno;
    Eigen::FullPivLU<MatrixXd> lu(ctc);
    if(lu.rank() < ctc.rows()){
        //LOGGER.e(0, "Can't invert the covariate matrix. It might be conlinear problem.\nIf you have confirmed there are no problem, add --ignore-check to go through.");
        LOGGER.w(0, "The covariate matrix has some similar or little information columns, however the results will not be affected mostly.");
    }
    covar = X;
    
    //MatrixXd interm = X * lu.solve(MatrixXd::Identity(X.cols(), X.cols())) * X.transpose();
    //H = X * lu.solve(MatrixXd::Identity(X.cols(), X.cols())); 
    H = lu.solve(X.transpose()); 
}


void FastFAM::conditionCovarReg(VectorXd &y, VectorXd &condPheno){
    if(covarFlag){
        condPheno.noalias() = H * y;
    }else{
        condPheno = y;
    }
}

void FastFAM::conditionCovarBinReg(Eigen::Ref<VectorXd> y){

    //Eigen::SparseVector<double> xvec_sp = y.sparseView(1e-6);
    double *Hy = new double[num_indi];
    const char nT = 'N';
    const double a1 = 1.0;
    const double a2 = -1.0;
    const double b1 = 0;
    const double b2 = 1.0;
    const int incr = 1;
    dgemv(&nT, &num_covar, &numi_indi, &a1, H.data(), &num_covar, y.data(), &incr, &b1, Hy, &incr);
    //VectorXd Hy = H * xvec_sp;
    dgemv(&nT, &numi_indi, &num_covar, &a2, covar.data(), &numi_indi, Hy, &incr, &b2, y.data(), &incr);
    //dgemv(&nT, &numi_indi, &num_covar, &a2, covar.data(), &numi_indi, Hy.data(), &incr, &b2, y.data(), &incr);
    //y = y - covar * (H * y);
    delete[] Hy;
}

void FastFAM::conditionCovarReg(Eigen::Ref<VectorXd> y){
 
    /*Eigen::Ref<VectorXd>
        FILE * out = fopen((options["out"] + ".covar.bin1").c_str(), "wb");
        fwrite(covar.data(), sizeof(double), covar.rows() * covar.cols(), out);
        fclose(out);

        FILE * pout = fopen((options["out"] + ".phenob1").c_str(), "wb");
        fwrite(pheno.data(), sizeof(double), pheno.size(), pout);
        fclose(pout);
        */


    //MatrixXd t_covar = covar.transpose();
    //VectorXd beta = lu.solve(cty);
    //LOGGER.i(0, "DEBUG: condition betas:");
    //LOGGER << beta << std::endl;
    //VectorXd beta = covar.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(pheno);
    if(covarFlag){
        //y = y - H * (covar.transpose() * y);
        y = y - covar * (H * y);
    }
    /*
        FILE * pout1 = fopen((options["out"] + ".phenoa0").c_str(), "wb");
        fwrite(pheno.data(), sizeof(double), pheno.size(), pout1);
        fclose(pout1);
        FILE * pout2 = fopen((options["out"] + ".beta").c_str(), "wb");
        fwrite(beta.data(), sizeof(double), beta.size(), pout2);
        fclose(pout2);
        */



    //double pheno_mean = pheno.mean();
    //pheno -= (VectorXd::Ones(pheno.size())) * pheno_mean;
}

void generateRandom(Ref<MatrixXd> mat){
    uint64_t row = mat.rows();
    uint64_t col = mat.cols();
    
    boost::mt19937 rng( std::chrono::system_clock::now().time_since_epoch().count() + 77 );
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
        randn(rng, boost::normal_distribution<>(0.0, 1.0));
    for(int i = 0; i < col; i++){
        for(int j = 0; j < row; j++){
            mat(j, i) = randn();
        }
    }
}

void FastFAM::genRandY(uint64_t *buf, int num_marker){
    MatrixXd X(num_indi, num_marker);
    double * ptr = X.data();
    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < num_marker; i++){
        geno->makeMarkerX(buf, i, ptr + i * num_indi, true, true);
    }

    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < mcTrails; i++){
        rand_y.col(i) += X * rand_beta.block(finished_rand_marker, i, num_marker, 1);
    }

    finished_rand_marker += num_marker;
}

void FastFAM::estBeta(uint64_t *buf, int num_marker){
    MatrixXd X(num_indi, num_marker);
    double * ptr = X.data();
    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < num_marker; i++){
        geno->makeMarkerX(buf, i, ptr + i * num_indi, true, true);
    }

    phenoEstBeta.block(finished_rand_marker, 0, num_marker, 1) = X.transpose() * VinvY;

    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < mcTrails; i++){
        randEstBeta.block(finished_rand_marker, i, num_marker, 1) = X.transpose() * randVinvY.col(i);
    }

    finished_rand_marker += num_marker;
    
}
//MCREML
double FastFAM::fitREML(double logdet, const SpMat &fam, const VectorXd &pheno){
    SpMat eye(fam.rows(), fam.cols());
    eye.setIdentity();

    uint32_t n = fam.cols();
    double det = exp(logdet);
    double sqrt_det = sqrt(det);

    //solve the V
    SpMat V = fam + det * eye;
    Eigen::SimplicialLDLT<SpMat> solver;
    solver.compute(V);
    if(solver.info() != Eigen::Success){
        LOGGER.e(0, "can't inverse the V");
    }
    
    vector<uint32_t> marker_index = marker->get_extract_index();
    uint32_t m = marker_index.size();

    finished_rand_marker = 0;
    vector<function<void (uint64_t *, int)>> callBacks;
    if(!geno->filterMAF()){
        callBacks.push_back(std::bind(&Geno::freq64, geno, _1, _2));
    }
    callBacks.push_back(std::bind(&FastFAM::genRandY, this, _1, _2));
    geno->loop_64block(marker_index, callBacks, true);
    //geno->resetLoop();

    VectorXd sum_e2_rand(mcTrails); // sum e rand ^ 2
    double det2 = det * det;
    #pragma omp parallel for
    for(int i = 0; i < mcTrails; i++){
        rand_y.col(i) = rand_y.col(i) + sqrt_det * rand_e.col(i);
        randVinvY.col(i) = solver.solve(rand_y.col(i));
        sum_e2_rand(i) = det2 * (randVinvY.col(i).squaredNorm());
    }

    LOGGER << "solve VinvY" << std::endl;
    VinvY = solver.solve(pheno);
    double sum_e2 = det2 * (VinvY.squaredNorm()); 

    finished_rand_marker = 0;
    vector<function<void (uint64_t *, int)>> callBacks2;
    if(!geno->filterMAF()){
        callBacks2.push_back(std::bind(&Geno::freq64, geno, _1, _2));
    }
    callBacks2.push_back(std::bind(&FastFAM::estBeta, this, _1, _2));
    geno->loop_64block(marker_index, callBacks2, true);
    //geno->resetLoop();

    LOGGER << "sum beta" << std::endl;

    VectorXd sum_beta2_rand(mcTrails);
    double invm2 = 1.0 / m / m;
    #pragma omp parallel for
    for(int i = 0; i < mcTrails; i++){
        sum_beta2_rand(i) = invm2 * randEstBeta.col(i).squaredNorm();
    }

    //sum beta 2
    LOGGER << "Get the sum beta" << std::endl;
    double sum_beta2 = invm2 * phenoEstBeta.squaredNorm();

    double f = log( (sum_beta2 / sum_e2) / (sum_beta2_rand.sum() / sum_e2_rand.sum()));
    return f;
}


double FastFAM::MCREML(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig){
    uint64_t n = pheno.size();
    uint64_t m = marker->get_extract_index().size();
    //get random vectors
    mcTrails = 4e9 / n;
    mcTrails = mcTrails > 15 ? 15 : mcTrails;
    mcTrails = mcTrails < 3 ? 3 : mcTrails;

    rand_y.resize(n, mcTrails);
    rand_y.setZero();

    rand_beta.resize(m, mcTrails);
    rand_e.resize(n, mcTrails); 
    generateRandom(rand_beta);
    rand_beta = rand_beta.array() * sqrt(1.0/m);

    generateRandom(rand_e);

    randVinvY.resize(n, mcTrails);
    VinvY.resize(n);

    phenoEstBeta.resize(m, mcTrails);
    randEstBeta.resize(m, mcTrails);

    double hsq1 = 0.25;
    double logdet1 = log((1-hsq1)/hsq1);
    LOGGER << "Fitting MCREML, logdet: " << logdet1 << " h2: " << hsq1 << std::endl;
    double f1 =  fitREML(logdet1, fam, pheno);

    double hsq2;
    if(f1 < 0){
        hsq2 = 0.125;
    }else{
        hsq2 = 0.50;
    }
    double logdet2 = log((1-hsq2)/hsq2);
    LOGGER << "Fitting MCREML, logdet: " << logdet2 << " h2: " << hsq2 << std::endl;
    double f2 = fitREML(logdet2, fam, pheno);

    double prev2_logdet = logdet1, prev_logdet = logdet2;
    double prev2_f = f1, prev_f = f2;
    double f, logdet;
    for(int s = 0; s < 5; s++){
        logdet = (prev2_logdet * prev_f - prev_logdet * prev2_f) / (prev_f - prev2_f);
        if(std::abs(logdet - prev_logdet) < 0.01){
            break;
        }
        LOGGER << "Fitting MCREML, logdet: " << logdet << std::endl;
        f = fitREML(logdet, fam, pheno);
        prev2_logdet = prev_logdet;
        prev_logdet = logdet;
        prev2_f = prev_f;
        prev_f = f;
    }

    double h2 = pheno.dot(VinvY) / n;
    LOGGER << "  Vg = " << h2 << std::endl;
    isSig = true;
    return h2;
}

double FastFAM::HEreg(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig){
    int num_covar = 1;
    int num_component = 1;
    int col_X = num_covar + num_component;
    MatrixXd XtX = MatrixXd::Zero(col_X, col_X);
    VectorXd XtY = VectorXd::Zero(col_X);
    double SSy = 0;

    uint64_t size = fam.cols() * fam.rows();
    XtX(0, 0) = size;

    int n = fam.cols();
    VectorXd SSys(n);
    VectorXd XtY0s(n);
    VectorXd XtY1s(n);
    VectorXd XtX01s(n);
    VectorXd XtX11s(n);
    #pragma omp parallel for
    for(int i = 1; i < n; i++){
        double temp_pheno = pheno[i];
        auto fam_block = fam.block(0, i, i, 1);
        auto pheno_block = pheno.head(i) * temp_pheno;
        SSys[i] = pheno_block.dot(pheno_block);
        XtY0s[i] = pheno_block.sum();
        XtY1s[i] = (pheno_block.transpose() * fam_block)[0];
        XtX01s[i] = fam_block.sum();
        XtX11s[i] = (fam_block.cwiseProduct(fam_block)).sum();
    }
    SSy = SSys.sum();
    XtY[0] = XtY0s.sum();
    XtY[1] = XtY1s.sum();
    XtX(0, 1) = XtX01s.sum();
    XtX(1, 1) = XtX11s.sum();

    //MatrixXd XtXi = XtX.selfadjointView<Eigen::Upper>().inverse();
    XtX(1,0) = XtX(0,1);
    LOGGER << "XtX:" << endl;
    LOGGER << XtX << endl;

    Eigen::FullPivLU<MatrixXd> lu(XtX);
    if(lu.rank() < XtX.rows()){
        LOGGER.w(0, "the XtX matrix is invertable.");
        isSig = false;
        return std::numeric_limits<double>::quiet_NaN();
    }

    MatrixXd XtXi = lu.inverse();

    VectorXd betas = XtXi * XtY;
    LOGGER << "beta:" << endl;
    LOGGER << betas.transpose() << endl;

    double sse = (SSy - betas.dot(XtY)) / (size - col_X);
    LOGGER << "SSE: " << sse << endl;

    VectorXd SDs = sse * XtXi.diagonal();
    LOGGER << "SD: " << SDs.transpose() << endl;

    double hsq = betas[betas.size() - 1];
    double SD = SDs[SDs.size() - 1];
    double Zsq = hsq * hsq / SD;
    double p = StatLib::pchisqd1(Zsq);

    double Vpheno = pheno.array().square().sum() / (pheno.size() - 1);
    //LOGGER.i(2, "Vg = " + to_string(hsq) + ", se = " + to_string(sqrt(SD)) +  ", P = " + to_string(p));
    LOGGER << "\nSource\tVariance\tSE" << std::endl;
    LOGGER << "Vg\t" << hsq << "\t" << sqrt(SD) << std::endl;
    LOGGER << "Ve\t" << Vpheno - hsq << std::endl;
    LOGGER << "Vp\t" << Vpheno << std::endl;
    
    LOGGER << "\nHeritability = " << hsq / Vpheno << " (Pval = " << p << ")" << std::endl;

    if(p > 0.05){
        isSig = false;
    }else{
        isSig = true;
    }
    return hsq;
}

void FastFAM::logLREML(const Ref<const VectorXd> pheno, vector<double> &varcomp, double &logL, double *Hinv){
    int n_comp = varcomp.size();
    int n = A[0].cols();

    SpMat V(n, n);
    for(int j = 0; j < n_comp; j++){
        V += varcomp[j] * A[j];
    }
    //LOGGER << "Non zeros V: " << V.nonZeros() << std::endl;

    Eigen::SimplicialLDLT<SpMat> solverV;
    solverV.compute(V);
    if(solverV.info() != Eigen::Success){
        LOGGER.e(0, "can't inverse the V");
    }

    VectorXd d = solverV.vectorD();
    double logdet_V = d.array().log().sum(); 
    //double logdet_V = solverV.logAbsDeterminant();

    MatrixXd ViX = solverV.solve(covar); // n*c

    MatrixXd XtViX = covar.transpose() * ViX; // c*c
    INVmethod method = INV_FQR;
    double logdet_XtViX;
    int rank;
    if(!SquareMatrixInverse(XtViX, logdet_XtViX, rank, method)){
        LOGGER.e(0, "can't inverse XtViX");
    }

    MatrixXd b_proj_t = ViX * XtViX; //n*c
    VectorXd b2 = b_proj_t.transpose() * pheno; // c*1
    VectorXd resi = pheno - covar * b2; // n*1
    VectorXd Py = solverV.solve(resi);

    logL = -0.5 * (logdet_V + logdet_XtViX + pheno.dot(Py));

    if(Hinv){
        MatrixXd APy(n, n_comp);
        for(int i = 0; i < n_comp; i++){
            APy.col(i) = A[i] * Py;
        }

        MatrixXd Hi(n_comp, n_comp);
        for(int i = 0; i < n_comp; i++){
            VectorXd cur_APy = APy.col(i);
            VectorXd cvec = solverV.solve(cur_APy) - ViX * (b_proj_t.transpose() * cur_APy);
            Hi(i, i) = 0.5 * cur_APy.dot(cvec);
            for(int j = 0; j < i; j++){
                Hi(j, i) = 0.5 * APy.col(j).dot(cvec);
                Hi(i, j) = Hi(j, i);
            }
        }
        INVmethod mtd_hi = INV_FQR;
        double logdet_hi;
        int rank;
        if(!SquareMatrixInverse(Hi, logdet_hi, rank, mtd_hi)){
            LOGGER.e(0, "Hi can't be inverted.");
        }
        memcpy(Hinv, Hi.data(), n_comp * n_comp * sizeof(double));
    }
}


double FastFAM::spREML(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig){
    int n_comp = 2;

    int n = pheno.size();
    int n_covar = covar.cols();

    A.resize(2);
    A[0] = fam;
    A[1].resize(n, n); 
    A[1].setIdentity();

    double Vp = pheno.array().square().sum() / (n - 1);

    //vector<int> trails = {101, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11};
    vector<int> trails = {101, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16};
    double MAX_Vg = Vp * options_d["h2_limit"];
    double start = 0;
    double end = MAX_Vg;
    double end_logL = 0;
    vector<bool> bEndNAN(trails.size());
    vector<double> starts(trails.size());
    vector<double> ends(trails.size());
    vector<double> steps(trails.size());

    bool b_reml = false;
    std::ofstream out;
    if(options.find("reml_detail") != options.end()){
        b_reml = true;
        out.open(options["out"] + ".reml");
        if(!out){
            LOGGER.e(0, "Error to open the " + options["out"] + ".reml");
        }

        out << std::setprecision( std::numeric_limits<double>::digits10);
        out << "index\tlogL\tVg\tVe" << std::endl;

    }

    for(int iter = 0; iter < trails.size(); iter++){
        int n_trails = trails[iter];
        double step = (end - start) / (n_trails - 1);
        vector<double> varcomps(n_trails);
        vector<double> logLs(n_trails);

        #pragma omp parallel for
        for(int i = 0; i < n_trails; i++){
            vector<double> varcomp(n_comp);
            varcomp[0] = start + i * step;
            varcomp[1] = Vp - varcomp[0];
            double logL;
            logLREML(pheno, varcomp, logL);
            varcomps[i] = varcomp[0];
            logLs[i] = logL;
        }

        double max_logL = -1e300;
        double max_varcomp = 0;

        int max_index = 0;
        for(int i = 0; i < n_trails; i++){
            double logL = logLs[i];
            if(logL > max_logL){
                max_logL = logL;
                max_varcomp = varcomps[i];
                max_index = i;
            }
        }
        if(max_index == 0){
            start = varcomps[max_index];
            end = varcomps[max_index + 1];
            end_logL = logLs[max_index + 1];
        }else if(max_index == n_trails - 1){
            start = varcomps[max_index - 1];
            end = varcomps[max_index];
            end_logL = logLs[max_index];
        }else{
            start = varcomps[max_index - 1];
            end = varcomps[max_index + 1];
            end_logL = logLs[max_index + 1];
        }

        if(b_reml){
            for(int i = 0; i < logLs.size(); i++){
                out << i << "\t" << logLs[i] << "\t" << varcomps[i] << "\t" << Vp - varcomps[i] << std::endl;
            }
            out << "----------------------------" << std::endl;
        }

        if(!std::isfinite(end_logL)){
            bEndNAN[iter] = true;
        }else{
            bEndNAN[iter] = false;
        }
        starts[iter] = start;
        ends[iter] = end;
        steps[iter] = step;


        LOGGER << "Iteration " << iter + 1 << ", step size: " << step << ", logL: " << max_logL
               << ". Vg: " << max_varcomp
               << ", searching range: " << start << " to " << end << std::endl;
    }

    if(b_reml){
        out.close();
    }

    double Vg = (start + end) / 2.0;
    double Ve = Vp - Vg;

    if(b_reml){
        LOGGER << "\n" << "Up boundary detail: " << std::endl;
        LOGGER << "  " << "iter\tstep\tstart\tend\tendNAN" << std::endl;
        LOGGER << std::boolalpha;
        for(int iter = 0; iter < trails.size(); iter++){
            LOGGER << "  " << iter << "\t" << steps[iter] << "\t" << starts[iter] << "\t" << ends[iter] << "\t" << bEndNAN[iter] << std::endl; 
        }
        LOGGER << std::endl;
    }

    if(bEndNAN[0] && Ve < 0){
        std::streamsize ss = std::cout.precision();
        LOGGER << "Best guess Vg " 
            << std::setprecision( std::numeric_limits<double>::digits10) 
            << "range: " << start << " to " << end 
            << ", Vp: " << Vp << std::setprecision(ss) << std::endl;
        LOGGER.e(0, "fastGWA-REML can't converge.");
    }

    if(MAX_Vg - end <= 1e-10){
        std::streamsize ss = std::cout.precision();
        LOGGER << "Best guess Vg " 
            << std::setprecision( std::numeric_limits<double>::digits10) 
            << "range: " << start << " to " << end 
            << ", Vp: " << Vp << std::setprecision(ss) << std::endl;

        LOGGER.e(0, "fastGWA-REML can't converge, hit upper limit!");
    }

    LOGGER << "fastGWA-REML converged." << std::endl;

    //std::streamsize ss = std::cout.precision();
    //LOGGER << "Estimated Vg: " << std::setprecision( std::numeric_limits<double>::digits10) << Vg << ", range: " << start << " to " << end << std::setprecision(ss) << std::endl;

    double logL;
    MatrixXd Hinv(n_comp, n_comp);
    vector<double> varcomp(n_comp);
    varcomp[0] = Vg;
    varcomp[1] = Ve;
    logLREML(pheno, varcomp, logL, Hinv.data());

    LOGGER << "logL: " << logL << std::endl;

    LOGGER <<"Sampling variance/covariance of the estimates of Vg and Ve:" << std::endl; 
    LOGGER << Hinv << std::endl;

    LOGGER << "\nSource\tVariance\tSE" << std::endl;
    LOGGER << "Vg" << "\t" << varcomp[0] << "\t" << sqrt(Hinv(0, 0)) << std::endl;
    LOGGER << "Ve" << "\t" << varcomp[1] << "\t" << sqrt(Hinv(1, 1)) << std::endl;
    LOGGER << "Vp" << "\t" << Vp << std::endl;

    for(int i = 0; i < n_comp; i++){
        A[i].resize(0, 0);
    }
    A.resize(0);

    double zsq = Vg * Vg / Hinv(0, 0);
    double p = StatLib::pchisqd1(zsq);

    LOGGER << "\nHeritability = " << Vg / Vp << " (Pval = " << p << ")" << std::endl;
    
    if(p>0.05){
        isSig = false;
    }else{
        isSig = true;
    }

    return Vg;
}

double FastFAM::HEreg(vector<double> &Zij, vector<double> &Aij, bool &isSig){
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
        isSig = false;
    }else{
        isSig = true;
    }

    return hsq;
}
    

void FastFAM::readFAM(string filename, SpMat& fam, const vector<string> &ids, vector<uint32_t> &remain_index){
    LOGGER.i(0, "Reading the sparse GRM file from [" + filename + "]...");
    uint32_t num_indi = ids.size();
    vector<string> sublist = Pheno::read_sublist(filename + ".grm.id");
    vector<uint32_t> fam_index;
    vector_commonIndex(sublist, ids, fam_index, remain_index);
    //LOGGER.i(0, "DEBUG: " + to_string(fam_index.size()) + " subjects remained");

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

    vector<uint32_t> num_elements(remain_index.size(), 0);

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
    /*
    std::ofstream out("test.txt");
    for(int i = 0; i < id1.size(); i++){
        out << id1[i] << "\t" <<  id2[i] << "\t" << grm[i] << std::endl;
    }
    out.close();
    */

    auto sorted_index = sort_indexes(id2, id1);

    fam.resize(ordered_fam_index.size(), ordered_fam_index.size());
    fam.reserve(num_elements);

    for(auto index : sorted_index){
        fam.insertBackUncompressed(id1[index], id2[index]) = grm[index];
    }
    fam.finalize();
    fam.makeCompressed();

}

void FastFAM::grammar_func(uintptr_t *genobuf, const vector<uint32_t> &markerIndex){
    int nMarker = markerIndex.size();
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nMarker; i++){
        int index_cur_marker = num_grammar_markers + i;
        GenoBufItem item;
        item.extractedMarkerIndex = markerIndex[i];
        geno->getGenoDouble(genobuf, i, &item);
        bValids[index_cur_marker] = item.valid;
        if(item.valid){
            Map< VectorXd > curGeno(item.geno.data(), num_indi);
            conditionCovarReg(curGeno);
            VectorXd Vg = solver.solve(curGeno);
            double gt_Vg = curGeno.dot(Vg);
            double g_Vi_y = curGeno.dot(Vi_y);
            double temp_chisq = g_Vi_y * g_Vi_y / gt_Vg;

            v_chisq[index_cur_marker] = temp_chisq;

            if(temp_chisq < 5){
                double gt_g = curGeno.dot(curGeno);
                double tmp_cinf = gt_Vg / gt_g;
                v_c_infs[index_cur_marker] = tmp_cinf;
            }else{
                bValids[index_cur_marker] = false;
            }
        }
    }
    num_grammar_markers += nMarker;
}


void FastFAM::grammar(SpMat& fam, double VG, double VR){
    int num_marker_rand = 1000;

    LOGGER.i(0, "\nTuning parameters using " + to_string(num_marker_rand) + " null SNPs...");

    SpMat eye(fam.rows(), fam.cols());
    eye.setIdentity();

    fam *= VG;
    fam += eye * VR;

    //LOGGER.i(0, "Estimating conjugate gradient...");
    LOGGER.ts("tuning");
    solver.compute(fam);
    if(solver.info() != Eigen::Success){
        LOGGER.e(0, "can't invert the V matrix");
    }
    //LOGGER << "TCG compute time: " << LOGGER.tp("TCG") << std::endl;

    //LOGGER.i(0, "Solving Vi * y via conjugate gradient...");
    //LOGGER.ts("vi_y");
    Vi_y = solver.solve(phenoVec);
    //LOGGER << "  time: " << LOGGER.tp("vi_y") << std::endl;


    // get 1000 random SNPs
    auto total_markers_index = marker->get_extract_index_autosome();
    if(total_markers_index.size() < num_marker_rand){
        LOGGER.e(0, "can't read " + to_string(num_marker_rand) + " SNPs from autosome for tuning.");
    }
    std::default_random_engine generator;
    //generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    uint32_t seed;
    if(options_d["seed"] == 0){
        seed = pheno->getSeed();
    }else{
        seed = options_d["seed"];
        LOGGER << "Use seed: " << seed << std::endl;
    }
    //LOGGER << "Setting seed: " << seed << std::endl;
    generator.seed(seed);
    std::uniform_int_distribution<uint32_t> distribution(0, total_markers_index.size() - 1);
    vector<uint32_t> seq_marker_index(num_marker_rand);

    for(int i = 0; i < num_marker_rand; i++){
        seq_marker_index[i] = distribution(generator);
    }
    std::sort(seq_marker_index.begin(), seq_marker_index.end());
    auto last = std::unique(seq_marker_index.begin(), seq_marker_index.end());
    seq_marker_index.erase(last, seq_marker_index.end());

    num_marker_rand = seq_marker_index.size();
    vector<uint32_t> marker_index(num_marker_rand);
    std::transform(seq_marker_index.begin(), seq_marker_index.end(), marker_index.begin(), [&total_markers_index](size_t pos){return total_markers_index[pos];});

    int nMarker = 100;

    //get previous threshold
    double preAF = geno->getMAF();
    double preInfo = geno->getFilterInfo();
    double preMiss = geno->getFilterMiss();
    if(options.find("c-inf-no-filter") != options.end()){
        if(preAF < 0.01){
            geno->setMAF(0.01);
        }
        if(preInfo < 0.3 && hasInfo){
            geno->setFilterInfo(0.3);
        }
        if(preMiss < 0.9){
            geno->setFilterMiss(0.9);
        }
    }

    v_chisq.resize(num_marker_rand);
    v_c_infs.resize(num_marker_rand);
    bValids.resize(num_marker_rand);
    num_grammar_markers = 0;

    LOGGER << "  reading genotypes..." << std::endl; 
    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    callBacks.push_back(bind(&FastFAM::grammar_func, this, _1, _2));
    geno->loopDouble(marker_index, nMarker, true, true, false, false, callBacks);

    if(num_marker_rand != num_grammar_markers){
        LOGGER.e(0, "some SNPs didn't read successfully!");
    }

    //reset back
    geno->setMAF(preAF);
    geno->setFilterInfo(preInfo);
    geno->setFilterMiss(preMiss);

    
    double tmp_cinf = 0;
    int n_valid_null = 0;

    std::ofstream o_inf;
    bool out_inf = false;
    if(options.find("c-inf") != options.end()){
        o_inf.open(options["out"] + ".cinf");
        o_inf << "CHR\tSNP\tPOS\tA1\tA2\tcinf\tchisq" << std::endl;
        out_inf = true;
    }
    for(int i=0; i<num_marker_rand; i++){
        if(bValids[i]){
            tmp_cinf += v_c_infs[i];
            n_valid_null++;
            if(out_inf){
                o_inf << marker->getMarkerStrExtract(marker_index[i]) << "\t" << v_c_infs[i] << "\t" << v_chisq[i] << std::endl;
            }
        }
    }
    if(out_inf){
        o_inf.close();
    }
    //LOGGER.i(0, "Got " + to_string(n_valid_null) + " null SNPs");

    if(n_valid_null < 100){
        LOGGER.e(0, "Not enough null SNPs (<100).");
    }

    c_inf = tmp_cinf / n_valid_null;
    LOGGER.i(0, "Mean gamma = " + to_string(c_inf));
    LOGGER << "Tuning of gamma finished " << LOGGER.tp("tuning") << " seconds." << std::endl;
    Vi_y_cinf = Vi_y.array() / c_inf;
    v_chisq.resize(0);
    v_c_infs.resize(0);
    bValids.resize(0);

    if(Vi_y_cinf.size() != num_indi){
        LOGGER.e(0, "Inconsistent sample size, there may be some unknown bugs.");
    }
    if(options.find("save_resi") != options.end()){
        std::ofstream oresi((options["out"] + ".fastGWA.residual").c_str());
        std::ofstream ogam((options["out"] + ".fastGWA.gamma").c_str());
        if((!ogam) || (!oresi)){
            LOGGER.e(0, "can't write fastGWA MLM parameters to the file.");
        }
        ogam << c_inf << std::endl;
        ogam.close();

        for(int i = 0; i < num_indi; i++){
            oresi << (pheno->get_id(i,i, "\t")[0]) << "\t" << Vi_y_cinf[i] << "\n";
        }
        oresi.close();
        LOGGER << "Saved grammar gamma residual to [" << options["out"] << ".fastGWA(.residual, .gamma)" << "]." << std::endl;
    }
    /*
    //calculate the c_inf
    VectorXd c_infs(num_marker_tune);
    // shuffle the vector
    vector<uint32_t> est_index(num_marker_rand);
    std::iota(est_index.begin(), est_index.end(), 0);
    std::shuffle(est_index.begin(), est_index.end(), generator);

    int valid_null = 0;
    LOGGER.ts("CG_g");
    while(!est_index.empty() && valid_null < num_marker_tune){
        uint32_t cur_index = est_index.back();
        est_index.pop_back();
        VectorXd curGeno = randGenoBuffer.col(cur_index);

        VectorXd Vg = solver.solve(curGeno);

        double gt_Vg = curGeno.dot(Vg);
        double g_Vi_y = curGeno.dot(Vi_y);
        double temp_chisq = g_Vi_y * g_Vi_y / gt_Vg;

        if(temp_chisq <= 5){
            double gt_g = curGeno.dot(curGeno);
            double tmp_cinf = gt_Vg / gt_g;
            c_infs[valid_null] = tmp_cinf;
            valid_null++;
        }
    }

    LOGGER << "CG genotype time: " << LOGGER.tp("CG_g") << std::endl;

    if(valid_null != num_marker_tune){
        LOGGER.e(0, "Got only " + to_string(valid_null) + " null SNPs. This might due to the huge number of significat SNPs. You can run --fastGWA without --grm-sparse to check");
    }
    c_inf = c_infs.mean();
    LOGGER << "c_infs:" << std::endl;
    LOGGER << c_infs << std::endl;
    LOGGER.i(0, "Mean c_inf = " + to_string(c_inf));
    */
}
/*
void FastFAM::readGenoSample(uint64_t *buf, int num_marker){

    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < num_marker; i++){
        VectorXd curGeno(num_indi);
        geno->makeMarkerX(buf, i, curGeno.data(), true, false);
        conditionCovarReg(curGeno);

        VectorXd Vg = solver.solve(curGeno);

        double gt_Vg = curGeno.dot(Vg);
        double g_Vi_y = curGeno.dot(Vi_y);
        double temp_chisq = g_Vi_y * g_Vi_y / gt_Vg;

        int index_cur_marker = finished_rand_marker + i;
        v_chisq[index_cur_marker] = temp_chisq;

        if(temp_chisq <= 5){
            double gt_g = curGeno.dot(curGeno);
            double tmp_cinf = gt_Vg / gt_g;
            v_c_infs[index_cur_marker] = tmp_cinf;
        }
 
    }

    finished_rand_marker += num_marker;
}
*/

 

void FastFAM::inverseFAM(SpMat& fam, double VG, double VR){
    LOGGER.i(0, "\nInverting the variance-covarinace matrix via " + options["inv_method"] + " (This may take a long time)...");
    //LOGGER.i(0, "DEUBG: Inverse Threads " + to_string(Eigen::nbThreads()));
    LOGGER.ts("INVERSE_FAM");
    SpMat eye(fam.rows(), fam.cols());
    //LOGGER.i(0, "FAM " + to_string(fam.rows()) + " * " + to_string(fam.cols()));
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

        V_inverse = solver.solve(eye); //extreme slow for large matrix
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
    }else if(options["inv_method"] == "inv-t1"){
        /*
           LOGGER.ts("lu");
           Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<long long>> lu_solver;
           lu_solver.analyzePattern(fam);
           lu_solver.factorize(fam);
           LOGGER << " LU logdet: " << lu_solver.logAbsDeterminant() << std::endl;
           LOGGER << " LU time: " << LOGGER.tp("lu") << std::endl; // LU time: 0.620486
           LOGGER.ts("lu_sol");
           LOGGER << lu_solver.solve(phenoVec).sum() << std::endl;
           LOGGER << " LU solve time: " << LOGGER.tp("lu_sol"); // LU solve time: 0.0302136
           */

        //SpMat ori(fam);
        LOGGER.ts("ldlt");
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.compute(fam);
 
        VectorXd d = solver.vectorD();
        //double logdet = 0;
        //for (int i = 0; i < d.size(); i++) logdet += log(d[i]);
        double logdet2 = d.array().log().sum(); 
        double invtrace = (1.0 / d.array()).sum();
        //LOGGER << " logdet = " << logdet << std::endl;
        LOGGER << " logdet2 = " << logdet2 << std::endl;
        LOGGER << " inv trace = " << invtrace << std::endl;

        LOGGER << "  compute LDLT time: " << LOGGER.tp("ldlt") << std::endl; // compute LDLT time: 0.244541
        if(solver.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the FAM");
        }

        LOGGER.ts("ldlt_sol");
        LOGGER << solver.solve(phenoVec).sum() << std::endl;
        LOGGER << "  LDLT solve y time: " << LOGGER.tp("ldlt_sol") << std::endl;//LDLT solve time: 0.0224732

        LOGGER.ts("ldlt_solA");
        //SpMat solvedA = solver.solve(ori);
        //LOGGER << solvedA.diagonal().sum() << std::endl;
        //LOGGER << "  ori non zeros: " << ori.nonZeros() << std::endl;
        //LOGGER << "  solved A non zeros: " << solvedA.nonZeros() << std::endl;
        LOGGER << "  LDLT solve A time: " << LOGGER.tp("ldlt_solA") << std::endl;

        LOGGER << "Spectra" << std::endl;
        LOGGER.ts("sp_solve");
        Spectra::SparseSymMatProd<double, Eigen::Lower, 0, long long> op(fam);
        LOGGER << "1" << std::endl;
        Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE,  Spectra::SparseSymMatProd<double, Eigen::Lower, 0, long long> > eigs(&op, 1000, 2000);
        LOGGER << "2" << std::endl;
        eigs.init();
        LOGGER << "3" << std::endl;
        int nconv = eigs.compute();
        LOGGER << " Eigens compute time: " << LOGGER.tp("sp_solve") << std::endl;

        LOGGER.ts("sp_get");
        Eigen::VectorXd evalues;
        if(eigs.info() == Spectra::SUCCESSFUL)
            evalues = eigs.eigenvalues();
        LOGGER << "logdet from eigen: " << evalues.array().log().sum() << std::endl;
        LOGGER << "trace of fam-1: " <<  (1.0/evalues.array().sum()) << std::endl;
        LOGGER << " eigen get time: " << LOGGER.tp("sp_get") << std::endl;

        LOGGER.e(0, "test finished");

        LOGGER.ts("sol");
        V_inverse = solver.solve(eye); //extreme slow for large matrix
        LOGGER << " inv trace true value = " << V_inverse.diagonal().sum() << std::endl;
        LOGGER << "  inverse time: " << LOGGER.tp("sol") << std::endl;
    }else{
        LOGGER.e(0, "Unknown inverse methods");
    }


    //solver.setTolerance(1e-3);
    ///solver.setMaxIterations(10);;

   //LOGGER.i(0, "# iteations: " + to_string(solver.iterations()));
    //LOGGER.i(0, "# error: " + to_string(solver.error()));

    LOGGER.i(0, "Inverted in " + to_string(LOGGER.tp("INVERSE_FAM")) + " sec.");
}

void FastFAM::calculate_gwa(uintptr_t * genobuf, const vector<uint32_t> &markerIndex){
/*
        std::ofstream pheno_w2(options["out"] + "_gwa_phen.txt");
        pheno_w2 << phenoVec << std::endl;
        pheno_w2.close();
        LOGGER << "LOOP: " << num_finished_marker << std::endl;
*/
        //MatrixXd dealGeno(num_indi, num_marker);
    //const double nan = std::numeric_limits<double>::quiet_NaN();
    //const double MAF_L_THRESH = 0.00001;
    //const double MAF_U_THRESH = 0.99999;
 
    static double iN = 1.0 /(num_indi - (covarFlag ? covar.cols() : 1.0) - 1.0);
    static double SSy = phenoVec.dot(phenoVec);

    int num_marker = markerIndex.size();
    vector<uint8_t> isValids(num_marker);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;

        geno->getGenoDouble(genobuf, i, &item);

        isValids[i] = item.valid;
        if(!item.valid) {
            continue;
        }

        Map< VectorXd > xMat(item.geno.data(), num_indi);
        //MatrixXd XMat_V;

        conditionCovarReg(xMat);
        //dealGeno.col(cur_marker) = xMat;

        //MatrixXd tMat_V = xMat.transpose();

        double xMat_V_x = 1.0 / xMat.dot(xMat);
        double xMat_V_p = xMat.dot(phenoVec);

        double temp_beta =  xMat_V_x * xMat_V_p;
        //uint32_t temp_N = geno->countMarkers[cur_raw_marker] >> 1;
        // slow: VectorXd res = phenoVec - xMat * temp_beta;
        //double sse = res.dot(res) / N;
        double sse = (SSy - temp_beta * xMat_V_p) * iN;
        double temp_se = sqrt(sse * xMat_V_x);


        //double temp_se = sqrt(xMat_V_x);
        double temp_z = temp_beta / temp_se;

        beta[i] = (float)temp_beta; //* geno->RDev[cur_raw_marker]; 
        se[i] = (float)temp_se;
        p[i] = StatLib::pchisqd1(temp_z * temp_z); 

        af[i] = (float)item.af;
        countMarkers[i] = item.nValidN;
        info[i] = item.info; 
        //}else{
            //beta[index] = nan;
            //se[index] = nan;
            //p[index] = nan;
        //}
    }

    output_res(isValids, markerIndex);
   /*
        std::ofstream o_geno(options["out"] + "_geno.txt");
        o_geno << dealGeno << std::endl;
        o_geno.close();
    */
}

void FastFAM::output_res_spa(const vector<uint8_t> &isValids, const vector<uint32_t> markerIndex){
    int numKept = 0;
    int num_marker = markerIndex.size();
    if(bSaveBin){
        int num_write = 1;
        for(int i = 0; i != num_marker; i++){
            if(isValids[i]){
                numKept++;
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\n";
                if(fwrite(&af[i], sizeof(float), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write allele frequency to [" + sFileName + ".bin].");
                }

                if(fwrite(&beta[i], sizeof(float), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write beta to [" + sFileName + ".bin].");
                }
                if(fwrite(&se[i], sizeof(float), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write se to [" + sFileName + ".bin].");
                }
                if(fwrite(&p[i], sizeof(double), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write p to [" + sFileName + ".bin].");
                }
                if(fwrite(&padj[i], sizeof(double), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write Padj to [" + sFileName + ".bin].");
                }

                if(fwrite(&countMarkers[i], sizeof(uint32_t), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write N to [" + sFileName + ".bin].");
                }
                if(hasInfo){
                    if(fwrite(&info[i], sizeof(float), num_write, bOut) != num_write){
                        LOGGER.e(0, "can't write INFO score to [" + sFileName + ".bin].");
                    }
                }


            }
        }
    }else{
        for(int i = 0; i != num_marker; i++){
            if(isValids[i]){
                numKept++;
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t" << countMarkers[i]
                    << "\t" << af[i] << "\t" << Tscore[i] << "\t" << Tse[i] << "\t" << p[i] 
                    << "\t" << beta[i] << "\t" << se[i] << "\t" << padj[i] << "\t" << ((int)rConverge[i]);
                if(hasInfo){
                    osOut << "\t" << info[i];
                }
                osOut << "\n";
            }else{
                if(bOutResAll){
                    numKept++;
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t" << countMarkers[i]
                    << "\t" << af[i] << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                if(hasInfo){
                    osOut << "\t" << info[i];
                }
                osOut << "\n";
 
                }
            }
        }
    }

    numMarkerOutput += numKept;
}
 


void FastFAM::output_res(const vector<uint8_t> &isValids, const vector<uint32_t> markerIndex){
    int numKept = 0;
    int num_marker = markerIndex.size();
    if(bSaveBin){
        int num_write = 1;
        for(int i = 0; i != num_marker; i++){
            if(isValids[i]){
                numKept++;
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\n";
                if(fwrite(&af[i], sizeof(float), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write allele frequency to [" + sFileName + ".bin].");
                }

                if(fwrite(&beta[i], sizeof(float), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write beta to [" + sFileName + ".bin].");
                }
                if(fwrite(&se[i], sizeof(float), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write se to [" + sFileName + ".bin].");
                }
                if(fwrite(&p[i], sizeof(double), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write p to [" + sFileName + ".bin].");
                }

                if(fwrite(&countMarkers[i], sizeof(uint32_t), num_write, bOut) != num_write){
                    LOGGER.e(0, "can't write N to [" + sFileName + ".bin].");
                }
                if(hasInfo){
                    if(fwrite(&info[i], sizeof(float), num_write, bOut) != num_write){
                        LOGGER.e(0, "can't write INFO score to [" + sFileName + ".bin].");
                    }
                }


            }
        }
    }else{
        for(int i = 0; i != num_marker; i++){
            if(isValids[i]){
                numKept++;
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t" << countMarkers[i]
                    << "\t" << af[i] << "\t" << beta[i] << "\t" << se[i] << "\t" << p[i];
                if(hasInfo){
                    osOut << "\t" << info[i];
                }
                osOut << "\n";
            }else{
                if(bOutResAll){
                    numKept++;
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t" << countMarkers[i]
                    << "\t" << af[i] << "\tNA\tNA\tNA";
                if(hasInfo){
                    osOut << "\t" << info[i];
                }
                osOut << "\n";
 
                }
            }
        }
    }

    numMarkerOutput += numKept;
}
 


void FastFAM::calculate_fam(uintptr_t *genobuf, const vector<uint32_t> &markerIndex){

    // Memory fam_size * 2 * 4 + (N * 8 * 2 ) * thread_num + M * 3 * 8  B
    /*
    LOGGER <<"new2" << std::endl;
    FILE *phen = fopen((options["out"] + ".phena").c_str(), "wb");
    fwrite(phenoVec.data(), sizeof(double), phenoVec.size(), phen);
    fclose(phen);
    MatrixXd dealGeno(num_indi, num_marker);
    */
    int num_marker = markerIndex.size();
    vector<uint8_t> isValids(num_marker);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;

        geno->getGenoDouble(genobuf, i, &item);

        isValids[i] = item.valid;
        if(!item.valid) {
            continue;
        }
 
        Map< VectorXd > xMat(item.geno.data(), num_indi);
        //MatrixXd XMat_V;

        conditionCovarReg(xMat);
        //dealGeno.col(cur_marker) = xMat;

        //norm
        /*
        MatrixXd xMat_V = xMat.transpose() * V_inverse;
        double xMat_V_x = 1.0 / (xMat_V * xMat)(0, 0); 
        double xMat_V_p = (xMat_V * phenoVec)(0, 0);
        */
        //end norm
        // Xt * V-1
        //VectorXd xMat_V = V_inverse.selfadjointView<Eigen::Lower>() * xMat;
        VectorXd xMat_V = V_inverse * xMat;
        double xMat_V_x = 1.0 / xMat_V.dot(xMat);
        double xMat_V_p = xMat_V.dot(phenoVec);

        double temp_beta =  xMat_V_x * xMat_V_p;
        double temp_se = sqrt(xMat_V_x);
        double temp_z = temp_beta / temp_se;


        beta[i] = (float)temp_beta; //* geno->RDev[cur_raw_marker]; 
        se[i] = (float)temp_se;
        p[i] = StatLib::pchisqd1(temp_z * temp_z); 

        af[i] = (float)item.af;
        countMarkers[i] = item.nValidN;
        info[i] = item.info;
    }
    output_res(isValids, markerIndex);
}

void FastFAM::calculate_grammar(uintptr_t *genobuf, const vector<uint32_t> &markerIndex){
    int num_marker = markerIndex.size();
    vector<uint8_t> isValids(num_marker);
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;
        geno->getGenoDouble(genobuf, i, &item);
        
        isValids[i] = item.valid;
        if(!item.valid){
            continue;
        }
        
        Map< VectorXd > xvec(item.geno.data(), num_indi);

        conditionCovarReg(xvec);

        double gtg = xvec.dot(xvec);
        //double gt_Vi_y = xvec.dot(Vi_y);
        double gt_Vi_y = xvec.dot(Vi_y_cinf);

        //double temp_beta = gt_Vi_y / gtg / c_inf;
        double temp_beta = gt_Vi_y / gtg;
        //double temp_chisq = temp_beta * gt_Vi_y;
        double temp_chisq = temp_beta * gt_Vi_y * c_inf;
        double temp_se = sqrt(temp_beta * temp_beta / temp_chisq);

        beta[i] = (float)temp_beta; //* geno->RDev[cur_raw_marker]; 
        se[i] = (float)temp_se;
        p[i] = StatLib::pchisqd1(temp_chisq); 
        af[i] = (float)item.af;
        countMarkers[i] = item.nValidN;
        info[i] = item.info;
    }

    output_res(isValids, markerIndex);
}


/*
void FastFAM::calculate_grammar_check(GenoBuf *buf, const vector<uint32_t> &markerIndex){
    
    static uint32_t n_sample = buf->n_sample;
    int num_marker = markerIndex.size();

    #pragma omp parallel for schedule(dynamic)
    for(int cur_marker = 0; cur_marker < num_marker; cur_marker++){
        if(buf->usedIndex[cur_marker]){
            Map< VectorXd > xvec(buf->geno.data() + cur_marker * n_sample, n_sample);

            conditionCovarReg(xvec);

            double gtg = xvec.dot(xvec);
            double gt_Vi_y = xvec.dot(Vi_y);

            double temp_beta = gt_Vi_y / gtg / c_inf;
            double temp_chisq = temp_beta * gt_Vi_y;
            double temp_se = sqrt(temp_beta * temp_beta / temp_chisq);

            uint32_t cur_raw_marker = num_finished_marker + cur_marker;

            beta[cur_raw_marker] = temp_beta; / geno->RDev[cur_raw_marker]; 
            se[cur_raw_marker] = temp_se;
            p[cur_raw_marker] = StatLib::pchisqd1(temp_chisq); 
            af[cur_marker] = buf->af[cur_marker];
            countMarkers[cur_marker] = buf->nValidN[cur_marker];
        }
    }
    if(bSaveBin){
        for(int i = 0; i != num_marker; i++){
            osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\n";
        }
        
        if(fwrite(af, sizeof(float), num_marker, bOut) != num_marker){
            LOGGER.e(0, "can't write allele frequency to [" + sFileName + ".bin].");
        }

        if(fwrite(beta, sizeof(float), num_marker, bOut) != num_marker){
            LOGGER.e(0, "can't write beta to [" + sFileName + ".bin].");
        }
         if(fwrite(se, sizeof(float), num_marker, bOut) != num_marker){
            LOGGER.e(0, "can't write se to [" + sFileName + ".bin].");
        }
        if(fwrite(p, sizeof(double), num_marker, bOut) != num_marker){
            LOGGER.e(0, "can't write p to [" + sFileName + ".bin].");
        }

        if(fwrite(countMarkers, sizeof(uint32_t), num_marker, bOut) != num_marker){
            LOGGER.e(0, "can't write N to [" + sFileName + ".bin].");
        }
    }else{
        for(int i = 0; i != num_marker; i++){
            if(buf->usedIndex[i]){
                osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t" << countMarkers[i]
                    << "\t" << af[i] << "\t" << beta[i] << "\t" << se[i] << "\t" << p[i] << "\n";
            }
        }
    }
}
*/

void FastFAM::processFAM(vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks){
    sFileName = options["out"]; 
    int buf_size = 23068672;
    osBuf.resize(buf_size);
    osOut.rdbuf()->pubsetbuf(&osBuf[0], buf_size);
    if(options.find("save_bin") == options.end()){
        bSaveBin = false;
        LOGGER << "fastGWA results will be saved in text format to [" << sFileName << "]." << std::endl;
        osOut.open(sFileName.c_str());
        vector<string> header = {"CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P"};
        if(hasInfo)header.push_back("INFO");
        if(bBinary){
            header = {"CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "T", "SE_T", "P_noSPA", "BETA", "SE",  "P",  "CONVERGE"};
            if(hasInfo)header.push_back("INFO");
        }
        string header_string = boost::algorithm::join(header, "\t");
        if(osOut.bad()){
            LOGGER.e(0, "can't open [" + sFileName + "] to write.");
        }
        osOut << header_string << std::endl;
    }else{
        bSaveBin = true;
        LOGGER << "fastGWA results will be saved in binary format to [" << sFileName << "(.snpinfo, .bin)]" << std::endl;
        osOut.open((sFileName + ".snpinfo").c_str());
        if(osOut.bad()){
            LOGGER.e(0, "can't open [" + sFileName + ".snpinfo] to write.");
        }
        osOut << "CHR\tSNP\tPOS\tA1\tA2" << std::endl;
        bOut = fopen((sFileName + ".bin").c_str(), "wb");
        if(bOut == NULL){
            LOGGER.e(0, "can't open [" + sFileName + ".bin] to write.");
        }
    }


    if(options.find("no_filter") == options.end()){
        bOutResAll = false;
        double preAF = geno->getMAF();
        double preInfo = geno->getFilterInfo();
        double preMiss = geno->getFilterMiss();
        if(preAF < 1e-10){
            geno->setMAF(0.0001);
            LOGGER << "  Filtering out variants with MAF < 0.0001, customise it with --maf flag." << std::endl; 
        }
        /*
        if(preInfo < 1e-10 && hasInfo){
            geno->setFilterInfo(0.3);
            LOGGER << "  Filtering out variants with imputation INFO score < 0.30, customise it with --info flag." << std::endl;
        }
        */
        if(preMiss < 1e-10){
            geno->setFilterMiss(0.9);
            LOGGER << "  Filtering out variants with missingness rate > 0.10, customise it with --geno flag." << std::endl;
        }
    }else{
        bOutResAll = true;
    }

    vector<uint32_t> extractIndex(marker->count_extract());
    std::iota(extractIndex.begin(), extractIndex.end(), 0);

    int nMarker = 1024;
 
    beta = new float[nMarker];
    se = new float[nMarker];
    p = new double[nMarker];
   
    countMarkers = new uint32_t[nMarker];
    af = new float[nMarker];
    info = new float[nMarker];

    numMarkerOutput = 0;

    bool bCenter = true;
    if(bBinary){
        Tscore = new float[nMarker];
        Tse = new float[nMarker];
        padj = new double[nMarker];
        rConverge = new uint8_t[nMarker];
    }
    geno->loopDouble(extractIndex, nMarker, true, bCenter, false, false, callBacks);

    osOut.flush();
    osOut.close();
    if(bOut){
        fflush(bOut);
        fclose(bOut);
    }
    LOGGER << "Saved " << numMarkerOutput << " SNPs." << std::endl;

    delete[] beta;
    delete[] se;
    delete[] p;
    delete[] countMarkers;
    delete[] af;
    delete[] info;
    if(bBinary){
        delete[] padj;
        delete[] rConverge;
        delete[] Tscore;
        delete[] Tse;
    }

}

/*
void FastFAM::processFAM(vector<function<void (GenoBuf *, const vector<uint32_t> &)>> callBacks){
    sFileName = options["out"]; 
    int buf_size = 23068672;
    osBuf.resize(buf_size);
    osOut.rdbuf()->pubsetbuf(&osBuf[0], buf_size);
    if(options.find("save_bin") == options.end()){
        bSaveBin = false;
        osOut.open(sFileName.c_str());
        vector<string> header{"CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P"};
        string header_string = boost::algorithm::join(header, "\t");
        if(osOut.bad()){
            LOGGER.e(0, "can't open [" + sFileName + "] to write.");
        }
        osOut << header_string << std::endl;
    }else{
        bSaveBin = true;
        osOut.open((sFileName + ".snpinfo").c_str());
        if(osOut.bad()){
            LOGGER.e(0, "can't open [" + sFileName + ".snpinfo] to write.");
        }
        osOut << "CHR\tSNP\tPOS\tA1\tA2" << std::endl;
        bOut = fopen((sFileName + ".bin").c_str(), "wb");
        if(bOut == NULL){
            LOGGER.e(0, "can't open [" + sFileName + ".bin] to write.");
        }
    }

    numMarkerOutput = 0;
    vector<uint32_t> extractIndex(marker->count_extract());
    std::iota(extractIndex.begin(), extractIndex.end(), 0);
    GenoBuf gbuf;
    gbuf.center = true;
    gbuf.std = false;
    gbuf.saveGeno = true;
    gbuf.saveMiss = false;
    int nMarker = 128;

    vector<uint32_t> rawIndex = marker->get_extract_index();

    geno->preProcess(&gbuf, nMarker, &rawIndex);
 
    double preAF = geno->getMAF();
    double preInfo = geno->getFilterInfo();
    double preMiss = geno->getFilterMiss();
    if(preAF < 1e-10){
        geno->setMAF(0.0001);
        LOGGER << "  Filtering out variants with MAF < 0.0001, customise it with --maf flag." << std::endl; 
    }
    if(preInfo < 1e-10 && gbuf.hasInfo){
        geno->setFilterInfo(0.3);
        LOGGER << "  Filtering out variants with imputation INFO score < 0.30, customise it with --info flag." << std::endl;
    }
    if(preMiss < 1e-10){
        geno->setFilterMiss(0.9);
        LOGGER << "  Filtering out variants with missingness rate > 0.10, customise it with --geno flag." << std::endl;
    }

    beta = new float[nMarker];
    se = new float[nMarker];
    p = new double[nMarker];
    countMarkers = new uint32_t[nMarker];
    af = new float[nMarker];
    */

    /* use no extract no in code
    int curMarker = 0;
    uint32_t nTMarker = extractIndex.size();
    LOGGER.ts("FAM");
    while(curMarker < nTMarker){
        int endMarker = curMarker + nMarker;
        endMarker = endMarker > nTMarker ? nTMarker : endMarker;
        vector<uint32_t> curMarkers(extractIndex.begin() + curMarker, extractIndex.begin() + endMarker);
        geno->getGenoArrayExtract(curMarkers, &gbuf);

        for(auto callback : callBacks){
            callback(&gbuf, curMarkers);
        }
        curMarker = endMarker;
    }
    LOGGER << "Time: " << LOGGER.tp("FAM") << " sec"<< std::endl;
    */
 
/*
    geno->loopDoublePre(extractIndex, &gbuf, callBacks);
    osOut.flush();
    osOut.close();
    if(bOut){
        fflush(bOut);
        fclose(bOut);
    }
    LOGGER << "Saved " << numMarkerOutput << " SNPs." << std::endl;
    geno->endProcess();

    delete[] beta;
    delete[] se;
    delete[] p;
    delete[] countMarkers;
    delete[] af;
}
*/


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

        beta[cur_raw_marker] = temp_beta; // geno->RDev[cur_raw_marker]; 
        se[cur_raw_marker] = temp_se;
        {
            std::lock_guard<std::mutex> lock(chisq_lock);
            p[cur_raw_marker] = StatLib::pchisqd1(temp_z * temp_z); 
        } 
    }
    delete[] w_buf;
}
*/

/*
void FastFAM::output(string filename){
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double MAF_L_THRESH = 0.00001;
    const double MAF_U_THRESH = 0.99999;
    //const double qchisq05 = 0.4549364;
    if(options.find("save_bin") == options.end()){
        std::ofstream out(filename.c_str());
        if(options.find("no_marker") == options.end()){
            vector<string> header{"CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P"};
            //std::copy(header.begin(), header.end(), std::ostream_iterator<string>(out, "\t"));
            string header_string = boost::algorithm::join(header, "\t");
            out << header_string << std::endl;
            for(int index = 0; index != num_marker; index++){
                double af = geno->AFA1[index];
                uint32_t real_index = marker->getRawIndex(index);
                uint32_t countSample = geno->countMarkers[index] >> 1;
                if(af > MAF_L_THRESH && af < MAF_U_THRESH){
                    out << marker->get_marker(real_index) << "\t" << countSample << "\t" << 
                        af << "\t" << beta[index] << "\t" << se[index] << "\t" << p[index] << std::endl;
                }else{
                    out << marker->get_marker(marker->getRawIndex(index)) << "\t" << countSample << "\t"<<
                        af << "\t" << nan << "\t" << nan << "\t" << nan << std::endl;
                }
            }
        }else{
            vector<string> header{"N", "AF1", "BETA", "SE", "P"};
            //std::copy(header.begin(), header.end(), std::ostream_iterator<string>(out, "\t"));
            string header_string = boost::algorithm::join(header, "\t");
            out << header_string << std::endl;
            for(int index = 0; index != num_marker; index++){
                double af = geno->AFA1[index];
                uint32_t countSample = geno->countMarkers[index] >> 1;
 
                if(af > MAF_L_THRESH && af < MAF_U_THRESH){
                    out << countSample << "\t" << af << "\t" << beta[index] << "\t" << se[index] << "\t" << p[index] << std::endl;
                }else{
                    out << countSample << "\t" << af << "\t" << nan << "\t" << nan << "\t" << nan << std::endl;
                }
            }
            LOGGER.i(0, "No SNP information saved, " + to_string(num_marker) + " SNPs saved");
        }
        out.close();
        LOGGER.i(0, "The association results have been saved in [" + filename +"].");
    }else{
        if(options.find("no_marker") == options.end()){
            std::ofstream out((filename + ".snp").c_str());
            for(int index = 0; index != num_marker; index++){
                out << marker->get_marker(marker->getRawIndex(index)) << std::endl;
            }
            out.close();
            LOGGER.i(0, "The SNP inf of association results has been saved to [" + filename + ".snp].");
        }else{
            LOGGER.i(0, "No SNP information saved, " + to_string(num_marker) + " SNPs saved");
        }

        float * afa1 = new float[num_marker];
        uint32_t * countSamples = new uint32_t[num_marker];
        for(int index = 0; index != num_marker; index++){
            double af = geno->AFA1[index];
            afa1[index] = af;
            if(af < MAF_L_THRESH || af > MAF_U_THRESH){
                beta[index] = nan;
                se[index] = nan;
                p[index] = nan;
            }
            countSamples[index] = geno->countMarkers[index] >> 1;
        }

        FILE * h_out = fopen((filename + ".bin").c_str(), "wb");
        if(!h_out){LOGGER.e(0,  "can't open [" + filename + ".bin] to write.");}
        if(fwrite(afa1, sizeof(float), num_marker, h_out) != num_marker){
            LOGGER.e(0, "can't write AF to [" + filename + ".bin].");
        }
        delete[] afa1;

        if(fwrite(beta, sizeof(float), num_marker, h_out) != num_marker){
            LOGGER.e(0, "can't write beta to [" + filename + ".bin].");
        }
         if(fwrite(se, sizeof(float), num_marker, h_out) != num_marker){
            LOGGER.e(0, "can't write se to [" + filename + ".bin].");
        }
        if(fwrite(p, sizeof(double), num_marker, h_out) != num_marker){
            LOGGER.e(0, "can't write p to [" + filename + ".bin].");
        }

        if(fwrite(countSamples, sizeof(uint32_t), num_marker, h_out) != num_marker){
            LOGGER.e(0, "can't write N to [" + filename + ".bin].");
        }
        delete[] countSamples;
        fclose(h_out);
        LOGGER.i(0, "The association results have been saved to [" + filename + ".bin] in binary format.");
    }
        
}
*/

int FastFAM::registerOption(map<string, vector<string>>& options_in){
    int returnValue = 0;
    //DEBUG: change to .fastFAM
    options["out"] = options_in["out"][0] + ".fastGWA";

    string curFlag = "--fastGWA";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        returnValue++;
        options_in.erase(curFlag);
        LOGGER.e(0, "Obsoleted flag, try --fastGWA-mlm or --fastGWA-mlm-exact");
    }

    curFlag = "--fastGWA-gram";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        options["grammar"] = "true";
        returnValue++;
        options_in.erase(curFlag);
        LOGGER.e(0, "Obsoleted flag, try --fastGWA-mlm or --fastGWA-mlm-exact");
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

    curFlag = "--fastGWA-mlm";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        options["grammar"] = "true";
        if(options.find("grmsparse_file") == options.end()){
            LOGGER.e(0, "--fastGWA-mlm must run with --grm-sparse");
        }
        returnValue++;
        options_in.erase(curFlag);
    }

    curFlag = "--fastGWA-mlm-binary";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        options["binary"] = "true";
        if(options.find("grmsparse_file") == options.end()){
            LOGGER.e(0, "--fastGWA-mlm-binary must run with --grm-sparse");
        }
        returnValue++;
        options_in.erase(curFlag);
    }
    
    curFlag = "--burden";
    if(options_in.find(curFlag) != options_in.end()){
        if(options.find("binary") != options.end()){
            LOGGER.e(0, "burden test can only work in binary trait only currently!");
        }else{
            if(options_in[curFlag].size() != 1){
                options["genset"] = options_in[curFlag][0];
            }
        }
    }


    curFlag = "--save-fastGWA-mlm-residual";
    if(options_in.find(curFlag) != options_in.end()){
        if(options.find("grammar") == options.end()){
            LOGGER.e(0, curFlag + " can only work with --fastGWA-mlm");
        }
        options["save_resi"] = "true";
        options_in.erase(curFlag);
    }

    curFlag = "--fastGWA-mlm-exact";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        if(options.find("grmsparse_file") == options.end()){
            LOGGER.e(0, "--fastGWA-mlm-exact must run with --grm-sparse");
        }
        returnValue++;
        options_in.erase(curFlag);
    }

    curFlag = "--fastGWA-lr";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("fast_fam");
        if(options.find("grmsparse_file") != options.end()){
            LOGGER.e(0, "--fastGWA-lr can't run with --grm-sparse currently");
        }
        returnValue++;
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
    vector<string> flags = {"--cg", "--ldlt", "--llt", "--pardiso", "--tcg", "--lscg", "--inv-t1"};
    for(auto curFlag : flags){
        if(options_in.find(curFlag) != options_in.end()){
            boost::erase_all(curFlag, "--");
            options["inv_method"] = curFlag;
            options_in.erase(curFlag);
        }
    }

    options["VgEstMethod"] = "REML";
    curFlag = "--est-vg";
    if(options_in.find(curFlag) != options_in.end()){
        auto cur_option = options_in[curFlag];
        if(cur_option.size() >=1){
            options["VgEstMethod"] = cur_option[0];
        }
    }

    curFlag = "--save-inv";
    if(options_in.find(curFlag) != options_in.end()){
        options["save_inv"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--save-bin";
    if(options_in.find(curFlag) != options_in.end()){
        options["save_bin"] = "yes";
        //options_in.erase(curFlag);
    }

    curFlag = "--no-marker";
    if(options_in.find(curFlag) != options_in.end()){
        options["no_marker"] = "yes";
        //options_in.erase(curFlag);
    }

    curFlag = "--nofilter";
    if(options_in.find(curFlag) != options_in.end()){
        options["no_filter"] = "yes";
        //options_in.erase(curFlag);
    }


    curFlag = "--c-inf";
    if(options_in.find(curFlag) != options_in.end()){
        options["c-inf"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--c-inf-no-filter";
    if(options_in.find(curFlag) != options_in.end()){
        options["c-inf-no-filter"] = "yes";
        options_in.erase(curFlag);
    }



    curFlag = "--load-inv";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["inv_file"] = options_in[curFlag][0];
        }else{
            LOGGER.e(0, "can't load multiple --load-inv files");
        }
        options_in.erase(curFlag);
    }

    curFlag = "--model-only";
    bool model_only = false;
    if(options_in.find(curFlag) != options_in.end()){
        options["model_only"] = "yes";
        model_only = true;
        options_in.erase(curFlag);
    }

    curFlag = "--load-model";
    options["model_file"] = "";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            options["model_file"] = options_in[curFlag][0];
            if(!checkFileReadable(options["model_file"] + ".mdl.id")){
                LOGGER.e(0, "can't read " + options["model_file"] + ".mdl.id");
            }
            if(checkFileReadable(options["model_file"] + ".mdl.bin")){
                processFunctions.push_back("fast_fam");
            }else if(checkFileReadable(options["model_file"] + ".mdl.bin2")){
                processFunctions.push_back("fast_fam");
                options["binary"] = "yes";
            }else{
                LOGGER.e(0, "can't read the model binary file.");
            }
            returnValue++;
        }else{
            LOGGER.e(0, "can't load multiple --load-model files");
        }
       options_in.erase(curFlag);
    }

    if(model_only && options["model_file"] != ""){
        LOGGER.e(0, "can't model and load model at the same time");
    }


    curFlag = "--save-pheno";
    if(options_in.find(curFlag) != options_in.end()){
        options["save_pheno"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--joint-covar";
    if(options_in.find(curFlag) != options_in.end()){
        options["adj_covar"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--force-gwa";
    if(options_in.find(curFlag) != options_in.end()){
        options["force_gwa"] = "yes";
        options_in.erase(curFlag);
    }

    // no contrain by default
    options["no_constrain"] = "yes";
    curFlag = "--gwa-no-constrain";
    if(options_in.find(curFlag) != options_in.end()){
        options["no_constrain"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--reml-detail";
    if(options_in.find(curFlag) != options_in.end()){
        options["reml_detail"] = "yes";
        options_in.erase(curFlag);
    }

    curFlag = "--h2-limit";
    if(options_in.find(curFlag) != options_in.end()){
        auto cur_option = options_in[curFlag];
        if(cur_option.size() >=1){
            options_d["h2_limit"] = stod(cur_option[0]);
        }
    }else{
        //default value 1.6
        options_d["h2_limit"] = 1.6;
    }

    //random seed
    curFlag = "--seed";
    if(options_in.find(curFlag) != options_in.end()){
        auto cur_option = options_in[curFlag];
        if(cur_option.size() >=1){
            options_d["seed"] = stod(cur_option[0]);
        }
    }else{
        options_d["seed"] = 0;
    }


    curFlag = "--rel-only";
    if(options_in.find(curFlag) != options_in.end()){
        options["rel_only"] = "yes";
        options_in.erase(curFlag);
    }else{
        options["rel_only"] = "no";
    }

    curFlag = "--cv-threshold";
    if(options_in.find(curFlag) != options_in.end()){
        auto cur_option = options_in[curFlag];
        if(cur_option.size() >=1){
            options_d["cv_threshold"] = stod(cur_option[0]);
        }
    }else{
        options_d["cv_threshold"] = 0.1;
    }

    curFlag = "--tao-start";
    if(options_in.find(curFlag) != options_in.end()){
        auto cur_option = options_in[curFlag];
        if(cur_option.size() >=1){
            options_d["tao_start"] = stod(cur_option[0]);
        }
    }else{
        options_d["tao_start"] = 0.0;
    }


    return returnValue;
}


void FastFAM::processMain(){
    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    //THREADS.JoinAll();
    for(auto &process_function : processFunctions){
        if(process_function == "fast_fam"){
            FastFAM ffam;
            if(options.find("save_inv") != options.end()){
                LOGGER.i(0, "Use --load-inv to load the inversed file for fastGWA");
                return;
            }
            if(options.find("model_only") != options.end()){
                LOGGER.i(0, "Use \"--load-model " + options["out"] + "\" to load the model for further analysis.");
                LOGGER << "Note: the sample ID (FID IID) in model has to exsit in genotype for further association step." << std::endl;
                return;
            }
            bool bBinary = false;
            if(options.find("binary") != options.end()){
                bBinary = true;
            }
            //Eigen::setNbThreads(1);
            if(options.find("grmsparse_file") != options.end() && ffam.fam_flag){
                if(bBinary){
                    LOGGER.i(0, "\nPerforming fastGWA generalized linear mixed model association analysis...");
                    callBacks.push_back(bind(&FastFAM::calculate_spa, &ffam, _1, _2));
                }else{
                    if(options.find("grammar") == options.end()){
                        LOGGER.i(0, "\nPerforming fastGWA mixed model association analysis (extact test)...");
                        callBacks.push_back(bind(&FastFAM::calculate_fam, &ffam, _1, _2));
                    }else{
                        LOGGER.i(0, "\nPerforming fastGWA mixed model association analysis...");
                        callBacks.push_back(bind(&FastFAM::calculate_grammar, &ffam, _1, _2));
                    }
                }
            }else{
                LOGGER.i(0, "\nPerforming fastGWA linear regression analysis...");
                callBacks.push_back(bind(&FastFAM::calculate_gwa, &ffam, _1, _2));
            }
            ffam.processFAM(callBacks);
            callBacks.clear();
        }
    }
}

bool FastFAM::covarGLM(const VectorXd& phenoVec, const MatrixXd& covar, Ref<VectorXd> est_beta, int maxIter, double thresh){
    bool bConverge = false;
    int curIter = 0;
    /*
    std::ofstream out("test.txt");
    out << covar << std::endl;
    out.close();
    std::ofstream out2("pheno.txt");
    out2 << phenoVec << std::endl;
    out2.close();
    */
    while(curIter < maxIter){
        VectorXd theta = covar * est_beta;
        VectorXd theta_exp = theta.array().exp();
        
        VectorXd mu = theta_exp.array() / (theta_exp.array() + 1);
        VectorXd diags = mu.array() * (-mu.array() + 1);

        // covar.transpose().array().colwise() * diags.array()
        MatrixXd XtWX =  covar.transpose() * diags.asDiagonal() * covar;
        VectorXd comp2 = covar.transpose() * (phenoVec - mu);
        auto XtWX_solver = XtWX.colPivHouseholderQr();
        if(!XtWX_solver.isInvertible()){
            LOGGER.e(0, "XtWX is not invertable!");
        }
        VectorXd comp1_2 = XtWX_solver.solve(comp2);

        VectorXd new_beta = est_beta + comp1_2;
        VectorXd mChange = (est_beta - new_beta).array().abs() / (est_beta + new_beta).array().abs();
        double cur_thresh = mChange.maxCoeff();
        est_beta = new_beta;
        if(cur_thresh < thresh){
            bConverge = true;
            break;
        }
    }
    return bConverge;
}

double varVector(const Ref<VectorXd> vec){
    VectorXd vec_mean = vec.array() - vec.mean();
    return(vec_mean.dot(vec_mean) / (vec.size() - 1));
}

bool checkNAN(const MatrixXd &mat){
    return (!mat.array().isFinite()).any();
}


double FastFAM::binLogL(double cur_tao, const SpMat& fam, const SpMat& W, const Ref<VectorXd> Y, const Ref<MatrixXd> X){
    SpMat V = W + cur_tao * fam;
    Eigen::SimplicialLDLT<SpMat> solverV;
    solverV.compute(V);
    if(solverV.info() != Eigen::Success){
        LOGGER.e(0, "can't inverse the V matrix!");
    }

    MatrixXd ViX = solverV.solve(X); // n*c
    if(solverV.info() != Eigen::Success){
        LOGGER.e(0, "can't get the ViX matrix!");
    }

    VectorXd d = solverV.vectorD();
    double logdet_V = d.array().log().sum(); 

    auto XtVX_solver = (X.transpose() * ViX).colPivHouseholderQr();
    if(!XtVX_solver.isInvertible()){
        LOGGER.e(0, "XtViX is not invertable!");
    }
    double logdet_XtVX = XtVX_solver.logAbsDeterminant();

    MatrixXd inv_XtVX_ViX = XtVX_solver.solve(ViX.transpose()); // c*n
    VectorXd PY = solverV.solve(Y) - ViX * (inv_XtVX_ViX * Y);

    return(-0.5 * (logdet_V + logdet_XtVX + Y.dot(PY)));

}

bool FastFAM::binGridREML(const SpMat& fam, Ref<VectorXd> est_a, int maxIter, double threshold){
    LOGGER << "Performing GLM to get the initial guess of beta..." << std::endl;
    if(!covarGLM(phenoVec, covar, est_a, maxIter, threshold)){
        LOGGER.w(0, "GLM didn't reach convergence");
    }
    LOGGER << "GLM finished, fixed effects: " << est_a.transpose() << std::endl;


    SpMat W(num_indi, num_indi);
    W.setIdentity();

   //init first run
    VectorXd Xa_b = covar * est_a;
    VectorXd Xa_b_exp = Xa_b.array().exp();
    mu = Xa_b_exp.array() / (Xa_b_exp.array() + 1);
    VectorXd var_mu_i = 1.0 / (mu.array() * (-mu.array() + 1));
    // working vector
    VectorXd Y = Xa_b.array() + var_mu_i.array() * (phenoVec - mu).array();
    LOGGER << "Init Var(Y): " << varVector(Y) << std::endl;
    
    double cur_tao = options_d["tao_start"] * varVector(Y); 

    {
        // test the tao val
    
        W.diagonal() = var_mu_i;
        SpMat V = W + cur_tao * fam;
        solverV.compute(V);

        if(solverV.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the V matrix!");
        }

        ViX = solverV.solve(covar); // n*c

        auto XtVX_solver = (covar.transpose() * ViX).colPivHouseholderQr();
        if(!XtVX_solver.isInvertible()){
            LOGGER.e(0, "XtViX is not invertable!");
        }

        inv_XtVX_ViX = XtVX_solver.solve(ViX.transpose()); // c*n
        VectorXd ViY = solverV.solve(Y);
  
        est_a = inv_XtVX_ViX * Y;
        Xa_b = Y.array() - var_mu_i.array() * (ViY - ViX * est_a).array();
        Xa_b_exp = Xa_b.array().exp();
        mu = Xa_b_exp.array() / (Xa_b_exp.array() + 1);
        var_mu_i = 1.0 / (mu.array() * (-mu.array() + 1));
        Y = Xa_b.array() + var_mu_i.array() * (phenoVec - mu).array();

        LOGGER << "Init with tao: " << cur_tao << ", Var(Y): " << varVector(Y) << cur_tao << ", fixed effects: " << est_a.transpose() << "." << std::endl;
    }
 

    bool b_reml = false;
    std::ofstream out;
    if(options.find("reml_detail") != options.end()){
        b_reml = true;
        out.open(options["out"] + ".reml");
        if(!out){
            LOGGER.e(0, "Error to open the " + options["out"] + ".reml");
        }

        out << std::setprecision( std::numeric_limits<double>::digits10);
        out << "index\tlogL\tVg\tVe" << std::endl;

    }

    //double cur_tao = -100;
    int curIter = 0;
    bool bConverge = false;
    double tol = 5e-5;
    double abTol = 1e-6;
    int numAbTol = 0;
    vector<int> trailsFull = {801, 16, 16, 16, 16, 16, 16};
    vector<int> trailsHalf = {201, 16, 16, 16, 16, 16, 16};
    vector<int> trailsFine = {51, 16, 16, 16, 16, 16, 16};
    double cutThresh = 0.1;
    int numFullGrid = 3;
    vector<double> hisTao(maxIter);
    double startTao, endTao;
    vector<int> trails;
    while(curIter < maxIter && (!bConverge) && numAbTol < 5){
        LOGGER << "------------------------------------" << std::endl;
        double pre_tao = cur_tao;
        VectorXd pre_est_a = est_a;
        W.diagonal() = var_mu_i;

        double varY = varVector(Y);
        //double MAX_Vg = Vp * options_d["h2_limit"];
        if(curIter == 0){
            trails = trailsFull;
            startTao = 0;
            endTao = varY * 1.01;
            LOGGER << "Fine tuning within " << startTao << " ~ " << endTao << " with " << trails[0] << " steps." << std::endl;
        } else if(curIter > 0 & curIter < numFullGrid){
            trails = trailsHalf;
            startTao = 0;
            endTao = pre_tao * 10;
            if(endTao > varY){
                endTao = varY * 1.01;
            }
            // make sure don't sacrafice the precision
            int temp_step = (int)((endTao / varY) * trailsFull[0]);
            if(temp_step > trails[0]){
                trails[0] = temp_step;
            }
            LOGGER << "Fine tuning within " << startTao << " ~ " << endTao << " with " << trails[0] << " steps." << std::endl;
        }else if(curIter == numFullGrid){
            trails = trailsFine;
            double sum = 0;
            double maxTao = 0, minTao = 10000000;
            for(int i = 0; i < numFullGrid; i++){
                double curTAO = hisTao[i];
                sum += curTAO;
                if(maxTao < curTAO) maxTao = curTAO;
                if(minTao > curTAO) minTao = curTAO;
            }
            double meanTao = sum / numFullGrid; 
            if(meanTao > cutThresh){
                startTao = 0.8 * minTao;
                endTao = 1.2 * maxTao;
            }else{
                startTao = 0;
                endTao = 1.2 * maxTao;
            }
            LOGGER << "Mean tao in first 3 iterations: " << meanTao << ". Fine tuning within " << startTao << " ~ " << endTao << " ." << std::endl;
        }else if(curIter > 2 * numFullGrid){
            double sum = 0;
            double maxTao = 0, minTao = 10000000;
            for(int i = curIter - numFullGrid; i < curIter; i++){
                double curTAO = hisTao[i];
                sum += curTAO;
                if(maxTao < curTAO) maxTao = curTAO;
                if(minTao > curTAO) minTao = curTAO;
            }
            double meanTao = sum / numFullGrid; 
            if(meanTao > cutThresh){
                startTao = 0.8 * minTao;
                endTao = 1.2 * maxTao;
            }else{
                startTao = 0;
                endTao = 1.2 * maxTao;
            }
            LOGGER << "Mean tao in past 3 iterations: " << meanTao << ". Fine tuning within " << startTao << " ~ " << endTao << " .\n" << std::endl;
 
        }
        double end_logL = 0;
        vector<bool> bEndNAN(trails.size());
        vector<double> starts(trails.size());
        vector<double> ends(trails.size());
        vector<double> steps(trails.size());
        
        double start = startTao;
        double end = endTao;

       for(int iter = 0; iter < trails.size(); iter++){
            int n_trails = trails[iter];
            double step = (end - start) / (n_trails - 1);
            vector<double> varcomps(n_trails);
            vector<double> logLs(n_trails);

            #pragma omp parallel for
            for(int i = 0; i < n_trails; i++){
                double temp_tao = start + i * step;
                varcomps[i] = temp_tao;
                logLs[i] = binLogL(temp_tao, fam, W, Y, covar);
            }

            double max_logL = -1e300;
            double max_varcomp = 0;

            int max_index = 0;
            for(int i = 0; i < n_trails; i++){
                double logL = logLs[i];
                if(logL > max_logL){
                    max_logL = logL;
                    max_varcomp = varcomps[i];
                    max_index = i;
                }
            }
            if(max_index == 0){
                start = varcomps[max_index];
                end = varcomps[max_index + 1];
                end_logL = logLs[max_index + 1];
            }else if(max_index == n_trails - 1){
                start = varcomps[max_index - 1];
                end = varcomps[max_index];
                end_logL = logLs[max_index];
            }else{
                start = varcomps[max_index - 1];
                end = varcomps[max_index + 1];
                end_logL = logLs[max_index + 1];
            }

            if(b_reml){
                for(int i = 0; i < logLs.size(); i++){
                    out << i << "\t" << logLs[i] << "\t" << varcomps[i] << std::endl;
                }
                out << "----------------------------" << std::endl;
            }

            if(!std::isfinite(end_logL)){
                bEndNAN[iter] = true;
            }else{
                bEndNAN[iter] = false;
            }
            starts[iter] = start;
            ends[iter] = end;
            steps[iter] = step;


            LOGGER << "Iteration " << iter + 1 << ", step size: " << step << ", logL: " << max_logL
                << ". Tao: " << max_varcomp
                << ", searching range: " << start << " to " << end << std::endl;
            if(b_reml){
                LOGGER << "\n" << "Up boundary detail: " << std::endl;
                LOGGER << "  " << "iter\tstep\tstart\tend\tendNAN" << std::endl;
                LOGGER << std::boolalpha;
                for(int iter = 0; iter < trails.size(); iter++){
                    LOGGER << "  " << iter << "\t" << steps[iter] << "\t" << starts[iter] << "\t" << ends[iter] << "\t" << bEndNAN[iter] << std::endl; 
                }
                LOGGER << std::endl;
            }


        }
        cur_tao = (start + end) / 2.0;
        hisTao[curIter] = cur_tao;
        if(cur_tao < abTol){
            numAbTol++;
        }else{
            numAbTol = 0;
        }
        SpMat V = W + cur_tao * fam;
        solverV.compute(V);

        if(solverV.info() != Eigen::Success){
            LOGGER.e(0, "can't inverse the V matrix!");
        }

        ViX = solverV.solve(covar); // n*c

        auto XtVX_solver = (covar.transpose() * ViX).colPivHouseholderQr();
        if(!XtVX_solver.isInvertible()){
            LOGGER.e(0, "XtViX is not invertable!");
        }

        inv_XtVX_ViX = XtVX_solver.solve(ViX.transpose()); // c*n
        VectorXd ViY = solverV.solve(Y);
  
        est_a = inv_XtVX_ViX * Y;
        Xa_b = Y.array() - var_mu_i.array() * (ViY - ViX * est_a).array();
        Xa_b_exp = Xa_b.array().exp();
        mu = Xa_b_exp.array() / (Xa_b_exp.array() + 1);
        var_mu_i = 1.0 / (mu.array() * (-mu.array() + 1));
        Y = Xa_b.array() + var_mu_i.array() * (phenoVec - mu).array();
        //LOGGER << "Tao value: " << cur_tao << ", Var(Y): " << varVector(Y) << std::endl;
        
        //his_tao[curIter] = cur_tao;
        //his_a[curIter] = est_a;
        LOGGER << "Iter " << curIter << ", tao: " << cur_tao << ", Var(Y): " << varVector(Y) << ", fixed effects: " << est_a.transpose() << "." << std::endl;

        double thresh = ((est_a - pre_est_a).array().abs() / (est_a.array().abs() + pre_est_a.array().abs() + tol)).maxCoeff();
        thresh = std::max(thresh, std::abs((cur_tao - pre_tao) / (cur_tao + pre_tao + tol)) );
        if(thresh  < tol){
            bConverge = true;
            break;
        }


        curIter++;
    }

    if(b_reml){
        out.close();
    }

    if(!bConverge){
        if(numAbTol >=5){
            LOGGER << "fastGWA-GLM-REML stopped at a very small tao value." << std::endl;
        }else{
            LOGGER.w(0, "fastGWA-GLM-REML didn't converge, pay attention to the results, they mostly work OK.");
        }
    }else{
        LOGGER << "fastGWA-GLM-REML converged." << std::endl;
    }
    // for analysis
    W.diagonal() = 1.0 / var_mu_i.array();
    MatrixXd XtW = covar.transpose() * W;
    MatrixXd XtWX = XtW * covar;
    auto t_solver = XtWX.colPivHouseholderQr();
    if(!t_solver.isInvertible()){
        LOGGER.e(0, "XtWX is not invertable!");
    }
    H = t_solver.solve(XtW);

    // for 
    W.diagonal() = var_mu_i;
    SpMat V = W + cur_tao * fam;
    solverV.compute(V);

    taoVal = cur_tao;

    if(solverV.info() != Eigen::Success){
        LOGGER.e(0, "can't inverse the V matrix!");
    }

    ViX = solverV.solve(covar); // n*c

    auto XtVX_solver = (covar.transpose() * ViX).colPivHouseholderQr();
    if(!XtVX_solver.isInvertible()){
        LOGGER.e(0, "XtViX is not invertable!");
    }

    inv_XtVX_ViX = XtVX_solver.solve(ViX.transpose()); // c*n
    //VectorXd PY = ViY - ViX * est_a;
    //template of P
    //VectorXd PY = solverV.solve(Y) - ViX * inv_XtVX_ViX * Y;

    return bConverge;


}

void FastFAM::initVar(){
    numi_indi = num_indi;
    phenoVecMu = phenoVec - mu;
    SPA::setMu(mu);
    dWp = mu.array() * (-mu.array() + 1);

    if(spaCutOff < 0.1){
        spaCutOff = 0.1;
    }
}

void FastFAM::initBinary(const SpMat &fam){

    num_covar = covar.cols();

    VectorXd est_beta = VectorXd::Zero(covar.cols());
    //est_beta.setZero();
    double thresh = 1e-6;
    int maxIter = 200;

    LOGGER.ts("binREML");
    binGridREML(fam, est_beta, maxIter, thresh);
    LOGGER << "Grid REML finished in " << LOGGER.tp("binREML") << " seconds." << std::endl;

    initVar();

    if(!bGLMMExact){
        estBinGamma();
    }
    
    if(options.find("model_only") != options.end()){
        LOGGER << "Saving fastGWA GLM model information..." << std::endl;
        string id_file = options["out"] + ".mdl.id";
        std::ofstream inv_id(id_file.c_str());
        if(!inv_id) LOGGER.e(0, "failed to write " + options["out"]+".grm.id");
        vector<string> ids = pheno->get_id(0, num_indi - 1, "\t");
        for(auto & t_id : ids){
            inv_id << t_id << "\n";
        }
        inv_id.close();
        LOGGER << "Sample information has been saved to [" << id_file << "]." << std::endl;

        char magic[5];
        magic[0] = 'f';
        magic[1] = 'G';
        magic[2] = 'L';
        magic[3] = 'M';
        magic[4] = '\0';
 
        string bin_file = options["out"] + ".mdl.bin2";
        FILE * pFile = fopen(bin_file.c_str(), "wb");
        // TODO check the output
        fwrite(magic, 5, sizeof(char), pFile);
        fwrite(&num_indi, 1, sizeof(uint32_t), pFile);
        uint32_t col_covar = covar.cols();
        fwrite(&col_covar, 1, sizeof(uint32_t), pFile);
        fwrite(&bPreciseCovar, 1, sizeof(bool), pFile);

        fwrite(&taoVal, 1, sizeof(double), pFile);
        fwrite(&c_inf, 1, sizeof(double), pFile);
        fwrite(mu.data(), mu.size(), sizeof(double), pFile);
        fwrite(phenoVec.data(), phenoVec.size(), sizeof(double), pFile);
        fwrite(covar.data(), covar.rows() * covar.cols(), sizeof(double), pFile);
        fwrite(H.data(), H.rows() * H.cols(), sizeof(double), pFile);
        fclose(pFile);
        LOGGER << "Model has been saved to [" << bin_file << "]." << std::endl;
    }
    ViX.resize(0,0);
    inv_XtVX_ViX.resize(0,0);
    //solverV.~SimplicialLDLT<SpMat>();
}

void FastFAM::binGrammar_func(uintptr_t *genobuf, const vector<uint32_t> &markerIndex){
    int nMarker = markerIndex.size();
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nMarker; i++){
        int index_cur_marker = num_grammar_markers + i;
        GenoBufItem item;
        item.extractedMarkerIndex = markerIndex[i];
        geno->getGenoDouble(genobuf, i, &item);
        bValids[index_cur_marker] = item.valid;
        if(item.valid){
            Map< VectorXd > curGeno(item.geno.data(), num_indi);
            if(bPreciseCovar) conditionCovarBinReg(curGeno);
            //VectorXd PY = solverV.solve(Y) - ViX * inv_XtVX_ViX * Y;
            VectorXd PG = solverV.solve(curGeno) - ViX * (inv_XtVX_ViX * curGeno);
            //VectorXd centerGeno = curGeno.array() - curGeno.mean();
            //VectorXd WG = dWp.cwiseProduct(centerGeno);
            //double temp_gamma = curGeno.dot(PG) / centerGeno.dot(WG);
            VectorXd WG = dWp.cwiseProduct(curGeno);
            double temp_gamma = curGeno.dot(PG) / curGeno.dot(WG);

            v_c_infs[index_cur_marker] = temp_gamma;
        }
    }
    num_grammar_markers += nMarker;
}


void FastFAM::estBinGamma(){
    int num_marker_rand = 2000;
    LOGGER.i(0, "\nTuning parameters using null SNPs...");  
    auto total_markers_index = marker->get_extract_index_autosome();
    if(total_markers_index.size() < num_marker_rand){
        LOGGER.e(0, "can't read " + to_string(num_marker_rand) + " SNPs from autosome for tuning.");
    }
    std::random_device rd;
    std::mt19937 gen(seed);
    std::shuffle(total_markers_index.begin(), total_markers_index.end(), gen);
    
   //get previous threshold
    double preAF = geno->getMAF();
    double preInfo = geno->getFilterInfo();
    double preMiss = geno->getFilterMiss();
    if(options.find("c-inf-no-filter") != options.end()){
        double mac20 = 10.0 / num_indi;
        if(preAF < mac20){
            // mac 20
            geno->setMAF(mac20);
        }
        if(preInfo < 0.3 && hasInfo){
            geno->setFilterInfo(0.3);
        }
        if(preMiss < 0.9){
            geno->setFilterMiss(0.9);
        }
    }

    num_grammar_markers = 0;

    int nModelSNP = 200;
    int nMarker = 100;
    int curNumModel = 0;

    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    callBacks.push_back(bind(&FastFAM::binGrammar_func, this, _1, _2));
    double cvThreshold = options_d["cv_threshold"];

    LOGGER.ts("tuning");
    while(true){
        curNumModel += nModelSNP;
        if(curNumModel > num_marker_rand){
            LOGGER.e(0, "fail to estimate gamma, too many signals may exist! Or use --cv-threshold to lower the threshold");
        }

        LOGGER << "  reading " << nModelSNP << " SNP..." << std::endl; 
        v_c_infs.resize(curNumModel);
        bValids.resize(curNumModel);

        vector<uint32_t> marker_index(total_markers_index.begin() + curNumModel - nModelSNP, total_markers_index.begin() + curNumModel);
        std::sort(marker_index.begin(), marker_index.end());
       geno->loopDouble(marker_index, nMarker, true, true, false, false, callBacks);

        if(curNumModel != num_grammar_markers){
            LOGGER.e(0, "some SNPs didn't read successfully!");
        }

        double tmp_cinf = 0;
        double temp_cinf2 = 0;
        int n_valid_null = 0;

        for(int i=0; i < curNumModel; i++){
            if(bValids[i]){
                double curInf = v_c_infs[i];
                tmp_cinf += curInf;
                temp_cinf2 += curInf * curInf;
                n_valid_null++;
            }
        }
        if(n_valid_null < 100){
            LOGGER << "Not enough null SNP, rerun." << std::endl;
            continue;
        }

        c_inf = tmp_cinf / n_valid_null;

        double sumsq = temp_cinf2 - c_inf * c_inf * n_valid_null;
        if(sumsq < 0 && sumsq > -1e-8){
            sumsq = 0;
        }

        double sd_c_inf = std::sqrt(sumsq / (n_valid_null - 1));
        double cv = sd_c_inf / c_inf;
        if(cv < cvThreshold){
            break;
        }else{
            LOGGER << "  gamma = " << c_inf << ", CV = " << cv << ", rerun " << std::endl;
        }

    }

    //reset back
    geno->setMAF(preAF);
    geno->setFilterInfo(preInfo);
    geno->setFilterMiss(preMiss);


    LOGGER.i(0, "Mean gamma = " + to_string(c_inf));
    LOGGER << "Tuning time: " << LOGGER.tp("tuning") << " seconds." << std::endl;
    std::ofstream o_inf;
    bool out_inf = false;
    if(options.find("c-inf") != options.end()){
        o_inf.open(options["out"] + ".cinf");
        o_inf << "CHR\tSNP\tPOS\tA1\tA2\tcinf" << std::endl;
        out_inf = true;
    }
    for(int i=0; i<num_marker_rand; i++){
        if(bValids[i]){
            if(out_inf){
                //o_inf << marker->getMarkerStrExtract(marker_index[i]) << "\t" << v_c_infs[i] << std::endl;
            }
        }
    }
    if(out_inf){
        o_inf.close();
    }
    //LOGGER.i(0, "Got " + to_string(n_valid_null) + " null SNPs");
    v_c_infs.resize(0);
    bValids.resize(0);
}

void FastFAM::calculate_gene(uintptr_t *genobuf, const vector<uint32_t> &markerIndex){

}

void FastFAM::calculate_spa(uintptr_t *genobuf, const vector<uint32_t> &markerIndex){
    int num_marker = markerIndex.size();
    vector<uint8_t> isValids(num_marker);
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;
        geno->getGenoDouble(genobuf, i, &item);
        
        isValids[i] = item.valid;
        if(!item.valid){
            continue;
        }

        Map< VectorXd > xvec(item.geno.data(), num_indi);
        VectorXd xvec2 = xvec;

        if(bPreciseCovar) conditionCovarBinReg(xvec);

        double varSNP = std::sqrt(xvec.dot(dWp.cwiseProduct(xvec)) * c_inf);
 
        SPARes res;
        res.score = xvec.dot(phenoVecMu);
        double chisq = std::abs(res.score) / varSNP;

        res.p = StatLib::pchisqd1(chisq * chisq);

        res.bConverge = true;
        if( chisq < spaCutOff){
            res.p_adj = res.p;
        }else{
            vector<uint32_t> index0;
            index0.reserve(num_indi);
            double thresh = -item.mean + 1e-6;
            //double thresh = 1e-6;
            for(uint32_t i = 0; i < num_indi; i++){
                if(xvec2[i] > thresh){
                    index0.push_back(i);
                }
            }

            if(!bPreciseCovar) conditionCovarBinReg(xvec);

            double q = xvec.dot(phenoVec);
            double qinv = q - res.score - res.score;
            SPA spa(q, qinv, xvec, index0);
            spa.saddleProb(&res);
        }

        Tscore[i] = (float)res.score; //* geno->RDev[cur_raw_marker]; 
        Tse[i] = (float)varSNP;
        p[i] = res.p; 
        padj[i] = res.p_adj;
        rConverge[i] = res.bConverge;
        af[i] = (float)item.af;
        countMarkers[i] = item.nValidN;
        info[i] = item.info;
        double temp_beta = res.score / (varSNP *varSNP);
        beta[i] = (float) temp_beta;
        se[i] = std::abs(temp_beta) / sqrt(StatLib::qchisqd1(res.p_adj));
    }

    output_res_spa(isValids, markerIndex);

}
