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

#ifndef GCTA2_FASTFAM_H
#define GCTA2_FASTFAM_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include "Logger.h"
#include "Geno.h"
#include "Pheno.h"
#include "Marker.h" 
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <mutex>
#include <omp.h>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Dynamic;
using Eigen::Ref;
using std::vector;

typedef SparseMatrix<double, Eigen::ColMajor, long long> SpMat;

class FastFAM {
public:
    FastFAM();
    ~FastFAM();

    void calculate_fam(uintptr_t *buf, const vector<uint32_t> &markerIndex);
    void calculate_grammar(uintptr_t *buf, const vector<uint32_t> &markerIndex);
    void calculate_gwa(uintptr_t * geno, const vector<uint32_t> &markerIndex);
    void output_res(const vector<uint8_t> &isValids, const vector<uint32_t> markerIndex);

    //void readGenoSample(uint64_t *buf, int num_marker);
    void genRandY(uint64_t *buf, int num_marker);
    void estBeta(uint64_t *buf, int num_marker);

    //void output(string filename);
    //void initMarkerVars();

    static void readFAM(string filename, SpMat& fam, const vector<string> &ids, vector<uint32_t> &remain_index);
    static double HEreg(vector<double> &Zij, vector<double> &Aij, bool &isSig);
    static double HEreg(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig);
    double MCREML(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig);
    double spREML(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig);

    void conditionCovarReg(Eigen::Ref<VectorXd> pheno);
    void conditionCovarReg(VectorXd &pheno, VectorXd &condPheno);
    
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();
    void processFAM(vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks);
    //void processFAM();


private:
    Pheno *pheno;
    Marker *marker;
    Geno *geno;
    uint32_t num_indi;
    uint32_t num_marker;
    uint32_t num_finished_marker = 0;
    float *beta = NULL;
    float *se = NULL;
    double *p = NULL;
    uint32_t * countMarkers = NULL;
    float *af = NULL;
    float *info = NULL;
    bool bOutResAll = false;
    
    bool fam_flag;
    MatrixXd H, covar;
    bool covarFlag = false;

    bool bGrammar = false;
    int num_grammar_markers;
    VectorXd Vi_y;
    VectorXd Vi_y_cinf;
    double c_inf;
    uint64_t finished_rand_marker = 0;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> solver;
    void grammar_func(uintptr_t *genobuf, const vector<uint32_t> &markerIndex);
    vector<double> v_chisq;
    vector<double> v_c_infs;
    vector<uint8_t> bValids;
 

    //std::mutex chisq_lock;
    //
    //REML
    vector<SpMat> A;
    void logLREML(const Ref<const VectorXd> pheno, vector<double> &varcomp, double &logL, double *Hinv=NULL);

    void loadModel();

    SpMat V_inverse;
    vector<double> phenos;
    VectorXd phenoVec;
    VectorXd rawPhenoVec;

    void inverseFAM(SpMat& fam, double VG, double VR);
    void makeIH(MatrixXd &X);
    void grammar(SpMat& fam, double VG, double VR);

    static map<string, string> options;
    static map<string, double> options_d;
    static vector<string> processFunctions;
    double fitREML(double logdet, const SpMat &fam, const VectorXd &pheno);

    //for BOLT reml
    MatrixXd rand_y;
    MatrixXd rand_beta;
    MatrixXd rand_e;

    MatrixXd randVinvY;
    VectorXd VinvY;

    MatrixXd randEstBeta;
    VectorXd phenoEstBeta;
    uint64_t finished_boltreml_marker = 0;
    
    int mcTrails;


    //void reg_thread(uint8_t *buf, int from_marker, int to_marker);
    bool bSaveBin = false;
    bool bNoSaveMarker = false;
    bool hasInfo = false;
    string sFileName;

    std::ofstream osOut;
    FILE * bOut = NULL;
    vector<char> osBuf;
    uint32_t numMarkerOutput = 0;

};


#endif


