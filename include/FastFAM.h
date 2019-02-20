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
    FastFAM(Geno *geno);
    ~FastFAM();

    void calculate_fam(uint64_t *buf, int num_marker);
    void calculate_gwa(uint64_t *buf, int num_marker);
    void output(string filename);
    void initMarkerVars();

    static void readFAM(string filename, SpMat& fam, const vector<string> &ids, vector<uint32_t> &remain_index);
    static double HEreg(vector<double> &Zij, vector<double> &Aij, bool &isSig);
    static double HEreg(const Ref<const SpMat> fam, const Ref<const VectorXd> pheno, bool &isSig);
    static void conditionCovarReg(VectorXd &pheno, MatrixXd &covar);
    
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();

private:
    Geno *geno;
    uint32_t num_indi;
    uint32_t num_marker;
    uint32_t num_finished_marker = 0;
    double *beta = NULL;
    double *se = NULL;
    double *p = NULL;
    bool fam_flag;

    std::mutex chisq_lock;
    

    SpMat V_inverse;
    vector<double> phenos;
    VectorXd phenoVec;

    void inverseFAM(SpMat& fam, double VG, double VR);

    static map<string, string> options;
    static map<string, double> options_d;
    static vector<string> processFunctions;

    //void reg_thread(uint8_t *buf, int from_marker, int to_marker);

};


#endif


