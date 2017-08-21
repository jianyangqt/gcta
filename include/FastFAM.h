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

#include "Logger.h"
#include "Geno.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Dynamic;

typedef SparseMatrix<double, Eigen::ColMajor, long long> SpMat;

class FastFAM {
public:
    FastFAM(Geno *geno);
    ~FastFAM(){
        delete[] beta;
        delete[] se;
        delete[] p;
    }

    void calculate_fam(uint8_t *buf, int num_marker);
    void output(string filename);

    static void readFAM(string filename, SpMat& fam);
    
    static int registerOption(map<string, vector<string>>& options_in);
    static void processMain();

private:
    Geno *geno;
    uint32_t num_indi;
    uint32_t num_marker;
    uint32_t num_finished_marker = 0;
    double *beta;
    double *se;
    double *p;

    SpMat V_inverse;
    VectorXd phenoVec;

    void inverseFAM(SpMat& fam, double VG, double VR);

};


#endif


