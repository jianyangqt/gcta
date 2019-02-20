/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Some statistic functions

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/
#include "StatLib.h"
#include <cmath>
#include <cstdio>
#include <limits>
#include <boost/math/distributions/chi_squared.hpp>
#include <mkl.h>

using namespace boost::math;

namespace StatLib{
    chi_squared dist1(1);

    double pchisqd1(double x){
        if(x > 0){
            return cdf(complement(dist1, x));
        }else{
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    //n rank size,  Z shall be n * n double memory.
    // return true, success; false, unsuccess
    bool rankContrast(int n, double * Z){
        double mean = 1.0 - (1.0 + n) / 2.0;

        //outer op
        double* X = new double[n * n];
        for(int i = 0; i < n; i++){
            int base_index = i * n;
            for(int j = 0; j < n; j++){
                X[base_index + j] = pow(j + mean, i);
            }
        }

        double* tau = new double[n];
        double* work = new double[n];
        int info = 0;
        int lda = n;
        int lwork = n;
        dgeqrf(&n, &n, X, &lda, tau, work, &lwork, &info);
        if(info != 0){
            return false;
        }

        double *c = Z;
        for(int i = 0; i < n; i++){
            int base_index = i * n;
            for(int j = 0; j < n; j++){
                int index = base_index + j;
                if(i != j){
                    c[index] = 0.0;
                }else{
                    c[index] = X[index];
                }
            }
        }

        char side = 'L';
        char t = 'N';

        dormqr(&side, &t, &n, &n, &n, X, &lda, tau, c, 
                &lda, work, &lwork, &info);

        if(info != 0){
            return false;
        }

        for(int i = 0; i < n; i++){
            double rdiv = 0;
            int base_index = i * n;
            for(int j = 0; j < n; j++){
                rdiv += c[base_index +j] * c[base_index + j];
            }
            rdiv = 1.0 / sqrt(rdiv);
            for(int j = 0; j < n; j++){
                c[base_index + j] = c[base_index + j] * rdiv;
            }
        }

        return true;
    }


}
