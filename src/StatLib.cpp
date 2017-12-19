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

namespace StatLib{

    //pchisq(chisq value, degree of freedom)
    static long double log_igf(long double S, long double Z);
    static long double KM(long double S, long double Z);

    double pchisq(double x, double df){
        if(x < 0 || df < 0){
            return std::numeric_limits<double>::quiet_NaN();
        }

        double K = ((double)df) * 0.5;
        double X = x * 0.5;
        if(df == 2){
            return exp(-1.0 * X);
        }
        long double PValue, Gam;
        long double ln_PV;
        ln_PV = log_igf(K, X);

        Gam = lgammal(K);

        ln_PV -= Gam;
        PValue = 1.0 - expl(ln_PV);

        return (double)PValue;
    }


    static long double log_igf(long double S, long double Z){
        if(Z < 0.0){
            return 0.0;
        }
        long double Sc, K;
        Sc = (logl(Z) * S) - Z - logl(S);

        K = KM(S, Z);

        return logl(K) + Sc;
    }

    // 1000 interation
    static long double KM(long double S, long double Z){
        long double Sum = 1.0;
        long double Nom = 1.0;
        long double Denom = 1.0;

        for(int I = 0; I < 1000; I++){
            Nom *= Z;
            S++;
            Denom *= S;
            Sum += (Nom / Denom);
        }

        return Sum;
    }
    // End chisq
}
