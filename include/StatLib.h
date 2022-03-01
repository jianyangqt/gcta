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

#ifndef STAT_LIB_H
#define STAT_LIB_H
#include <Eigen/Eigen>
using Eigen::VectorXd;

namespace StatLib{
    double pchisqd1(double x);
    double pchisqd2(double x); 
    double qchisqd1(double x);

    double pnorm(double x, bool bLowerTail=false);

    bool rankContrast(int n, double *Z);

    VectorXd weightBetaMAF(const VectorXd& MAF, double weight_alpha, double weight_beta);
}
#endif //STAT_LIB_H
