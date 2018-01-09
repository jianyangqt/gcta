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

using namespace boost::math;

namespace StatLib{
    chi_squared dist1(1);

    double pchisqd1(double x){
        return cdf(complement(dist1, x));
    }
}
