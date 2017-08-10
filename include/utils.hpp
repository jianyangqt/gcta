/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Some utils in header files

   Singleton pattern to use in the whole program.

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef UTILS_H
#define UTILS_H
#if defined(_MSC_VER)
#include <intrin.h>
#define popcount __popcnt64
#elif defined(__GNUC__) || defined(__GNUG__)
#define popcount __builtin_popcountll 
#else
#error "Compiler is not supported"
#endif

#include <string>

std::string getHostName();
std::string getLocalTime();

#endif
