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
#include <vector>
#include <numeric>
#include <algorithm>

std::string getHostName();
std::string getLocalTime();

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
    std::vector<size_t> index(v.size());
    std::iota(index.begin(), index.end(), 0);

    std::sort(index.begin(), index.end(),
            [&v](size_t item1, size_t item2) {
                return v[item1] < v[item2];
            });

    return index;
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v1, const std::vector<T> &v2) {
    if(v1.size() != v2.size())return std::vector<size_t>();

    std::vector<size_t> index(v1.size());
    std::iota(index.begin(), index.end(), 0);

    std::sort(index.begin(), index.end(),
            [&v1, &v2](size_t item1, size_t item2) {
                if(v1[item1] < v1[item2]){
                    return true;
                }else if(v1[item1] == v1[item2]){
                    return v2[item1] < v2[item2];
                }else{
                    return false;
                }
            }
    );

    return index;
}

#endif
