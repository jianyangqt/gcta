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
#include <sstream>

std::string getHostName();
std::string getLocalTime();
std::string getFileName(const std::string & path);
std::string getPathName(const std::string & path);
std::string joinPath(const std::string & dir, const std::string & path);
std::string getSSEvar();
std::string getOSName();

template <typename T>
bool hasVectorDuplicate(const std::vector<T> &v){
    std::vector<T> t = v;
    std::sort(t.begin(), t.end());
    auto last = std::unique(t.begin(), t.end());
    return last != t.end();
}

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

template <typename T, typename P>
void vector_commonIndex(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<P>& k1, std::vector<P>& k2){
	k1.clear();
	k2.clear();
	if(v1 == v2){
		k1.resize(v1.size());
		std::iota(k1.begin(), k1.end(), 0);
		
		k2.resize(v2.size());
		std::iota(k2.begin(), k2.end(), 0);
		return;
	}
        std::vector<size_t> v1_index = sort_indexes(v1);
        std::vector<size_t> v2_index = sort_indexes(v2);
	
        std::vector<T> sorted_v1(v1.size());
        std::vector<T> sorted_v2(v2.size());
	std::transform(v1_index.begin(), v1_index.end(), sorted_v1.begin(), [&v1](size_t pos){return v1[pos];});
	std::transform(v2_index.begin(), v2_index.end(), sorted_v2.begin(), [&v2](size_t pos){return v2[pos];});
	
	auto v1_begin = sorted_v1.begin();
	
	for(uint32_t v2_ind = 0; v2_ind != v2.size(); v2_ind++){
		T v2_temp = sorted_v2[v2_ind];
		auto v1_iter = std::lower_bound(v1_begin, sorted_v1.end(), v2_temp);
		
		if( v1_iter != sorted_v1.end() && !(v2_temp < *v1_iter)){
			k1.push_back(v1_index[v1_iter - sorted_v1.begin()]);
			k2.push_back(v2_index[v2_ind]);
			v1_begin = v1_iter;
		}
	}
        // order sometimes matters	
	//std::sort(k1.begin(), k1.end());
	//std::sort(k2.begin(), k2.end());
}

template <typename T, typename P>
void vector_commonIndex_sorted1(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<P>& k1, std::vector<P>& k2){
    vector_commonIndex(v1, v2, k1, k2);
    std::vector<size_t> k1_index = sort_indexes(k1);
    std::vector<P> k1_sorted(k1.size()), k2_sorted(k2.size());
    std::transform(k1_index.begin(), k1_index.end(), k1_sorted.begin(), [&k1](size_t pos){return k1[pos];});
    std::transform(k1_index.begin(), k1_index.end(), k2_sorted.begin(), [&k2](size_t pos){return k2[pos];});
    k1 = k1_sorted;
    k2 = k2_sorted;
}

/* Permute vector elements to all the combinations
 * @param elements: vector of each possible values
 * @param combines: out: combinations of these vectors;
 * @param temp and col:  only for internal use, don't change these values;
 * @example
 *   #include <vector>
     vector<vector<int>> elements = {
           {0, 1, 2},
           {10, 11, 12},
           {20, 21, 22}
     };
     vector<vector<int>> combines;
     permute_vector(elements, combines);
 * @depends vector c++11
 * @Bug report: zhili<zhilizheng@outlook.com>
 */
template <typename T>
void permute_vector(const std::vector<std::vector<T>> &elements, std::vector<std::vector<T>> &combines, std::vector<T> temp = {}, size_t col = 0){
    size_t e_size = elements.size();
    temp.resize(e_size);
    if(col < e_size){
        for(size_t i = 0; i < elements[col].size(); i++){
            temp[col] = elements[col][i];
            permute_vector(elements, combines, temp, col + 1);
        }
    }else{
        combines.emplace_back(temp);
   }
}

template <typename T>
std::string to_string_precision(const T value, const int n = 0){
    std::ostringstream out;
    if(n > 0){
        out.precision(n);
    }
    out << std::fixed << value;
    return out.str();
} 
#endif
