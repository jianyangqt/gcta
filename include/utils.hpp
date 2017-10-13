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

template <typename T, typename P>
void vector_commonIndex(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<P>& k1, std::vector<P>& k2){
	k1.clear();
	k2.clear();
	if(v1 == v2){
		k1.resize(v1.size());
		std::iota(k1.begin(), k1.end(), 0);
		
		k2.resize(v2.size());
		std::iota(k2.begin(), k2.begin(), 0);
		return;
	}
        std::vector<size_t> v1_index = sort_indexes(v1);
        std::vector<size_t> v2_index = sort_indexes(v2);
	
        std::vector<T> sorted_v1(v1.size(), 0);
        std::vector<T> sorted_v2(v2.size(), 0);
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


#endif
