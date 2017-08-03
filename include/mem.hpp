/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Mocks the same memory allocation api cross platform

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GCTA2_MEM_HPP
#define GCTA2_MEM_HPP

#ifdef _WIN32
#include <malloc.h>
#define posix_memalign(p, a, s) ( ((*(p)) = _aligned_malloc((s), (a))), *(p) ? 0 : errno )
#define posix_mem_free _aligned_free
#else
#include <stdlib.h>
#define posix_mem_free free
#endif

#endif //GCTA2_MEM_HPP
