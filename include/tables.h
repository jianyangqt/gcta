/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Lookup table for genotype, get the allele counts

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TABLES_H
#define TABLES_H
#define GENO_TABLE_COUNT 65536
#include <stdint.h>
#include <cstddef>

class GBitCountTable{
    public:
        GBitCountTable();
        uint16_t get(const uint16_t& index, int col);
        
    private:
        uint16_t geno_count8[256][4];
        uint16_t geno_count_table[GENO_TABLE_COUNT][4];
        void count_geno_init(uint16_t geno8[][4], size_t length);
        void count_geno(uint16_t geno_array[][4], size_t length, uint16_t geno8[][4]);
};

#endif
