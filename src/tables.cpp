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
#include "tables.h"

GBitCountTable::GBitCountTable(){
    count_geno_init(geno_count8, 256);
    count_geno(geno_count_table, GENO_TABLE_COUNT, geno_count8);
}

uint16_t GBitCountTable::get(const uint16_t& index, int col){
    return geno_count_table[index][col];
}

void GBitCountTable::count_geno_init(uint16_t geno8[][4], size_t length){
    uint16_t mask = 0b11;
    for(uint16_t index = 0; index < length; index++){
        for(int bit_shift = 0; bit_shift < 8; bit_shift += 2){
            uint8_t cur_bit = (index >> bit_shift) & mask;
            switch(cur_bit){
                case 0b00:
                    geno8[index][0]++;
                    break;
                case 0b10:
                    geno8[index][1]++;
                    break;
                case 0b11:
                    geno8[index][2]++;
                    break;
                default:
                    geno8[index][3]++;
            }
        }
    }
}

void GBitCountTable::count_geno(uint16_t geno_array[][4], size_t length, uint16_t geno8[][4]){
    for (size_t index = 0; index < length; index++){
        uint16_t lower = index & 0xFF;
        uint16_t upper = index >> 8;
        geno_array[index][0] = geno8[lower][0] + geno8[upper][0];
        geno_array[index][1] = geno8[lower][1] + geno8[upper][1];
        geno_array[index][2] = geno8[lower][2] + geno8[upper][2];
        geno_array[index][3] = geno8[lower][3] + geno8[upper][3];
    }
}

