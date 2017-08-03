/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Asynchronous tri-buffer for parallel loading and processing.

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef GCTA2_ASYNCBUFFER_H
#define GCTA2_ASYNCBUFFER_H
#include <mutex>
#include <tuple>
using std::mutex;
using std::lock_guard;
using std::tuple;
using std::tie;

struct BufferStat{
    uint8_t write_count[3] = {0, 0, 0};
    uint8_t read_count[3] = {0, 0, 0};
    uint8_t nextBufferRead = 0;
    uint8_t nextBufferWrite = 0;
    bool eof[3] = {false, false, false};
    bool accessed[3] = {true, true, true};
};

template <typename T>
class AsyncBuffer {
public:
    AsyncBuffer(uint64_t bufferSize){
        buffer[0] = new T[bufferSize];
        buffer[1] = new T[bufferSize];
        buffer[2] = new T[bufferSize];
    }

    ~AsyncBuffer(){
        delete[] buffer[0];
        delete[] buffer[1];
        delete[] buffer[2];
    }


    T* start_write(){
        //std::lock_guard<std::mutex> lock(buf_mutex);
        uint8_t curBufferWrite = m_stat.nextBufferWrite;
        if(m_stat.accessed[curBufferWrite] &&
                m_stat.write_count[curBufferWrite] == 0 &&
                m_stat.read_count[curBufferWrite] == 0 &&
                (!m_stat.eof[curBufferWrite])){
            m_stat.write_count[curBufferWrite] += 1;
            m_stat.accessed[curBufferWrite] = false;
            return buffer[curBufferWrite];
        }else{
            return NULL;
        }
    }

    /*set current buffer to EOF
     * Please don't call this if the stream to read is not end;
    */
    void setEOF(){
        //std::lock_guard<std::mutex> lock(buf_mutex);
        m_stat.eof[m_stat.nextBufferWrite] = true;
    }

    void end_write(){
        //std::lock_guard<std::mutex> lock(buf_mutex);
        uint8_t curBufferWrite = m_stat.nextBufferWrite;
        m_stat.write_count[curBufferWrite] -= 1;
        if(!m_stat.eof[curBufferWrite]){
            m_stat.nextBufferWrite = (curBufferWrite + 1) % 3;
        }
    }

    tuple<T*, bool> start_read(){
        //std::lock_guard<std::mutex> lock(buf_mutex);
        uint8_t curBufferRead = m_stat.nextBufferRead;
        if(m_stat.write_count[curBufferRead] == 0 && !m_stat.accessed[curBufferRead]){
            m_stat.read_count[curBufferRead] += 1;
            return tuple<T*, bool>{buffer[curBufferRead], m_stat.eof[curBufferRead]};
        }else{
            return tuple<T*, bool>{NULL, m_stat.eof[curBufferRead]};
        }
    }

    void end_read(){
        //std::lock_guard<std::mutex> lock(buf_mutex);
        uint8_t curBufferRead = m_stat.nextBufferRead;
        m_stat.read_count[curBufferRead] -= 1;
        if((m_stat.read_count[curBufferRead] == 0) && (!m_stat.eof[curBufferRead])){
            m_stat.nextBufferRead = (curBufferRead + 1) % 3;
            m_stat.accessed[curBufferRead] = true;
        }
    }

private:
    //mutex buf_mutex;
    BufferStat m_stat;
    T* buffer[3];
};
#endif //GCTA2_ASYNCBUFFER_H
