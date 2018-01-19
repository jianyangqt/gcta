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
#include <condition_variable>
#include <tuple>
using std::mutex;
using std::lock_guard;
using std::condition_variable;
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
        uint8_t curBufferWrite = m_stat.nextBufferWrite;
        if(m_stat.accessed[curBufferWrite] &&
                m_stat.write_count[curBufferWrite] == 0 &&
                m_stat.read_count[curBufferWrite] == 0 &&
                (!m_stat.eof[curBufferWrite])){
            m_stat.write_count[curBufferWrite] += 1;
            m_stat.accessed[curBufferWrite] = false;
            return buffer[curBufferWrite];
        }else{
            std::unique_lock<std::mutex> lock(rmut);
            wcv.wait(lock);
            return start_write();
            //return NULL;
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
        //std::unique_lock<std::mutex> lock(wmut);
        uint8_t curBufferWrite = m_stat.nextBufferWrite;
        m_stat.write_count[curBufferWrite] -= 1;
        if(!m_stat.eof[curBufferWrite]){
            m_stat.nextBufferWrite = (curBufferWrite + 1) % 3;
        }
        //lock.unlock();
        rcv.notify_one();
    }

    tuple<T*, bool> start_read(){
        uint8_t curBufferRead = m_stat.nextBufferRead;
        if(m_stat.write_count[curBufferRead] == 0 && !m_stat.accessed[curBufferRead]){
            m_stat.read_count[curBufferRead] += 1;
            return tuple<T*, bool>{buffer[curBufferRead], m_stat.eof[curBufferRead]};
        }else{
            std::unique_lock<std::mutex> lock(rmut);
            rcv.wait(lock);
            return start_read();
            //return tuple<T*, bool>{NULL, m_stat.eof[curBufferRead]};
        }
    }

    void end_read(){
        //std::unique_lock<std::mutex> lock(rmut);
        uint8_t curBufferRead = m_stat.nextBufferRead;
        m_stat.read_count[curBufferRead] -= 1;
        if((m_stat.read_count[curBufferRead] == 0) && (!m_stat.eof[curBufferRead])){
            m_stat.nextBufferRead = (curBufferRead + 1) % 3;
            m_stat.accessed[curBufferRead] = true;
        }
        //lock.unlock();
        wcv.notify_one();
    }

private:
    mutex rmut;
    BufferStat m_stat;
    condition_variable rcv, wcv;
    T* buffer[3];
};
#endif //GCTA2_ASYNCBUFFER_H
