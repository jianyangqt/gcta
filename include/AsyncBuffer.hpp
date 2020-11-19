/*
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
#include <chrono>
#include "mem.hpp"

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
        uint64_t bufferRawSize = bufferSize * sizeof(T);
        int ret1 = posix_memalign((void **) &(buffer[0]), 64, bufferRawSize);
        int ret2 = posix_memalign((void **) &(buffer[1]), 64, bufferRawSize);
        int ret3 = posix_memalign((void **) &(buffer[2]), 64, bufferRawSize);
        if(ret1 != 0 || ret2 != 0 || ret3 != 0){
            initStatus = false;
        }else{
            initStatus = true;
        }

    }

    ~AsyncBuffer(){
        posix_mem_free(buffer[0]);
        posix_mem_free(buffer[1]);
        posix_mem_free(buffer[2]);
    }

    bool init_status(){
        return initStatus;
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
            std::unique_lock<std::mutex> lock(wmut);
            wcv.wait_for(lock, std::chrono::milliseconds(20000));
            lock.unlock();

            return start_write();
        }
    }

    /*set current buffer to EOF
     * Please don't call this if the stream to read is not end;
    */
    void setEOF(){
        m_stat.eof[m_stat.nextBufferWrite] = true;
    }

    void end_write(){
        uint8_t curBufferWrite = m_stat.nextBufferWrite;
        m_stat.write_count[curBufferWrite] -= 1;
        if(!m_stat.eof[curBufferWrite]){
            m_stat.nextBufferWrite = (curBufferWrite + 1) % 3;
        }
        rcv.notify_one();
    }

    tuple<T*, bool> start_read(){
        uint8_t curBufferRead = m_stat.nextBufferRead;
        if(m_stat.write_count[curBufferRead] == 0 && !m_stat.accessed[curBufferRead]){
            m_stat.read_count[curBufferRead] += 1;
            return tuple<T*, bool>{buffer[curBufferRead], m_stat.eof[curBufferRead]};
        }else{
            std::unique_lock<std::mutex> lock(rmut);
            rcv.wait_for(lock, std::chrono::milliseconds(5000));
            lock.unlock();
            
            return start_read();
        }
    }

    void end_read(){
        uint8_t curBufferRead = m_stat.nextBufferRead;
        m_stat.read_count[curBufferRead] -= 1;
        if((m_stat.read_count[curBufferRead] == 0) && (!m_stat.eof[curBufferRead])){
            m_stat.nextBufferRead = (curBufferRead + 1) % 3;
            m_stat.accessed[curBufferRead] = true;
        }
        wcv.notify_one();
    }

private:
    mutex rmut, wmut;
    BufferStat m_stat;
    condition_variable rcv, wcv;
    T* buffer[3];
    bool initStatus;
};
#endif //GCTA2_ASYNCBUFFER_H
