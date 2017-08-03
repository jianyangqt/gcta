/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Robust thread pool with C++11 only, cross platform.

   Singleton pattern to use in the whole program.

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GCTA2_THREADPOOL_H
#define GCTA2_THREADPOOL_H

#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#define THREADS (*ThreadPool::GetPool())
#define THREADS_P ThreadPool::GetPool()

class ThreadPool {
public:
    static ThreadPool *GetPool(int threadCount = 5);
    void AddJob(std::function<void(void)> job);
    void JoinAll();
    void WaitAll();
    int GetRemainCount();
    int getThreadCount();
    ~ThreadPool(){JoinAll();};

private:
    explicit ThreadPool(int threadCount);
    ThreadPool(const ThreadPool&){};
    ThreadPool& operator=(const ThreadPool&){return *this;};

    static ThreadPool *m_pThis;

    void MainLoop();

    std::vector<std::thread> threads;

    std::queue<std::function<void(void)>> thread_queue;
    std::mutex mutex_queue;

    int num_job_left = 0;
    std::mutex mutex_job_left;
    bool is_exit_mainloop = false;
    std::condition_variable cond_jobAvail;
    std::condition_variable cond_wait;

};


#endif //GCTA2_THREADPOOL_H
