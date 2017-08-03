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

#include "ThreadPool.h"

ThreadPool* ThreadPool::m_pThis = NULL;

ThreadPool* ThreadPool::GetPool(int threadCount){
    if(m_pThis == NULL) {
        m_pThis = new ThreadPool(threadCount);
    }
    return m_pThis;
}

ThreadPool::ThreadPool(int threadCount) {
    num_job_left = 0;
    is_exit_mainloop = false;
    threads = std::vector<std::thread>(threadCount);

    for(int i = 0; i != threadCount; i++){
        threads[i] = std::thread([this]{
            this->MainLoop();
        });
    }
}

void ThreadPool::AddJob(std::function<void(void)> job_function){
    {
        std::lock_guard<std::mutex> lock(mutex_queue);
        thread_queue.emplace(job_function);
    }
    {
        std::lock_guard<std::mutex> lock(mutex_job_left);
        num_job_left++;
    }
    cond_jobAvail.notify_one();
}

void ThreadPool::JoinAll(){
    {
        std::lock_guard<std::mutex> lock(mutex_queue);
        if(is_exit_mainloop){ return; }
    }

    cond_jobAvail.notify_all();

    for(auto &thread : threads){
        if(thread.joinable()){
            thread.join();
        }
    }
}

void ThreadPool::WaitAll() {
    std::unique_lock<std::mutex> lock(mutex_job_left);
    if(num_job_left > 0){
        cond_wait.wait(lock, [this]{
            return num_job_left == 0;
        });
    }
}

int ThreadPool::GetRemainCount() {
    std::unique_lock<std::mutex> lock(mutex_job_left);
    return num_job_left;
}

int ThreadPool::getThreadCount(){
    return threads.size();
}

void ThreadPool::MainLoop() {
    std::function<void(void)> job_function;
    while(true){
        {
            std::unique_lock<std::mutex> lock(mutex_queue);
            if(is_exit_mainloop){ return; }

            cond_jobAvail.wait(lock, [this]{
                return thread_queue.size() > 0 || is_exit_mainloop;
            });

            if(is_exit_mainloop){ return; }

            job_function = thread_queue.front();
            thread_queue.pop();
        }

        job_function();

        {
            std::lock_guard<std::mutex> lock(mutex_job_left);
            num_job_left--;
        }

        cond_wait.notify_one();

    }
}
