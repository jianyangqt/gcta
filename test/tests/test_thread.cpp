//
// Created by Zhili Zheng on 9/7/17.
//

#include "gtest/gtest.h"
#include "Logger.h"
#include <chrono>
#include <thread>
#include "test_config.h"
#include "ThreadPool.h"

class DemoArray{
public:
    DemoArray(int cache_size){
        buffer = new int[cache_size]();
        this->cache_size = cache_size;
    }
    void add_value(int pos, int value){
        buffer[pos] += value;

    }
    void print(){
        for(int index = 0; index != cache_size; index++){
            std::cout << index << ": " << buffer[index] << std::endl;
        }
    }
private:
    int * buffer;
    int cache_size;
};

TEST(ThreadTest, Init){
    ThreadPool *threadPool = ThreadPool::GetPool(10);
    int JOB_COUNT = 30;

    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < JOB_COUNT; ++i)
        THREADS.AddJob([]() {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            //std::lock_guard<std::mutex> lock(m_print);
            std::cout << "Out " << std::endl;
        });

    std::cout << "Out block 1 finished" << std::endl;
    std::cout << "Sleep begin" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(10));
    std::cout << "Sleep end" << std::endl;

    THREADS.WaitAll();
    auto end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;  //15.002

    for (int i = 0; i < JOB_COUNT; ++i)
        THREADS.AddJob([]() {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            //std::lock_guard<std::mutex> lock(m_print);
            std::cout << "Out2 " << std::endl;
        });

    THREADS.WaitAll();
    end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;  //15.002
    std::cout << "block 2 ends." << std::endl;

    DemoArray test(20);
    for(int i = 0; i < 20; i++){
        THREADS.AddJob(std::bind(&DemoArray::add_value, &test, i, i));
    }
    THREADS.WaitAll();
    test.print();
    std::cout << "finished" << std::endl;


}