//
// Created by Zhili Zheng on 28/3/17.
//

#include "gtest/gtest.h"
#include "Logger.h"
#include <chrono>
#include <thread>
#include "test_config.h"

TEST(LoggerTest, Init){
    Logger* log = Logger::GetLogger();
    log->open(CUR_OUT_DIR + "/test_log.log");
    log->i(0, "Test message");
    log->i(1,"Green","test message2");
    log->w(2,"warning message");
    log->d(0,"debug message");
    log->p(3,"Progress 1");
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    log->p(3, "Progress 2");
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    log->i(5, "Success","progress finished");
    LOGGER.i(5, "Success", "logger from new macro");
    try{
        log->e(0,"Error message");
    }catch(char const* e){
        cout << "Error captured: " << e << endl;
        //cout << "Error captured" << endl;
    }
    LOGGER.m(5, "Only here to prompt, but not to file");
}
