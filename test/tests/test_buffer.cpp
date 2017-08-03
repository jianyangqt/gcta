//
// Created by Zhili Zheng on 31/3/17.
//
#include "gtest/gtest.h"
#include "test_config.h"
#include "Logger.h"
#include <thread>
#include <chrono>
#include "AsyncBuffer.hpp"
#include <tuple>
using std::thread;

TEST(buffer_test, single_test){
    AsyncBuffer<uint8_t> buf(3);
    thread writer([&](){
        LOGGER.m(0, "Write start");
        int i = 0;
        uint8_t *w_buf = NULL;
        while(true){
            w_buf = NULL;
            while(true){
                w_buf = buf.start_write();
                if(!w_buf){
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }else{
                    break;
                }
            }
            *w_buf = i;
            *(w_buf+1) = i + 10;
            *(w_buf+2) = i + 100;
            LOGGER.m(0, "writing " + to_string(*w_buf) + " " + to_string(*(w_buf+1)) + " " + to_string(*(w_buf+2)));
            if(i > 50){
                buf.setEOF();
                buf.end_write();
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            buf.end_write();
            LOGGER.m(2, "written");
            //std::this_thread::sleep_for(std::chrono::milliseconds(500));
            i++;
        }
    });
/*
    thread reader2([&](){
        LOGGER.m(0,"Read 2 start");
        bool c_flag = true;
        uint16_t *r_buf = NULL;
        bool isEOF = false;
        while(c_flag){
            r_buf = NULL;
            isEOF = false;
            while(true){
                std::tie(r_buf, isEOF) = buf.start_read();
                if(!r_buf){
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }else{
                    break;
                }
            }
            LOGGER.m(0, "read " + to_string(*r_buf) + " " + to_string(*(r_buf+1)) + " " + to_string(*(r_buf+2)));
            buf.end_read();
            if(isEOF) {
                c_flag = false;
            }
        }
    });
    */

    thread reader([&](){
        LOGGER.m(0,"Read 1 start");
        bool c_flag = true;
        uint8_t *r_buf = NULL;
        bool isEOF = false;
        while(c_flag){
            r_buf = NULL;
            isEOF = false;
            while(true){
                std::tie(r_buf, isEOF) = buf.start_read();
                if(!r_buf){
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }else{
                    break;
                }
            }
            LOGGER.m(0, "reading " + to_string(*r_buf) + " " + to_string(*(r_buf+1)) + " " + to_string(*(r_buf+2)));
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            buf.end_read();
            LOGGER.m(2, "end read");
            if(isEOF) {
                c_flag = false;
            }
        }
    });

    writer.join();
    reader.join();
    //reader2.join();
}

