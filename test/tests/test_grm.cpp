//
// Created by Zhili Zheng on 16/4/17.
//

#include "gtest/gtest.h"
#include "Logger.h"
#include "test_config.h"
#include "GRM.h"
#include <functional>
#include "ThreadPool.h"
using std::bind;
using std::function;
using std::placeholders::_1;
using std::placeholders::_2;

TEST(test_grm, test_grm){

    int num_threads = 4;

    ThreadPool *threadPool = ThreadPool::GetPool(num_threads);

    LOGGER.open(CUR_OUT_DIR + "/test_geno.log");
    Marker marker(CUR_SRC_DIR + "/data/test.bim");
    Pheno pheno(CUR_SRC_DIR + "/data/test.fam");
    Geno geno(CUR_SRC_DIR + "/data/test.bed", &pheno, &marker);
    GRM grm(&geno, 1, 1, num_threads, "./test");
    vector<function<void (uint8_t *, int)>> CALLBacks;
    CALLBacks.push_back(bind(&GRM::calculate_GRM, &grm, _1, _2));
    geno.freq(CALLBacks);
    geno.out_freq("./test");
    grm.deduce_GRM();

}
