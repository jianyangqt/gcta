//
// Created by Zhili Zheng on 31/3/17.
//
#include "gtest/gtest.h"
#include "test_config.h"
#include "Geno.h"
#include "Logger.h"

TEST(test_geno, valid_geno){
    EXPECT_EQ(108, 0x6c);
    EXPECT_EQ(27, 0x1b);
    EXPECT_EQ(1, 0x01);
    LOGGER.open(CUR_OUT_DIR + "/test_geno.log");
    Marker marker(CUR_SRC_DIR + "/data/test.bim");
    Pheno pheno(CUR_SRC_DIR + "/data/test.fam");
    Geno geno(CUR_SRC_DIR + "/data/test.bed", &pheno, &marker);
}
