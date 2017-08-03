//
// Created by Zhili Zheng on 29/3/17.
//

#include "gtest/gtest.h"
#include "Pheno.h"
#include "Logger.h"
#include "test_config.h"
#include <bitset>

TEST(pheno_test, init_test){
    LOGGER.open(CUR_OUT_DIR + "/test_pheno.log");
    Pheno test(CUR_SRC_DIR + "/data/test.fam");
    EXPECT_EQ(3925, test.count_raw());
    /*
    EXPECT_EQ(0, test.get_mask(10000).mask);
    EXPECT_EQ(test.get_mask(0).mask, 3);
    EXPECT_EQ(test.get_mask(0).shift, 0);
    EXPECT_EQ(test.get_mask(3925).mask, 0);
    EXPECT_EQ(test.get_mask(3925).shift, 0);
    EXPECT_EQ(test.get_mask(3924).mask, 3);
    EXPECT_EQ(test.get_mask(3924).shift, 0);
    EXPECT_EQ(test.get_mask(3923).mask, 192);
    EXPECT_EQ(test.get_mask(3923).shift, 6);
     */
    cout << bitset<8>(192) << endl;

    std::vector<std::string> ids = test.get_id(1,3);
    cout << ids[0] << ", " << ids[1] << ", " << ids.size() << endl;
}


