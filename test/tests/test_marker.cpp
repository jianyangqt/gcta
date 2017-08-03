//
// Created by Zhili Zheng on 29/3/17.
//

#include "gtest/gtest.h"
#include "test_config.h"
#include "Marker.h"
#include "Logger.h"

TEST(test_marker, init){
    LOGGER.open(CUR_OUT_DIR + "/test_marker.log");
    Marker marker(CUR_SRC_DIR + "/data/test.bim");
    EXPECT_EQ(1000, marker.count_raw());
}
