#include <gtest/gtest.h>
#include "tables.h"
#include <string>
#include <iostream>
using std::string;
using std::cout;

extern int test_argc;
extern char** test_argv;

TEST(Table, Correct0){
    GBitCountTable g_table;
    EXPECT_EQ(8, g_table.get(0,0));
    EXPECT_EQ(0,g_table.get(0,1));
    EXPECT_EQ(0,g_table.get(0,2));
    EXPECT_EQ(0, g_table.get(0,3));
    //std::cout << test_argv[0] << std::endl;
}

