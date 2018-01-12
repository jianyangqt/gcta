#include <gtest/gtest.h>
#include "StatLib.h"
#include <string>
#include <iostream>
using std::string;
using std::cout;

extern int test_argc;
extern char** test_argv;

TEST(StatLib, validStat){
    EXPECT_TRUE(false) << "test: " << std::scientific << StatLib::pchisqd1(300) << std::endl;
}

