#include <gtest/gtest.h> // googletest header file
#include <string>
#include <iostream>
using std::string;
using std::cout;

extern int test_argc;
extern char** test_argv;

const char *actualValTrue  = "hello gtest";
const char *actualValFalse = "hello world";
const char *expectVal      = "hello gtest";

TEST(StrCompare, CStrEqual) {
    EXPECT_STREQ(expectVal, actualValTrue);
    cout << "ARGC " << test_argc << test_argv[0] << std::endl;
}

TEST(StrCompare, CStrNotEqual) {
    EXPECT_STRNE(expectVal, actualValFalse);
}



