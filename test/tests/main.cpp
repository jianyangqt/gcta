#include <gtest/gtest.h>
int test_argc;
char** test_argv;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv); 
    test_argc = argc;
    test_argv = argv;
    return RUN_ALL_TESTS();
}
