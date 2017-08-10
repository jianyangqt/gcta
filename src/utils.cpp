#include <utils.hpp>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

std::string getHostName(){
    char *temp = NULL;
    std::string computerName;

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
    temp = getenv("COMPUTERNAME");
    if (temp != NULL) {
        computerName = temp;
    }
#else
    temp = getenv("HOSTNAME");
    if (temp != NULL) {
        computerName = temp;
    } else {
//#include <unistd.h>
        temp = new char[512];
//        if (gethostname(temp, 512) == 0) {
//            computerName = temp;
//        }
        delete []temp;
    }
    return computerName;
#endif
}

std::string getLocalTime(){
    std::time_t t = std::time(nullptr);
    auto tm = std::put_time(std::localtime(&t), "%c %Z");
    std::ostringstream oss;
    oss << tm;
    return oss.str();
}
