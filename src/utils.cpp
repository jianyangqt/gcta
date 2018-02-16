#include <utils.hpp>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

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
        temp = new char[512];
        if (gethostname(temp, 512) == 0) {
            computerName = temp;
        }
        delete []temp;
    }
#endif
    return computerName;
}

std::string getLocalTime(){
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    auto tm = std::put_time(std::localtime(&now_c), "%T %Z on %a %b %d %Y");
    std::ostringstream oss;
    oss << tm;
    return oss.str();
}
