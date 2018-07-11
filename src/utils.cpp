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

std::string getFileName(const std::string & path){
    auto index = path.find_last_of("\\/");
    if(index == std::string::npos){
        return path;
    }else{
        return path.substr(index + 1);
    }
}

// get the path without the last slash
std::string getPathName(const std::string & path){
    auto index = path.find_last_of("\\/");
    if(index == std::string::npos){
        return "";
    }else{
        return path.substr(0, index);
    }
}

std::string joinPath(const std::string & dir, const std::string & path){
    #ifdef WIN32
    std::string sep = "\\";
    #else
    std::string sep="/";
    #endif
    if(dir.empty()){
        return path;
    }else{
        return dir + path;
    }
}
