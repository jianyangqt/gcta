/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Global Singleton logger system with robust functions and thread safe.

      * open a logger file first, otherwise no output into log
      * e(int level, string message):  prompt error, and exit the program, level is the number of indent spaces
      * i:  prompt information
      * w:  prompt warning message
      * d:  debug message, only seen in the debug mode
      * m:  message that only show on the terminal that not log into logger file
      * l:  log into logger file only;
      * p:  progress, that always show in one line in the terminal but no output into logger file.
      * << :  use like std::cout

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#include "Logger.h"
#include <iostream>

#ifdef _WIN32
#include <io.h>
#include <stdio.h>
#define ISATTY(x) _isatty(_fileno(x))
#else 
#include <unistd.h>
#define ISATTY(x) isatty(fileno(x))
#endif

using std::cout;

std::mutex Logger::log_mutex;
std::string Logger::empty = {};
std::map<string, std::chrono::time_point<std::chrono::steady_clock>> Logger::time_map;
#ifdef _WIN32 
std::map<Logger::Type, string> Logger::style_map = {{Logger::INFO, ""}, {Logger::PROMPT, ""},
                                                     {Logger::PROGRESS, "\r"}, {Logger::WARN, ""},
                                                     {Logger::ERROR, ""}, {Logger::DEBUG, ""}};
#else
std::map<Logger::Type, string> Logger::style_map = {{Logger::INFO, "\033[0m"}, {Logger::PROMPT, "\033[0;32m"},
                                                     {Logger::PROGRESS, "\r"}, {Logger::WARN, "\033[0;33m"},
                                                     {Logger::ERROR, "\033[0;31m"}, {Logger::DEBUG, "\033[0;34m"}};
#endif
Logger* Logger::m_pThis = NULL;
Logger::Type Logger::m_stat = Logger::INFO;
std::ofstream Logger::m_logFile;
string Logger::m_FileName = "";

Logger::Logger(){
    if(!ISATTY(stdout)){
        Logger::style_map =  {{Logger::INFO, ""}, {Logger::PROMPT, ""},
                              {Logger::PROGRESS, "\n"}, {Logger::WARN, ""},
                              {Logger::ERROR, ""}, {Logger::DEBUG, ""}};
    }

}

Logger* Logger::GetLogger(){
    if(m_pThis == NULL) {
        m_pThis = new Logger();
    }
    return m_pThis;
}

bool Logger::check(){
    if(m_pThis == NULL){
        return false;
    }
    if(!m_logFile.is_open()){
        return false;
    }
    if(m_FileName.empty()){
        return false;
    }
    return true;
}

void Logger::open(string ofile){
    if(check()){
        cout << "Logger has been set, not support to set another time" << endl;
    }else{
        m_logFile.open(ofile, std::ios::out);
        m_FileName = ofile;
        if(!check()){
            m_pThis->e(0, "can't write to log file [" + ofile + "].\nPlease check file/folder permission or disk quota.");
        }
    }
}

void Logger::close(){
    m_logFile.close();
}

void Logger::flush(){
    m_logFile.flush();
}

void Logger::Log(int level, Type type, const string& prompt, const string& message){
    string spaces(level * 2, ' ');
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << spaces << type << prompt << INFO << message << endl;
}

void Logger::e(int level, const string& message, const string& title){
    string head = title.empty() ? "Error: " : (title+" ");
    m_pThis->Log(level, ERROR, head, message);
    m_pThis->Log(level, INFO, "", "An error occurs, please check the options or data");
    exit(1);
}

int Logger::precision(int p){
    cout.precision(p);
    m_logFile.precision(p);
    return cout.precision();
}

string Logger::setprecision(int p){
    cout.precision(p);
    m_logFile.precision(p);
    return("");
}

int Logger::precision(){
    return cout.precision();
}

void Logger::i(int level, const string& message, const string& title){
    string head = title.empty() ? "" : (title + " ");
    m_pThis->Log(level, PROMPT, head, message);
}

void Logger::w(int level, const string &message, const string& title) {
    string head = title.empty() ? "Warning: " : (title + " ");
    m_pThis->Log(level, WARN, head, message);
}

void Logger::d(int level, const string &message, const string& title) {
    #ifndef NDEBUG
    string head = title.empty() ? "Debug: " : (title + " ");
    m_pThis->Log(level, DEBUG, head, message);
    #endif
}

void Logger::p(int level, const string &message, const string& title) {
    string head = title.empty() ? "" : (title + " ");
    string spaces(level * 2, ' ');
    m_stat = PROGRESS;
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << spaces << head << message << PROGRESS << std::flush;
    m_stat = INFO;
}

void Logger::m(int level, const string &message, const string& title){
    string head = title.empty() ? "" : (title + " ");
    string spaces(level * 2, ' ');
    m_stat = PROGRESS;
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << spaces << PROMPT << head << INFO << message << endl;
    m_stat = INFO;
}

void Logger::l(int level, const string &message, const string &title){
    string head = title.empty() ? "" : (title + " ");
    string spaces(level * 2, ' ');
    m_logFile << spaces << head << message << endl;
}

void Logger::ts(string key){
    time_map[key] = std::chrono::steady_clock::now();
}

float Logger::tp(string key){
    auto end = std::chrono::steady_clock::now();
    float secs = 0.0;
    if(time_map.find(key) != time_map.end()){
        auto duration = end - time_map[key];
        secs = std::chrono::duration_cast<std::chrono::duration<float>>(duration).count();
    }else{
        m_pThis->w(2, "can't find key of time start");
    }
    return secs;
}

Logger& Logger::operator<<(Type type){
    m_stat = type;
    cout << style_map[type];
    return *m_pThis;
}

Logger& Logger::operator<<(std::ostream& (*op)(std::ostream&)){
    (*op)(cout);
    if(m_stat !=PROGRESS){
        (*op)(m_logFile);
        //m_logFile.flush();
    }
    return *m_pThis;
}

Logger& Logger::operator<<(std::ios& (*pf)(std::ios&)){
    cout << pf;
    m_logFile << pf;
    return *m_pThis;
}

Logger& Logger::operator<<(std::ios_base& (*pf)(std::ios_base&)){
    cout << pf;
    m_logFile << pf;
    return *m_pThis;
}


