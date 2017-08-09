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
using std::cout;

std::mutex Logger::log_mutex;
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

Logger::Logger(){}

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
            throw("Error: can't write to log file [" + ofile + "]");
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
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << level << type << prompt << INFO << message << endl;
}

void Logger::e(int level, const string& message){
    m_pThis->Log(level, ERROR, "Error: ", message);
    m_pThis->Log(level, INFO, "", "An error occurs, please check the data");
    exit(1);
}

void Logger::i(int level, const string& message){
    m_pThis->Log(level, INFO, "", message);
}

void Logger::i(int level, const string& prompt, const string& message){
    m_pThis->Log(level, PROMPT, prompt+" ", message);
}

void Logger::w(int level, const string &message) {
    m_pThis->Log(level, WARN, "Warn: ", message);
}

void Logger::d(int level, const string &message) {
    #ifndef NDEBUG
    m_pThis->Log(level, DEBUG, "Debug: ", message);
    #endif
}

void Logger::p(int level, const string &message) {
    m_stat = PROGRESS;
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << level << message << PROGRESS << std::flush;
}

void Logger::m(int level, const string &message){
    m_stat = PROGRESS;
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << level << message << endl;
}

void Logger::l(int level, const string &message){
    string spaces(level * 2, ' ');
    m_logFile << spaces << message << endl;
}

Logger& Logger::operator<<(const string& message){
    cout << message;
    if(m_stat != PROGRESS){
        m_logFile << message;
        //m_logFile.flush();
    }
    return *m_pThis;
}

Logger& Logger::operator<<(int level){
    string spaces(level * 2, ' ');
    cout << spaces;
    if(m_stat != PROGRESS) {
        m_logFile << spaces;
        //m_logFile.flush();
    }
    return *m_pThis;
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


