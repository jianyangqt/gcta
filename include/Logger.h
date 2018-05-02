/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Global Singleton logger system with robust functions and thread safe.

      * open a logger file first, otherwise no output into log
      * e(int level, string message):  prompt error, and exit the program, level is the number of indent spaces
      * i:  prompt information
      * w:  prompt warning message
      * d:  debug message, only seen in the debug mode
      * m:  message that only show on the terminal that not log into logger file
      * l:  log into logger file only
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

#ifndef GCTA2_LOGGER_H
#define GCTA2_LOGGER_H
#if !defined(_DEBUG) && defined(_WIN32)
#define NDEBUG
#endif
#include <mutex>
#include <string>
#include <map>
#include <fstream>
#include <ios>
#include <chrono>
#include <iostream>
#define LOGGER (*Logger::GetLogger())
#define LOGGER_P Logger::GetLogger()
using std::string;
using std::endl;

class Logger {
public:
    enum Type{INFO, PROMPT, PROGRESS, WARN, ERROR, DEBUG};

    void Log(int level, Type type, const string& prompt, const string& message);
    void e(int level, const string& message, const string& title = empty);
    void i(int level, const string& message, const string& title = empty);
    void d(int level, const string& message, const string& title = empty);
    void w(int level, const string& message, const string& title = empty);
    void p(int level, const string& message, const string& title = empty);
    void m(int level, const string& message, const string& title = empty);
    void l(int level, const string& message, const string& title = empty);
    void ts(string marker);
    float tp(string marker);
    Logger& operator<<(Type type);
    Logger& operator<<(std::ostream& (*manip)(std::ostream&));
    Logger& operator<<(std::ios& (*pf)(std::ios&));
    Logger& operator<<(std::ios_base& (*pf)(std::ios_base&));

    template<typename T>
    Logger& operator<<(const T& t){
        std::cout << t;
        if(m_stat != PROGRESS){
            m_logFile << t;
        }
        return *m_pThis;
    };

    int precision(int p);
    int precision();
    string setprecision(int p);

    void open(string ofile);
    void close();
    void flush();
    static Logger* GetLogger();

private:
    static std::map<Logger::Type, string> style_map;
    static std::map<string, std::chrono::time_point<std::chrono::steady_clock>> time_map;
    static string empty;

    static std::mutex log_mutex;

    Logger();
    Logger(const Logger&){};
    Logger& operator=(const Logger&){return *this;};
    bool check();
    static Type m_stat;
    static string m_FileName;
    static Logger* m_pThis;
    static std::ofstream m_logFile;
};

#endif //GCTA2_LOGGER_H
