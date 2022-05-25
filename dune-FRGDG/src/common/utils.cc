// Standard Library
#include <algorithm>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <chrono>
#include <iomanip>
#define GetCurrentDir getcwd

// Own includes
#include <dune/FRGDG/common/utils.hh>

namespace utils
{
  std::string getCurrentWorkingDir() 
  {
    char* curr_working_dir = get_current_dir_name();
    return std::string(curr_working_dir);
  }
  /*
   * a set of utlities mainly for IO and file management of Simulations.
   */
  
  /*
   * pads a string to size pad, using the character paddingChar
   */
  std::string padTo(std::string str, const size_t pad, const char paddingChar)
  {
    if(pad > str.size())
      str.insert(0, pad - str.size(), paddingChar);
    return str;
  }
  
  bool pathExists(const std::string path)
  {    
    struct stat info;
    
    if(stat( path.c_str(), &info ) != 0)
      return 0;
    else if(info.st_mode & S_IFDIR)
      return 1;
    else
      return 0;
  }
  
  /*
   * add a trailing '/' to a string, in order for it to be in standard form of a path
   */
  std::string makePath(std::string path)
  {
    if(path.back() != '/')
      return path + '/';
    else
      return path;
  }
  
  /*
   * returns seperately every directory up to the highest level of path
   */
  std::vector<std::string> pathIncrements(std::string path)
  {
    std::vector<std::string> segments;
    std::istringstream f(path);
    std::string s;
    while (std::getline(f, s, '/'))
      segments.push_back(s);
    
    std::vector<std::string> increments(segments.size());
    for(unsigned int n = 0; n<segments.size(); ++n)
      std::for_each(increments.begin()+n, increments.end(), [segments,n](std::string& path) {path += segments[n] + "/";});
    
    return increments;
  }
  
  /*
   * creates the directory path, even if its parent directories should not exist
   */
  bool createPathsRecursive(std::string path)
  {
    bool success = true;
    
    path = makePath(path);
    std::vector<std::string> increments  = pathIncrements(path);

    for(auto pathIncrement : increments)
      if(!pathExists(pathIncrement))
        success &= (system((std::string("mkdir ") + pathIncrement + " >nul 2>nul").c_str()) == 0);
      
    return success;
  }
  
  /*
   * creates all given directory paths, even if their parent directories should not exist
   */
  bool createPathsRecursive(std::vector<std::string> paths)
  {
    bool ret = true;
    std::for_each(paths.begin(), paths.end(), [&ret](std::string path){ret &= createPathsRecursive(path);});
    return ret;
  }
  
  /*
   * creates a set of directories as subdirectories of outputPath. A number of steps directories are created,
   * with names starting at start and incrementing with stepSize
   */
  bool createIncrementalStructure(std::string outputPath, int start, int steps, int stepSize, bool check, const size_t pad)
  {
    std::vector<std::string> pathList(steps+1);
    int count = start;
    bool existing = false;
    std::for_each(pathList.begin(), pathList.end(), 
                  [outputPath, stepSize, &count, &existing, pad](std::string& path){ 
                    path += outputPath + padTo(std::to_string(count), pad) + "/";
                    count += stepSize;
                    existing |= pathExists(path);
                  });
    
    if(check && existing)
      return false;
    
    return createPathsRecursive(pathList);
  }
  
  std::string timeFormat(size_t time)
  {
    if(time < 120)
      return std::to_string(time) + "s";
    else if(time/60. < 60.)
      return std::to_string((int)std::round(time/60.)) + "min";
    else if(time/60./60. < 24.)
      return std::to_string((int)(time/60./60.)) + "h" + std::to_string((int)std::round(time/60.)-(int)(time/60./60.)*60) + "min";
    return std::to_string((int)(time/60./60./24.)) + "d" + std::to_string((int)(time/60./60.)-(int)(time/60./60./24.)*24) + "h" + std::to_string((int)std::round(time/60.)-(int)(time/60./60.)*60) + "min";
  }

  template<typename T>
  std::string time_in_HH_MM_SS_MMM(T now)
  {
    using namespace std::chrono;

    // get number of milliseconds for the current second
    // (remainder after division into seconds)
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    // convert to std::time_t in order to convert to std::tm (broken time)
    auto timer = system_clock::to_time_t(now);

    // convert to broken time
    std::tm bt = *std::localtime(&timer);

    std::ostringstream oss;

    oss << std::put_time(&bt, "%H:%M:%S"); // HH:MM:SS
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();

    return oss.str();
  }

  constexpr char _Timeformat[] = "%Y-%m-%d";

  std::string time_date(std::chrono::system_clock::time_point t)
  {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), _Timeformat);
    ss << " ";
    ss << time_in_HH_MM_SS_MMM(now);
    return ss.str();
  }

  bool after(std::string t1_s, std::string t2_s)
  {
    std::tm t1 = {};
    std::tm t2 = {};
    std::istringstream s1(t1_s);
    std::istringstream s2(t2_s);
    s1 >> std::get_time(&t1, _Timeformat);
    s2 >> std::get_time(&t2, _Timeformat);
    if( (t1.tm_year > t2.tm_year) &&
        (t1.tm_yday > t2.tm_yday) &&
        (t1.tm_hour > t2.tm_hour) &&
        (t1.tm_min > t2.tm_min) &&
        (t1.tm_sec > t2.tm_sec) )
      return true;
    return false;
  }

  void touchFile(std::string file)
  {
    std::ofstream ofs(file, std::ios::out);
  }
}
