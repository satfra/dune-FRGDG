#include <fstream>
#include <chrono>
#include <iomanip>
#include <iostream>

#include <dune/FRGDG/common/logger.hh>
#include <dune/FRGDG/common/utils.hh>

Logger& discard()
{
  static Logger d("","");
  return d;
}

void cleanString(std::string& str)
{
  //str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
}

Logger::Logger(std::string path_, std::string filename_) : path(utils::makePath(path_)), filename(filename_) 
{
  utils::createPathsRecursive(path);
}

void Logger::log(std::string msg)
{
  // discard the message if necessary
  if(path == "" || filename == "")
    return;
  cleanString(msg);
  std::string log = "<" + utils::time_date() + ">    " + msg + "\n";
  SaveLog(log);
}

void Logger::SaveLog(std::string& log)
{
  mutex.lock();
  std::fstream file(path + filename, file.out | file.app);
  if(!file)
    throw std::runtime_error("Failed to open file " + path + filename);
  file << log;
  file.close();
  mutex.unlock();
}
