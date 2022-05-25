#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <fstream>

class Logger
{
public:
  Logger(std::string outputPath_, std::string filename_);
  static Logger& discard();

  void log(std::string msg);

protected:
  std::mutex mutex;
  void SaveLog(std::string& log);
  std::string path;
  std::string filename;
};
