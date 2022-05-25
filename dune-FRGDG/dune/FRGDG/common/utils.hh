#pragma once

// Standard Library
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <tuple>
#include <sstream>
#include <iostream>
#include <iomanip>

// Dune
#include <dune/common/fvector.hh>

namespace utils
{
  /*
   * a simple power function, use like powr<3>(4.23). Can be evaluated at compile time.
   */
  template <int n, typename RF>
    constexpr RF powr(const RF& x)
    {
      if constexpr(n == 0)
        return RF(1);
      else if constexpr(n > 1)
        return x*powr<n-1, RF>(x);
      else if constexpr(n < 1)
        return 1./powr<-n, RF>(x);
      else
        return x;
    }

  /*
   * a set of utlities mainly for IO and file management of Simulations.
   */

  /*
   * pads a string to size pad, using the character paddingChar
   */
  std::string padTo(std::string str, const size_t pad, const char paddingChar = '0');

  bool pathExists(const std::string path);

  std::string getCurrentWorkingDir();

  /*
   * add a trailing '/' to a string, in order for it to be in standard form of a path
   */
  std::string makePath(std::string path);

  /*
   * returns seperately every directory up to the highest level of path
   */
  std::vector<std::string> pathIncrements(std::string path);


  /*
   * creates the directory path, even if its parent directories should not exist
   */
  bool createPathsRecursive(std::string path);

  /*
   * creates all given directory paths, even if their parent directories should not exist
   */
  bool createPathsRecursive(std::vector<std::string> paths);

  /*
   * creates a set of directories as subdirectories of outputPath. A number of steps directories are created,
   * with names starting at start and incrementing with stepSize
   */
  bool createIncrementalStructure(std::string outputPath, int start, int steps, int stepSize, bool check, const size_t pad=6);

  /* 
   * formats time [s] as a string in in [m], [h], <[d]
   */
  std::string timeFormat(size_t time);

  bool after(std::string t1_s, std::string t2_s);

  std::string time_date(std::chrono::system_clock::time_point t = std::chrono::system_clock::now());

  void touchFile(std::string filename);

  template<typename T>
    bool isEqual(T a, T b, T eps_ = std::numeric_limits<T>::epsilon())
    {
      T diff = std::fabs(a - b);
      if (diff < std::fmax(std::fabs(a), std::fabs(b)) * eps_)
        return true;
      return false;
    }
  
  /*
   * Return number with fixed precision after the decimal point
   */
  template<typename T>
  std::string getWithPrecision(int precision, T number)
    {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(precision) << number;
      return stream.str();
    }
  
  /*
   * Return number with fixed significant digits
   */
  template<typename T>
  std::string getWithSignificantDigits(int digits, T number)
    {
      const int d = int(std::ceil(std::log10(number < 0 ? -number : number)));
      std::stringstream stream;
      stream << std::fixed << std::setprecision(std::max(digits-d,0)) << number;
      return stream.str();
    }
  
  /*
   * Returns amax divided by the smallest edge length of cell geometry geo. Only going to work for cubic cells.
   */
  template<typename RF, typename GEO>
  RF divideBySmallestEdge(RF& amax, GEO& geo)
    {
      RF dx(1e100);
      const auto lower = geo.global(Dune::FieldVector<RF,geo.mydimension>(0.));
      const auto upper = geo.global(Dune::FieldVector<RF,geo.mydimension>(1.));
      for (int i=0;i<geo.mydimension;i++)
        if (upper[i]-lower[i]<dx)
          dx=upper[i]-lower[i];
      
      return amax/dx;
    }
  template<typename GEO>
  double returnSmallestEdge(GEO& geo)
    {
      double dx(1e100);
      const auto lower = geo.global(Dune::FieldVector<double,geo.mydimension>(0.));
      const auto upper = geo.global(Dune::FieldVector<double,geo.mydimension>(1.));
      for (int i=0;i<geo.mydimension;i++)
        if (upper[i]-lower[i]<dx)
          dx=upper[i]-lower[i];
      
      return dx;
    }

  template<typename RF, typename F>
    RF integrate(F function, RF min, RF max, unsigned subdiv)
    {
      const RF len = max - min;
      const RF subdiv_len = len / RF(subdiv);
      RF result = 0.;
      for(unsigned i = 0; i < subdiv; ++i)
      {
        const RF pos = (RF(i)+0.5) * subdiv_len;
        result += function(pos);
      }
      return result;
    }
}
