#pragma once

// this module
#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/common/logger.hh>

#include <dune/common/parametertree.hh>

// std
#include <string>

// In the Simulation Set, an alias like "using Limiters = Set_Limiters<SimSet, ...>" should be defined to 
// allow the solvers to access limiter functions. If no limiter should be compiled into the runtime,
// use simply "using Limiters = Set_Limiters<SimSet>"
// This class is recursively defined; each hierarchy level adds its own limiter to the vector of all limiters,
// and then lets its base class add its, too.
template<typename SimSet, template<typename, typename> class ...Limiter_>
class Set_Limiters
{
  public:
    Set_Limiters() = delete;

  template <unsigned idx, typename GFS>
  static std::vector<typename SimSet::Traits::template LimiterPtr<idx,GFS>> 
    getSlopeLimiters(const GFS& gfs, 
        typename SimSet::Traits::template Numflux<idx>& numflux, 
        Dune::ParameterTree ptree, Logger& log)
  {
    using LimiterPtr = typename SimSet::Traits::template LimiterPtr<idx,GFS>;
    std::vector<LimiterPtr> limiters;
    return limiters;
  }

  template <unsigned idx, typename GFS>
  static void 
    insertLimiter(const GFS& gfs, 
        typename SimSet::Traits::template Numflux<idx>& numflux, 
        Dune::ParameterTree ptree, Logger& log, 
        std::vector<typename SimSet::Traits::template LimiterPtr<idx,GFS>>& limiters)
    {
    }
};

template<typename SimSet, template<typename, typename> class curLimiter_, 
  template<typename, typename> class ...Limiter_>
class Set_Limiters<SimSet, curLimiter_, Limiter_...> : protected Set_Limiters<SimSet, Limiter_...>
{
  protected:
    using Base = Set_Limiters<SimSet, Limiter_...>;
  public:
    Set_Limiters() = delete;

  template <unsigned idx, typename GFS>
  static std::vector<typename SimSet::Traits::template LimiterPtr<idx,GFS>> 
    getSlopeLimiters(const GFS& gfs, 
        typename SimSet::Traits::template Numflux<idx>& numflux, 
        Dune::ParameterTree ptree, Logger& log)
  {
    using LimiterPtr = typename SimSet::Traits::template LimiterPtr<idx,GFS>;
    std::vector<LimiterPtr> limiters;
    insertLimiter<idx,GFS>(gfs, numflux, ptree, log, limiters);

    return limiters;
  }

  template <unsigned idx, typename GFS>
  static void 
    insertLimiter(const GFS& gfs, 
        typename SimSet::Traits::template Numflux<idx>& numflux, 
        Dune::ParameterTree ptree, Logger& log, 
        std::vector<typename SimSet::Traits::template LimiterPtr<idx,GFS>>& limiters)
    {
      using LimiterPtr = typename SimSet::Traits::template LimiterPtr<idx,GFS>;
      limiters.push_back(
          LimiterPtr(new curLimiter_<GFS, typename SimSet::Traits::template Numflux<idx>>(gfs, numflux, ptree, log)));
      Base::template insertLimiter<idx,GFS>(gfs, numflux, ptree, log, limiters);
    }
};
