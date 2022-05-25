#pragma once

#include <dune/common/parametertree.hh>

#include <gsl/gsl_integration.h> // GSL numerical integration

namespace utils
{
  template<typename RF>
  class Integrator
  {
    public:
      Integrator(Dune::ParameterTree ptree_)
        : ptree(ptree_)
      {
        gslRule = ptree.get("LDG.gslRule" , 1); // rule for gsl, just use 1
        gslLimit = ptree.get("LDG.gslLimit" , 100); // number of bisections before gsl throws an error
        gslRelativeError = ptree.get("LDG.gslRelativeError" , 1e-7); // relative error goal
        gslWorkspace = gsl_integration_workspace_alloc(gslLimit);
      }

      template<typename FUN>
        RF integrate(FUN f, RF lower, RF upper)
        {
          gsl_function_pp<FUN> Fp(f);
          gsl_function *F = static_cast<gsl_function*>(&Fp);
          RF result, error;
          gsl_integration_qag(F, // function
              lower, //lower boundary
              upper, // upper boundary
              0, // absolute error goal, 0 just gets ignored
              gslRelativeError, // relative error goal
              gslLimit, // section number limit
              gslRule, // which rule to use
              gslWorkspace, // workspace
              &result, // result
              &error //absolute error
              );
          return result;
        }

    private:
      template<typename F>  
        class gsl_function_pp : public gsl_function {
          public:
            gsl_function_pp(const F& func) : _func(func) {
              function = &gsl_function_pp::invoke;
              params=this;
            }
          private:
            const F& _func;
            static double invoke(double x, void *params) {
              return static_cast<gsl_function_pp*>(params)->_func(x);
            }
        };

      Dune::ParameterTree ptree;

      RF gslRelativeError; // relative error goal
      int gslLimit; // max number of bisections
      int gslRule; // gsl rule
      gsl_integration_workspace * gslWorkspace; // workspace
  };
}
