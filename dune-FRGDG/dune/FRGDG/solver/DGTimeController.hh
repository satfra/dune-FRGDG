#pragma once

#include <cmath>
#include <dune/pdelab/instationary/explicitonestep.hh>

namespace Dune {
  namespace PDELab {
    template<typename RF, class IGOS, class LOP>
      class DGTimeController : public TimeControllerInterface<RF>
    {
      private:
        static constexpr RF EPSILON = 1e-14;

        RF dx;
        RF polynomialDegree;
        RF methodOrder;
        RF gridDistance;
        const IGOS& igos;
        int verbosity;
        bool& atGridPoint;
        RF maxTimeStep;
        RF alphaExp;

      public:
        DGTimeController (RF dx_, RF polynomialDegree_, RF methodOrder_, RF gridDistance_, const IGOS& igos_, bool& atGridPoint_, RF alphaExp_ = 1., int verbosity_ = 0) 
          : dx(dx_), polynomialDegree(polynomialDegree_), methodOrder(methodOrder_), gridDistance(gridDistance_), igos(igos_), verbosity(verbosity_), atGridPoint(atGridPoint_), maxTimeStep(1e100), alphaExp(alphaExp_)
        {}

        virtual RF suggestTimestep (RF time, RF givendt) override
        {
          atGridPoint = false;
          // in our code, givendt is a prefactor, scaling the timestep; 
          // therefore we just apply the formula for a DG hyperbolic problem's timestep and multiply it by the prefactor
          RF dt = std::pow(givendt,alphaExp) / (2*this->polynomialDegree + 1);
          // divide the timestep by the maximum wavenumber - which is done by the spatial grid operator
          dt = igos.suggestTimestep(dt);

          if(dt > maxTimeStep)
            dt = maxTimeStep;

          // we want the code to snap to a time grid
          const RF mod = std::fmod(time,gridDistance);
          if(mod + dt >= gridDistance)
          {
            dt = gridDistance - mod + EPSILON;
            atGridPoint = true;
          }

          if(verbosity >= 1)
            std::cout << "Using timestep dt = " << dt << std::endl;

          return dt;
        }

        void setVerbosityLevel(const int& verbosity_)
        {
          verbosity = verbosity_;
        }

        void setMaxTimestep(const RF& max_)
        {
          maxTimeStep = max_;
        }
    };
  }
}
