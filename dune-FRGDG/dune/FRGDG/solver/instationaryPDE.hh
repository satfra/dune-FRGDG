#pragma once

#include <exception>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/solver/newton.hh>
#include <dune/pdelab/constraints/p0.hh>

#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/common/typeFactory.hh>
#include <dune/FRGDG/common/logger.hh>

#include <dune/FRGDG/solver/DGTimeController.hh>
#include <dune/FRGDG/solver/explicitonestepLDG.hh>

template<typename SET, int order, typename LOP, typename TLOP, typename NUMFLUX>
class InstationaryPDE_ex
{
  public:
    static constexpr int dim = SET::Traits::dim;
    static constexpr int m = NUMFLUX::m; 

    using RF = typename SET::Traits::RF;
    using GV = typename SET::Traits::GridConstructor::GV;
    using GridConstructor = typename SET::Traits::GridConstructor;
    using FEMDG = typename SET::Traits::template FEM<order>;

    using CON = Dune::PDELab::P0ParallelConstraints;
    using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;
    using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;

    // Power Multi-Component Function space for use by the operators
    using GFS = typename utils::template MakeGFSType<GV, FEMDG, m>::Type;

    // Matrix Backend to use
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    using C = typename GFS::template ConstraintsContainer<RF>::Type;
    using GO0 = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, C, C>;
    using GO1 = Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, RF, RF, RF, C, C>;
    using IGO = Dune::PDELab::OneStepGridOperator<GO0, GO1, false>;

    using TimeController = Dune::PDELab::DGTimeController<RF, IGO,LOP>;

    //using LS = Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS>;
    using LS = typename SET::Traits::template LS_instat<GO0,C>; //Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO>;

    using Domain = typename IGO::Traits::Domain;

    using OneStepMethod = Dune::PDELab::ExplicitOneStepMethodLDG<RF, IGO, LS, Domain, Domain, TimeController>;

  protected:
    // ptree parameters
    RF orig_prefactor;
    RF cur_prefactor;
    RF checkPointDistance;
    int verbosity;

    Dune::FieldVector<RF, dim> gridLength;
    std::array<int, dim> gridCells;
    RF polynomialDegree;
    int timesteppingOrder;
    RF timeGridDistance;
    RF limitTime;

    RF minTimeStep;
    RF orig_maxTimeStep;
    RF cur_maxTimeStep;
    int maxLoops;
    RF alphaExp;

    bool atTimeGridPoint;

    // non-config stuff
    NUMFLUX& numflux;
    LOP& lop;
    TLOP& tlop;
    RF alpha;

    Dune::ParameterTree& ptree;
    Logger& log;

    const GFS& gfs;

    // Create grid operators for spatial and temporal residuals, and construct a 
    // one-step grid operator for integrating the equations
    MBE mbe; 
    C cg;
    GO0 go0;
    GO1 go1;
    IGO igo;
    LS ls;

    Dune::PDELab::ExplicitEulerParameter<RF> method1;
    Dune::PDELab::HeunParameter<RF> method2;
    Dune::PDELab::Shu3Parameter<RF> method3;
    Dune::PDELab::RK4Parameter<RF> method4;

    Dune::PDELab::TimeSteppingParameterInterface<RF> *method;
    
    // For Explicit Methods, the timeController sets the correct time-step size
    std::shared_ptr<TimeController> timeController;
    // Initialize with maximal number of nonzeroes per row

    std::shared_ptr<OneStepMethod> osm;
  public:

    InstationaryPDE_ex(const GFS& gfs_, NUMFLUX& numflux_, LOP& lop_, TLOP& tlop_, Dune::ParameterTree &ptree_, Logger& log_ = Logger::discard()) : 
      numflux(numflux_), lop(lop_), tlop(tlop_), ptree(ptree_), log(log_),
      gfs(gfs_), mbe(2*dim +1), cg(), go0(gfs, cg, gfs, cg, lop, mbe), go1(gfs, cg, gfs, cg, tlop, mbe), 
      igo(go0, go1), ls(gfs)
    {
      gridLength = ptree.get<Dune::FieldVector<RF, dim>>("grid.L");
      gridCells = ptree.get<std::array<int, dim>>("grid.N");

      polynomialDegree = ptree.get("fem.degree", RF(1.));
      timesteppingOrder = ptree.get("fem.torder", int(1));

      timeGridDistance = ptree.get("time.timeGrid", RF(1e15));
      orig_prefactor = ptree.get("time.safetyfactor", RF(1.));
      cur_prefactor = orig_prefactor;
      checkPointDistance = ptree.get("time.checkPointDist", RF(0.25));
      orig_maxTimeStep = ptree.get("time.maxTimeStep", RF(1.0));
      cur_maxTimeStep = orig_maxTimeStep;
      atTimeGridPoint = false;
      minTimeStep = ptree.get("time.minTimeStep", RF(1e-14));

      verbosity = ptree.get<int>("solver.verbosity", int(0));

      Setup();
    }
    
    template<class Xpass>
    RF Step(RF& time, Domain& xold, Domain& x, RF prefactor, Xpass& xpass)
    {
      timeController->setMaxTimestep(cur_maxTimeStep);
      return osm->apply_LDG(time, prefactor, xold, x, xpass);
    }

    RF Step(RF& time, Domain& xold, Domain& x, RF prefactor)
    {
      timeController->setMaxTimestep(cur_maxTimeStep);
      return osm->apply(time, prefactor, xold, x);
    }
   
    bool IsAtTimeGridPoint()
    {
      return atTimeGridPoint;
    }

  protected:
    void Setup()
    {
      // get variables for cfl condition
      RF dx = GridConstructor::getMinimalCellEdge(ptree);

      timeController = std::shared_ptr<TimeController>(new TimeController(dx, polynomialDegree, timesteppingOrder, timeGridDistance, igo, atTimeGridPoint, 1.));
      timeController->setVerbosityLevel(verbosity);
      timeController->setMaxTimestep(orig_maxTimeStep);

      // select and prepare time-stepping scheme
      switch(timesteppingOrder)
      {
        case 1:
          method = &method1; // explicit Euler
          break;
        case 2:
          method = &method2; // Heun
          break;
        case 3:
          method = &method3; // Shu 3
          break;
        case 4:
          method = &method4; // RK4
          break;
        default:
          throw(std::runtime_error("fem.torder must be in [1,4]!"));
          break;
      }
      igo.setMethod(*method);

      // the One-Step Method will perform the integration of the equations towards the next time step
      osm = std::shared_ptr<OneStepMethod>( new OneStepMethod(*method, igo, ls, *timeController));
      osm->setVerbosityLevel(verbosity);
    }
};
