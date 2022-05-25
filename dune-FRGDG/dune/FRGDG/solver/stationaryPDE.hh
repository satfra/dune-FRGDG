#pragma once

#include <exception>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/solver/newton.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/solver/newton.hh>

#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/common/typeFactory.hh>
#include <dune/FRGDG/common/logger.hh>

template<typename SET, unsigned order, typename LOP_, typename NUMFLUX>
class StationaryPDE
{
  public:
    static constexpr unsigned dim = SET::Traits::dim;
    static constexpr unsigned m = NUMFLUX::m; 

    using RF = typename SET::Traits::RF;
    using GV = typename SET::Traits::GridConstructor::GV;
    using FEMDG = typename SET::Traits::template FEM<order>;

    using CON = Dune::PDELab::P0ParallelConstraints;
    using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;
    using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;

    // Power Multi-Component Function space for use by the operators
    using GFS = typename utils::template MakeGFSType<GV, FEMDG, m>::Type;

    static constexpr unsigned idxModel = NUMFLUX::Model::idx;

    // Matrix Backend to use
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    using C = typename GFS::template ConstraintsContainer<RF>::Type;
    using LOP = LOP_;
    using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, C, C>;

    using Domain = typename GO::Traits::Domain;

    using LS = typename SET::Traits::template LS_stat<GO,C>; //Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO>;
    using SLP = Dune::PDELab::StationaryLinearProblemSolver<GO,LS,Domain>;
    using NLP = Dune::PDELab::NewtonMethod<GO, LS>;

  protected:
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

    Dune::ParameterTree& ptree;
    Logger& log;
    const GFS& gfs;
    
    // Initialize with maximal number of nonzeroes per row
    MBE mbe; 

    // Create grid operators for spatial and temporal residuals, and construct a 
    // one-step grid operator for integrating the equations
    C cg;
    GO go;
    std::shared_ptr<LS> ls;
    std::shared_ptr<SLP> slp;
    std::shared_ptr<NLP> nlp;

    bool firstRun;

    int LSIter;
    RF min_defect;
    RF reduction;
    unsigned idx;

  public:
    StationaryPDE(const GFS& gfs_, NUMFLUX& numflux_, LOP& lop_, Dune::ParameterTree &ptree_, Logger& log_ = Logger::discard(), unsigned idx_ = 1) : 
      numflux(numflux_), lop(lop_), ptree(ptree_), log(log_),
      gfs(gfs_), mbe(2*dim +1), go(gfs, cg, gfs, cg, lop, mbe),
      firstRun(true), idx(idx_)
    {
      gridLength = ptree.get<Dune::FieldVector<RF, dim>>("grid.L");
      gridCells = ptree.get<std::array<int, dim>>("grid.N");

      polynomialDegree = ptree.get("fem.degree", RF(1.));
      timesteppingOrder = ptree.get("fem.torder", int(1));

      timeGridDistance = ptree.get("time.timeGrid", RF(1e15));
      orig_prefactor = ptree.get("time.safetyfactor", RF(1.));
      cur_prefactor = orig_prefactor;
      checkPointDistance = ptree.get("time.checkPointDist", RF(0.25));
      alphaExp = ptree.get("time.alphaExp", RF(1.));
      orig_maxTimeStep = ptree.get("time.maxTimeStep", RF(1.0));
      cur_maxTimeStep = orig_maxTimeStep;
      atTimeGridPoint = false;
      minTimeStep = ptree.get("time.minTimeStep", RF(1e-14));

      verbosity = ptree.get<int>("solver.verbosity", int(0));

      LSIter = ptree.get("solver.SLSIter", int(500));
      min_defect = ptree.get("solver.min_defect", RF(1e-99));
      reduction = ptree.get("solver.reduction", RF(1e-10));

      setup();
    }
    
    void solve(Domain& x)
    {
      if constexpr ( SET::Traits::linear ) 
      {
        if(verbosity <= 0)
        {
          std::cout.setstate(std::ios_base::failbit);
          slp->apply(x, !firstRun);
          std::cout.clear();
        }
        else
          slp->apply(x, !firstRun);
      }
      else
      {
        if(verbosity <= 0)
        {
          std::cout.setstate(std::ios_base::failbit);
          nlp->apply(x);
          std::cout.clear();
        }
        else
          nlp->apply(x);
      }

      firstRun = false;
    }

  protected:
    void setup()
    {
      // Select a linear solver backend
      ls = std::shared_ptr<LS>(new LS(gfs));

      // Assemble the linear problem
      slp = std::shared_ptr<SLP>(new SLP(go,*ls,reduction,min_defect,verbosity));

      // Assemble the nonlinear problem
      nlp = std::shared_ptr<NLP>(new NLP(go,*ls));
    }
};

