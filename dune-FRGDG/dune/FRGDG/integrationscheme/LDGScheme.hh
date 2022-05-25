#pragma once

#include <vector>
#include <mpi.h>

#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/grid/gridDataAccessor.hh>

#include <dune/FRGDG/solver/rejectStep.hh>
#include <dune/FRGDG/solver/instationaryPDE.hh>
#include <dune/FRGDG/solver/stationaryPDE.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

template<typename GFS0_, typename GFS1_>
class LDGDataSet
{
  public:
    using GFS0 = GFS0_;
    using GV0 = typename GFS0::Traits::GridViewType;
    static constexpr int m0 = GFS0::Traits::noChild;
    static constexpr int dim0 = GV0::dimension;
    using RF0 = typename GV0::Grid::ctype;
    using ElementType0 = typename GV0::Traits::template Codim<0>::Entity;
    using Domain0 = Dune::PDELab::Backend::Vector<GFS0,RF0>;
    using GDA0 = typename Dune::PDELab::GridDataAccessor<GFS0>;

    using GFS1 = GFS1_;
    using GV1 = typename GFS1::Traits::GridViewType;
    static constexpr int m1 = GFS1::Traits::noChild;
    static constexpr int dim1 = GV1::dimension;
    using RF1 = typename GV1::Grid::ctype;
    using ElementType1 = typename GV1::Traits::template Codim<0>::Entity;
    using Domain1 = Dune::PDELab::Backend::Vector<GFS1, RF1>;
    using GDA1 = typename Dune::PDELab::GridDataAccessor<GFS1>;

    static constexpr unsigned int size  = 2;
    using type_0 = Domain0;
    using type_1 = Domain1;

    LDGDataSet(const GFS0& gfs0_, const GFS1& gfs1_)
      : gda0(gfs0_), gfs0(gfs0_), x0(gfs0, 0.0),
      gda1(gfs1_), gfs1(gfs1_), x1(gfs1, 0.0) {}

    LDGDataSet& operator=(LDGDataSet& rhs)
    {
      this->x0 = rhs.x0;
      this->x1 = rhs.x1;
      return *this;
    }

    void empty()
    {
      x0 = Domain0(gfs0, 0.0);
      x1 = Domain1(gfs1, 0.0);
    }

    template<int idx>
      typename std::enable_if<(idx == 0) || (idx == 1), utils::static_switch<idx,type_0&,type_1&>>::type
      getData()
      {
        if constexpr(idx == 0)
          return x0;
        else if constexpr(idx == 1)
          return x1;
      }

    template<int idx>
      typename std::enable_if<(idx == 0) || (idx == 1), typename utils::static_switch<idx,
               typename GDA0::Caches&,
               typename GDA1::Caches&
                 >>::type
                 getCaches()
                 {
                   if constexpr(idx == 0)
                     return gda0.getCaches();
                   else if constexpr(idx == 1)
                     return gda1.getCaches();
                 }

    template<int idx>
      typename std::enable_if<(idx == 0) || (idx == 1), void>::type
      setData(utils::static_switch<idx,const type_0&,const type_1&> x)
      {
        if constexpr(idx == 0)
          x0 = x;
        else if constexpr(idx == 1)
          x1 = x;
      }

    template<int idx>
      typename std::enable_if<(idx == 0) || (idx == 1), std::vector<utils::static_switch<idx,RF0,RF1>>>::type
      viewData(utils::static_switch<idx, const ElementType0&, const ElementType1&> e)
      {
        if constexpr(idx == 0)
          return gda0.viewData(x0, e);
        else if constexpr(idx == 1)
          return gda1.viewData(x1, e);
      }

    template<int idx>
      typename std::enable_if<(idx == 0) || (idx == 1), void>::type
      writeData(utils::static_switch<idx, const ElementType0&, const ElementType1&> e, 
          const std::vector<utils::static_switch<idx, RF0, RF1>>& xl)
      {
        if constexpr(idx == 0)
          gda0.writeData(x0, e, xl);
        else if constexpr(idx == 1)
          gda1.writeData(x1, e, xl);
      }

    template<int idx>
      typename std::enable_if<(idx == 0) || (idx == 1), utils::static_switch<idx,const GFS0&,const GFS1&>>::type
      getGFS()
      {
        if constexpr(idx == 0)
          return gfs0;
        else if constexpr(idx == 1)
          return gfs1;
      }

  private:
    GDA0 gda0;
    const GFS0& gfs0;
    Domain0 x0;

    GDA1 gda1;
    const GFS1& gfs1;
    Domain1 x1;
};

template<typename SET, int order>
class LDGScheme
{
  public:
    static constexpr unsigned int dim = SET::Traits::dim;

    using RF = typename SET::Traits::RF;
    using GV = typename SET::Traits::GV;
    using FEMDG = typename SET::Traits::template FEM<order>;
    using CommType = Dune::Communication<MPI_Comm>;

    using MODEL0 = typename SET::Traits::template Model<0>;
    using MODEL1 = typename SET::Traits::template Model<1>;

    using NUMFLUX0 = typename SET::Traits::template Numflux<0>;
    using NUMFLUX1 = typename SET::Traits::template Numflux<1>;

    static constexpr unsigned int m0 = MODEL0::m;
    static constexpr unsigned int m1 = MODEL1::m;

    using GFS0 = typename utils::template MakeGFSType<GV, FEMDG, m0>::Type;
    using GFSDG0 = typename utils::template MakeGFSType<GV, FEMDG, m0>::Child;
    using GFS1 = typename utils::template MakeGFSType<GV, FEMDG, m1>::Type;
    using GFSDG1 = typename utils::template MakeGFSType<GV, FEMDG, m1>::Child;

    using DataStorage = LDGDataSet<GFS0, GFS1>;

    using LOP0 = typename SET::Traits::template LOP<0, order, DataStorage>;
    using LOP1 = typename SET::Traits::template LOP<1, order, DataStorage>;

    using TLOP = typename SET::Traits::template TLOP<order>;

    using PDESolver_ex = InstationaryPDE_ex<SET, order, LOP0, TLOP, NUMFLUX0>;
    using StatSolver = StationaryPDE<SET, order, LOP1, NUMFLUX1>;

    using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
    using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;

    using REJECTOR = Dune::PDELab::RejectStep<GFS0>;

  private:
    const GV& gv;
    const FEMDG& femdg;
    Dune::ParameterTree ptree;

    CommType mpicomm;
    Logger& log;

    GFSDG0 gfsdg0;
    GFS0 gfs0;
    GFSDG1 gfsdg1;
    GFS1 gfs1;

    DataStorage data_prev;
    DataStorage data;
    DataStorage data_next;

    MODEL0 model0;
    NUMFLUX0 numflux0;
    MODEL1 model1;
    NUMFLUX1 numflux1;

    LOP0 lop0;
    LOP1 lop1;
    TLOP tlop;

    PDESolver_ex pdeSolver_ex;
    StatSolver statSolver;

    REJECTOR rejector;

  public: 
    LDGScheme(const GV& gv_, const FEMDG &femdg_, Dune::ParameterTree &_ptree, CommType mpicomm_ = CommType(), Logger& log_ = Logger::discard())
      : gv(gv_), femdg(femdg_), ptree(_ptree), mpicomm(mpicomm_), log(log_), 
        gfsdg0(gv_, femdg_), gfs0(gfsdg0), gfsdg1(gv_, femdg_), gfs1(gfsdg0), 
        data_prev(gfs0, gfs1), data(gfs0, gfs1), data_next(gfs0, gfs1),
        model0(ptree), numflux0(model0),
        model1(ptree), numflux1(model1),
        lop0(numflux0,numflux1,data,2*order), lop1(numflux0,numflux1,data,2*order), tlop(numflux0), 
        pdeSolver_ex(gfs0, numflux0, lop0, tlop, ptree, log),
        statSolver(gfs1, numflux1, lop1, ptree, log),
        rejector(gfs0)
    {
      // create an initial condition and interpolate u0 by xold
      auto u0lambda = [&](const auto &i, const auto &x) { return numflux0.model().u0(i, x); };
      auto u0 = Dune::PDELab::makeGridFunctionFromCallable(gv, u0lambda);
      Dune::PDELab::interpolate(u0, gfs0, data_prev.template getData<0>());
      data = data_prev;
    }

    RF getAlpha()
    {
      return lop0.getAlpha();
    }

    RF step(RF time, RF prefactor) 
    {
      auto ldgcon = [&](auto& x) {
        data.template setData<0>(x);
        statSolver.solve(data.template getData<1>());
      };

      //RF Step(RF& time, Domain& xold, Domain& x, RF prefactor, Xpass& xpass)
      const RF dt = pdeSolver_ex.Step(time, data_prev.template getData<0>(), data_next.template getData<0>(), prefactor, ldgcon);
      return dt;
    }

    void accept()
    {
      data_prev.template setData<0>(data_next.template getData<0>());
      data_prev.template setData<1>(data.template getData<1>());
      data = data_prev;
      data_next.empty();
    }

    bool isAtTimeGridPoint()
    {
      return pdeSolver_ex.IsAtTimeGridPoint();
    }

    bool checkConsistency()
    {
      // check whether the step done previously is valid
      int reject = rejector.CheckStep(data_prev.template getData<0>());
      mpicomm.max<int>(&reject, 1);
      return !bool(reject);
    }

    VTKSEQUENCEWRITER constructVTKWriter(std::string name, unsigned int subsampling)
    {
      // prepare VTK writer and write first file
      utils::createPathsRecursive(name);

      VTKWRITER vtkwriter(gv, Dune::refinementIntervals(subsampling));
      VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter), name, name, "");
      Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter, gfs0, data_prev.template getData<0>(), Dune::PDELab::vtk::DefaultFunctionNameGenerator("u"));
      Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter, gfs1, data_prev.template getData<1>(), Dune::PDELab::vtk::DefaultFunctionNameGenerator("g"));

      return vtkSequenceWriter;
    }

    DataStorage& getDataStorage()
    {
      return data;
    }
};
