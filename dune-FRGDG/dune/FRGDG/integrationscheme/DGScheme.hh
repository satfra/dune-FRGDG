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

template<typename GFS0_>
class DGDataSet
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

    static constexpr unsigned int size  = 1;
    using type_0 = Domain0;

    DGDataSet(const GFS0& gfs0_)
      : gda0(gfs0_), gfs0(gfs0_), x0(gfs0, 0.0) {}

    DGDataSet& operator=(DGDataSet& rhs)
    {
      this->x0 = rhs.x0;
      return *this;
    }

    void empty()
    {
      x0 = Domain0(gfs0, 0.0);
    }

    template<int idx>
      typename std::enable_if<(idx == 0), typename GDA0::Caches&>::type
      getCaches()
      {
        return gda0.getCaches();
      }

    template<int idx>
      typename std::enable_if<(idx == 0), type_0&>::type
      getData()
      {
        return x0;
      }

    template<int idx>
      typename std::enable_if<(idx == 0), void>::type
      setData(const type_0& x)
      {
        x0 = x;
      }

    template<int idx>
      typename std::enable_if<(idx == 0), std::vector<RF0>>::type
      viewData(const ElementType0& e)
      {
        return gda0.viewData(x0, e);
      }

    template<int idx>
      typename std::enable_if<(idx == 0), void>::type
      writeData(const ElementType0& e, 
          const std::vector<RF0>& xl)
      {
        gda0.writeData(x0, e, xl);
      }

    template<int idx>
      typename std::enable_if<(idx == 0), const GFS0&>::type
      getGFS()
      {
        return gfs0;
      }

  private:
    GDA0 gda0;
    const GFS0& gfs0;
    Domain0 x0;
};

template<typename SET, int order>
class DGScheme
{
  public:
    static constexpr unsigned int dim = SET::Traits::dim;

    using RF = typename SET::Traits::RF;
    using GV = typename SET::Traits::GV;
    using FEMDG = typename SET::Traits::template FEM<order>;
    using CommType = Dune::Communication<MPI_Comm>;

    using MODEL0 = typename SET::Traits::template Model<0>;

    using NUMFLUX0 = typename SET::Traits::template Numflux<0>;

    static constexpr unsigned int m0 = MODEL0::m;

    using GFS0 = typename utils::template MakeGFSType<GV, FEMDG, m0>::Type;
    using GFSDG0 = typename utils::template MakeGFSType<GV, FEMDG, m0>::Child;

    using DataStorage = DGDataSet<GFS0>;

    using LOP0 = typename SET::Traits::template LOP<0, order, DataStorage>;

    using TLOP = typename SET::Traits::template TLOP<order>;

    using PDESolver_ex = InstationaryPDE_ex<SET, order, LOP0, TLOP, NUMFLUX0>;

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

    DataStorage data_prev;
    DataStorage data;
    DataStorage data_next;

    MODEL0 model0;
    NUMFLUX0 numflux0;

    LOP0 lop0;
    TLOP tlop;

    PDESolver_ex pdeSolver_ex;

    REJECTOR rejector;

  public: 
    DGScheme(const GV& gv_, const FEMDG &femdg_, Dune::ParameterTree &_ptree, CommType mpicomm_ = CommType(), Logger& log_ = Logger::discard())
      : gv(gv_), femdg(femdg_), ptree(_ptree), mpicomm(mpicomm_), log(log_), 
        gfsdg0(gv_, femdg_), gfs0(gfsdg0),
        data_prev(gfs0), data(gfs0), data_next(gfs0),
        model0(ptree), numflux0(model0),
        lop0(numflux0,data,2*order), tlop(numflux0), 
        pdeSolver_ex(gfs0, numflux0, lop0, tlop, ptree, log),
        rejector(gfs0)
    {
      // create an initial condition and interpolate u0 by xold
      auto u0lambda = [&](const auto &i, const auto &x) { return numflux0.model().u0(i, x); };
      auto u0 = Dune::PDELab::makeGridFunctionFromCallable(gv, u0lambda);
      Dune::PDELab::interpolate(u0, gfs0, data_prev.template getData<0>());
    }

    RF getAlpha()
    {
      return lop0.getAlpha();
    }

    RF step(RF time, RF prefactor) 
    {
      //RF Step(RF& time, Domain& xold, Domain& x, RF prefactor, Xpass& xpass)
      const RF dt = pdeSolver_ex.Step(time, data_prev.template getData<0>(), data_next.template getData<0>(), prefactor);
      return dt;
    }

    void accept()
    {
      data_prev = data_next;
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

      return vtkSequenceWriter;
    }

    DataStorage& getDataStorage()
    {
      return data;
    }
};
