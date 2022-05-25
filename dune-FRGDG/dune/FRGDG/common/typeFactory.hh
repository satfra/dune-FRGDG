#pragma once

// dune
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/solver/newton.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

namespace utils
{
  template<typename GV, typename FEMDG, int m>
    class MakeGFSType
    {
      private:
        using RF = typename GV::ctype;

        using CON = Dune::PDELab::P0ParallelConstraints;
        using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;
        using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;

        // Single-Component Function
        using GFSDG = Dune::PDELab::GridFunctionSpace<GV, FEMDG, CON, VBE0>;
        // Power Multi-Component Function space for use by the operators
        using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
        using GFS = Dune::PDELab::PowerGridFunctionSpace<GFSDG, m, VBE, OrderingTag>;
        using X = Dune::PDELab::Backend::Vector<GFS,RF>;
      public:
        using Child = GFSDG;
        using Type = GFS;
        using Domain = X;
    };

  template<std::size_t N, typename... T>
    using static_switch = typename std::tuple_element<N, std::tuple<T...> >::type;

  template<std::size_t N, unsigned ...T>
    constexpr unsigned get_unsigned()
    {
      return std::get<N>(std::make_tuple(std::forward<int>(T)...));
    }
}
