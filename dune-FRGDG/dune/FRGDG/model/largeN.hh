#pragma once

#include <cmath>

// this module
#include <dune/FRGDG/math/thermodynamicEquations.hh>
#include <dune/FRGDG/grid/gridConstructor.hh>
#include <dune/FRGDG/common/utils.hh>

#include <dune/FRGDG/numericalflux/LLFFlux.hh>

#include <dune/FRGDG/localoperator/DG_LOP_ConEq.hh>
#include <dune/FRGDG/localoperator/DG_TLOP_ConEq.hh>

#include <dune/FRGDG/model/modelInterface.hh>

#include <dune/FRGDG/solver/ovlpistlsolverbackend_extension.hh>

#include <dune/FRGDG/integrationscheme/DGScheme.hh>

#include <dune/pdelab/finiteelementmap/qkdg.hh>

namespace largeN
{
  using namespace ThermodynamicEquations;

  template <typename GV>
  class iModel : public ModelInterfaceConEq<GV, 1>
  {
    protected:
      using MI = ModelInterfaceConEq<GV, 1>;

    public:
      using RF = typename MI::RF;
      using Range = typename MI::Range;

      using MI::dim, MI::m;

    protected:
      using MI::k2, MI::k5;
      using MI::ptree;

      RF l, hs;
      RF T, mu;

      static constexpr RF nf = 2.;
      
    public:
      iModel(Dune::ParameterTree ptree_)
        : MI(ptree_) 
      {
        T = ptree.get("param.T", RF(1e-3));
        mu = ptree.get("param.mu", RF(0.));

        l = ptree.get("param.l" , RF(71.6));
        hs = ptree.get("param.hs" , RF(3.6));
      }

      static std::string getName(const Dune::ParameterTree& ptree)
      {
        return "_mu=" + ptree.template get<std::string>("param.mu") + "_T=" + ptree.template get<std::string>("param.T");
      }

      template <typename E, typename X>
      Range u0(const E &e, const X &x) const
      {
        const X xg = e.geometry().global(x);
        return xg[0]*l/2.;
      }

      template <typename E, typename X>
      std::vector<Dune::FieldMatrix<RF,m,m>> Jacobian(const E &cell, 
          const X &x,
          const Range &u) const
      {
        std::vector<Dune::FieldMatrix<RF,m,dim>> res(dim);
        const RF ep = std::sqrt(u[0] + k2);
        const RF pion  = k5/(24.*pi2) * ( dcothS(ep,T) - cothS(ep,T)/ep ) / (ep*ep);  
        res[0][0][0] = (nf*nf)*pion;
        return res;
      }

      template <typename E, typename X, typename RF>
      void max_eigenvalue(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const Range &u_s, const Range &u_n,
          Dune::FieldMatrix<RF, m, dim> &alpha, RF &alpha_t) const
      {
        const auto Jacobian_s = Jacobian(inside, x_inside, u_s);
        const auto Jacobian_n = Jacobian(outside, x_outside, u_n);

        alpha[0][0] = std::max(std::abs(Jacobian_s[0][0][0]),std::abs(Jacobian_n[0][0][0]));
        alpha_t = alpha[0][0];
      }

      RF fluxPion(const Dune::FieldVector<RF, m> &u) const
      {
        const RF ep2 = u[0] + k2;
        const RF ep = std::sqrt(ep2);
        return k5/(12.*pi2) * cothS(ep,T)/ep;
      }

      template<typename X>
      RF fluxQuark(const X& xg) const
      {
        const RF eq = std::sqrt(2.*xg[0] * hs*hs + k2);
        return k5/(3.*pi2) * (-nf) * 0.5*(tanhS(eq+mu,T) + tanhS(eq-mu,T))/eq;
      }

      template <typename E, typename X, typename RF>
      void flux(const E &e, const X &x,
          const Range &u,
          Dune::FieldMatrix<RF, m, dim> &F) const
      {
        const X xg = e.geometry().global(x);
        F[0][0] = (nf*nf)*fluxPion(u) + 3.*fluxQuark(xg);
      }
  };

  class SimSet
  {
    public:
      class Traits
      {
        public:
          Traits() = delete;

          using RF = double;
          static constexpr unsigned dim = 1;
          using GridConstructor = CubicYaspGridConstructor<RF, dim>;
          using GV = typename GridConstructor::GV;

          template<int order>
            using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<RF, double, order, dim, Dune::PDELab::QkDGBasisPolynomial::legendre>;
          static constexpr bool DIAGONAL_FEM = true;
          template<typename GO, typename CC>
            using LS_stat = Dune::PDELab::ISTLBackend_OVLP_Diagonal<typename GO::Traits::TrialGridFunctionSpace>;
          template<typename GO, typename CC>
            using LS_instat = Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<typename GO::Traits::TrialGridFunctionSpace>;
          static constexpr bool linear = true;

          template<int order> using Scheme = DGScheme<SimSet, order>;

          template<int order> 
            using TLOP = Dune::PDELab::DGConservationEqTemporalOperator<FEM<order>, SimSet>;
          template<unsigned idx, unsigned order, typename DATA>
            using LOP = utils::static_switch<idx,Dune::PDELab::DGConservationEqSpatialOperator<FEM<order>, SimSet, DATA>>;

          template<unsigned idx>
            using Model = utils::static_switch<idx,iModel<GV>>;
          template<unsigned idx>
            using Numflux = utils::static_switch<idx,LLFflux<Model<idx>>>;

          static constexpr std::array<int,4> orders {1,2,3,4};
      };
  };
}
