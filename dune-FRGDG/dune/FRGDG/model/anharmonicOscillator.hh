#pragma once

#include <cmath>

// this module
#include <dune/FRGDG/grid/gridConstructor.hh>
#include <dune/FRGDG/common/utils.hh>

#include <dune/FRGDG/numericalflux/LLFFlux.hh>
#include <dune/FRGDG/numericalflux/centralFlux.hh>

#include <dune/FRGDG/localoperator/LDG_LOP_ConEq.hh>
#include <dune/FRGDG/localoperator/DG_TLOP_ConEq.hh>

#include <dune/FRGDG/model/modelInterfaceLDG.hh>

#include <dune/FRGDG/solver/ovlpistlsolverbackend_extension.hh>

#include <dune/FRGDG/integrationscheme/LDGScheme.hh>

#include <dune/pdelab/finiteelementmap/qkdg.hh>

namespace anharmonicOscillator
{
  using utils::powr;

  template <typename GV>
    class iModel : public ModelInterfaceLDGinstat<GV, 1, 1>
  {
    protected:
      using MI = ModelInterfaceLDGinstat<GV, 1, 1>;

    public:
      using RF = typename MI::RF;

      using MI::dim, MI::m;
      static constexpr unsigned m0 = MI::ms[0];
      static constexpr unsigned m1 = MI::ms[1];
      using Range0 = Dune::FieldVector<RF,m0>;
      using Range1 = Dune::FieldVector<RF,m1>;

    protected:
      using MI::k, MI::k2;
      using MI::ptree;

      RF l, ml, v_0, A;
      bool riemann;

    public:
      iModel(Dune::ParameterTree ptree_)
        : MI(ptree_) 
      {				
        l  = ptree.get("param.l" , RF(1));
        ml = ptree.get("param.m" , RF(0));
        A  = ptree.get("param.A" , RF(1));
        v_0 = ptree.get("param.v_0" , RF(0));
        riemann = ptree.get("param.riemann", true);
      }

      static std::string getName(const Dune::ParameterTree& ptree)
      {
        return "_l=" + ptree.template get<std::string>("param.l") + "_ml=" + ptree.template get<std::string>("param.ml");
      }

      template <typename E, typename X>
        Range0 u0(const E &e, const X &x) const
        {
          const X xg = e.geometry().global(x);
          Range0 u(0.0);

          if(riemann)
            u[0] = xg[0] >= 5 ? A : v_0;
          else
            u[0] = 3. * powr<2>(xg[0]) * l/4. + ml;
          return u;
        }

      template <typename E, typename X>
        void numericalDiffFlux(const E &inside, const X &x_inside,
            const E &outside, const X &x_outside,
            const Range0 &u_s, const Range0 &u_n,
            const Range1 &p_s, const Range1 &p_n,
            Dune::FieldMatrix<RF, m, dim> &A,
            Dune::FieldMatrix<RF, m, dim> &beta, RF &alpha_t) const
        {
          const RF diffusion = k * 0.5 *(1./(k2 + u_n[0]) + 1./(k2 + u_s[0]));

          if (!utils::isEqual(u_s[0], u_n[0]))
            A[0][0] = - k * std::log((k2 + u_s[0]) / (k2 +u_n[0])) / (u_s[0] - u_n[0]) * 0.5 * (p_s[0]+p_n[0]);
          else
            A[0][0] = - diffusion * 0.5*(p_s[0] + p_n[0]);

          beta[0][0] = diffusion * (p_n[0] - p_s[0]);
          alpha_t = std::max(alpha_t, powr<2>(diffusion));
        }

      template <typename E, typename X>
        void diffFlux(const E &e, const X &x,
            const Range0 &u, const Range1 &p,
            Dune::FieldMatrix<RF, m, dim> &A) const
        {
          A[0][0] = - k / (k2 + u[0]) * p[0];
        }
  };

  template <typename GV>
    class sModel : public ModelInterfaceLDGstat<GV, 1, 1>
  {
    protected:
      using MI = ModelInterfaceLDGstat<GV, 1, 1>;

    public:
      using RF = typename MI::RF;

      using MI::dim, MI::m;
      using MI::depm, MI::curm;

      static constexpr unsigned m0 = depm;
      static constexpr unsigned m1 = curm;

      using Range0 = Dune::FieldVector<RF,m0>;
      using Range1 = Dune::FieldVector<RF,m1>;

    protected:
      using MI::k, MI::k2;

    public:
      sModel(Dune::ParameterTree ptree_)
        : MI(ptree_) 
      {
      }

      //Flux function
      template <typename E, typename X>
        void flux(const E &e, const X &x,
            const Range0 &u,
            Dune::FieldMatrix<RF, m1, dim> &F) const
        {
          if (u[0] > 0)
            F[0][0] = - k * std::log(k2 + u[0]);
          else 
            F[0][0] = k * std::log(1./(k2 + u[0]));
        }

      template <typename E, typename X>
        void numericalDiffFlux(const E &inside, const X &x_inside,
            const E &outside, const X &x_outside,
            const Range0 &u_s, const Range0 &u_n,
            Dune::FieldMatrix<RF, m, dim> &beta) const
        {
          const RF v_s = k / (k2 + u_s[0]);
          const RF v_n = k / (k2 + u_n[0]);
          const RF diffusion = 0.5 *( v_s + v_n );

          beta[0][0] = - diffusion * (u_n[0] - u_s[0]);
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
          static constexpr std::array<int,4> orders {1,2,3,4};

          template<int order> 
            using TLOP = Dune::PDELab::DGConservationEqTemporalOperator<FEM<order>, SimSet>;
          template<unsigned idx, unsigned order, typename DATA>
            using LOP = utils::static_switch<idx, Dune::PDELab::LDGConservationEqSpatialOperator<FEM<order>, SimSet, DATA>, Dune::PDELab::LDGConservationEqSpatialOperator_stat<FEM<order>, SimSet, DATA>>;

          template<unsigned idx>
            using Model = utils::static_switch<idx,iModel<GV>,sModel<GV>>;
          template<unsigned idx>
            using Numflux = utils::static_switch<idx,LLFfluxLDG<Model<idx>>,LLFfluxLDG<Model<idx>>>;

          template<typename GO, typename CC>
            using LS_stat = Dune::PDELab::ISTLBackend_OVLP_Diagonal<typename GO::Traits::TrialGridFunctionSpace>;
          template<typename GO, typename CC>
            using LS_instat = Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<typename GO::Traits::TrialGridFunctionSpace>;
          static constexpr bool linear = true;

          template<int order> using Scheme = LDGScheme<SimSet, order>;

      };
  };
}
