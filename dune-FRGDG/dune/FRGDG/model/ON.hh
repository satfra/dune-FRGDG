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

namespace ON
{
  using utils::powr;

  template <typename GV>
  class iModel : public ModelInterfaceLDGinstat<GV, 2, 1>
  {
    protected:
      using MI = ModelInterfaceLDGinstat<GV, 2, 1>;

    public:
      using RF = typename MI::RF;

      using MI::dim, MI::m;
      static constexpr unsigned m0 = MI::ms[0];
      static constexpr unsigned m1 = MI::ms[1];
      using Range0 = Dune::FieldVector<RF,m0>;
      using Range1 = Dune::FieldVector<RF,m1>;

    protected:
      static constexpr int N = 4;
      static constexpr int d = 4;

      using MI::k, MI::k2;
      using MI::ptree;

      RF l2, l4, l8;

      static constexpr unsigned factorial(unsigned i)
      {
        return i == 0 ? 1 : i * factorial(i-1);
      }

      RF A_d() const
      {
        const RF v_d = 2. * std::pow(M_PI, d/2.) / ( d * std::tgamma(d/2.) ); 
        return v_d * powr<d>(k) / ( 2. * powr<d>(2.*M_PI) );
      }

      template<int i, int deriv = 0>
        RF f(const RF& x) const
        {
          return A_d() * RF(i) / powr<deriv+1>(1. + x / k2) / (powr<deriv>(-k2) * factorial(deriv));
        }

    public:
      iModel(Dune::ParameterTree ptree_)
        : MI(ptree_) 
      {				
        l2 = ptree.get("param.l2", RF(0.));
        l4 = ptree.get("param.l4" , RF(71.6));
        l8 = ptree.get("param.l8" , RF(71.6));
      }

      static std::string getName(const Dune::ParameterTree& ptree)
      {
        return "_l2=" + ptree.template get<std::string>("param.l2") + "_l4=" + ptree.template get<std::string>("param.l4") + "_l8=" + ptree.template get<std::string>("param.l8");
      }

      template <typename E, typename X>
        Range0 u0(const E &e, const X &x) const
        {
          const X xg = e.geometry().global(x);
          Range0 u(0.0);
          u[0] = l2 + xg[0]*l4 + powr<3>(xg[0])*l8;
          u[1] = l4 + 3.*powr<2>(xg[0])*l8;
          return u;
        }

      template <typename E, typename X>
      std::vector<Dune::FieldMatrix<RF,m,m>> Jacobian(const E &cell, 
          const X &x,
          const Range0 &u) const
      {
        const X xg = cell.geometry().global(x);

        const RF duF = f<N-1,1>(u[0]) + f<1,1>(u[0] + 2.*xg[0]*u[1]);
        const RF dvF = f<1,1>(u[0] + 2.*xg[0]*u[1])*2.*xg[0];
        const RF duG = u[1] * f<N-1,2>(u[0]);
        const RF dvG = f<N-1,1>(u[0]);

        std::vector<Dune::FieldMatrix<RF,m,m>> res(dim);
        res[0][0][0] = duF;
        res[0][0][1] = dvF;
        res[0][1][0] = duG;
        res[0][1][1] = dvG;
        return res;
      }

      template <typename E, typename X, typename RF>
      void max_eigenvalue(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const Range0 &u_s, const Range0 &u_n,
          const Range1 &p_s, const Range1 &p_n,
          Dune::FieldMatrix<RF, m, dim> &alpha, RF &alpha_t) const
      {
        const auto Jacobian_s = Jacobian(inside, x_inside, u_s);
        const auto Jacobian_n = Jacobian(outside, x_outside, u_n);

        alpha[0][0] = std::max(std::abs(Jacobian_s[0][0][0]), std::abs(Jacobian_n[0][0][0]));
        alpha[1][0] = std::max(std::abs(Jacobian_s[0][1][1]), std::abs(Jacobian_n[0][1][1]));
        alpha_t = std::max(alpha[0][0], alpha[1][0]);
      }

      //Flux function
      template <typename E, typename X, typename RF>
      void flux(const E &e, const X &x,
          const Range0 &u,
          const Range1 &p,
          Dune::FieldMatrix<RF, m, dim> &Flux) const
      {
        const X xg = e.geometry().global(x);

        const RF F = f<N-1>(u[0]) + f<1>(u[0] + 2. * xg[0] *u[1]);
        const RF G = u[1] * f<N-1,1>(u[0]);

        Flux[0][0] = F;
        Flux[1][0] = G;
      }

      template <typename E, typename X>
        void numericalDiffFlux(const E &inside, const X &x_inside,
            const E &outside, const X &x_outside,
            const Range0 &u_s, const Range0 &u_n,
            const Range1 &p_s, const Range1 &p_n,
            Dune::FieldMatrix<RF, m, dim> &A,
            Dune::FieldMatrix<RF, m, dim> &beta, RF &alpha_t) const
        { 
          const X xg_s = inside.geometry().global(x_inside);
          const X xg_n = outside.geometry().global(x_outside);

          const RF s_s = u_s[0] + 2.*xg_s[0] * u_s[1];
          const RF s_n = u_n[0] + 2.*xg_n[0] * u_n[1];

          const RF diffusion = 0.5*(std::sqrt(A_d())*k/(k2 + s_n) 
              + std::sqrt(A_d())*k/(k2 + s_s));

          if (!utils::isEqual(s_s, s_n))
						A[1][0] = - std::sqrt(A_d())*k * std::log((k2 + s_s)/(k2 + s_n))/(s_s - s_n) * 0.5*(p_s[0] + p_n[0]);
          else
            A[1][0] = - diffusion * 0.5 * (p_s[0] + p_n[0]);
          beta[1][0] = diffusion*(p_n[0] - p_s[0]);
        }

      template <typename E, typename X>
        void diffFlux(const E &e, const X &x,
            const Range0 &u, const Range1 &p,
            Dune::FieldMatrix<RF, m, dim> &A) const
        {
          const X xg = e.geometry().global(x);
          const RF s = u[0] + 2.*xg[0]*u[1];

          A[1][0] = - std::sqrt(A_d()) * k / (k2 + s) * p[0];
        }
  };

  template <typename GV>
    class sModel : public ModelInterfaceLDGstat<GV, 2, 1>
  {
    protected:
      using MI = ModelInterfaceLDGstat<GV, 2, 1>;

    public:
      using RF = typename MI::RF;

      using MI::dim, MI::m;
      using MI::depm, MI::curm;
      static constexpr unsigned m0 = depm;
      static constexpr unsigned m1 = curm;
      using Range0 = Dune::FieldVector<RF,m0>;
      using Range1 = Dune::FieldVector<RF,m1>;

    protected:
      using MI::k2, MI::k;

      static constexpr int d = 4;

      RF A_d() const
      {
        const RF v_d = 2. * std::pow(M_PI, d/2.) / ( d * std::tgamma(d/2.) ); 
        return v_d * powr<d>(k) / ( 2. * powr<d>(2.*M_PI) );
      }

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
          const X xg = e.geometry().global(x);
          const RF s = u[0] + 2.*xg[0]*u[1];

          if (s > 0)
						F[0][0] = - std::sqrt(A_d())*k * std::log(k2 + s);
					else 
						F[0][0] =  std::sqrt(A_d())*k * std::log(1./(k2 + s));
        }

      template <typename E, typename X>
        void numericalDiffFlux(const E &inside, const X &x_inside,
            const E &outside, const X &x_outside,
            const Range0 &u_s, const Range0 &u_n,
            Dune::FieldMatrix<RF, m, dim> &beta) const
        {
          const X xg_s = inside.geometry().global(x_inside);
          const X xg_n = outside.geometry().global(x_outside);

          const RF s_s = u_s[0] + 2.*xg_s[0]*u_s[1];
          const RF s_n = u_n[0] + 2.*xg_n[0]*u_n[1];

          const RF diffusion = std::sqrt(A_d())*k * 0.5*(1/(k2 + s_n) + 1/(k2 + s_s));

          beta[0][0] = - diffusion * (u_n[1] - u_s[1]);
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
