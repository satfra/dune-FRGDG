#pragma once

#include<vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

namespace Dune {
  namespace PDELab {
    /** a local operator for the mass operator of a vector valued lfs (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    template<typename FEM, typename SET>
      class DGConservationEqTemporalOperator :
        public NumericalJacobianApplyVolume<DGConservationEqTemporalOperator<FEM,SET> >,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename SET::Traits::RF>
    {
      private:
        using RF = typename SET::Traits::RF;
        using NUMFLUX = typename SET::Traits::template Numflux<0>;

        static constexpr unsigned int dim = NUMFLUX::dim;
        static constexpr unsigned int m = NUMFLUX::m;

        NUMFLUX& numflux;
        int overintegration;

        using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
        using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
        std::vector<Cache> cache;

      public:
        // pattern assembly flags
        enum { doPatternVolume = true };

        // residual assembly flags
        enum { doAlphaVolume = true };

        DGConservationEqTemporalOperator (NUMFLUX& numflux_, unsigned int overintegration_=0)
          : numflux(numflux_), overintegration(overintegration_), cache(20)
        {}

        // define sparsity pattern of operator representation
        template<typename LFSU, typename LFSV, typename LocalPattern>
          void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
              LocalPattern& pattern) const
          {
            for (unsigned int k=0; k<TypeTree::degree(lfsv); k++)
              for (unsigned int i=0; i<lfsv.child(k).size(); ++i)
                for (unsigned int j=0; j<lfsu.child(k).size(); ++j)
                  pattern.addLink(lfsv.child(k),i,lfsu.child(k),j);
          }

        // volume integral depending on test and ansatz functions
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            // get types
            using namespace Indices;

            // get local function space that is identical for all components
            const auto& dgspace = child(lfsv,_0);

            // get geometry
            const auto& geo = eg.geometry();

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m> u(0.0);

            // loop over quadrature points
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+2*order;
            for (const auto& ip : quadratureRule(geo,intorder))
            {
              // evaluate basis functions
              const auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

              // evaluate u
              u = 0.0;
              for (unsigned int k=0; k<m; k++) // for all components
                for (unsigned int j=0; j<dgspace.size(); j++) // for all basis functions
                  u[k] += x(lfsv.child(k),j)*phi[j];

              // integrate
              const auto factor = ip.weight() * geo.integrationElement(ip.position());
              for (unsigned int k=0; k<m; k++) // for all components
                for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component
                  r.accumulate(lfsv.child(k),i, u[k]*phi[i]*factor);
            }
          }

        // jacobian of volume term
        template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
          void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
              M & mat) const
          {
            // get types
            using namespace Indices;

            // get local function space that is identical for all components
            const auto& dgspace = child(lfsv,_0);

            // get geometry
            const auto& geo = eg.geometry();

            // loop over quadrature points
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+2*order;
            for (const auto& ip : quadratureRule(geo,intorder))
            {
              // evaluate basis functions
              const auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

              // integrate
              const auto factor = ip.weight() * geo.integrationElement(ip.position());
              if constexpr(SET::Traits::DIAGONAL_FEM)
                for (unsigned int k=0; k<m; k++) // for all components
                  for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component
                    mat.accumulate(lfsv.child(k),i,lfsu.child(k),i, phi[i]*phi[i]*factor);
              else
                for (unsigned int k=0; k<m; k++) // for all components
                  for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component
                    for (unsigned int j=0; j<dgspace.size(); j++) // for all ansatz functions of this component
                      mat.accumulate(lfsv.child(k),i,lfsu.child(k),j, phi[j]*phi[i]*factor);
            }
          }

        static constexpr RF REALBIGNUMBER = 1e50;

        RF suggestTimestep (RF dt) const
        {
          return REALBIGNUMBER;
        }
    };
  }
}
