#pragma once

#include <vector>

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
#include <dune/pdelab/backend/interface.hh>

#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/grid/gridDataAccessor.hh>

namespace Dune {
  namespace PDELab {
    /** Spatial local operator for local discontinuous Galerkin method
      for the system of stationary equations:

      \nabla \cdot \{ G(u) \}  = g in \Omega

      Where g = (g1,...,gm) is the solution with m components
      and u = (u1,...,un) is another vector with n components.

      - Assumes that the local function space is a power space
      with m identical components.
      - Assumes Galerkin method, i.e. U=V

      \tparam FEM Finite Element Map needed to select the cache
      \tparam SET Simulation Set class from which we pull some types
      \tparam DATA Data class needed for some schemes
      */
    template<typename FEM, typename SET, typename DATA, unsigned i0 = 0, unsigned i1 = 1>
      class LDGConservationEqSpatialOperator_stat :
        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename SET::Traits::RF>
    {
      private:
        using RF = typename SET::Traits::RF;
        using GV = typename SET::Traits::GV;

        using NUMFLUX0 = typename SET::Traits::template Numflux<i0>;
        using NUMFLUX1 = typename SET::Traits::template Numflux<i1>;

        static constexpr unsigned int dim = NUMFLUX1::dim;
        static constexpr unsigned int m0 = NUMFLUX0::m; //number of components
        static constexpr unsigned int m1 = NUMFLUX1::m; //number of components

        NUMFLUX0& numflux0;
        NUMFLUX1& numflux1;

        unsigned int overintegration;

        using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
        using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
        std::vector<Cache> cache;

        DATA& data;

      public:
        // pattern assembly flags
        enum { doPatternVolume = true };
        enum { doPatternSkeleton = true };

        // residual assembly flags
        enum { doAlphaVolume  = true };
        enum { doLambdaVolume  = true };
        enum { doLambdaSkeleton  = true };
        enum { doLambdaBoundary  = true };

        LDGConservationEqSpatialOperator_stat (NUMFLUX0& numflux0_, NUMFLUX1& numflux1_, DATA& data_, unsigned int overintegration_=0)
          : numflux0(numflux0_), numflux1(numflux1_), overintegration(overintegration_), cache(20), data(data_) {}

        // volume integral depending on test and ansatz functions
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            using namespace Indices;

            // get local function space that is identical for all components
            const auto& dgspace = child(lfsv,_0);

            // Get geometry
            const auto geo = eg.geometry();

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m1> p(0.0);

            // loop over quadrature points
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+2*order;

            for (const auto& ip : quadratureRule(geo,intorder))
            {
              const auto qp = ip.position();

              // evaluate basis functions
              auto& phi = cache[order].evaluateFunction(qp,dgspace.finiteElement().localBasis());

              // evaluate g
              p = 0.0;
              for (unsigned int k=0; k<m1; k++) // for all components
                for (unsigned int j=0; j<dgspace.size(); j++) // for all basis functions
                  p[k] += x(lfsu.child(k),j)*phi[j];

              // Integrate
              const auto factor = ip.weight() * geo.integrationElement(qp);

              for (unsigned int k=0; k<dgspace.size(); k++)
                for (unsigned int i=0; i<m1; i++)
                  r.accumulate(lfsv.child(i),k,p[i]*phi[k]*factor);
            }
          }

        // apply jacobian of volume term
        template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
          void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
          {
            alpha_volume(eg,lfsu,x,lfsv,y);
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
              {
                for (unsigned int k=0; k<m1; k++) // for all components
                  for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component
                    mat.accumulate(lfsv.child(k),i,lfsu.child(k),i, phi[i]*phi[i]*factor);
              }
              else
              {
                for (unsigned int k=0; k<m1; k++) // for all components
                  for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component
                    for (unsigned int j=0; j<dgspace.size(); j++) // for all ansatz functions of this component
                      mat.accumulate(lfsv.child(k),i,lfsu.child(k),j, phi[j]*phi[i]*factor);
              }
            }
          }

        // volume integral depending on test functions
        template<typename EG, typename LFSV, typename R>
          void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
          {
            using namespace Indices;

            // get local function space that is identical for all components
            const auto& dgspace = child(lfsv,_0);

            // Get geometry
            auto geo = eg.geometry();
            const auto& cell = eg.entity();

            // Transformation
            typename EG::Geometry::JacobianInverseTransposed jac;

            const auto xl0 = data.template viewData<i0>(cell);

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m0> u(0.0);
            std::vector<Dune::FieldVector<RF,dim>> gradphi(dgspace.size());
            Dune::FieldMatrix<RF,m1,dim> F;

            // loop over quadrature points
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+2*order;

            for (const auto& ip : quadratureRule(geo,intorder))
            {
              const auto qp = ip.position();

              // evaluate basis functions
              const auto& phi = cache[order].evaluateFunction(qp,dgspace.finiteElement().localBasis());

              // evaluate u
              u = 0.0;
              for (unsigned int k=0; k<m0; k++) // for all components
                for (unsigned int j=0; j<dgspace.size(); j++) // for all basis functions
                  u[k] += xl0[k*dgspace.size() + j]*phi[j];

              // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
              const auto& js = cache[order].evaluateJacobian(qp,dgspace.finiteElement().localBasis());

              // compute global gradients
              jac = geo.jacobianInverseTransposed(qp);
              for (size_t i=0; i<dgspace.size(); i++)
                jac.mv(js[i][0],gradphi[i]);

              // flux term
              numflux1.model().flux(eg,qp,u,F);
              const auto q = numflux0.model().q(cell,qp,u);

              // Integrate
              const auto factor = ip.weight() * geo.integrationElement(qp);

              for (unsigned int k=0; k<dgspace.size(); k++)
                for (unsigned int i=0; i<m1; i++)
                  for (unsigned int j=0; j<dim; j++)
                    r.accumulate(lfsv.child(i),k,-F[i][j]*gradphi[k][j]*factor);

              for (unsigned int k=0; k<m1; k++) // for all components
                for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component 
                  r.accumulate(lfsv.child(k),i, -q[k]*phi[i]*factor);
            }
          }

        // skeleton integral depending on test functions (for communication at boundaries between single elements)
        // each face is only visited ONCE!
        template<typename IG, typename LFSV, typename R>
          void lambda_skeleton (const IG& ig, const LFSV& lfsv_s,  const LFSV& lfsv_n,
              R& r_s, R& r_n) const
          {
            const int face = ig.indexInInside();
            const int axis = face / 2;
            // Get types
            using namespace Indices;

            // Get local function space that is identical for all components
            const auto& dgspace_s = child(lfsv_s,_0);
            const auto& dgspace_n = child(lfsv_n,_0);

            // References to inside and outside cells
            const auto& cell_inside = ig.inside();
            const auto& cell_outside = ig.outside();

            // Get geometries
            const auto& geo = ig.geometry();
            //const auto& geo_inside = cell_inside.geometry();
            //const auto& geo_outside = cell_outside.geometry();
            const auto& geo_in_inside = ig.geometryInInside();
            const auto& geo_in_outside = ig.geometryInOutside();

            const auto xl0_s = data.template viewData<i0>(cell_inside);
            const auto xl0_n = data.template viewData<i0>(cell_outside);

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m0> u_s(0.0);
            Dune::FieldVector<RF,m0> u_n(0.0);
            Dune::FieldVector<RF,m1> f(0.0);

            // Loop over quadrature points
            const unsigned int order_s = dgspace_s.finiteElement().localBasis().order();
            const unsigned int order_n = dgspace_n.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+1+2*std::max(order_s,order_n);

            for (const auto& ip : quadratureRule(geo,intorder))
            {
              const auto qp = ip.position();
              // Position of quadrature point in local coordinates of elements
              const auto iplocal_s = geo_in_inside.global(qp);
              const auto iplocal_n = geo_in_outside.global(qp);

              // Evaluate basis functions
              const auto& phi_s = cache[order_s].evaluateFunction(
                  iplocal_s,dgspace_s.finiteElement().localBasis());
              const auto& phi_n = cache[order_n].evaluateFunction(
                  iplocal_n,dgspace_n.finiteElement().localBasis());

              // Evaluate u from inside and outside
              u_s = 0.0;
              u_n = 0.0;
              for (unsigned int i=0; i<m0; i++) // for all components
              {
                for (unsigned int k=0; k<dgspace_s.size(); k++)
                {
                  u_s[i] += xl0_s[i*dgspace_s.size() + k]*phi_s[k];
                  u_n[i] += xl0_n[i*dgspace_n.size() + k]*phi_n[k];
                }
              }

              // Compute numerical flux at the integration point
              numflux1.numericalFlux(cell_inside, iplocal_s,
                  cell_outside, iplocal_n,
                  ig.unitOuterNormal(qp), u_s, u_n, axis,
                  f);

              // Integrate
              const auto factor = ip.weight() * geo.integrationElement(qp);
              // loop over all vector-valued basis functions
              for (unsigned int k=0; k<dgspace_s.size(); k++)
                for (unsigned int i=0; i<m1; i++) // loop over all components
                  r_s.accumulate(lfsv_s.child(i),k, f[i]*phi_s[k]*factor);
              // loop over all vector-valued basis functions
              for (unsigned int k=0; k<dgspace_n.size(); k++)
                for (unsigned int i=0; i<m1; i++) // loop over all components
                  r_n.accumulate(lfsv_n.child(i),k, -f[i]*phi_n[k]*factor);
            }
          }

        //Boundary conditions
        template<typename IG, typename LFSV, typename R>
          void lambda_boundary (const IG& ig, const LFSV& lfsv_s, R& r_s) const
          {
            const int face = ig.indexInInside();
            const int axis = face / 2;
            // Get types
            using namespace Indices;

            // Get local function space that is identical for all components
            const auto& dgspace_s = child(lfsv_s,_0);

            // Reference to inside cell
            const auto& cell_inside = ig.inside();

            // Get geometries
            const auto& geo = ig.geometry();
            //const auto& geo_inside = cell_inside.geometry();
            const auto& geo_in_inside = ig.geometryInInside();
            const auto xl0_s = data.template viewData<i0>(cell_inside);

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m0> u_s(0.0);
            Dune::FieldVector<RF,m1> f(0.0);

            // Loop over quadrature points
            const unsigned int order_s = dgspace_s.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+1+2*order_s;
            for (const auto& ip : quadratureRule(geo,intorder))
            {
              const auto qp = ip.position();
              //const auto nf_local = ig.unitOuterNormal(qp);
              // Position of quadrature point in local coordinates of elements
              const auto iplocal_s = geo_in_inside.global(qp);
              // Evaluate basis functions
              const auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,dgspace_s.finiteElement().localBasis());

              // Evaluate u from inside
              u_s = 0.0;
              for (unsigned int i=0; i<m0; i++) // for all components
                for (unsigned int k=0; k<dgspace_s.size(); k++) // for all basis functions
                  u_s[i] += xl0_s[i*dgspace_s.size() + k]*phi_s[k];
              // Evaluate p from inside

              // Evaluate boundary condition
              Dune::FieldVector<RF, m0> u_n(numflux0.model().g(cell_inside, iplocal_s, u_s));

              // Compute numerical flux at integration point
              numflux1.numericalFlux(cell_inside, iplocal_s,
                  cell_inside, iplocal_s,
                  ig.unitOuterNormal(qp),u_s,u_n,axis,
                  f);

              // Evaluate Flux boundary condition
              Dune::FieldVector<RF, m1> j(numflux1.model().j(ig.intersection(), qp, u_s, f));

              const auto factor = ip.weight() * geo.integrationElement(qp);

              for (unsigned int k=0; k<dgspace_s.size(); k++)
                for (unsigned int i=0; i<m1; i++) // loop over all components
                  r_s.accumulate(lfsv_s.child(i),k,j[i]*phi_s[k]*factor);

            }
          }
    };

    template<typename FEM, typename SET, typename DATA>
      using LDGConservationEqSpatialOperator_stat1 = LDGConservationEqSpatialOperator_stat<FEM,SET,DATA,0,1>;
    template<typename FEM, typename SET, typename DATA>
      using LDGConservationEqSpatialOperator_stat2 = LDGConservationEqSpatialOperator_stat<FEM,SET,DATA,1,2>;
  }
}
