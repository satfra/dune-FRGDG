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

#include <dune/FRGDG/common/utils.hh>

namespace Dune {
  namespace PDELab {
    /** Spatial local operator for discontinuous Galerkin method
      for the system of hyperbolic conservation laws:

      \nabla \cdot \{ F(u) \}  = 0 in \Omega

      Where u = (u1,...,um) is the solution with m components

      - Assumes that the local function space is a power space
      with m identical components.
      - Assumes Galerkin method, i.e. U=V

      \tparam FEM Finite Element Map needed to select the cache
      \tparam SET Simulation Set class from which we pull some types
      \tparam DATA Data class needed for some schemes
      */

    template<typename FEM, typename SET, typename DATA>
      class DGConservationEqSpatialOperator :
        public NumericalJacobianApplyVolume<DGConservationEqSpatialOperator<FEM,SET,DATA> >,
        public NumericalJacobianVolume<DGConservationEqSpatialOperator<FEM,SET,DATA> >,
        public NumericalJacobianApplySkeleton<DGConservationEqSpatialOperator<FEM,SET,DATA> >,
        public NumericalJacobianSkeleton<DGConservationEqSpatialOperator<FEM,SET,DATA> >,
        public NumericalJacobianApplyBoundary<DGConservationEqSpatialOperator<FEM,SET,DATA> >,
        public NumericalJacobianBoundary<DGConservationEqSpatialOperator<FEM,SET,DATA> >,
        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename SET::Traits::RF>
    {
      private:
        using RF = typename SET::Traits::RF;
        using GV = typename SET::Traits::GV;

        using NUMFLUX = typename SET::Traits::template Numflux<0>;

        static constexpr unsigned int dim = NUMFLUX::dim;
        static constexpr unsigned int m = NUMFLUX::m; //number of components

        NUMFLUX& numflux;

        unsigned int overintegration;

        using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
        using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
        std::vector<Cache> cache;

        mutable RF alpha_overall;

        RF cur_dt;
        RF cur_time;
        RF post_stage_time;
        int stages;
        int cur_stage;
        
        // dummy thing. we don't use it
        DATA& data;

      public:
        // pattern assembly flags
        enum { doPatternVolume = true };
        enum { doPatternSkeleton = true };

        // residual assembly flags
        enum { doAlphaVolume  = true };
        enum { doAlphaSkeleton  = true };
        enum { doAlphaBoundary  = true };
        enum { doLambdaVolume  = false };

        DGConservationEqSpatialOperator (NUMFLUX& numflux_, DATA& data_, unsigned int overintegration_=0)
          : numflux(numflux_), overintegration(overintegration_), cache(20), alpha_overall(1.), data(data_) {}

        RF getAlpha()
        {
          return alpha_overall;
        }

        void preStep (RF time_, RF dt_, int stages_)
        {
          alpha_overall = 0.;

          cur_time = time_;
          cur_dt = dt_;
          stages = stages_;

          numflux.model().setTime(time_);
        }

        void postStep()
        {
          numflux.model().setTime(cur_time + cur_dt);
        }

        void preStage (RF time_, int r_)
        {
          cur_stage = r_;
          post_stage_time = time_;
        }

        void postStage()
        {
          numflux.model().setTime(post_stage_time);
        }

        // volume integral depending on test and ansatz functions
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            using namespace Indices;

            // get local function space that is identical for all components
            const auto& dgspace = child(lfsv,_0);

            // Get geometry
            auto geo = eg.geometry();
            const auto& cell = eg.entity();

            // Transformation
            typename EG::Geometry::JacobianInverseTransposed jac;

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m> u(0.0);
            std::vector<Dune::FieldVector<RF,dim>> gradphi(dgspace.size());

            // loop over quadrature points
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const unsigned int intorder = overintegration+2*order;

            for (const auto& ip : quadratureRule(geo,intorder))
            {
              const auto qp = ip.position();

              // evaluate basis functions
              auto& phi = cache[order].evaluateFunction(qp,dgspace.finiteElement().localBasis());

              // evaluate u
              u = 0.0;
              for (unsigned int k=0; k<m; k++) // for all components
                for (unsigned int j=0; j<dgspace.size(); j++) // for all basis functions
                  u[k] += x(lfsv.child(k),j)*phi[j];

              // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
              auto& js = cache[order].evaluateJacobian(qp,dgspace.finiteElement().localBasis());

              // compute global gradients
              jac = geo.jacobianInverseTransposed(qp);
              for (size_t i=0; i<dgspace.size(); i++)
                jac.mv(js[i][0],gradphi[i]);

              // flux term
              Dune::FieldMatrix<RF,m,dim> F(0.);
              numflux.model().flux(eg,qp,u,F);

              // source term 
              auto q = numflux.model().q(cell,qp, u);

              // Integrate
              auto factor = ip.weight() * geo.integrationElement(qp);

              for (unsigned int k=0; k<m; k++) // for all components
                for (unsigned int i=0; i<dgspace.size(); i++) // for all test functions of this component 
                  r.accumulate(lfsv.child(k),i, -q[k]*phi[i]*factor);
              // - F(u) \grad phi

              for (unsigned int k=0; k<dgspace.size(); k++)
                for (unsigned int i=0; i<m; i++)
                  for (unsigned int j=0; j<dim; j++)
                    r.accumulate(lfsv.child(i),k,-F[i][j]*gradphi[k][j]*factor);
            }
          }

        // skeleton integral depending on test and ansatz functions (for communication at boundaries between single elements)
        // each face is only visited ONCE!
        template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_skeleton (const IG& ig,
              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
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
            const auto& geo_inside = cell_inside.geometry();
            //const auto& geo_outside = cell_outside.geometry();
            const auto& geo_in_inside = ig.geometryInInside();
            const auto& geo_in_outside = ig.geometryInOutside();

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m> u_s(0.0);
            Dune::FieldVector<RF,m> u_n(0.0);

            RF amax(0.0);

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
              for (unsigned int i=0; i<m; i++) // for all components
                for (unsigned int k=0; k<dgspace_s.size(); k++)
                  u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
              u_n = 0.0;
              for (unsigned int i=0; i<m; i++) // for all components
                for (unsigned int k=0; k<dgspace_n.size(); k++)
                  u_n[i] += x_n(lfsv_n.child(i),k)*phi_n[k];

              // Compute numerical flux at the integration point
            Dune::FieldVector<RF,m> f(0.0);
              numflux.numericalFlux(cell_inside, iplocal_s,
                  cell_outside, iplocal_n,
                  ig.unitOuterNormal(qp),u_s,u_n,axis,
                  f,amax);

              // Integrate
              const auto factor = ip.weight() * geo.integrationElement(qp);
              // loop over all vector-valued basis functions
              for (unsigned int k=0; k<dgspace_s.size(); k++)
                for (unsigned int i=0; i<m; i++) // loop over all components
                  r_s.accumulate(lfsv_s.child(i),k, f[i]*phi_s[k]*factor);
              // loop over all vector-valued basis functions
              for (unsigned int k=0; k<dgspace_n.size(); k++)
                for (unsigned int i=0; i<m; i++) // loop over all components
                  r_n.accumulate(lfsv_n.child(i),k, -f[i]*phi_n[k]*factor);
            }

            alpha_overall = std::max(alpha_overall, utils::divideBySmallestEdge(amax, geo_inside));
          }

        // Skeleton integral depending on test and ansatz functions

        //Boundary conditions
        template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_boundary (const IG& ig,
              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
              R& r_s) const
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
            const auto& geo_inside = cell_inside.geometry();
            const auto& geo_in_inside = ig.geometryInInside();

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m> u_s(0.0);

            RF amax(0.0);

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

              // Evaluate u from inside and outside
              u_s = 0.0;
              for (unsigned int i=0; i<m; i++) // for all components
                for (unsigned int k=0; k<dgspace_s.size(); k++) // for all basis functions
                  u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];

              // Evaluate boundary condition
              Dune::FieldVector<RF, m> u_n(numflux.model().g(cell_inside, iplocal_s, u_s));

              // Compute numerical flux at integration point
              Dune::FieldVector<RF,m> f(0.0);
              numflux.numericalFlux(cell_inside, iplocal_s,
                  cell_inside, iplocal_s,
                  ig.unitOuterNormal(qp),u_s,u_n,axis,
                  f,amax);

              // Evaluate Flux boundary condition
              Dune::FieldVector<RF, m> j(numflux.model().j(ig.intersection(), qp, u_s, f));

              const auto factor = ip.weight() * geo.integrationElement(qp);

              for (unsigned int k=0; k<dgspace_s.size(); k++)
                for (unsigned int i=0; i<m; i++) // loop over all components
                  r_s.accumulate(lfsv_s.child(i),k,j[i]*phi_s[k]*factor);

            }
            alpha_overall = std::max(alpha_overall, utils::divideBySmallestEdge(amax, geo_inside));
          }

        RF suggestTimestep (RF dt) const
        {
          return (dt / alpha_overall);
        }
    };
  }
}
