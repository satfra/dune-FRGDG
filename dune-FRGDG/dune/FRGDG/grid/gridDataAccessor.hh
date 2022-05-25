#pragma once

#include <dune/pdelab/finiteelementmap/monomfem.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include <dune/FRGDG/common/utils.hh>

namespace Dune {
  namespace PDELab {
    template <typename GFS>
      class GridDataAccessor
      {
        public:
          using FEM = typename GFS::ChildType::Traits::FiniteElementMapType;
          using GV = typename GFS::Traits::GridViewType;
          static constexpr unsigned m = GFS::Traits::noChilds;
          using LeafIdxSet = typename GV::Grid::LeafIndexSet;
          using LeafIdx = typename GV::Grid::LeafIndexSet::IndexType;
          static constexpr unsigned dim = GV::dimension;
          using RF =typename GV::Grid::ctype;
          using ElementType = typename GV::Traits::template Codim<0>::Entity;
          using VertexType = typename GV::Traits::template Codim<dim>::Entity;
          using IntersectionIterator = typename GV::IntersectionIterator;
          using Intersection = typename IntersectionIterator::Intersection;
          using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
          using JacobianType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
          using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
          using FE = typename FEM::Traits::FiniteElementType;
          using X = Dune::PDELab::Backend::Vector<GFS,RF>;
          using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
          using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
          using Caches = std::vector<Cache>;

        protected:
          template<unsigned d>
            static std::array<unsigned, d> multiindex (unsigned i, unsigned k)
            {
              std::array<unsigned, d> alpha;
              for (unsigned j=0; j<d; j++)
              {
                alpha[j] = i % (k+1);
                i = i/(k+1);
              }
              return alpha;
            }
          const GFS& gfs;
          mutable LFS lfs;
          mutable LFSCache lfs_cache;
          std::vector<Cache> cache;

          FE fe;

          static GeometryType getCube() {
            return GeometryType(GeometryTypes::cube(dim));
          }

        public:
          GridDataAccessor(const GFS& gfs_) :
            gfs(gfs_), lfs(gfs), lfs_cache(lfs), cache(20) {}

          std::vector<Cache>& getCaches()
          {
            return cache;
          }
          
          template<typename ElementType>
          std::vector<RF> viewData(X& x, const ElementType& e) const
          {
            typename X::template LocalView<LFSCache> x_view(x);
            lfs.bind(e);
            lfs_cache.update();
            std::vector<RF> xl(this->lfs.size());
            x_view.bind(this->lfs_cache);
            x_view.read(xl);
            x_view.unbind();
            return xl;
          }

          const LFS& getLFS()
          {
            return lfs;
          }

          template<typename ElementType>
          void writeData(X& x, const ElementType& e, const std::vector<RF>& xl)
          {
            typename X::template LocalView<LFSCache> x_view(x);
            lfs.bind(e);
            lfs_cache.update();
            x_view.bind(lfs_cache);
            x_view.write(xl);
            x_view.unbind();
          }
         
          Dune::FieldVector<RF,m> smoothnessIndicator(std::vector<RF> xl, const ElementType& e) const
          {
            const auto& geo = e.geometry();
            const auto x_lc = geo.center();
            const auto jac = geo.jacobianTransposed(x_lc);
            const auto& dx = jac.diagonal();

            using JacInvTr = typename ElementType::Geometry::JacobianInverseTransposed;
            static_assert(std::is_same<JacInvTr, DiagonalMatrix<RF,dim>>::value, "Only implemented for axis-aligned grids");
            const JacInvTr jacinv = geo.jacobianInverseTransposed(x_lc);

            // Initialize vectors outside for loop
            Dune::FieldVector<RF,m> beta(0.0);

            using namespace Indices;
            const auto& dgspace = child(this->lfs,_0);

            // Loop over quadrature points
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const unsigned int intorder = 2*dgspace.finiteElement().localBasis().order();
            for (const auto& ip : quadratureRule(geo,intorder))
            {
              const auto qp = ip.position();
              const auto intfactor = ip.weight() * geo.integrationElement(qp);

              const auto totIndices = utils::powr<dim>(order);
              for (unsigned int i = 1; i < totIndices; ++i)
              {
                std::array<unsigned int,dim> derivOrder = multiindex<dim>(i,order);
                std::vector<Dune::FieldVector<RF,1>> dnphi;
                dgspace.finiteElement().localBasis().partial(derivOrder,qp, dnphi); 

                RF derivtransf = 1.;
                RF fac = 1.;
                for(unsigned int d = 0; d < dim; ++d)
                {
                  derivtransf *= std::pow(jacinv.diagonal(d), RF(derivOrder[d]));
                  fac *= std::pow(dx[d], RF(2*derivOrder[d]-1));
                }

                for(unsigned int c = 0; c < m; ++c)
                  for(unsigned int j = 0; j < dgspace.size(); ++j)
                    beta[c] += utils::powr<2>(derivtransf*dnphi[j][0]*xl[c*dgspace.size() + j]) * fac * intfactor;
              }
            }

            return beta;
          }
      };
  }
}
