#pragma once

#include <dune/localfunctions/monomial.hh>

#include <dune/pdelab/finiteelementmap/monomfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/common/logger.hh>
#include <dune/FRGDG/grid/gridDataAccessor.hh>

namespace Dune {
  namespace PDELab {
    template<unsigned dim>
      class TaylorStuff
      {
        public:
          static constexpr unsigned int numOrderK(unsigned int d, unsigned int K)
          {
            if (K == 0)
              return 1;
            else 
              if (d > 1)
              {
                if (K > 1)
                {
                  unsigned int val = 0;
                  for(unsigned int l = 1; l <= K; ++l)
                    val += numOrderK(d-1,K-l);
                  return val;
                }
                else
                  return 1;
              }
              else
                return K;
          }

          static constexpr unsigned int numTotal(unsigned int d, unsigned int K, unsigned int Kmin = 0)
          {
            unsigned int val = 0;
            for(unsigned int l = Kmin; l <= K; ++l)
              val += numOrderK(d,l);
            return val;
          }

          static constexpr unsigned int pow(unsigned int a, unsigned int b)
          {
            if (b > 1)
              return a*pow(a, b-1);
            else 
              return a;
          }

          static constexpr std::array<unsigned int, dim> polySig(unsigned int i, unsigned int K)
          {
            std::array<unsigned int, dim> value;
            // 0 0
            // 1 0 
            // 0 1 
            // 1 1
            // 2 0 
            // 0 2

            if (K == 2 && dim == 2)
            {
              switch(i)
              {
                case 0: 
                  value[0] = 0.; value[1] = 0.; break;
                case 1: 
                  value[0] = 1.; value[1] = 0.; break;
                case 2: 
                  value[0] = 0.; value[1] = 1.; break;
                case 3: 
                  value[0] = 1.; value[1] = 1.; break;
                case 4: 
                  value[0] = 2.; value[1] = 0.; break;
                case 5: 
                  value[0] = 0.; value[1] = 2.; break;
              }
            }
            else
            {
              for(unsigned int k = 1; k < K; ++k)
              {
                unsigned int cN = numTotal(dim, k-1, 0);
                if(i > cN)
                {
                  value[0];
                  break;
                }
              }
            }
            return value;
          }
      };
    template<typename GFS>
      class LinearInterpolationArgument
      {
        public:
          using GV = typename GFS::Traits::GridViewType;
          using RF = typename GV::ctype;
          using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
          using ElementType = typename GV::Traits::template Codim<0>::Entity;
          enum { dim = ElementType::Geometry::coorddimension };
          struct Traits
          {
            using RangeType = RF;
          };

          LinearInterpolationArgument (const ElementType& e_, const Dune::FieldVector<RF,dim+1>& coeff_) 
            : e(e_), coeff(coeff_)
          {
            center = e_.geometry().center();
          }

          template<typename DT>
            RF operator () (const DT& x) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              RF y = coeff[0];
              for (unsigned int i=0; i<dim; i++) 
                y += coeff[i+1]*(global[i]-center[i]);
              return y;
            }
/*
          template<typename DT, typename RT>
            inline void evaluate (const DT& x, RT& y) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              y = coeff[0];
              for (unsigned int i=0; i<dim; i++) 
                y += coeff[i+1]*(global[i]-center[i]);
            }
*/
        private:
          const ElementType& e;
          const Dune::FieldVector<RF,dim+1>& coeff;
          Dune::FieldVector<RF,dim> center;
      };

    template<typename GFS>
      class PointInterpolationArgument1d
      {
        public:
          using GV = typename GFS::Traits::GridViewType;
          using RF = typename GV::ctype;
          using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
          using ElementType = typename GV::Traits::template Codim<0>::Entity;
          enum { dim = ElementType::Geometry::coorddimension };
          struct Traits
          {
            using RangeType = RF;
          };
          using DataType = std::pair<Dune::FieldVector<RF,dim>,RF>;

          PointInterpolationArgument1d(const ElementType& e_, const std::vector<DataType>& coeff_) : e(e_), coeff(coeff_)
        {
        }

          RF findNearest(const Dune::FieldVector<RF,dim>& pos) const
          {
            RF ret = coeff[0].second;
            RF dist = Dune::FieldVector<RF,dim>(coeff[0].first-pos).two_norm();
            for(const auto& p : coeff)
            {
              RF curdist = Dune::FieldVector<RF,dim>(p.first-pos).two_norm();
              if(curdist < dist)
              {
                dist = curdist;
                ret = p.second;
              }
            }
            return ret;
          }

          template<typename DT>
            RF operator () (const DT& x) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              RF y = RF(findNearest(global));
              return y;
            }
/*
          template<typename DT, typename RT>
            inline void evaluate (const DT& x, RT& y) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              y = RT(findNearest(global));
            }
*/
        private:
          const ElementType& e;
          const std::vector<DataType>& coeff;
          Dune::FieldVector<RF,dim> center;
      };

    template<typename GFS>
      class QuadraticInterpolationArgument_1d
      {
        public:
          using GV = typename GFS::Traits::GridViewType;
          using RF = typename GV::ctype;
          using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
          using ElementType = typename GV::Traits::template Codim<0>::Entity;
          enum { dim = ElementType::Geometry::coorddimension };
          struct Traits
          {
            using RangeType = RF;
          };

          QuadraticInterpolationArgument_1d (const ElementType& e_, const Dune::FieldVector<RF,3>& coeff_) : e(e_), coeff(coeff_)
        {
          center = e.geometry().center();
        }

          template<typename DT>
            RF operator () (const DT& x) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              RF y = coeff[0];
              y += coeff[1]*(global[0]-center[0]);
              y += coeff[2]*utils::powr<2>(global[0]-center[0]);
              return y;
            }
/*
          template<typename DT, typename RT>
            inline void evaluate (const DT& x, RT& y) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              y = coeff[0];
              y += coeff[1]*(global[0]-center[0]);
              y += coeff[2]*utils::powr<2>(global[0]-center[0]);
            }
*/
        private:
          const ElementType& e;
          const Dune::FieldVector<RF,3>& coeff;
          Dune::FieldVector<RF,dim> center;
      };

    template<typename GFS>
      class QuadraticInterpolationArgument_2d
      {
        public:
          using GV = typename GFS::Traits::GridViewType;
          using RF = typename GV::ctype;
          using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
          using ElementType = typename GV::Traits::template Codim<0>::Entity;
          enum { dim = ElementType::Geometry::coorddimension };
          struct Traits
          {
            using RangeType = RF;
          };

          QuadraticInterpolationArgument_2d (const ElementType& e_, const Dune::FieldVector<RF,6>& coeff_) : e(e_), coeff(coeff_)
        {
          center = e.geometry().center();
        }

          template<typename DT>
            RF operator () (const DT& x) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              RF y = coeff[0];
              // linear coefficients
              for(unsigned int d = 0; d < dim; ++d)
                y += coeff[1 + d] * (global[d]-center[d]);
              // quadratic coefficients
              for(unsigned int d = 0; d < dim; ++d)
                y += coeff[1 + dim + d] * utils::powr<2>(global[d]-center[d]);
              // mixed coefficient
              y += coeff[1 + dim + dim] * (global[0]-center[0]) * (global[1]-center[1]);
              return y;
            }
/*
          template<typename DT, typename RT>
            inline void evaluate (const DT& x, RT& y) const
            {
              Dune::FieldVector<RF,dim> global = e.geometry().global(x);
              y = coeff[0];
              // linear coefficients
              for(unsigned int d = 0; d < dim; ++d)
                y += coeff[1 + d] * (global[d]-center[d]);
              // quadratic coefficients
              for(unsigned int d = 0; d < dim; ++d)
                y += coeff[1 + dim + d] * utils::powr<2>(global[d]-center[d]);
              // mixed coefficient
              y += coeff[1 + dim + dim] * (global[0]-center[0]) * (global[1]-center[1]);
            }
*/
        private:
          const ElementType& e;
          const Dune::FieldVector<RF,6>& coeff;
          Dune::FieldVector<RF,dim> center;
      };

    template<typename GFS>
      class AverageInterpolationArgument
      {
        public:
          using GV = typename GFS::Traits::GridViewType;
          using RF = typename GV::ctype;
          using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
          using FEM = typename GFS::ChildType::Traits::FiniteElementMapType;
          using ElementType = typename GV::Traits::template Codim<0>::Entity;
          enum { dim = ElementType::Geometry::coorddimension };
          struct Traits
          {
            using RangeType = RF;
          };
          using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
          using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;

          AverageInterpolationArgument (unsigned int c_, std::vector<Cache>& cache_, const LFS& lfs_, RF value_, RF lambda_, const std::vector<RF>& xl_)
            : c(c_), cache(cache_), lfs(lfs_), value(value_), lambda(lambda_), xl(xl_)
          {
          }

          template<typename DT>
            RF operator () (const DT& x) const
            {
              using namespace Indices;
              const auto& dgspace = child(this->lfs,_0);
              //Dune::FieldVector<RF,dim> xglobal = it->geometry().global(x);
              const unsigned int order = dgspace.finiteElement().localBasis().order();
              auto& phi = cache[order].evaluateFunction(x,dgspace.finiteElement().localBasis());
              RF y = 0;
              for (unsigned int i=0; i<dgspace.size(); i++)
                y += xl[i + dgspace.size()*c]*phi[i];
              y = (1.-lambda)*value + lambda*y;
              return y;
            }
/*
          template<typename DT, typename RT>
            inline void evaluate (const DT& x, RT& y) const
            {
              using namespace Indices;
              const auto& dgspace = child(this->lfs,_0);
              //Dune::FieldVector<RF,dim> xglobal = it->geometry().global(x);
              const unsigned int order = dgspace.finiteElement().localBasis().order();
              auto& phi = cache[order].evaluateFunction(x,dgspace.finiteElement().localBasis());
              y = 0;
              for (unsigned int i=0; i<dgspace.size(); i++)
                y += xl[i + dgspace.size()*c]*phi[i];
              y = (1.-lambda)*value + lambda*y;
            }
*/
        private:
          unsigned int c;
          std::vector<Cache>& cache;
          const LFS& lfs;
          RF value;
          RF lambda;
          const std::vector<RF>& xl;
      };
  }
}
