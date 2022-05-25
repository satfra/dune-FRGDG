#pragma once

#include <mpi.h>

#include <dune/common/fvector.hh>

#include <dune/common/parallel/mpicommunication.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/FRGDG/common/utils.hh>

namespace Dune {
  namespace PDELab {

    template <typename DATA>
    class IntersectionFinder
    {
      static constexpr double EPSILON = 1e-14;
      private:
        using GV0 = typename DATA::GV0;
        using IndexSet = typename GV0::IndexSet;
        using ElementType = typename GV0::Traits::template Codim<0>::Entity;
        using ElementIterator = typename GV0::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator;
        using IntersectionIterator = typename GV0::IntersectionIterator;
        using Intersection = typename IntersectionIterator::Intersection;
        using GFS0 = typename DATA::GFS0;
        using LFS0 = Dune::PDELab::LocalFunctionSpace<GFS0>;
        using LFSCache = Dune::PDELab::LFSIndexCache<LFS0>;
        using size_type = typename LFS0::Traits::SizeType;
        using FEM = typename GFS0::ChildType::Traits::FiniteElementMap;
        using LB = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
        using RangeFieldType = typename LB::Traits::RangeFieldType;
        using RF = typename LB::Traits::RangeFieldType;
        using X = Dune::PDELab::Backend::Vector<GFS0,RangeFieldType>;
        using RangeType = typename LB::Traits::RangeType;
        using CommType = Dune::Communication<MPI_Comm>;

        enum { dim = GV0::dimension };
        static constexpr int m = dim;
        DATA& data;
        const GFS0& gfs0;
        LFS0 lfs0;

        CommType mpicomm;

      public:
        IntersectionFinder(DATA& data_, CommType mpicomm_ = CommType()) :
          data(data_), gfs0(data_.template getGFS<0>()), lfs0(gfs0), mpicomm(mpicomm_)
      {}
        //! update the construction values for a given timestep
        Dune::FieldVector<RangeFieldType, dim> getIntersection(Dune::FieldVector<RangeFieldType,m> val)
        {
          if constexpr (dim == 1)
            return getIntersection_1d(val);
          else
            return getIntersection_nd(val);
        }

        //! update the construction values for a given timestep
        Dune::FieldVector<RangeFieldType, dim> getIntersection_1d(Dune::FieldVector<RangeFieldType,m> val)
        {
          Dune::FieldVector<Dune::FieldVector<RangeFieldType, m>, 2*dim+1> values;
          Dune::FieldVector<Dune::FieldVector<RangeFieldType, dim>, 2*dim+1> xglobal;
          Dune::FieldVector<RangeFieldType, dim> slopes(0.0);
          std::array<std::vector<Dune::FieldVector<RangeFieldType,dim>>,m> cross_verts;

          auto& cache = data.template getCaches<0>();

          RF _tol = 100000000.;
          RF mindx = 1000000000.;

          // for each cell
          ElementIterator it = gfs0.gridView().template begin<0,Dune::Interior_Partition>();
          ElementIterator endit = gfs0.gridView().template end<0,Dune::Interior_Partition>();
          for (; it != endit; ++it)
          {
            const auto xl = data.template viewData<0>(*it);
            const auto geo = it->geometry();
            const auto center = it->geometry().center();

            Dune::FieldMatrix<RangeFieldType,dim,dim> jac = geo.jacobianTransposed(center);
            RF dx = jac[0][0];
            mindx = std::min(mindx,dx);
            _tol = std::min(dx/2., _tol);

            lfs0.bind(*it);
            using namespace Indices;
            auto dgspace = child(lfs0,_0);

            // evaluate basis functions
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            auto getMaxMin = [&](Dune::FieldVector<RF,dim>& min, Dune::FieldVector<RF,dim>& max, bool& noMin, const unsigned& c)
            {
              noMin = false;
              const auto mid = 0.5*(max + min);
              const auto phi_min = cache[order].evaluateFunction(geo.local(min),dgspace.finiteElement().localBasis());
              const auto phi_mid = cache[order].evaluateFunction(geo.local(mid),dgspace.finiteElement().localBasis());
              const auto phi_max = cache[order].evaluateFunction(geo.local(max),dgspace.finiteElement().localBasis());

              RF val_min = -val[c];
              RF val_mid = -val[c];
              RF val_max = -val[c];
              for (unsigned i=0; i<dgspace.size(); i++) 
              {
                val_min += xl[i + dgspace.size()*c]*phi_min[i];
                val_mid += xl[i + dgspace.size()*c]*phi_mid[i];
                val_max += xl[i + dgspace.size()*c]*phi_max[i];
              }

              if(val_min * val_mid < 0)
                max = mid;
              else if(val_max * val_mid < 0)
                min = mid;
              else
                noMin = true;
            };

            auto checkIntersection = [&](const auto& iit, const unsigned& c)
            {
              std::vector<RF> xl_n = data.template viewData<0>(iit->outside());

              const auto& itgeo = iit->geometry();
              const auto& ig = *iit;
              const auto& geo_in_inside = ig.inside().geometry();
              const auto& geo_in_outside = ig.outside().geometry();
              // Position of quadrature point in local coordinates of elements
              const auto pos_s = geo_in_inside.local(itgeo.center());
              const auto pos_n = geo_in_outside.local(itgeo.center());

              const auto phi_s = cache[order].evaluateFunction(pos_s,dgspace.finiteElement().localBasis());
              const auto phi_n = cache[order].evaluateFunction(pos_n,dgspace.finiteElement().localBasis());

              RF val_s = -val[c];
              RF val_n = -val[c];
              for (unsigned i=0; i<dgspace.size(); i++) 
              {
                val_s += xl[i + dgspace.size()*c]*phi_s[i];
                val_n += xl[i + dgspace.size()*c]*phi_n[i];
              }

              if(val_s * val_n < 0)
                return true;
              else 
                return false;
            };

            //check if there is a zero crossing
            for(unsigned c = 0; c < m; c++)
            {
              auto min = geo.corner(0);
              auto max = geo.corner(1);
              bool noMin_ = false;
              do
              {
                getMaxMin(min, max, noMin_, c);
              }
              while(std::abs(max[0]-min[0]) > dx/100. && !noMin_);

              if(!noMin_)
              {
                auto pos = 0.5*(min+max);
                cross_verts[c].push_back(pos);
              }
              else
              {
                IntersectionIterator endit = gfs0.gridView().iend(*it);
                IntersectionIterator iit = gfs0.gridView().ibegin(*it);
                for (; iit!=endit; ++iit)
                {
                  if(iit->neighbor())
                  {
                    if(checkIntersection(iit,c))
                    {
                      auto pos = iit->geometry().center();
                      cross_verts[c].push_back(pos);
                    }
                  }
                }
              }
            }
          }
          
          std::array<bool, m> has_intersection;
          std::vector<int> relevant_components;
          std::vector<int> irrelevant_components;
          std::vector<int> irrelevant_components_g;
          for(int c = 0; c < m; c++)
          {
            has_intersection[c] = (cross_verts[c].size() != 0);
            if(has_intersection[c])
              relevant_components.push_back(c);
            else
              irrelevant_components.push_back(c);
            bool has_intersection_g = mpicomm.max(int(has_intersection[c]));
            if(!has_intersection_g)
              irrelevant_components_g.push_back(c);
          }

          // first get rid of double countings
          for(const auto& c : relevant_components)
          {
            std::vector<Dune::FieldVector<RangeFieldType,dim>> new_cross;
            std::vector<bool> has_partner(cross_verts[c].size());
            for(unsigned int i = 0; i < cross_verts[c].size(); ++i)
              has_partner[i] = false;

            for(unsigned int i = 0; i < cross_verts[c].size()-1; ++i)
            {
              const auto& pos1 = cross_verts[c][i];
              for(unsigned int j = i+1; j < cross_verts[c].size(); ++j)
              {
                const auto& pos2 = cross_verts[c][j];
                const auto diff = pos2 - pos1;
                if(diff.two_norm() < _tol)
                {
                  has_partner[i] = true;
                  has_partner[j] = true;
                  new_cross.push_back(pos1 + diff/2.);
                }
              }
            }
            for(unsigned int i = 0; i < cross_verts[c].size(); ++i)
              if(!has_partner[i])
                new_cross.push_back(cross_verts[c][i]);

            cross_verts[c].clear();
            cross_verts[c].assign(new_cross.begin(), new_cross.end());
          }
          
          std::vector<Dune::FieldVector<RangeFieldType,dim>> crossings;
          std::vector<Dune::FieldVector<bool,m>> crossings_components;
          std::vector<RangeFieldType> crossings_dist;
          if(relevant_components.size() != 0)
          {
            crossings.assign(cross_verts[relevant_components.back()].begin(), cross_verts[relevant_components.back()].end());
            Dune::FieldVector<bool,m> comp(false);
            comp[relevant_components.back()] = true;
            for(unsigned int i = 0; i < crossings.size(); ++i)
              crossings_components.push_back(comp);
            crossings_dist.resize(crossings.size());
            relevant_components.pop_back();
          }

          // now see where the EOM are fulfilled at the same time in different components
          for(const auto& c : relevant_components)
          {
            std::vector<Dune::FieldVector<RangeFieldType,dim>> new_crossings;
            std::vector<RangeFieldType> new_crossings_dist;
            std::vector<Dune::FieldVector<bool,m>> new_crossings_components;

            for(unsigned int i = 0; i < crossings.size(); ++i)
            {
              const auto& pos1 = crossings[i];
              const auto& dist = crossings_dist[i];
              const auto& comp = crossings_components[i];

              for(unsigned int j = 0; j < cross_verts[c].size(); ++j)
              {
                const auto& pos2 = cross_verts[c][j];
                const auto diff = pos2 - pos1;
                if(diff.two_norm() < mindx/1.4)
                {
                  double newdist = dist + std::pow(diff.two_norm(),2);
                  new_crossings.push_back(pos1 + diff/2.);
                  new_crossings_dist.push_back(newdist);

                  Dune::FieldVector<bool,m> newcomp = comp;
                  newcomp[c] = true;
                  new_crossings_components.push_back(newcomp);
                }
              }
            }

            crossings.clear();
            crossings.assign(new_crossings.begin(), new_crossings.end());
            crossings_dist.clear();
            crossings_dist.assign(new_crossings_dist.begin(), new_crossings_dist.end());
            crossings_components.clear();
            crossings_components.assign(new_crossings_components.begin(), new_crossings_components.end());
          }

          for(unsigned int c = 0; c < m; ++c)
          {
            for(unsigned int i = 0; i < crossings.size(); ++i)
            {
              if(!crossings_components[i][c])
              {
                if(crossings[i][c] < mindx)
                  crossings[i][c] = 0.;
                else
                  crossings_dist[i] += std::pow(crossings[i][c],2);
              }
            }
          }

          Dune::FieldVector<RangeFieldType, dim> ret(0.);
          RangeFieldType dist = 1./EPSILON;
        
          unsigned int crossing_exists = crossings.size();
          crossing_exists = mpicomm.max(crossing_exists);
          if(crossing_exists == 0)
            return ret;
          
          unsigned int resulting = 0;
          for(unsigned int i = 1; i < crossings.size(); ++i)
            if(crossings_dist[i] < crossings_dist[resulting])
              resulting = i;
          if(crossings.size() != 0)
          {
            ret = crossings[resulting];
            dist = crossings_dist[resulting];
          }
          
          RangeFieldType distMin = mpicomm.min(dist);
          if(!utils::isEqual(dist,distMin,EPSILON))
            for(unsigned int c = 0; c < m; c++)
              ret[c] = 0.;

          auto retget = mpicomm.sum(ret);
          int fail = 0;
          if(utils::isEqual(dist,distMin,EPSILON))
            for(unsigned int c = 0; c < m; ++c)
              if(!utils::isEqual(retget[c],ret[c], EPSILON))
                fail = 1;
          fail = mpicomm.max(fail);
          if(fail == 1)
          {
            //std::cout << "-1: " << ret << " " << retget << "\n";
            return Dune::FieldVector<RangeFieldType, dim>(-1.);
          }

          return retget;
        }

        //! update the construction values for a given timestep
        Dune::FieldVector<RangeFieldType, dim> getIntersection_nd(Dune::FieldVector<RangeFieldType,m> val)
        {
          Dune::FieldVector<Dune::FieldVector<RangeFieldType, m>, 2*dim+1> values;
          Dune::FieldVector<Dune::FieldVector<RangeFieldType, dim>, 2*dim+1> xglobal;
          Dune::FieldVector<RangeFieldType, dim> slopes(0.0);
          std::array<std::vector<Dune::FieldVector<RangeFieldType,dim>>,m> cross_verts;

          auto& cache = data.template getCaches<0>();

          // for each cell
          ElementIterator it = gfs0.gridView().template begin<0,Dune::Interior_Partition>();
          ElementIterator endit = gfs0.gridView().template end<0,Dune::Interior_Partition>();
          Dune::FieldMatrix<RangeFieldType,dim,dim> jac = it->geometry().jacobianTransposed(it->geometry().center());
          double dx = jac[0][0];
          RangeFieldType _tol = dx/2.;
          for (; it != endit; ++it)
          {
            Dune::FieldMatrix<RangeFieldType,dim,dim> jac = it->geometry().jacobianTransposed(it->geometry().center());

            const auto xl = data.template viewData<0>(*it);

            lfs0.bind(*it);
            using namespace Indices;
            auto dgspace = child(lfs0,_0);

            const auto center = it->geometry().center();

            // evaluate basis functions
            const unsigned int order = dgspace.finiteElement().localBasis().order();
            const auto& phi = cache[order].evaluateFunction(center,dgspace.finiteElement().localBasis());
            // evaluate u the potential

            // compute cell value
            for (size_type c=0; c<m; c++)
              values[0][c] = -val[c];
            for (size_type i=0; i<dgspace.size(); i++)
              for (size_type c=0; c<m; c++)
                values[0][c] += xl[i + dgspace.size()*c]*phi[i + dgspace.size()*c];
            xglobal[0] = (*it).geometry().center();

            // compute neighbor values and check if its zero
            IntersectionIterator iit = gfs0.gridView().ibegin(*it);
            IntersectionIterator endiit = gfs0.gridView().iend(*it);
            for (; iit !=endiit; ++iit)
            {
              int face = iit->indexInInside();
              int axis = face/2;
              int dir = face%2;
              // evaluate only the intersection if there is another element owning the intersection
              if (iit->neighbor())
              {
                const auto xlnb = data.template viewData<0>(iit->outside());

                const auto& dgspace_nb = child(lfs0,_0);

                const auto center = iit->outside().geometry().center();

                const auto& phi = cache[order].evaluateFunction(center,dgspace_nb.finiteElement().localBasis());

                // evaluate u
                for (unsigned int c=0; c<m; c++)
                  values[1+axis*2+dir][c] = -val[c];
                for (unsigned int i=0; i<dgspace_nb.size(); i++) 
                  for (unsigned int c=0; c<m; c++)
                    values[1+axis*2+dir][c] += xlnb[i + dgspace_nb.size()*c]*phi[i];
                xglobal[1+axis*2+dir] = iit->outside().geometry().center();
              }
              else
              {    
                values[1+axis*2+dir] = values[0];
                xglobal[1+axis*2+dir] = xglobal[0];
              }
            }

            //check if there is a zero crossing
            for(int c = 0; c < m; c++)
            {
              for (int a = 0; a<dim; a++)
              {
                //check for zero crossing
                if (values[1+2*a][c]*values[1+2*a+1][c]<0)
                {
                  if (values[1+2*a][c]*values[0][c]<0)
                  {
                    auto diff = xglobal[0] - xglobal[1+2*a];
                    slopes[a] = (values[0][c] - values[1+2*a][c])/diff.two_norm();
                  }
                  else if (values[0][c]*values[1+2*a+1][c]<0)
                  {
                    auto diff = xglobal[0] - xglobal[1+2*a+1];
                    slopes[a] = (values[0][c] - values[1+2*a+1][c])/diff.two_norm();
                  }
                  auto pos = xglobal[0];
                  pos[a] = xglobal[0][a] - values[0][c]/slopes[a];  
                  for(int d = 0; d < dim; ++d)
                    if(pos[d] < jac[d][d])
                      pos[d] = 0;
                  cross_verts[c].push_back(pos);
                }
              }
            }
          }
          
          std::array<bool, m> has_intersection;
          std::vector<int> relevant_components;
          std::vector<int> irrelevant_components;
          std::vector<int> irrelevant_components_g;
          for(int c = 0; c < m; c++)
          {
            has_intersection[c] = (cross_verts[c].size() != 0);
            if(has_intersection[c])
              relevant_components.push_back(c);
            else
              irrelevant_components.push_back(c);
            bool has_intersection_g = mpicomm.max(int(has_intersection[c]));
            if(!has_intersection_g)
              irrelevant_components_g.push_back(c);
          }

          // first get rid of double countings
          for(const auto& c : relevant_components)
          {
            std::vector<Dune::FieldVector<RangeFieldType,dim>> new_cross;
            std::vector<bool> has_partner(cross_verts[c].size());
            for(unsigned int i = 0; i < cross_verts[c].size(); ++i)
              has_partner[i] = false;

            for(unsigned int i = 0; i < cross_verts[c].size()-1; ++i)
            {
              const auto& pos1 = cross_verts[c][i];
              for(unsigned int j = i+1; j < cross_verts[c].size(); ++j)
              {
                const auto& pos2 = cross_verts[c][j];
                const auto diff = pos2 - pos1;
                if(diff.two_norm() < _tol)
                {
                  has_partner[i] = true;
                  has_partner[j] = true;
                  new_cross.push_back(pos1 + diff/2.);
                }
              }
            }
            for(unsigned int i = 0; i < cross_verts[c].size(); ++i)
              if(!has_partner[i])
                new_cross.push_back(cross_verts[c][i]);

            cross_verts[c].clear();
            cross_verts[c].assign(new_cross.begin(), new_cross.end());
          }
          
          std::vector<Dune::FieldVector<RangeFieldType,dim>> crossings;
          std::vector<Dune::FieldVector<bool,m>> crossings_components;
          std::vector<RangeFieldType> crossings_dist;
          if(relevant_components.size() != 0)
          {
            crossings.assign(cross_verts[relevant_components.back()].begin(), cross_verts[relevant_components.back()].end());
            Dune::FieldVector<bool,m> comp(false);
            comp[relevant_components.back()] = true;
            for(unsigned int i = 0; i < crossings.size(); ++i)
              crossings_components.push_back(comp);
            crossings_dist.resize(crossings.size());
            relevant_components.pop_back();
          }

          // now see where the EOM are fulfilled at the same time in different components
          for(const auto& c : relevant_components)
          {
            std::vector<Dune::FieldVector<RangeFieldType,dim>> new_crossings;
            std::vector<RangeFieldType> new_crossings_dist;
            std::vector<Dune::FieldVector<bool,m>> new_crossings_components;

            for(unsigned int i = 0; i < crossings.size(); ++i)
            {
              const auto& pos1 = crossings[i];
              const auto& dist = crossings_dist[i];
              const auto& comp = crossings_components[i];

              for(unsigned int j = 0; j < cross_verts[c].size(); ++j)
              {
                const auto& pos2 = cross_verts[c][j];
                const auto diff = pos2 - pos1;
                if(diff.two_norm() < dx/1.4)
                {
                  double newdist = dist + std::pow(diff.two_norm(),2);
                  new_crossings.push_back(pos1 + diff/2.);
                  new_crossings_dist.push_back(newdist);

                  Dune::FieldVector<bool,m> newcomp = comp;
                  newcomp[c] = true;
                  new_crossings_components.push_back(newcomp);
                }
              }
            }

            crossings.clear();
            crossings.assign(new_crossings.begin(), new_crossings.end());
            crossings_dist.clear();
            crossings_dist.assign(new_crossings_dist.begin(), new_crossings_dist.end());
            crossings_components.clear();
            crossings_components.assign(new_crossings_components.begin(), new_crossings_components.end());
          }

          for(unsigned int c = 0; c < m; ++c)
          {
            for(unsigned int i = 0; i < crossings.size(); ++i)
            {
              if(!crossings_components[i][c])
              {
                if(crossings[i][c] < dx)
                  crossings[i][c] = 0.;
                else
                  crossings_dist[i] += std::pow(crossings[i][c],2);
              }
            }
          }

          Dune::FieldVector<RangeFieldType, dim> ret(0.);
          RangeFieldType dist = 1./EPSILON;
        
          unsigned int crossing_exists = crossings.size();
          crossing_exists = mpicomm.max(crossing_exists);
          if(crossing_exists == 0)
            return ret;
          
          unsigned int resulting = 0;
          for(unsigned int i = 1; i < crossings.size(); ++i)
            if(crossings_dist[i] < crossings_dist[resulting])
              resulting = i;
          if(crossings.size() != 0)
          {
            ret = crossings[resulting];
            dist = crossings_dist[resulting];
          }
          
          RangeFieldType distMin = mpicomm.min(dist);
          if(!utils::isEqual(dist,distMin,EPSILON))
            for(unsigned int c = 0; c < m; c++)
              ret[c] = 0.;

          auto retget = mpicomm.sum(ret);
          int fail = 0;
          if(utils::isEqual(dist,distMin,EPSILON))
            for(unsigned int c = 0; c < m; ++c)
              if(!utils::isEqual(retget[c],ret[c], EPSILON))
                fail = 1;
          fail = mpicomm.max(fail);
          if(fail == 1)
          {
            //std::cout << "-1: " << ret << " " << retget << "\n";
            return Dune::FieldVector<RangeFieldType, dim>(-1.);
          }

          return retget;
        }
    };
  } // end namespace PDELab
} // end namespace Dune
