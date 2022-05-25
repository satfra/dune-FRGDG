#pragma once

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/FRGDG/grid/gridDataAccessor.hh>

namespace Dune {
  namespace PDELab {
    template <typename GFS>
      class RejectStep : GridDataAccessor<GFS>
      {
        private:
          using GDA = GridDataAccessor<GFS>;
          using GDA::gfs;
          using GDA::lfs;
          using GV = typename GFS::Traits::GridViewType;

          static constexpr unsigned int dim = GV::dimension;
          static constexpr unsigned int m = GFS::CHILDREN;

          using RF = typename GDA::RF;
          using X = typename GDA::X;

          int verbosity;

        public:
          RejectStep(const GFS& gfs_, int _verbosity = 0) : GDA(gfs_), verbosity(_verbosity) {}

          bool CheckStep(X& x)
          {
            bool rej = false;
            // for each cell
            auto it = gfs.gridView().template begin<0,Dune::InteriorBorder_Partition>();
            auto endit = gfs.gridView().template end<0,Dune::InteriorBorder_Partition>();
            for (; it != endit; ++it)
            {
              const auto xl = this->viewData(x, *it);

              using namespace Indices;
              const auto& dgspace = child(lfs,_0);

              // we simply check all coefficients for finiteness - if something went wrong, it will be nan there
              for (unsigned int j=0; j<m; j++) // for all components
                for (unsigned int i=0; i<dgspace.size(); i++) // for all basis functions
                  rej = rej || std::isnan(xl[i + dgspace.size()*j]) || std::isinf(xl[i + dgspace.size()*j]);

              if(rej)
              {
                auto geo = it->geometry();
                if(verbosity != 0)
                  std::cout << "failure at " << geo.center() << "\n";
                return rej;
              }
            }   
            return rej;
          }
      };
  } // end namespace PDELab
} // end namespace Dune
