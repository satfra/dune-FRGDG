#pragma once

#include<algorithm>

#include <mpi.h>

#include <dune/common/parametertree.hh>
#include <dune/common/parallel/mpicommunication.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
//#include <dune/grid/utility/parmetisgridpartitioner.hh>

#include <dune/pdelab/finiteelementmap/qkdg.hh>

#include <dune/alugrid/grid.hh>

/* 
 * Custom implementation for YaspGrid load balancing, prefers rectangles that are as close as possible to squares.
 */

template<int d>
  class YLoadBalanceCustom : public Dune::YLoadBalance<d>
  {
  public:
    typedef std::array<int, d> iTupel;
    virtual ~YLoadBalanceCustom() {}

    /** \brief Distribute a structured grid across a set of processors
     *
     * \param [in] size Number of elements in each coordinate direction, for the entire grid
     * \param [in] P Number of processors
     */
    virtual void loadbalance (const iTupel& size, int P, iTupel& dims) const
    {
      int opt=0;
      iTupel trydims;

      optimize_dims(d-1,size,P,dims,trydims,opt);
    }
  private:
    void optimize_dims (int i, const iTupel& size, int P, iTupel& dims, iTupel& trydims, int &opt ) const
    {
      if (i>0) // test all subdivisions recursively
      {
        for (int k=1; k<=P; k++)
          if (P%k==0)
          {
            // P divisible by k
            trydims[i] = k;
            optimize_dims(i-1,size,P/k,dims,trydims,opt);
          }
      }
      else
      {
        // found a possible combination
        trydims[0] = P;

        // check for optimality
        int m = 2147483647;

        for (int k=0; k<d; k++)
        {
          int mm=size[k]/trydims[k];
          if ( mm < m ) m = mm;
        }
        if (m>opt)
        {
          opt = m;
          dims = trydims;
        }
      }
    }
  };

/* 
 * Base class for constructing a YaspGrid.
 * Implements getGrid and getMinimalCellEdge, relies on a passed Constructor for getRefinementPerElement.
 */

template<typename RF, unsigned _dim, class Constructor>
class CubicYaspGridConstructorBase
{
  public:
    static constexpr unsigned dim = _dim;
    using Coordinates = Dune::TensorProductCoordinates<RF, dim>;
    using Grid = Dune::YaspGrid<dim,Coordinates>;
    using GV = typename Grid::LeafGridView;
    
    /* 
     * Wrapper to get the load balancer.
     */
    
    static const YLoadBalanceCustom<dim>* customLoadbalancer()
    {
      static YLoadBalanceCustom<dim> lb;
      return & lb;
    }
    
    /* 
     * Function to construct the final grid.
     * Relies on getRefinementPerElement for special refinement.
     * Implements global refinement on the cells from the special refinement.
     */

    static Dune::GridPtr<Grid> getGrid(const Dune::ParameterTree& ptree, const Dune::Communication<MPI_Comm>& mpicomm)
    {
      const Dune::FieldVector<RF, dim> origin = ptree.get<Dune::FieldVector<RF, dim>>("grid.origin");
      const Dune::FieldVector<RF, dim> length = ptree.get<Dune::FieldVector<RF, dim>>("grid.L");
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");

      const unsigned refinement = ptree.get<unsigned>("grid.refinement", 0); // dune refinement, applies to every cell after the other refinements

      // initialize, how often each of the N cells should be halved
      std::array<std::vector<unsigned>,dim> refinementPerElement = Constructor::getRefinementPerElement(ptree);
      
      // options for grid
      std::bitset<dim> periodic(ptree.get("grid.periodic",false));
      int overlap = 1;
      // create cell boundaries vectors, creates N cells, and divides every cell into 2*refinementPerElement[dim][cell] cells if refinementPerElement[dim][cell] is not zero
      std::array<std::vector<RF>, dim> coords;
      for(unsigned d = 0; d < dim; ++d)
      {
        const RF elLen = length[d]/RF(elements[d]);
        coords[d].push_back(origin[d]);
        for(unsigned i = 0; i < elements[d]; ++i)
        {
          for(unsigned j = 1; j < 2*refinementPerElement[d][i]; ++j)
            coords[d].push_back(origin[d] + i*elLen + j*elLen/RF(2*refinementPerElement[d][i]));
          coords[d].push_back(origin[d] + (i+1)*elLen);
        }
      }
      // Create grid from cell boundaries
      auto grid = Dune::GridPtr<Grid>(new Grid(coords, periodic, overlap, mpicomm, customLoadbalancer()));
      grid->globalRefine(refinement);
      return grid;
    }
    
    /* 
     * Function to get the minimal edge along all dimensions, relies on getRefinementPerElement.
     */
    
    static RF getMinimalCellEdge(const Dune::ParameterTree& ptree)
    {
      const Dune::FieldVector<RF, dim> origin = ptree.get<Dune::FieldVector<RF, dim>>("grid.origin");
      const Dune::FieldVector<RF, dim> length = ptree.get<Dune::FieldVector<RF, dim>>("grid.L");
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");
      const unsigned refinement = ptree.get<unsigned>("grid.refinement", 0);
      
      std::array<std::vector<unsigned>,dim> refinementPerElement = Constructor::getRefinementPerElement(ptree);
      
      std::vector<RF> dxs;
      RF dx;
      unsigned maxRefinement;
      for(unsigned d = 0; d < dim; ++d)
      {
        maxRefinement = *std::max_element(refinementPerElement[d].begin(), refinementPerElement[d].end());
        if (maxRefinement == 0) maxRefinement = 1;
        else                    maxRefinement = 2*maxRefinement;
        
        dx = (length[d]-origin[d])/elements[d]/maxRefinement/std::pow(2.,RF(refinement));
        dxs.push_back(dx);
      }
      
      return *std::min_element(dxs.begin(), dxs.end());
    }
};

/* 
 * Implements local refinement, inherits from CubicYaspGridConstructorBase.
 */

template<typename RF, unsigned _dim>
class CubicYaspGridConstructor : public CubicYaspGridConstructorBase<RF, _dim, CubicYaspGridConstructor<RF, _dim>>
{
  public:
    static constexpr unsigned dim = _dim;
    using Coordinates = Dune::TensorProductCoordinates<RF, dim>;
    using Grid = Dune::YaspGrid<dim,Coordinates>;
    using GV = typename Grid::LeafGridView;
    
    /* 
     * Refinement starting from the origin.
     */
    
    static std::array<std::vector<unsigned>,dim> getRefinementPerElement(const Dune::ParameterTree& ptree)
    {
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");
      const unsigned localRefinement = ptree.get<unsigned>("grid.localRefinement", 0);
      
      std::array<std::vector<unsigned>,dim> refinementPerElement;
      for(unsigned d = 0; d < dim; ++d)
        refinementPerElement[d].resize(elements[d]);
      for(unsigned r = 0; r < localRefinement; ++r)
      {
        const std::array<unsigned,dim> ref = ptree.get<std::array<unsigned,dim>>("grid.refine" + std::to_string(r));
        for(unsigned d = 0; d < dim; ++d)
          for(unsigned i = 0; i < elements[d]; ++i)
            refinementPerElement[d][i] += unsigned(i<ref[d]);
      }
      
      return refinementPerElement;
    }
};

/* 
 * Implements different version of local refinement where single cells can be refined, inherits from CubicYaspGridConstructorBase.
 */

template<typename RF, unsigned _dim>
class CubicYaspGridConstructorCellRefine : public CubicYaspGridConstructorBase<RF, _dim, CubicYaspGridConstructorCellRefine<RF, _dim>>
{
  public:
    static constexpr unsigned dim = _dim;
    using Coordinates = Dune::TensorProductCoordinates<RF, dim>;
    using Grid = Dune::YaspGrid<dim,Coordinates>;
    using GV = typename Grid::LeafGridView;
    
    /* 
     * Refinement starting from the origin (with a different form of the .ini parameters), and for single cells.
     */
    
    static std::array<std::vector<unsigned>,dim> getRefinementPerElement(const Dune::ParameterTree& ptree)
    {
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");
      // get local refinement 
      std::array<std::vector<unsigned>,dim> localRefinement;
      for(unsigned d = 0; d < dim; ++d)
        localRefinement[d] = ptree.get<std::vector<unsigned>>("grid.localRefinement" + std::to_string(d), {});
      // get local cell refinement 
      std::array<std::vector<unsigned>,dim> localCellRefinement;
      for(unsigned d = 0; d < dim; ++d)
        localCellRefinement[d] = ptree.get<std::vector<unsigned>>("grid.localCellRefinement" + std::to_string(d), {});
      
      // initialize, how often each of the N cells should be halved
      std::array<std::vector<unsigned>,dim> refinementPerElement;
      for(unsigned d = 0; d < dim; ++d)
        refinementPerElement[d].resize(elements[d]);
      // apply localRefinement
      for(unsigned d = 0; d < dim; ++d)
        for(unsigned l = 0; l < localRefinement[d].size(); l++)
          for(unsigned i = 0; i < elements[d]; ++i)
            refinementPerElement[d][i] += unsigned(i<localRefinement[d][l]);
      // apply localCellRefinement
      for(unsigned d = 0; d < dim; ++d)
        for(unsigned l = 0; l < localCellRefinement[d].size(); l++)
          for(unsigned i = 0; i < elements[d]; ++i)
            refinementPerElement[d][i] += unsigned(i==localCellRefinement[d][l]);
      
      return refinementPerElement;
    }
};


/* 
 * Base class for constructing an ALUGrid.
 * Implements getGrid and getMinimalCellEdge, relies on a passed Constructor for getRefinementPerElement.
 */

#define USE_ALUGRID_SFC_ORDERING 1

template<typename RF, unsigned _dim, class Constructor>
class CubicALUGridConstructorBase
{
  public:
    static constexpr unsigned dim = _dim;
    using Grid = Dune::ALUGrid< dim, dim, Dune::cube, Dune::nonconforming >;
    using GridFactory = Dune::StructuredGridFactory<Grid>;
    using GV = typename Grid::LeafGridView;
    
    /* 
     * Function to construct the final grid.
     * Relies on getRefinementPerElement for special refinement.
     * Implements global refinement on the cells from the special refinement.
     */

    static Dune::GridPtr<Grid> getGrid(const Dune::ParameterTree& ptree, const Dune::Communication<MPI_Comm>& mpicomm)
    {
      const Dune::FieldVector<RF, dim> origin = ptree.get<Dune::FieldVector<RF, dim>>("grid.origin");
      const Dune::FieldVector<RF, dim> length = ptree.get<Dune::FieldVector<RF, dim>>("grid.L");
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");

      const unsigned refinement = ptree.get<unsigned>("grid.refinement", 0); // dune refinement, applies to every cell after the other refinements

      // initialize, how often each of the N cells should be halved
      std::array<std::vector<unsigned>,dim> refinementPerElement = Constructor::getRefinementPerElement(ptree);
      
      Dune::FieldVector<RF, dim> upperRight;
      for(unsigned d = 0; d < dim; ++d)
        upperRight[d] = origin[d] + length[d];
      
      // create cell boundaries vectors, creates N cells, and divides every cell into 2*refinementPerElement[dim][cell] cells if refinementPerElement[dim][cell] is not zero
      std::array<std::vector<RF>, dim> coords;
      for(unsigned d = 0; d < dim; ++d)
      {
        const RF elLen = length[d]/RF(elements[d]);
        coords[d].push_back(origin[d]);
        for(unsigned i = 0; i < elements[d]; ++i)
        {
          for(unsigned j = 1; j < 2*refinementPerElement[d][i]; ++j)
            coords[d].push_back(origin[d] + i*elLen + j*elLen/RF(2*refinementPerElement[d][i]));
          coords[d].push_back(origin[d] + (i+1)*elLen);
        }
      }
      // Create grid from cell boundaries
      auto grid = Dune::GridPtr< Grid >(GridFactory::createCubeGrid(origin, upperRight, coords, mpicomm).release());
      grid->globalRefine(refinement);
      grid->loadBalance();
      return grid;
    }
    
    /* 
     * Function to get the minimal edge along all dimensions, relies on getRefinementPerElement.
     */
    
    static RF getMinimalCellEdge(const Dune::ParameterTree& ptree)
    {
      const Dune::FieldVector<RF, dim> origin = ptree.get<Dune::FieldVector<RF, dim>>("grid.origin");
      const Dune::FieldVector<RF, dim> length = ptree.get<Dune::FieldVector<RF, dim>>("grid.L");
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");
      const unsigned refinement = ptree.get<unsigned>("grid.refinement", 0);
      
      std::array<std::vector<unsigned>,dim> refinementPerElement = Constructor::getRefinementPerElement(ptree);
      
      std::vector<RF> dxs;
      RF dx;
      unsigned maxRefinement;
      for(unsigned d = 0; d < dim; ++d)
      {
        maxRefinement = *std::max_element(refinementPerElement[d].begin(), refinementPerElement[d].end());
        if (maxRefinement == 0) maxRefinement = 1;
        else                    maxRefinement = 2*maxRefinement;
        
        dx = (length[d]-origin[d])/elements[d]/maxRefinement/std::pow(2.,RF(refinement));
        dxs.push_back(dx);
      }
      
      return *std::min_element(dxs.begin(), dxs.end());
    }
};

/* 
 * Implements local refinement, inherits from CubicALUGridConstructorBase.
 */

template<typename RF, unsigned _dim>
class CubicALUGridConstructor : public CubicALUGridConstructorBase<RF, _dim, CubicALUGridConstructor<RF, _dim>>
{
  public:
    static constexpr unsigned dim = _dim;
    using Grid = Dune::ALUGrid< dim, dim, Dune::cube, Dune::nonconforming >;
    using GV = typename Grid::LeafGridView;
    
    /* 
     * Refinement starting from the origin.
     */
    
    static std::array<std::vector<unsigned>,dim> getRefinementPerElement(const Dune::ParameterTree& ptree)
    {
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");
      const unsigned localRefinement = ptree.get<unsigned>("grid.localRefinement", 0);
      
      std::array<std::vector<unsigned>,dim> refinementPerElement;
      for(unsigned d = 0; d < dim; ++d)
        refinementPerElement[d].resize(elements[d]);
      for(unsigned r = 0; r < localRefinement; ++r)
      {
        const std::array<unsigned,dim> ref = ptree.get<std::array<unsigned,dim>>("grid.refine" + std::to_string(r));
        for(unsigned d = 0; d < dim; ++d)
          for(unsigned i = 0; i < elements[d]; ++i)
            refinementPerElement[d][i] += unsigned(i<ref[d]);
      }
      
      return refinementPerElement;
    }
};

/* 
 * Implements different version of local refinement where single cells can be refined, inherits from CubicALUGridConstructorBase.
 */

template<typename RF, unsigned _dim>
class CubicALUGridConstructorCellRefine : public CubicALUGridConstructorBase<RF, _dim, CubicALUGridConstructorCellRefine<RF, _dim>>
{
  public:
    static constexpr unsigned dim = _dim;
    using Grid = Dune::ALUGrid< dim, dim, Dune::cube, Dune::nonconforming >;
    using GV = typename Grid::LeafGridView;
    
    /* 
     * Refinement starting from the origin (with a different form of the .ini parameters), and for single cells.
     */
    
    static std::array<std::vector<unsigned>,dim> getRefinementPerElement(const Dune::ParameterTree& ptree)
    {
      const std::array<unsigned, dim> elements = ptree.get<std::array<unsigned, dim>>("grid.N");
      // get local refinement 
      std::array<std::vector<unsigned>,dim> localRefinement;
      for(unsigned d = 0; d < dim; ++d)
        localRefinement[d] = ptree.get<std::vector<unsigned>>("grid.localRefinement" + std::to_string(d), {});
      // get local cell refinement 
      std::array<std::vector<unsigned>,dim> localCellRefinement;
      for(unsigned d = 0; d < dim; ++d)
        localCellRefinement[d] = ptree.get<std::vector<unsigned>>("grid.localCellRefinement" + std::to_string(d), {});
      
      // initialize, how often each of the N cells should be halved
      std::array<std::vector<unsigned>,dim> refinementPerElement;
      for(unsigned d = 0; d < dim; ++d)
        refinementPerElement[d].resize(elements[d]);
      // apply localRefinement
      for(unsigned d = 0; d < dim; ++d)
        for(unsigned l = 0; l < localRefinement[d].size(); l++)
          for(unsigned i = 0; i < elements[d]; ++i)
            refinementPerElement[d][i] += unsigned(i<localRefinement[d][l]);
      // apply localCellRefinement
      for(unsigned d = 0; d < dim; ++d)
        for(unsigned l = 0; l < localCellRefinement[d].size(); l++)
          for(unsigned i = 0; i < elements[d]; ++i)
            refinementPerElement[d][i] += unsigned(i==localCellRefinement[d][l]);
      
      return refinementPerElement;
    }
};