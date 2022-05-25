#ifndef DUNE_ALUGRID_STRUCTUREDGRIDFACTORY_HH
#define DUNE_ALUGRID_STRUCTUREDGRIDFACTORY_HH

#include <array>
#include <memory>
#include <vector>

#include <dune/common/version.hh>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/version.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/exceptions.hh>

#include <dune/alugrid/common/alugrid_assert.hh>
#include <dune/alugrid/common/declaration.hh>

#include <dune/alugrid/common/hsfc.hh>

// include DGF parser implementation for YaspGrid
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  class StructuredGridFactory;



  // StructuredGridFactory for ALUGrid
  // ---------------------------------

  template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refineType, class Comm >
  class StructuredGridFactory< ALUGrid< dim, dimworld, eltype, refineType, Comm > >
  {
  public:
    typedef ALUGrid< dim, dimworld, eltype, refineType, Comm > Grid;
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 6)
    typedef std::unique_ptr< Grid > SharedPtrType;
#else // #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 6)
    // mygrid_ptr in GridPtr is a derived from std::shared_ptr and implements a method release
    typedef typename Dune::GridPtr< Grid > :: mygrid_ptr  SharedPtrType;
#endif // #else // #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 6)

  protected:
    typedef StructuredGridFactory< Grid > This;

  private:
    // SimplePartitioner
    // -----------------
    template< class GV, PartitionIteratorType pitype, class IS = typename GV::IndexSet >
    class SimplePartitioner
    {
      typedef SimplePartitioner< GV, pitype, IS > This;

    public:
      typedef GV GridView;
      typedef typename GridView::Grid Grid;

      typedef IS IndexSet;

    protected:
      typedef typename IndexSet::IndexType IndexType;

      static const int dimension = Grid::dimension;

      typedef typename Grid::template Codim< 0 >::Entity Element;

      typedef typename Element::Geometry::GlobalCoordinate VertexType;

      // type of communicator
      typedef Dune :: CollectiveCommunication< typename MPIHelper :: MPICommunicator >
        CollectiveCommunication ;

#ifdef USE_ALUGRID_SFC_ORDERING
      typedef SpaceFillingCurveOrdering< VertexType >  SpaceFillingCurveOrderingType;
#endif

    public:
      SimplePartitioner ( const GridView &gridView, const CollectiveCommunication& comm,
                          const VertexType& lowerLeft, const VertexType& upperRight )
      : comm_( comm ),
        gridView_( gridView ),
        indexSet_( gridView_.indexSet() ),
        pSize_( comm_.size() ),
        elementCuts_( pSize_, -1 ),
#ifdef USE_ALUGRID_SFC_ORDERING
        sfc_( SpaceFillingCurveOrderingType::ZCurve, lowerLeft, upperRight, comm_ ),
#endif
        maxIndex_( -1.0 )
      {
#ifdef USE_ALUGRID_SFC_ORDERING
        const auto end = gridView_.template end< 0 > ();
        for( auto it = gridView_.template begin< 0 > (); it != end; ++it )
        {
          VertexType center = (*it).geometry().center();
          // get hilbert index in [0,1]
          const double hidx = sfc_.index( center );
          maxIndex_ = std::max( maxIndex_, hidx );
        }

        // adjust with number of elements
        maxIndex_ /= indexSet_.size( 0 );
#endif

        // compute decomposition of sfc
        calculateElementCuts();
      }

    public:
      template< class Entity >
      double index( const Entity &entity ) const
      {
        alugrid_assert ( Entity::codimension == 0 );
#ifdef USE_ALUGRID_SFC_ORDERING
        // get center of entity's geometry
        VertexType center = entity.geometry().center();
        // get hilbert index in [0,1]
        return sfc_.index( center );
#else
        return double(indexSet_.index( entity ));
#endif
      }

      template< class Entity >
      int rank( const Entity &entity ) const
      {
        alugrid_assert ( Entity::codimension == 0 );
#ifdef USE_ALUGRID_SFC_ORDERING
        // get center of entity's geometry
        VertexType center = entity.geometry().center();
        // get hilbert index in [0,1]
        const double hidx = sfc_.index( center );
        // transform to element index
        const long int index = (hidx / maxIndex_);
        //std::cout << "sfc index = " << hidx << " " << index << std::endl;
#else
        const long int index = indexSet_.index( entity );
#endif
        return rank( index );
      }

    protected:
      int rank( long int index ) const
      {
        if( index < elementCuts_[ 0 ] ) return 0;
        for( int p=1; p<pSize_; ++p )
        {
          if( index >= elementCuts_[ p-1 ] && index < elementCuts_[ p ] )
            return p;
        }
        return pSize_-1;
      }

      void calculateElementCuts()
      {
        const size_t nElements = indexSet_.size( 0 );

        // get number of MPI processes
        const int nRanks = pSize_;

        // get minimal number of entities per process
        const size_t minPerProc = (double(nElements) / double( nRanks ));
        size_t maxPerProc = minPerProc ;
        if( nElements % nRanks != 0 )
          ++ maxPerProc ;

        // calculate percentage of elements with larger number
        // of elements per process
        double percentage = (double(nElements) / double( nRanks ));
        percentage -= minPerProc ;
        percentage *= nRanks ;

        int rank = 0;
        size_t elementCount  = maxPerProc ;
        size_t elementNumber = 0;
        size_t localElementNumber = 0;
        const int lastRank = nRanks - 1;

        const size_t size = indexSet_.size( 0 );
        for( size_t i=0; i<size; ++i )
        {
          if( localElementNumber >= elementCount )
          {
            elementCuts_[ rank ] = i ;

            // increase rank
            if( rank < lastRank ) ++ rank;

            // reset local number
            localElementNumber = 0;

            // switch to smaller number if red line is crossed
            if( elementCount == maxPerProc && rank >= percentage )
              elementCount = minPerProc ;
          }

          // increase counters
          ++elementNumber;
          ++localElementNumber;
        }

        // set cut for last process
        elementCuts_[ lastRank ] = size ;

        //for( int p=0; p<pSize_; ++p )
        //  std::cout << "P[ " << p << " ] = " << elementCuts_[ p ] << std::endl;
      }

      const CollectiveCommunication& comm_;

      const GridView& gridView_;
      const IndexSet &indexSet_;

      const int pSize_;
      std::vector< long int > elementCuts_ ;

#ifdef USE_ALUGRID_SFC_ORDERING
      // get element to hilbert (or Z) index mapping
      SpaceFillingCurveOrdering< VertexType > sfc_;
#endif
      double maxIndex_ ;
    };

  public:
    typedef typename Grid::ctype ctype;
    typedef typename MPIHelper :: MPICommunicator MPICommunicatorType ;

    // type of communicator
    typedef Dune :: CollectiveCommunication< MPICommunicatorType >
        CollectiveCommunication ;

    static SharedPtrType
    createCubeGrid( const std::string& filename,
                    MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
    {
      std::ifstream file( filename.c_str() );
      if( ! file )
      {
        DUNE_THROW(InvalidStateException,"file not found " << filename );
      }
      return createCubeGrid( file, filename, mpiComm );
    }

    static SharedPtrType
    createCubeGrid( std::istream& input,
                    const std::string& name,
                    MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
    {
      CollectiveCommunication comm( MPIHelper :: getCommunicator() );
      static_assert( dim == dimworld, "YaspGrid is used for creation of the structured grid which only supports dim == dimworld");

      Dune::dgf::IntervalBlock intervalBlock( input );
      if( !intervalBlock.isactive() )
      {
        std::cerr << "No interval block found, using default DGF method to create ALUGrid!" << std::endl;
        return SharedPtrType( GridPtr< Grid > (input, mpiComm ).release());
      }

      if( intervalBlock.numIntervals() != 1 )
      {
        std::cerr << "ALUGrid creation from YaspGrid can only handle 1 interval block, using default DGF method to create ALUGrid!" << std::endl;
        return SharedPtrType( GridPtr< Grid > (input, mpiComm ).release());
      }

      if( intervalBlock.dimw() != dim )
      {
        std::cerr << "ALUGrid creation from YaspGrid only works for dim == dimworld, using default DGF method to create ALUGrid!" << std::endl;
        return SharedPtrType( GridPtr< Grid > (input, mpiComm ).release());
      }

      const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

      // only work for the new ALUGrid version
      // if creation of SGrid fails the DGF file does not contain a proper
      // IntervalBlock, and thus we cannot create the grid parallel,
      // we will use the standard technique
      std::array<int, dim> dims;
      FieldVector<ctype, dimworld> lowerLeft;
      FieldVector<ctype, dimworld> upperRight;
      for( int i=0; i<dim; ++i )
      {
        dims[ i ]       = interval.n[ i ] ;
        lowerLeft[ i ]  = interval.p[ 0 ][ i ];
        upperRight[ i ] = interval.p[ 1 ][ i ];
      }

      // broadcast array values
      comm.broadcast( &dims[ 0 ], dim, 0 );
      comm.broadcast( &lowerLeft [ 0 ], dim, 0 );
      comm.broadcast( &upperRight[ 0 ], dim, 0 );

      std::string nameYasp( name );
      nameYasp += " via YaspGrid";
      typedef StructuredGridFactory< Grid > SGF;
      return SGF :: createCubeGridImpl( lowerLeft, upperRight, dims, comm, nameYasp );
    }

    template < class int_t >
    static SharedPtrType
    createSimplexGrid ( const FieldVector<ctype,dimworld>& lowerLeft,
                        const FieldVector<ctype,dimworld>& upperRight,
                        const std::array< int_t, dim>& elements,
                        MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
    {
      // create DGF interval block and use DGF parser to create simplex grid
      std::stringstream dgfstream;
      dgfstream << "DGF" << std::endl;
      dgfstream << "Interval" << std::endl;
      dgfstream << lowerLeft  << std::endl;
      dgfstream << upperRight << std::endl;
      for( int i=0; i<dim; ++ i)
        dgfstream << elements[ i ] << " ";
      dgfstream << std::endl;
      dgfstream << "#" << std::endl;
      dgfstream << "Cube" << std::endl;
      dgfstream << "#" << std::endl;
      dgfstream << "Simplex" << std::endl;
      dgfstream << "#" << std::endl;

      std::cout << dgfstream.str() << std::endl;

      Dune::GridPtr< Grid > grid( dgfstream, mpiComm );
      return SharedPtrType( grid.release() );
    }

    static SharedPtrType
    createCubeGrid ( const FieldVector<ctype,dimworld>& lowerLeft,
                     const FieldVector<ctype,dimworld>& upperRight,
                     const std::array< int, dim>& elements,
                     MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
    {
      CollectiveCommunication comm( mpiComm );
      std::string name( "Cartesian ALUGrid via YaspGrid" );
      return createCubeGridImpl( lowerLeft, upperRight, elements, comm, name );
    }

    static SharedPtrType
    createCubeGrid ( const FieldVector<ctype,dimworld>& lowerLeft,
                     const FieldVector<ctype,dimworld>& upperRight,
                     std::array<std::vector<ctype>, dim>& coords,
                     MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
    {
      CollectiveCommunication comm( mpiComm );
      std::string name( "Cartesian ALUGrid via YaspGrid" );
      return createCubeGridImpl( lowerLeft, upperRight, coords, comm, name );
    }
    
  protected:
    template <int codim, class Entity>
    int subEntities ( const Entity& entity ) const
    {
      return entity.subEntities( codim );
    }

    static SharedPtrType
    createCubeGridImpl ( const FieldVector<ctype,dimworld>& lowerLeft,
                         const FieldVector<ctype,dimworld>& upperRight,
                         const std::array< int, dim>& elements,
                         const CollectiveCommunication& comm,
                         const std::string& name )
    {
      const int myrank = comm.rank();

      typedef YaspGrid< dimworld, EquidistantOffsetCoordinates<double,dimworld> > CartesianGridType ;
      std::array< int, dim > dims;
      for( int i=0; i<dim; ++i ) dims[ i ] = elements[ i ];

      CollectiveCommunication commSelf( MPIHelper :: getLocalCommunicator() );
      // create YaspGrid to partition and insert elements that belong to process directly
      CartesianGridType sgrid( lowerLeft, upperRight, dims, std::bitset<dim>(0ULL), 1, commSelf );

      typedef typename CartesianGridType :: LeafGridView GridView ;
      typedef typename GridView  :: IndexSet  IndexSet ;
      typedef typename IndexSet  :: IndexType IndexType ;
      typedef typename GridView  :: template Codim< 0 > :: Iterator ElementIterator ;
      typedef typename ElementIterator::Entity  Entity ;
      typedef typename GridView :: IntersectionIterator     IntersectionIterator ;
      typedef typename IntersectionIterator :: Intersection Intersection ;

      GridView gridView = sgrid.leafGridView();
      const IndexSet &indexSet = gridView.indexSet();

      // get decompostition of the marco grid
      SimplePartitioner< GridView, InteriorBorder_Partition > partitioner( gridView, comm, lowerLeft, upperRight );

      // create ALUGrid GridFactory
      GridFactory< Grid > factory;

      // map global vertex ids to local ones
      std::map< IndexType, unsigned int > vtxMap;
      std::map< double, const Entity > sortedElementList;

      const int numVertices = (1 << dim);
      std::vector< unsigned int > vertices( numVertices );

      const ElementIterator end = gridView.template end< 0 >();
      for( ElementIterator it = gridView.template begin< 0 >(); it != end; ++it )
      {
        const Entity &entity = *it;

        // if the element does not belong to our partition, continue
        if( partitioner.rank( entity ) != myrank )
          continue;

        const double elIndex = partitioner.index( entity );
        assert( sortedElementList.find( elIndex ) == sortedElementList.end() );
        sortedElementList.insert( std::make_pair( elIndex, entity ) );
      }

      int nextElementIndex = 0;
      const auto endSorted = sortedElementList.end();
      for( auto it = sortedElementList.begin(); it != endSorted; ++it )
      {
        const Entity &entity = (*it).second;

        // insert vertices and element
        const typename Entity::Geometry geo = entity.geometry();
        alugrid_assert( numVertices == geo.corners() );
        for( int i = 0; i < numVertices; ++i )
        {
          const IndexType vtxId = indexSet.subIndex( entity, i, dim );
          //auto result = vtxMap.insert( std::make_pair( vtxId, vtxMap.size() ) );
          std::pair< typename std::map< IndexType, unsigned int >::iterator, bool > result
            = vtxMap.insert( std::make_pair( vtxId, vtxMap.size() ) );
          if( result.second )
            factory.insertVertex( geo.corner( i ), vtxId );
          vertices[ i ] = result.first->second;
        }

        factory.insertElement( entity.type(), vertices );
        const int elementIndex = nextElementIndex++;

        //const auto iend = gridView.iend( entity );
        //for( auto iit = gridView.ibegin( entity ); iit != iend; ++iit )
        const IntersectionIterator iend = gridView.iend( entity );
        for( IntersectionIterator iit = gridView.ibegin( entity ); iit != iend; ++iit )
        {
          const Intersection &isec = *iit;
          const int faceNumber = isec.indexInInside();
          // insert boundary face in case of domain boundary
          if( isec.boundary() )
            factory.insertBoundary( elementIndex, faceNumber );
          // insert process boundary if the neighboring element has a different rank
          if( isec.neighbor() && (partitioner.rank( isec.outside() ) != myrank) )
            factory.insertProcessBorder( elementIndex, faceNumber );
        }
      }

      // for structured grids, do not mark longest edge
      // not necessary
      factory.setLongestEdgeFlag(false);

      // create shared grid pointer
      return SharedPtrType( factory.createGrid( true, true, name ) );
    }

    static SharedPtrType
    createCubeGridImpl ( const FieldVector<ctype,dimworld>& lowerLeft,
                         const FieldVector<ctype,dimworld>& upperRight,
                         std::array<std::vector<ctype>, dim>& coords,
                         const CollectiveCommunication& comm,
                         const std::string& name )
    {
      const int myrank = comm.rank();

      typedef YaspGrid< dimworld, TensorProductCoordinates<ctype, dim> > CartesianGridType ;

      CollectiveCommunication commSelf( MPIHelper :: getLocalCommunicator() );
      // create YaspGrid to partition and insert elements that belong to process directly
      CartesianGridType sgrid( coords, std::bitset<dim>(0ULL), 1, commSelf );

      typedef typename CartesianGridType :: LeafGridView GridView ;
      typedef typename GridView  :: IndexSet  IndexSet ;
      typedef typename IndexSet  :: IndexType IndexType ;
      typedef typename GridView  :: template Codim< 0 > :: Iterator ElementIterator ;
      typedef typename ElementIterator::Entity  Entity ;
      typedef typename GridView :: IntersectionIterator     IntersectionIterator ;
      typedef typename IntersectionIterator :: Intersection Intersection ;

      GridView gridView = sgrid.leafGridView();
      const IndexSet &indexSet = gridView.indexSet();

      // get decompostition of the marco grid
      SimplePartitioner< GridView, InteriorBorder_Partition > partitioner( gridView, comm, lowerLeft, upperRight );

      // create ALUGrid GridFactory
      GridFactory< Grid > factory;

      // map global vertex ids to local ones
      std::map< IndexType, unsigned int > vtxMap;
      std::map< double, const Entity > sortedElementList;

      const int numVertices = (1 << dim);
      std::vector< unsigned int > vertices( numVertices );

      const ElementIterator end = gridView.template end< 0 >();
      for( ElementIterator it = gridView.template begin< 0 >(); it != end; ++it )
      {
        const Entity &entity = *it;

        // if the element does not belong to our partition, continue
        if( partitioner.rank( entity ) != myrank )
          continue;

        const double elIndex = partitioner.index( entity );
        assert( sortedElementList.find( elIndex ) == sortedElementList.end() );
        sortedElementList.insert( std::make_pair( elIndex, entity ) );
      }

      int nextElementIndex = 0;
      const auto endSorted = sortedElementList.end();
      for( auto it = sortedElementList.begin(); it != endSorted; ++it )
      {
        const Entity &entity = (*it).second;

        // insert vertices and element
        const typename Entity::Geometry geo = entity.geometry();
        alugrid_assert( numVertices == geo.corners() );
        for( int i = 0; i < numVertices; ++i )
        {
          const IndexType vtxId = indexSet.subIndex( entity, i, dim );
          //auto result = vtxMap.insert( std::make_pair( vtxId, vtxMap.size() ) );
          std::pair< typename std::map< IndexType, unsigned int >::iterator, bool > result
            = vtxMap.insert( std::make_pair( vtxId, vtxMap.size() ) );
          if( result.second )
            factory.insertVertex( geo.corner( i ), vtxId );
          vertices[ i ] = result.first->second;
        }

        factory.insertElement( entity.type(), vertices );
        const int elementIndex = nextElementIndex++;

        //const auto iend = gridView.iend( entity );
        //for( auto iit = gridView.ibegin( entity ); iit != iend; ++iit )
        const IntersectionIterator iend = gridView.iend( entity );
        for( IntersectionIterator iit = gridView.ibegin( entity ); iit != iend; ++iit )
        {
          const Intersection &isec = *iit;
          const int faceNumber = isec.indexInInside();
          // insert boundary face in case of domain boundary
          if( isec.boundary() )
            factory.insertBoundary( elementIndex, faceNumber );
          // insert process boundary if the neighboring element has a different rank
          if( isec.neighbor() && (partitioner.rank( isec.outside() ) != myrank) )
            factory.insertProcessBorder( elementIndex, faceNumber );
        }
      }

      // for structured grids, do not mark longest edge
      // not necessary
      factory.setLongestEdgeFlag(false);

      // create shared grid pointer
      return SharedPtrType( factory.createGrid( true, true, name ) );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_ALUGRID_STRUCTUREDGRIDFACTORY_HH
