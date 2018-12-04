#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <dune/common/version.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/vtk/forward.hh>

namespace Dune {

template <class G, class Derived>
class GridCreatorInterface
{
public:
  using Grid = G;
  using GlobalCoordinate = typename Grid::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
  GridCreatorInterface (GridFactory<Grid>& factory)
    : factory_(&factory)
  {
#if DUNE_VERSION_LT(DUNE_GRID,2,7)
  // old GridFactory implementation does not provide access to its collective communication
  #if HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks_);
  #endif
#else
    rank_ = factory.comm().rank();
    numRanks_ = factory.comm().size();
#endif
  }

  /// Insert all points as vertices into the factory
  void insertVertices (std::vector<GlobalCoordinate> const& points,
                       std::vector<std::uint64_t> const& point_ids)
  {
    asDerived().insertVerticesImpl(points, point_ids);
  }

  /// Create elements based on type and connectivity description
  void insertElements (std::vector<std::uint8_t> const& types,
                       std::vector<std::int64_t> const& offsets,
                       std::vector<std::int64_t> const& connectivity)
  {
    asDerived().insertElementsImpl(types, offsets, connectivity);
  }

  /// Insert part of a grid stored in file into factory
  void insertPieces (std::vector<std::string> const& pieces)
  {
    asDerived().insertPiecesImpl(pieces);
  }

  /// Return the associated GridFactory
  GridFactory<Grid>& factory ()
  {
    return *factory_;
  }

  /// Return the MPI_Comm_rank of the factory, (or of MPI_COMM_WORLD)
  int rank () const
  {
    return rank_;
  }

  /// Return the MPI_Comm_size of the factory, (or of MPI_COMM_WORLD)
  int size () const
  {
    return numRanks_;
  }

protected: // cast to derived type

  Derived& asDerived ()
  {
    return static_cast<Derived&>(*this);
  }

  const Derived& asDerived () const
  {
    return static_cast<const Derived&>(*this);
  }

public: // default implementations

  void insertVerticesImpl (std::vector<GlobalCoordinate> const&,
                           std::vector<std::uint64_t> const&)
  {
    /* do nothing */
  }

  void insertElementsImpl (std::vector<std::uint8_t> const&,
                           std::vector<std::int64_t> const&,
                           std::vector<std::int64_t> const&)
  {
    /* do nothing */
  }

  void insertPiecesImpl (std::vector<std::string> const&)
  {
    /* do nothing */;
  }

protected:
  GridFactory<Grid>* factory_;
  int rank_ = 0;
  int numRanks_ = 1;
};

} // end namespace Dune
