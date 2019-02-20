#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <dune/common/version.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/vtk/forward.hh>

namespace Dune
{

  /// Base class for grid creators in a CRTP style.
  /**
   * Construct a grid from data read from VTK files.
   *
   * \tparam GridView   Model of Dune::GridView
   * \tparam Derived    Implementation of a concrete GridCreator.
   **/
  template <class G, class Derived>
  class GridCreatorInterface
  {
  public:
    using Grid = G;
    using GlobalCoordinate = typename Grid::template Codim<0>::Entity::Geometry::GlobalCoordinate;

  public:
    /// Constructor. Stores a reference to the passed GridFactory
    GridCreatorInterface (GridFactory<Grid>& factory)
      : factory_(&factory)
    {}

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

    /// Return the mpi collective communicator
    auto comm () const
    {
      return MPIHelper::getCollectiveCommunication();
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
  };

} // end namespace Dune
