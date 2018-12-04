#pragma once

#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/vtk/vtktypes.hh>
#include <dune/vtk/gridcreatorinterface.hh>

namespace Dune
{
  // Create a grid where the input points and connectivity is already
  // connected correctly.
  template <class Grid>
  struct ContinuousGridCreator
      : public GridCreatorInterface<Grid, ContinuousGridCreator<Grid>>
  {
    using Super = GridCreatorInterface<Grid, ContinuousGridCreator<Grid>>;
    using GlobalCoordinate = typename Super::GlobalCoordinate;

    ContinuousGridCreator (GridFactory<Grid>& factory)
      : Super(factory)
    {}

    using Super::factory;

    void insertVerticesImpl (std::vector<GlobalCoordinate> const& points,
                             std::vector<std::uint64_t> const& /*point_ids*/)
    {
      for (auto const& p : points)
        factory().insertVertex(p);
    }

    void insertElementsImpl (std::vector<std::uint8_t> const& types,
                             std::vector<std::int64_t> const& offsets,
                             std::vector<std::int64_t> const& connectivity)
    {
      std::size_t idx = 0;
      for (std::size_t i = 0; i < types.size(); ++i) {
        auto type = Vtk::to_geometry(types[i]);
        Vtk::CellType cellType{type};
        auto refElem = referenceElement<double,Grid::dimension>(type);

        int nNodes = offsets[i] - (i == 0 ? 0 : offsets[i-1]);
        assert(nNodes == refElem.size(Grid::dimension));
        std::vector<unsigned int> vtk_cell; vtk_cell.reserve(nNodes);
        for (int j = 0; j < nNodes; ++j)
          vtk_cell.push_back( connectivity[idx++] );

        if (cellType.noPermutation())
          factory().insertElement(type,vtk_cell);
        else {
          // apply index permutation
          std::vector<unsigned int> cell(nNodes);
          for (int j = 0; j < nNodes; ++j)
            cell[j] = vtk_cell[cellType.permutation(j)];

          factory().insertElement(type,cell);
        }
      }
    }
  };

} // end namespace Dune
