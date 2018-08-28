#pragma once

#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>

#include "vtktypes.hh"

namespace Dune
{
  // Create a grid where the input points and connectivity is already
  // connected correctly.
  struct DefaultGridCreator
  {
    template <class Grid, class Coord>
    static void create (GridFactory<Grid>& factory,
                        std::vector<Coord> const& points,
                        std::vector<std::uint8_t> const& types,
                        std::vector<std::int64_t> const& offsets,
                        std::vector<std::int64_t> const& connectivity)
    {
      for (auto const& p : points)
        factory.insertVertex(p);

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
          factory.insertElement(type,vtk_cell);
        else {
          // apply index permutation
          std::vector<unsigned int> cell(nNodes);
          for (int j = 0; j < nNodes; ++j)
            cell[j] = vtk_cell[cellType.permutation(j)];

          factory.insertElement(type,cell);
        }
      }
    }
  };


  // Create a grid where the input points are not connected and the connectivity
  // describes separated elements.
  struct ConnectedGridCreator
  {
    struct CoordLess
    {
      template <class T, int N>
      bool operator() (FieldVector<T,N> const& lhs, FieldVector<T,N> const& rhs) const
      {
        for (int i = 0; i < N; ++i) {
          if (std::abs(lhs[i] - rhs[i]) < std::numeric_limits<T>::epsilon())
            continue;
          return lhs[i] < rhs[i];
        }
        return false;
      }
    };

    template <class Grid, class Coord>
    static void create (GridFactory<Grid>& factory,
                        std::vector<Coord> const& points,
                        std::vector<std::uint8_t> const& types,
                        std::vector<std::int64_t> const& offsets,
                        std::vector<std::int64_t> const& connectivity)
    {
      std::size_t idx = 0;
      std::map<Coord, std::size_t, CoordLess> unique_points;
      for (auto const& p : points) {
        auto b = unique_points.emplace(std::make_pair(p,idx));
        if (b.second) {
          factory.insertVertex(p);
          ++idx;
        }
      }

      idx = 0;
      for (std::size_t i = 0; i < types.size(); ++i) {
        auto type = Vtk::to_geometry(types[i]);
        Vtk::CellType cellType{type};

        int nNodes = offsets[i] - (i == 0 ? 0 : offsets[i-1]);
        assert(nNodes > 0);
        std::vector<unsigned int> vtk_cell; vtk_cell.reserve(nNodes);
        for (int j = 0; j < nNodes; ++j) {
          std::size_t v_j = connectivity[idx++];
          std::size_t new_idx = unique_points[points[v_j]];
          vtk_cell.push_back(new_idx);
        }

        if (cellType.noPermutation())
          factory.insertElement(type,vtk_cell);
        else {
          // apply index permutation
          std::vector<unsigned int> cell(nNodes);
          for (int j = 0; j < nNodes; ++j)
            cell[j] = vtk_cell[cellType.permutation(j)];

          factory.insertElement(type,cell);
        }
      }
    }
  };

} // end namespace Dune
