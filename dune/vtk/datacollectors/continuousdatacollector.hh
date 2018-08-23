#pragma once

#include <dune/vtk/datacollector.hh>

namespace Dune { namespace experimental
{

/// Implementation of \ref DataCollector for linear cells, with continuous data.
template <class GridView>
class DefaultDataCollector
    : public DataCollectorInterface<GridView, DefaultDataCollector<GridView>>
{
  enum { dim = GridView::dimension };

  using Self = DefaultDataCollector;
  using Super = DataCollectorInterface<GridView, Self>;
  using Super::gridView_;

public:
  DefaultDataCollector (GridView const& gridView)
    : Super(gridView)
  {}

  /// Return number of grid vertices
  std::uint64_t numPointsImpl () const
  {
    return gridView_.size(dim);
  }

  /// Return the coordinates of all grid vertices in the order given by the indexSet
  template <class T>
  std::vector<T> pointsImpl () const
  {
    std::vector<T> data(this->numPoints() * 3);
    auto const& indexSet = gridView_.indexSet();
    for (auto const& vertex : vertices(gridView_)) {
      std::size_t idx = 3 * indexSet.index(vertex);
      auto v = vertex.geometry().center();
      for (std::size_t j = 0; j < v.size(); ++j)
        data[idx + j] = T(v[j]);
      for (std::size_t j = v.size(); j < 3u; ++j)
        data[idx + j] = T(0);
    }
    return data;
  }

  /// Return the types, offsets and connectivity of the cells, using the same connectivity as
  /// given by the grid.
  Cells cellsImpl () const
  {
    auto const& indexSet = gridView_.indexSet();
    auto types = indexSet.types(0);
    int maxVertices = std::accumulate(types.begin(), types.end(), 1, [](int m, GeometryType t) {
      auto refElem = referenceElement<double,dim>(t);
      return std::max(m, refElem.size(dim));
    });

    Cells cells;
    cells.connectivity.reserve(this->numCells() * maxVertices);
    cells.offsets.reserve(this->numCells());
    cells.types.reserve(this->numCells());

    std::int64_t old_o = 0;
    for (auto const& c : elements(gridView_)) {
      Vtk::CellType cellType(c.type());
      for (int j = 0; j < c.subEntities(dim); ++j)
        cells.connectivity.push_back( std::int64_t(indexSet.subIndex(c,cellType.permutation(j),dim)) );
      cells.offsets.push_back(old_o += c.subEntities(dim));
      cells.types.push_back(cellType.type());
    }

    return cells;
  }

  /// Evaluate the `fct` at the corners of the elements
  template <class T, class GlobalFunction>
  std::vector<T> pointDataImpl (GlobalFunction const& fct) const
  {
    std::vector<T> data(this->numPoints() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_)) {
      localFct.bind(e);
      Vtk::CellType cellType{e.type()};
      auto refElem = referenceElement(e.geometry());
      for (int j = 0; j < e.subEntities(dim); ++j) {
        std::size_t idx = fct.ncomps() * indexSet.subIndex(e,cellType.permutation(j),dim);
        for (int comp = 0; comp < fct.ncomps(); ++comp)
          data[idx + comp] = T(localFct.evaluate(comp, refElem.position(cellType.permutation(j),dim)));
      }
      localFct.unbind();
    }
    return data;
  }
};

}} // end namespace Dune::experimental
