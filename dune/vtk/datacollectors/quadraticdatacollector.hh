#pragma once

#include <dune/vtk/datacollector.hh>

namespace Dune { namespace experimental
{

/// Implementation of \ref DataCollector for quadratic cells, with continuous data.
template <class GridView>
class QuadraticDataCollector
    : public DataCollectorInterface<GridView, QuadraticDataCollector<GridView>>
{
  enum { dim = GridView::dimension };

  using Self = QuadraticDataCollector;
  using Super = DataCollectorInterface<GridView, Self>;
  using Super::gridView_;

public:
  QuadraticDataCollector (GridView const& gridView)
    : Super(gridView)
  {}

  /// Return number of vertices + number of edge
  std::uint64_t numPointsImpl () const
  {
    return gridView_.size(dim) + gridView_.size(dim-1);
  }

  /// Return a vector of point coordinates.
  /**
   * The vector of point coordinates is composed of vertex coordinates first and second
   * edge center coordinates.
   **/
  template <class T>
  std::vector<T> pointsImpl () const
  {
    std::vector<T> data(this->numPoints() * 3);
    auto const& indexSet = gridView_.indexSet();
    for (auto const& element : elements(gridView_, Partitions::interior)) {
      auto geometry = element.geometry();
      auto refElem = referenceElement<T,dim>(element.type());

      // vertices
      for (int i = 0; i < element.subEntities(dim); ++i) {
        std::size_t idx = 3 * indexSet.subIndex(element, i, dim);
        auto v = geometry.global(refElem.position(i,dim));
        for (std::size_t j = 0; j < v.size(); ++j)
          data[idx + j] = T(v[j]);
        for (std::size_t j = v.size(); j < 3u; ++j)
          data[idx + j] = T(0);
      }
      // edge centers
      for (int i = 0; i < element.subEntities(dim-1); ++i) {
        std::size_t idx = 3 * (indexSet.subIndex(element, i, dim-1) + gridView_.size(dim));
        auto v = geometry.global(refElem.position(i,dim-1));
        for (std::size_t j = 0; j < v.size(); ++j)
          data[idx + j] = T(v[j]);
        for (std::size_t j = v.size(); j < 3u; ++j)
          data[idx + j] = T(0);
      }
    }
    return data;
  }

  /// \brief Return cell types, offsets, and connectivity. \see Cells
  /**
   * The cell connectivity is composed of cell vertices first and second cell edges,
   * where the indices are grouped [vertex-indices..., (#vertices)+edge-indices...]
   **/
  Cells cellsImpl () const
  {
    Cells cells;
    cells.connectivity.reserve(this->numPoints());
    cells.offsets.reserve(this->numCells());
    cells.types.reserve(this->numCells());

    std::int64_t old_o = 0;
    auto const& indexSet = gridView_.indexSet();
    for (auto const& c : elements(gridView_, Partitions::interior)) {
      Vtk::CellType cellType(c.type(), Vtk::QUADRATIC);
      for (int j = 0; j < c.subEntities(dim); ++j) {
        int k = cellType.permutation(j);
        std::int64_t point_idx = indexSet.subIndex(c,k,dim);
        cells.connectivity.push_back(point_idx);
      }
      for (int j = 0; j < c.subEntities(dim-1); ++j) {
        int k = cellType.permutation(c.subEntities(dim) + j);
        std::int64_t point_idx = (indexSet.subIndex(c,k,dim-1) + gridView_.size(dim));
        cells.connectivity.push_back(point_idx);
      }
      cells.offsets.push_back(old_o += c.subEntities(dim)+c.subEntities(dim-1));
      cells.types.push_back(cellType.type());
    }

    return cells;
  }

  /// Evaluate the `fct` at element vertices and edge centers in the same order as the point coords.
  template <class T, class GlobalFunction>
  std::vector<T> pointDataImpl (GlobalFunction const& fct) const
  {
    std::vector<T> data(this->numPoints() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_, Partitions::interior)) {
      localFct.bind(e);
      Vtk::CellType cellType{e.type(), Vtk::QUADRATIC};
      auto refElem = referenceElement(e.geometry());
      for (int j = 0; j < e.subEntities(dim); ++j) {
        int k = cellType.permutation(j);
        std::size_t idx = fct.ncomps() * indexSet.subIndex(e, k, dim);
        for (int comp = 0; comp < fct.ncomps(); ++comp)
          data[idx + comp] = T(localFct.evaluate(comp, refElem.position(k, dim)));
      }
      for (int j = 0; j < e.subEntities(dim-1); ++j) {
        int k = cellType.permutation(e.subEntities(dim) + j);
        std::size_t idx = fct.ncomps() * (indexSet.subIndex(e, k, dim-1) + gridView_.size(dim));
        for (int comp = 0; comp < fct.ncomps(); ++comp)
          data[idx + comp] = T(localFct.evaluate(comp, refElem.position(k, dim-1)));
      }
      localFct.unbind();
    }
    return data;
  }
};

}} // end namespace Dune::extensions
