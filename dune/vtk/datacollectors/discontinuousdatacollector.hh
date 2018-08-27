#pragma once

#include <dune/vtk/datacollector.hh>

namespace Dune { namespace experimental
{

/// Implementation of \ref DataCollector for linear cells, with discontinuous data.
template <class GridView>
class DiscontinuousDataCollector
    : public DataCollectorInterface<GridView, DiscontinuousDataCollector<GridView>>
{
  enum { dim = GridView::dimension };

  using Self = DiscontinuousDataCollector;
  using Super = DataCollectorInterface<GridView, Self>;
  using Super::gridView_;

public:
  DiscontinuousDataCollector (GridView const& gridView)
    : Super(gridView)
  {
    this->update();
  }

  /// Create an index map the uniquely assignes an index to each pair (element,corner)
  void updateImpl ()
  {
    numPoints_ = 0;
    indexMap_.resize(gridView_.size(dim));
    std::int64_t vertex_idx = 0;
    auto const& indexSet = gridView_.indexSet();
    for (auto const& c : elements(gridView_, Partitions::interior)) {
      numPoints_ += c.subEntities(dim);
      for (unsigned int i = 0; i < c.subEntities(dim); ++i)
        indexMap_[indexSet.subIndex(c, i, dim)] = vertex_idx++;
    }
  }

  /// The number of pointsi approx. #cell * #corners-per-cell
  std::uint64_t numPointsImpl () const
  {
    return numPoints_;
  }

  /// Return the coordinates of the corners of all cells
  template <class T>
  std::vector<T> pointsImpl () const
  {
    std::vector<T> data(this->numPoints() * 3);
    auto const& indexSet = gridView_.indexSet();
    for (auto const& element : elements(gridView_, Partitions::interior)) {
      for (unsigned int i = 0; i < element.subEntities(dim); ++i) {
        std::size_t idx = 3 * indexMap_[indexSet.subIndex(element, i, dim)];
        auto v = element.geometry().corner(i);
        for (std::size_t j = 0; j < v.size(); ++j)
          data[idx + j] = T(v[j]);
        for (std::size_t j = v.size(); j < 3u; ++j)
          data[idx + j] = T(0);
      }
    }
    return data;
  }

  /// Connect the corners of each cell. The leads to a global discontinuous grid
  Cells cellsImpl () const
  {
    Cells cells;
    cells.connectivity.reserve(this->numPoints());
    cells.offsets.reserve(this->numCells());
    cells.types.reserve(this->numCells());

    std::int64_t old_o = 0;
    auto const& indexSet = gridView_.indexSet();
    for (auto const& c : elements(gridView_, Partitions::interior)) {
      Vtk::CellType cellType(c.type());
      for (unsigned int j = 0; j < c.subEntities(dim); ++j) {
        std::int64_t vertex_idx = indexMap_[indexSet.subIndex(c,cellType.permutation(j),dim)];
        cells.connectivity.push_back(vertex_idx);
      }
      cells.offsets.push_back(old_o += c.subEntities(dim));
      cells.types.push_back(cellType.type());
    }

    return cells;
  }

  /// Evaluate the `fct` in the corners of each cell
  template <class T, class GlobalFunction>
  std::vector<T> pointDataImpl (GlobalFunction const& fct) const
  {
    std::vector<T> data(this->numPoints() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_, Partitions::interior)) {
      localFct.bind(e);
      Vtk::CellType cellType{e.type()};
      auto refElem = referenceElement(e.geometry());
      for (unsigned int j = 0; j < e.subEntities(dim); ++j) {
        std::size_t idx = fct.ncomps() * indexMap_[indexSet.subIndex(e, cellType.permutation(j), dim)];
        for (int comp = 0; comp < fct.ncomps(); ++comp)
          data[idx + comp] = T(localFct.evaluate(comp, refElem.position(cellType.permutation(j),dim)));
      }
      localFct.unbind();
    }
    return data;
  }

private:
  std::uint64_t numPoints_ = 0;
  std::vector<std::int64_t> indexMap_;
};

}} // end namespace Dune::experimental
