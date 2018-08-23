#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "vtktypes.hh"

namespace Dune { namespace experimental {

struct Cells
{
  std::vector<std::uint8_t> types;
  std::vector<std::int64_t> offsets;
  std::vector<std::int64_t> connectivity;
};

template <class GridView, class Derived>
class DataCollectorInterface
{
public:

  DataCollectorInterface (GridView const& gridView)
    : gridView_(gridView)
  {}

  /// \brief Return the number of points in the grid
  std::uint64_t numPoints () const
  {
    return asDerived().numPointsImpl();
  }

  /// \brief Return the number of cells in the grid
  std::uint64_t numCells () const
  {
    return asDerived().numCellsImpl();
  }

  /// Update the DataCollector on the current GridView
  void update ()
  {
    asDerived().updateImpl();
  }

  /// \brief Return a flat vector of point coordinates
  /**
   * All coordinates are extended to 3 components and concatenated.
   * [p0_x, p0_y, p0_z, p1_x, p1_y, p1_z, ...]
   * If the GridView::dimensionworld < 3, the remaining components are
   * set to 0
   **/
  template <class T>
  std::vector<T> points () const
  {
    return asDerived().template pointsImpl<T>();
  }

  /// \brief Return cell types, offsets, and connectivity. \see Cells
  Cells cells () const
  {
    return asDerived().cellsImpl();
  }

  /// \brief Return a flat vector of function values evaluated at the grid points.
  /**
   * In case of a vector valued function, flat the vector entries:
   * [fct(p0)_0, fct(p0)_1, fct(p0)_2, fct(p1)_0, ...]
   * where the vector dimension must be 3 (possible extended by 0s)
   * In case of tensor valued function, flat the tensor row-wise:
   * [fct(p0)_00, fct(p0)_01, fct(p0)_02, fct(p0)_10, fct(p0)_11, fct(p0)_12, fct(p0)_20...]
   * where the tensor dimension must be 3x3 (possible extended by 0s)
   **/
  template <class T, class GlobalFunction>
  std::vector<T> pointData (GlobalFunction const& fct) const
  {
    return asDerived().template pointDataImpl<T>(fct);
  }

  /// \brief Return a flat vector of function values evaluated at the grid cells. \see pointData.
  template <class T, class GlobalFunction>
  std::vector<T> cellData (GlobalFunction const& fct) const
  {
    return asDerived().template cellDataImpl<T>(fct);
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


protected: // default implementations

  std::uint64_t numCellsImpl () const
  {
    return gridView_.size(0);
  }

  void updateImpl () { /* do nothing */ }

  // Evaluate `fct` in center of cell
  template <class T, class GlobalFunction>
  std::vector<T> cellDataImpl (GlobalFunction const& fct) const
  {
    std::vector<T> data(numCells() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_)) {
      localFct.bind(e);
      auto geometry = e.geometry();
      std::size_t idx = fct.ncomps() * indexSet.index(e);
      for (int comp = 0; comp < fct.ncomps(); ++comp)
        data[idx + comp] = T(localFct.evaluate(comp, geometry.center()));
      localFct.unbind();
    }
    return data;
  }


protected:

  GridView gridView_;
};


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
    for (auto const& c : elements(gridView_)) {
      numPoints_ += c.subEntities(dim);
      for (int i = 0; i < c.subEntities(dim); ++i)
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
    for (auto const& element : elements(gridView_)) {
      for (int i = 0; i < element.subEntities(dim); ++i) {
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
    for (auto const& c : elements(gridView_)) {
      Vtk::CellType cellType(c.type());
      for (int j = 0; j < c.subEntities(dim); ++j) {
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
    for (auto const& e : elements(gridView_)) {
      localFct.bind(e);
      Vtk::CellType cellType{e.type()};
      auto refElem = referenceElement(e.geometry());
      for (int j = 0; j < e.subEntities(dim); ++j) {
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
    for (auto const& element : elements(gridView_)) {
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
    for (auto const& c : elements(gridView_)) {
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
    for (auto const& e : elements(gridView_)) {
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

}} // end namespace Dune::experimental
