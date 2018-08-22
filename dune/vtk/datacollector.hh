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

class DataCollectorInterface
{
private:
  DataCollectorInterface () = default;

public:
  /// \brief Return the number of points in the grid
  std::uint64_t numPoints () const;

  /// \brief Return the number of cells in the grid
  std::uint64_t numCells () const;

  /// Update the DataCollector on the current GridView
  void update ();

  /// \brief Return a flat vector of point coordinates
  /**
   * All coordinates are extended to 3 components and concatenated.
   * [p0_x, p0_y, p0_z, p1_x, p1_y, p1_z, ...]
   * If the GridView::dimensionworld < 3, the remaining components are
   * set to 0
   **/
  template <class T>
  std::vector<T> points () const;

  /// \brief Return cell types, offsets, and connectivity. \see Cells
  Cells cells () const;

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
  std::vector<T> pointData (GlobalFunction const& fct) const;

  /// \brief Return a flat vector of function values evaluated at the grid cells. \see pointData.
  template <class T, class GlobalFunction>
  std::vector<T> cellData (GlobalFunction const& fct) const;
};


/// Implementation of \ref DataCollector for linear cells, with continuous data.
template <class GridView>
class DefaultDataCollector
{
  enum { dim = GridView::dimension };

public:
  DefaultDataCollector (GridView const& gridView)
    : gridView_(gridView)
  {}

  void update () {}

  std::uint64_t numPoints () const
  {
    return gridView_.size(dim);
  }

  std::uint64_t numCells () const
  {
    return gridView_.size(0);
  }

  template <class T>
  std::vector<T> points () const
  {
    std::vector<T> data(numPoints() * 3);
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

  Cells cells () const
  {
    auto const& indexSet = gridView_.indexSet();
    auto types = indexSet.types(0);
    int maxVertices = std::accumulate(types.begin(), types.end(), 1, [](int m, GeometryType t) {
      auto refElem = referenceElement<double,dim>(t);
      return std::max(m, refElem.size(dim));
    });

    Cells cells;
    cells.connectivity.reserve(numCells() * maxVertices);
    cells.offsets.reserve(numCells());
    cells.types.reserve(numCells());

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

  template <class T, class GlobalFunction>
  std::vector<T> pointData (GlobalFunction const& fct) const
  {
    std::vector<T> data(numPoints() * fct.ncomps());
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

  template <class T, class GlobalFunction>
  std::vector<T> cellData (GlobalFunction const& fct) const
  {
    std::vector<T> data(numCells() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_)) {
      localFct.bind(e);
      auto refElem = referenceElement(e.geometry());
      std::size_t idx = fct.ncomps() * indexSet.index(e);
      for (int comp = 0; comp < fct.ncomps(); ++comp)
        data[idx + comp] = T(localFct.evaluate(comp, refElem.position(0,0)));
      localFct.unbind();
    }
    return data;
  }

private:
  GridView gridView_;
};


/// Implementation of \ref DataCollector for linear cells, with discontinuous data.
template <class GridView>
class DiscontinuousDataCollector
{
  enum { dim = GridView::dimension };

public:
  DiscontinuousDataCollector (GridView const& gridView)
    : gridView_(gridView)
  {}

  void update ()
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

  std::uint64_t numPoints () const
  {
    return numPoints_;
  }

  std::uint64_t numCells () const
  {
    return gridView_.size(0);
  }

  template <class T>
  std::vector<T> points () const
  {
    std::vector<T> data(numPoints() * 3);
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

  Cells cells () const
  {
    Cells cells;
    cells.connectivity.reserve(numPoints());
    cells.offsets.reserve(numCells());
    cells.types.reserve(numCells());

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

  template <class T, class GlobalFunction>
  std::vector<T> pointData (GlobalFunction const& fct) const
  {
    std::vector<T> data(numPoints() * fct.ncomps());
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

  template <class T, class GlobalFunction>
  std::vector<T> cellData (GlobalFunction const& fct) const
  {
    std::vector<T> data(numCells() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_)) {
      localFct.bind(e);
      auto refElem = referenceElement(e.geometry());
      std::size_t idx = fct.ncomps() * indexSet.index(e);
      for (int comp = 0; comp < fct.ncomps(); ++comp)
        data[idx + comp] = T(localFct.evaluate(comp, refElem.position(0,0)));
      localFct.unbind();
    }
    return data;
  }

private:
  GridView gridView_;
  std::uint64_t numPoints_ = 0;
  std::vector<std::int64_t> indexMap_;
};


/// Implementation of \ref DataCollector for quadratic cells, with continuous data.
template <class GridView>
class QuadraticDataCollector
{
  enum { dim = GridView::dimension };

public:
  QuadraticDataCollector (GridView const& gridView)
    : gridView_(gridView)
  {}

  void update () {}

  std::uint64_t numPoints () const
  {
    return gridView_.size(dim) + gridView_.size(dim-1);
  }

  std::uint64_t numCells () const
  {
    return gridView_.size(0);
  }

  template <class T>
  std::vector<T> points () const
  {
    std::vector<T> data(numPoints() * 3);
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

  Cells cells () const
  {
    Cells cells;
    cells.connectivity.reserve(numPoints());
    cells.offsets.reserve(numCells());
    cells.types.reserve(numCells());

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

  template <class T, class GlobalFunction>
  std::vector<T> pointData (GlobalFunction const& fct) const
  {
    std::vector<T> data(numPoints() * fct.ncomps());
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

  template <class T, class GlobalFunction>
  std::vector<T> cellData (GlobalFunction const& fct) const
  {
    std::vector<T> data(numCells() * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_)) {
      localFct.bind(e);
      auto refElem = referenceElement(e.geometry());
      std::size_t idx = fct.ncomps() * indexSet.index(e);
      for (int comp = 0; comp < fct.ncomps(); ++comp)
        data[idx + comp] = T(localFct.evaluate(comp, refElem.position(0,0)));
      localFct.unbind();
    }
    return data;
  }

private:
  GridView gridView_;
};

}} // end namespace Dune::experimental
