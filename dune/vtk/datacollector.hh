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

}} // end namespace Dune::experimental
