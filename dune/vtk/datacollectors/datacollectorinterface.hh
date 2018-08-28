#pragma once

#include <dune/vtk/vtktypes.hh>

namespace Dune {

template <class GridView, class Derived>
class DataCollectorInterface
{
public:
  DataCollectorInterface (GridView const& gridView)
    : gridView_(gridView)
  {}

  /// Update the DataCollector on the current GridView
  void update ()
  {
    asDerived().updateImpl();
  }

  /// Return the number of overlapping elements
  int ghostLevel () const
  {
    return asDerived().ghostLevelImpl();
  }

  /// \brief Return the number of points in the grid
  std::uint64_t numPoints () const
  {
    return asDerived().numPointsImpl();
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

  /// \brief Return a flat vector of function values evaluated at the grid points.
  /**
   * In case of a vector valued function, flat the vector entries:
   * [fct(p0)_0, fct(p0)_1, fct(p0)_2, fct(p1)_0, ...]
   * where the vector dimension must be 3 (possible extended by 0s)
   * In case of tensor valued function, flat the tensor row-wise:
   * [fct(p0)_00, fct(p0)_01, fct(p0)_02, fct(p0)_10, fct(p0)_11, fct(p0)_12, fct(p0)_20...]
   * where the tensor dimension must be 3x3 (possible extended by 0s)
   **/
  template <class T, class VtkFunction>
  std::vector<T> pointData (VtkFunction const& fct) const
  {
    return asDerived().template pointDataImpl<T>(fct);
  }

  /// \brief Return a flat vector of function values evaluated at the grid cells. \see pointData.
  template <class T, class VtkFunction>
  std::vector<T> cellData (VtkFunction const& fct) const
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

public: // default implementations

  void updateImpl ()
  {
    /* do nothing */
  }

  int ghostLevelImpl () const
  {
    return gridView_.overlapSize(0);
  }

  // Evaluate `fct` in center of cell
  template <class T, class VtkFunction>
  std::vector<T> cellDataImpl (VtkFunction const& fct) const
  {
    std::vector<T> data(gridView_.size(0) * fct.ncomps());
    auto const& indexSet = gridView_.indexSet();
    auto localFct = localFunction(fct);
    for (auto const& e : elements(gridView_, Partitions::all)) {
      localFct.bind(e);
      auto refElem = referenceElement<T,GridView::dimension>(e.type());
      std::size_t idx = fct.ncomps() * indexSet.index(e);
      for (int comp = 0; comp < fct.ncomps(); ++comp)
        data[idx + comp] = T(localFct.evaluate(comp, refElem.position(0,0)));
      localFct.unbind();
    }
    return data;
  }

protected:
  GridView gridView_;
};

} // end namespace Dune
