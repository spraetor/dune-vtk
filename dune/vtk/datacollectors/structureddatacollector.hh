#pragma once

#include <dune/vtk/datacollectors/continuousdatacollector.hh>

namespace Dune { namespace experimental
{
template <class GridView, class Derived>
class StructuredDataCollectorInterface
    : public DataCollectorInterface<GridView, Derived>
{
protected:
  using Super = DataCollectorInterface<GridView, Derived>;
  using Super::gridView_;
  using ctype = typename GridView::ctype;

public:
  StructuredDataCollectorInterface (GridView const& gridView)
    : Super(gridView)
    , defaultDataCollector_(gridView)
  {}

  /// Return number of grid vertices
  std::uint64_t numPointsImpl () const
  {
    return gridView_.size(GridView::dimension);
  }

  void updateImpl ()
  {
    defaultDataCollector_.update();
  }

  std::array<int, 6> const& wholeExtent () const
  {
    return this->asDerived().wholeExtentImpl();
  }

  FieldVector<ctype, 3> const& origin () const
  {
    return this->asDerived().originImpl();
  }

  FieldVector<ctype, 3> const& spacing () const
  {
    return this->asDerived().spacingImpl();
  }

  template <class Writer>
  void writeLocalPiece (Writer const& writer) const
  {
    this->asDerived().writeLocalPieceImpl(writer);
  }

  template <class Writer>
  void writePieces (Writer const& writer) const
  {
    this->asDerived().writePiecesImpl(writer);
  }

  /// Return the coordinates of all grid vertices in the order given by the indexSet
  template <class T>
  std::vector<T> pointsImpl () const
  {
    return defaultDataCollector_.template points<T>();
  }

  /// Evaluate the `fct` at the corners of the elements
  template <class T, class GlobalFunction>
  std::vector<T> pointDataImpl (GlobalFunction const& fct) const
  {
    return defaultDataCollector_.template pointData<T>(fct);
  }

private:
  DefaultDataCollector<GridView> defaultDataCollector_;
};

}} // end namespace Dune::experimental
