#pragma once

#include <dune/vtk/datacollectors/continuousdatacollector.hh>

namespace Dune { namespace experimental
{

namespace Impl
{
  // Should be specialized for concrete structured grid
  template <class GridView, class Grid>
  struct StructuredDataCollector;
}

template <class GridView>
using StructuredDataCollector = typename Impl::StructuredDataCollector<GridView, typename GridView::Grid>::type;


/// The Interface for structured data-collectors
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
    , ghostLevel_(gridView.overlapSize(0))
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

  int ghostLevel () const
  {
    return this->asDerived().ghostLevelImpl();
  }

  /// Return the coordinates along the ordinates x, y, and z
  template <class T>
  std::array<std::vector<T>, 3> coordinates () const
  {
    return this->asDerived().template coordinatesImpl<T>();
  }


public:

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

  int ghostLevelImpl () const
  {
    return ghostLevel_;
  }

  template <class T>
  std::array<std::vector<T>, 3> coordinatesImpl () const
  {
    auto origin = this->origin();
    auto spacing = this->spacing();

    std::array<std::vector<T>, 3> ordinates{};
    writeLocalPiece([&ordinates,&origin,&spacing](auto const& extent) {
      for (std::size_t d = 0; d < GridView::dimension; ++d) {
        auto s = extent[2*d+1] - extent[2*d] + 1;
        ordinates[d].resize(s);
        for (std::size_t i = 0; i < s; ++i)
          ordinates[d][i] = origin[d] + (extent[2*d] + i)*spacing[d];
      }
    });

    for (std::size_t d = GridView::dimension; d < 3; ++d)
      ordinates[d].resize(1, T(0));

    return ordinates;
  }

private:
  DefaultDataCollector<GridView> defaultDataCollector_;
  int ghostLevel_;
};

}} // end namespace Dune::experimental
