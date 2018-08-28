#pragma once

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#include "structureddatacollector.hh"

namespace Dune
{
#if HAVE_DUNE_SPGRID

// Specialization for SPGrid
template <class GridView>
class SPDataCollector
    : public StructuredDataCollectorInterface<GridView, SPDataCollector<GridView>>
{
  enum { dim = GridView::dimension };

  using Self = SPDataCollector;
  using Super = StructuredDataCollectorInterface<GridView, Self>;
  using ctype = typename GridView::ctype;

public:
  SPDataCollector (GridView const& gridView)
    : Super(gridView)
    , wholeExtent_(filledArray<6,int>(0))
    , extent_(filledArray<6,int>(0))
    , origin_(0.0)
    , spacing_(0.0)
  {}

  std::array<int, 6> const& wholeExtentImpl () const
  {
    return wholeExtent_;
  }

  std::array<int, 6> const& extentImpl () const
  {
    return extent_;
  }

  auto const& originImpl () const
  {
    return origin_;
  }

  auto const& spacingImpl () const
  {
    return spacing_;
  }

  void updateImpl ()
  {
    Super::updateImpl();

    auto const& localMesh = gridView_.impl().gridLevel().localMesh();
    auto const& begin = gridView_.impl().gridLevel().globalMesh().begin();
    auto const& end = gridView_.impl().gridLevel().globalMesh().end();
    auto const& cube = gridView_.grid().domain().cube();

    for (int i = 0; i < dim; ++i) {
      wholeExtent_[2*i] = begin[i];
      wholeExtent_[2*i+1] = end[i];
      extent_[2*i] = localMesh.begin()[i];
      extent_[2*i+1] = localMesh.end()[i];
      spacing_[i] = cube.width()[i] / (end[i] - begin[i]);
      origin_[i] = cube.origin()[i];
    }
  }

protected:
  using Super::gridView_;
  std::array<int, 6> wholeExtent_;
  std::array<int, 6> extent_;
  FieldVector<ctype,3> spacing_;
  FieldVector<ctype,3> origin_;
  std::vector<std::size_t> indexMap_;
};

namespace Impl
{
  template <class GridView, class ct, int dim, template< int > class Ref, class Comm>
  struct StructuredDataCollectorImpl<GridView, SPGrid<ct,dim,Ref,Comm>>
  {
    using type = SPDataCollector<GridView>;
  };
}

#endif // HAVE_DUNE_SPGRID
} // end namespace Dune
