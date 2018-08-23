#pragma once

#include <dune/vtk/datacollector.hh>

namespace Dune { namespace experimental
{

/// Implementation of \ref DataCollector for linear cells, with continuous data.
template <class GridView>
class StructuredDataCollector
    : public DataCollectorInterface<GridView, StructuredDataCollector<GridView>>
{
  enum { dim = GridView::dimension };

  using Self = StructuredDataCollector;
  using Super = DataCollectorInterface<GridView, Self>;
  using Super::gridView_;
  using ctype = typename GridView::ctype;

  struct CoordLess
  {
    template <class T, int N>
    bool operator() (FieldVector<T,N> const& lhs, FieldVector<T,N> const& rhs) const
    {
      for (int i = N-1; i >= 0; --i) {
        if (std::abs(lhs[i] - rhs[i]) < std::numeric_limits<T>::epsilon())
          continue;
        return lhs[i] < rhs[i];
      }
      return false;
    }
  };

public:
  StructuredDataCollector (GridView const& gridView)
    : Super(gridView)
  {}

  /// Return number of grid vertices
  std::uint64_t numPointsImpl () const
  {
    return gridView_.size(dim);
  }

  std::array<int, 6> const& extent () const
  {
    return extent_;
  }

  auto const& origin () const
  {
    return origin_;
  }

  auto const& spacing () const
  {
    return spacing_;
  }

  void updateImpl ()
  {
    // TODO: extract this information from the grid
    int extent = GridView::dimension == 1 ? gridView_.size(dim) :
                 GridView::dimension == 2 ? isqrt(gridView_.size(dim)) :
                 GridView::dimension == 3 ? icbrt(gridView_.size(dim)) : 0;
    for (int i = 0; i < 3; ++i) {
      if (GridView::dimension > i) {
        extent_[2*i] = 0;
        extent_[2*i+1] = extent-1;
        spacing_[i] = gridView_.grid().domainSize()[i] / (extent-1);
      } else {
        extent_[2*i] = extent_[2*i+1] = 0;
        spacing_[i] = 0;
      }

      origin_[i] = 0;
    }
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
        int k = cellType.permutation(j);
        std::size_t idx = fct.ncomps() * indexSet.subIndex(e,k,dim);
        for (int comp = 0; comp < fct.ncomps(); ++comp)
          data[idx + comp] = T(localFct.evaluate(comp, refElem.position(k,dim)));
      }
      localFct.unbind();
    }
    return data;
  }

private:
  static constexpr std::uint32_t isqrt (std::uint64_t x)
  {
    std::uint32_t c = 0x8000;
    std::uint32_t g = 0x8000;

    while (true) {
      if (g*g > x)
        g ^= c;
      c >>= 1;
      if (c == 0)
        return g;
      g |= c;
    }
  }

  static constexpr std::uint32_t icbrt (std::uint64_t x)
  {
    std::uint32_t y = 0;
    for (int s = 63; s >= 0; s -= 3) {
      y += y;
      std::uint64_t b = 3*y*(std::uint64_t(y) + 1) + 1;
      if ((x >> s) >= b) {
        x -= b << s;
        y++;
      }
    }
    return y;
  }

private:
  std::array<int, 6> extent_;
  FieldVector<ctype,3> spacing_;
  FieldVector<ctype,3> origin_;
  std::vector<std::size_t> indexMap_;
};

}} // end namespace Dune::experimental
