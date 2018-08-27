#pragma once

#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/structureddatacollector.hh>

namespace Dune { namespace experimental
{
// Specialization for YaspGrid
template <class GridView>
class YaspDataCollector
    : public StructuredDataCollectorInterface<GridView, YaspDataCollector<GridView>>
{
  enum { dim = GridView::dimension };

  using Self = YaspDataCollector;
  using Super = StructuredDataCollectorInterface<GridView, Self>;
  using Super::gridView_;
  using ctype = typename GridView::ctype;

public:
  YaspDataCollector (GridView const& gridView)
    : Super(gridView)
    , wholeExtent_(filledArray<6,int>(0))
    , origin_(0.0)
    , spacing_(0.0)
    , level_(gridView.template begin<0,Interior_Partition>()->level())
  {}

  std::array<int, 6> const& wholeExtentImpl () const
  {
    return wholeExtent_;
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

    for (int i = 0; i < dim; ++i) {
      wholeExtent_[2*i] = 0;
      wholeExtent_[2*i+1] = gridView_.grid().levelSize(level_,i);
    }

    auto it = gridView_.grid().begin(level_);
    initGeometry(it->coords);
  }

  template <class Coords>
  void initGeometry (Coords const& coords) { DUNE_THROW(NotImplemented, "Coordinate-Type not implemented!"); }

  void initGeometry (EquidistantCoordinates<ctype,dim> const& coords)
  {
    for (int i = 0; i < dim; ++i) {
      spacing_[i] = coords.meshsize(i,0);
      origin_[i] = 0;
    }
  }

  void initGeometry (EquidistantOffsetCoordinates<ctype,dim> const& coords)
  {
    for (int i = 0; i < dim; ++i) {
      spacing_[i] = coords.meshsize(i,0);
      origin_[i] = coords.origin(i);
    }
  }

  void initGeometry (TensorProductCoordinates<ctype,dim> const& coords)
  {
    for (int i = 0; i < dim; ++i) {
      spacing_[i] = coords.meshsize(i,0); // is not constant, but also not used.
      origin_[i] = coords.coordinate(i,0);
    }
  }

  template <class Writer>
  void writeLocalPieceImpl (Writer const& writer) const
  {
    auto const& gl = *gridView_.grid().begin(level_);
    auto const& g = gl.interior[0];
    auto extent = filledArray<6,int>(0);
    auto const& gc = *g.dataBegin();
    for (int i = 0; i < dim; ++i) {
      extent[2*i] = gc.min(i);
      extent[2*i+1] = gc.max(i)+1;
    }
    writer(extent);
  }

  template <class Writer>
  void writePiecesImpl (Writer const& writer) const
  {
    auto extent = filledArray<6,int>(0);
    // for (auto const& part : gridView_.gridLevel().template partition<InteriorEntity>())
    // {
    //   for (int i = 0; i < dim; ++i) {
    //     extent[2*i] = part.begin()[i];
    //     extent[2*i+1] = part.end()[i]-1;
    //   }

    int num_ranks = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    for (int p = 0; p < num_ranks; ++p) {
      writer(p, extent, false);
    }
  }

  template <class T>
  std::array<std::vector<T>, 3> coordinatesImpl () const
  {
    auto it = gridView_.grid().begin(level_);
    auto const& coords = it->coords;

    std::array<std::vector<T>, 3> ordinates{};
    writeLocalPieceImpl([&ordinates,&coords](auto const& extent)
    {
      for (std::size_t d = 0; d < dim; ++d) {
        auto s = extent[2*d+1] - extent[2*d] + 1;
        ordinates[d].resize(s);
        for (std::size_t i = 0; i < s; ++i)
          ordinates[d][i] = coords.coordinate(d, extent[2*d] + i);
      }
    });

    for (std::size_t d = dim; d < 3; ++d)
      ordinates[d].resize(1, T(0));

    return ordinates;
  }

private:
  std::array<int, 6> wholeExtent_;
  FieldVector<ctype,3> spacing_;
  FieldVector<ctype,3> origin_;
  std::vector<std::size_t> indexMap_;
  int level_;
};

namespace Impl
{
  template<class GridView, int dim, class Coordinates>
  struct StructuredDataCollector<GridView, YaspGrid<dim,Coordinates>>
  {
    using type = YaspDataCollector<GridView>;
  };
}

}} // end namespace Dune::experimental
