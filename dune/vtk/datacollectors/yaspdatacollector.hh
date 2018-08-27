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
    , level_(gridView.template begin<0,All_Partition>()->level())
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
    localPieceCalled_ = false;

    for (int i = 0; i < dim; ++i) {
      wholeExtent_[2*i] = 0;
      wholeExtent_[2*i+1] = gridView_.grid().levelSize(level_,i);
    }

    auto it = gridView_.grid().begin(level_);
    initGeometry(it->coords);

#if HAVE_MPI
    int rank = -1;
    int num_ranks = -1;
    MPI_Comm_rank(gridView_.comm(), &rank);
    MPI_Comm_size(gridView_.comm(), &num_ranks);

    if (rank == 0) {
      extents_.resize(num_ranks);
      requests_.resize(num_ranks);
      for (int i = 1; i < num_ranks; ++i)
        MPI_Irecv(extents_[i].data(), extents_[i].size(), MPI_INT, i, /*tag=*/6, gridView_.comm(), &requests_[i]);
    }
#endif
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

#if HAVE_MPI
    if (!localPieceCalled_) {
      int rank = -1;
      MPI_Comm_rank(gridView_.comm(), &rank);
      if (rank != 0) {
        MPI_Isend(extent.data(), extent.size(), MPI_INT, 0, /*tag=*/6, gridView_.comm(), &sendRequest_);
      } else {
        extents_[0] = extent;
      }
      localPieceCalled_ = true;
    }
#endif

    writer(extent);
  }

  template <class Writer>
  void writePiecesImpl (Writer const& writer) const
  {
    assert(localPieceCalled_);
#if HAVE_MPI
    int num_ranks = -1;
    MPI_Comm_size(gridView_.comm(), &num_ranks);

    std::vector<MPI_Status> status(num_ranks-1);
    MPI_Waitall(num_ranks-1, requests_.data()+1, status.data());

    for (int p = 0; p < num_ranks; ++p) {
      writer(p, extents_[p], true);
    }
#endif
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

  mutable std::vector<std::array<int,6>> extents_;
  mutable std::vector<MPI_Request> requests_;
  mutable MPI_Request sendRequest_;
  mutable bool localPieceCalled_ = false;
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
