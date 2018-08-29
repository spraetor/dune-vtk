#pragma once

#include <cstdint>
#include <vector>

#include <dune/vtk/datacollectorinterface.hh>

namespace Dune {

struct Cells
{
  std::vector<std::uint8_t> types;
  std::vector<std::int64_t> offsets;
  std::vector<std::int64_t> connectivity;
};

template <class GridView, class Derived>
class UnstructuredDataCollectorInterface
    : public DataCollectorInterface<GridView, Derived>
{
  using Super = DataCollectorInterface<GridView, Derived>;

public:
  UnstructuredDataCollectorInterface (GridView const& gridView)
    : Super(gridView)
  {}

  /// \brief Return the number of cells in the grid
  std::uint64_t numCells () const
  {
    return this->asDerived().numCellsImpl();
  }

  /// \brief Return cell types, offsets, and connectivity. \see Cells
  Cells cells () const
  {
    return this->asDerived().cellsImpl();
  }

protected:
  using Super::gridView_;
};

} // end namespace Dune
