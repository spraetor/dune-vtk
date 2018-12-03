#pragma once

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

template <class GridView, class Derived>
  template <class T, class VtkFunction>
std::vector<T> DataCollectorInterface<GridView, Derived>
  ::cellDataImpl (VtkFunction const& fct) const
{
  std::vector<T> data(gridView_.size(0) * fct.ncomps());
  MultipleCodimMultipleGeomTypeMapper<GridView> mapper(gridView_, mcmgElementLayout());
  auto localFct = localFunction(fct);
  for (auto const& e : elements(gridView_, Partitions::all)) {
    localFct.bind(e);
    auto refElem = referenceElement<T,GridView::dimension>(e.type());
    std::size_t idx = fct.ncomps() * mapper.index(e);
    for (int comp = 0; comp < fct.ncomps(); ++comp)
      data[idx + comp] = T(localFct.evaluate(comp, refElem.position(0,0)));
    localFct.unbind();
  }
  return data;
}

} // end namespace Dune
