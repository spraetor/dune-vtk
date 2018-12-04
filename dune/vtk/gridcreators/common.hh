#pragma once

#include <type_traits>

namespace Dune
{
  template <class Factory, class... Args>
  using HasInsertVertex = decltype( std::declval<Factory>().insertVertex(std::declval<Args>()...) );

  namespace Impl
  {
    template <class GF, class = void>
    struct VertexIdType { using type = unsigned int; };

    template <class GF>
    struct VertexIdType<GF, typename GF::VertexId> { using type = typename GF::VertexId; };
  }

  template <class GF>
  using VertexId_t = typename Impl::VertexIdType<GF>::type;

} // end namespace Dune
