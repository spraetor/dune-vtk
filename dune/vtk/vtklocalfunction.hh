#pragma once

#include <memory>
#include <type_traits>

#include <dune/common/std/type_traits.hh>

#include "vtklocalfunctioninterface.hh"
#include "legacyvtkfunction.hh"
#include "defaultvtkfunction.hh"

namespace Dune { namespace experimental
{
  template <class GridView>
  class VtkLocalFunction
  {
    using Self = VtkLocalFunction;
    using Entity = typename GridView::template Codim<0>::Entity;
    using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

    template <class LF, class E>
    using HasBind = decltype(std::declval<LF>().bind(std::declval<E>()));

  public:
    template <class LF, disableCopyMove<Self, LF> = 0,
      std::enable_if_t<Std::is_detected<HasBind,LF,Entity>::value, int> = 0>
    VtkLocalFunction (LF&& lf)
      : localFct_(std::make_shared<LocalFunctionWrapper<GridView,LF>>(std::forward<LF>(lf)))
    {}

    VtkLocalFunction (std::shared_ptr<VTKFunction<GridView> const> const& lf)
      : localFct_(std::make_shared<VTKLocalFunctionWrapper<GridView>>(lf))
    {}

    VtkLocalFunction () = default;

    /// Bind the function to the grid entity
    void bind (Entity const& entity)
    {
      localFct_->bind(entity);
    }

    /// Unbind from the currently bound entity
    void unbind ()
    {
      localFct_->unbind();
    }

    /// Evaluate the `comp` component of the Range value at local coordinate `xi`
    double evaluate (int comp, LocalCoordinate const& xi) const
    {
      return localFct_->evaluate(comp, xi);
    }

  private:
    std::shared_ptr<VtkLocalFunctionInterface<GridView>> localFct_;
  };

}} // end namespace Dune::experimental
