#pragma once

#include <memory>

#include <dune/grid/io/file/vtk/function.hh>

#include "vtklocalfunctioninterface.hh"

namespace Dune { namespace experimental
{
  /// Type erasure for Legacy VTKFunction
  template <class GridView>
  class VTKLocalFunctionWrapper final
      : public VtkLocalFunctionInterface<GridView>
  {
    using Interface = VtkLocalFunctionInterface<GridView>;
    using Entity = typename Interface::Entity;
    using LocalCoordinate = typename Interface::LocalCoordinate;

  public:
    VTKLocalFunctionWrapper (std::shared_ptr<VTKFunction<GridView> const> const& fct)
      : fct_(fct)
    {}

    virtual void bind (Entity const& entity) override
    {
      entity_ = &entity;
    }

    virtual void unbind () override
    {
      entity_ = nullptr;
    }

    virtual double evaluate (int comp, LocalCoordinate const& xi) const override
    {
      return fct_->evaluate(comp, *entity_, xi);
    }

  private:
    std::shared_ptr<VTKFunction<GridView> const> fct_;
    Entity const* entity_;
  };

}} // end namespace Dune::experimental
