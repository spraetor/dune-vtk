#pragma once

#include "vtklocalfunctioninterface.hh"

namespace Dune
{
  /// Type erasure for dune-functions LocalFunction interface
  template <class GridView, class LocalFunction>
  class LocalFunctionWrapper final
      : public VtkLocalFunctionInterface<GridView>
  {
    using Self = LocalFunctionWrapper;
    using Interface = VtkLocalFunctionInterface<GridView>;
    using Entity = typename Interface::Entity;
    using LocalCoordinate = typename Interface::LocalCoordinate;

    template <class F, class D>
    using Range = std::decay_t<decltype(std::declval<F>()(std::declval<D>()))>;

    template <class F, class D>
    using VectorValued = decltype(std::declval<Range<F,D>>()[0u]);

  public:
    template <class LocalFct, disableCopyMove<Self, LocalFct> = 0>
    LocalFunctionWrapper (LocalFct&& localFct)
      : localFct_(std::forward<LocalFct>(localFct))
    {}

    virtual void bind (Entity const& entity) override
    {
      localFct_.bind(entity);
    }

    virtual void unbind () override
    {
      localFct_.unbind();
    }

    virtual double evaluate (int comp, LocalCoordinate const& xi) const override
    {
      return evaluateImpl(comp, xi, Std::is_detected<VectorValued,LocalFunction,LocalCoordinate>{});
    }

  private:
    // Evaluate a component of a vector valued data
    double evaluateImpl (int comp, LocalCoordinate const& xi, std::true_type) const
    {
      auto y = localFct_(xi);
      return comp < y.size() ? y[comp] : 0.0;
    }

    // Return the scalar values
    double evaluateImpl (int comp, LocalCoordinate const& xi, std::false_type) const
    {
      assert(comp == 0);
      return localFct_(xi);
    }

  private:
    LocalFunction localFct_;
  };

} // end namespace Dune
