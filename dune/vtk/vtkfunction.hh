#pragma once

#include <type_traits>

#include <dune/common/std/type_traits.hh>
#include <dune/functions/common/typeerasure.hh>

namespace Dune { namespace experimental
{
  /// An abstract base class for LocalFunctions
  template <class GridView>
  class VTKLocalFunctionInterface
  {
  public:
    using Entity = typename GridView::template Codim<0>::Entity;
    using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

    /// Bind the function to the grid entity
    virtual void bind (Entity const& entity) = 0;

    /// Unbind from the currently bound entity
    virtual void unbind () = 0;

    /// Evaluate single component comp in the entity at local coordinates xi
    virtual double evaluate (int comp, LocalCoordinate const& xi) const = 0;

    /// Virtual destructor
    virtual ~VTKLocalFunctionInterface () = default;
  };


  template <class GridView>
  struct VTKLocalFunctionImpl
  {
    template <class Wrapper>
    class Model : public Wrapper
    {
    public:
      using Wrapper::Wrapper;
      using Function = typename Wrapper::Wrapped;
      using Interface = VTKLocalFunctionInterface<GridView>;

      using Entity = typename Interface::Entity;
      using LocalCoordinate = typename Interface::LocalCoordinate;

      template <class F, class D>
      using Range = std::decay_t<decltype(std::declval<F>()(std::declval<D>()))>;

      template <class F, class D>
      using VectorValued = decltype(std::declval<Range<F,D>>()[0u]);

      virtual void bind (Entity const& entity) override
      {
        this->get().bind(entity);
      }

      virtual void unbind () override
      {
        this->get().unbind();
      }

      virtual double evaluate (int comp, LocalCoordinate const& xi) const override
      {
        return evaluateImpl(comp, xi, Std::is_detected<VectorValued,Function,LocalCoordinate>{});
      }

    private:
      // Evaluate a component of a vector valued data
      double evaluateImpl (int comp, LocalCoordinate const& xi, std::true_type) const
      {
        auto y = this->get()(xi);
        return comp < y.size() ? y[comp] : 0.0;
      }

      // Return the scalar values
      double evaluateImpl (int comp, LocalCoordinate const& xi, std::false_type) const
      {
        assert(comp == 0);
        return this->get()(xi);
      }
    };
  };


  template <class GridView>
  class VTKLocalFunction
      : public Functions::TypeErasureBase<VTKLocalFunctionInterface<GridView>,
                                          VTKLocalFunctionImpl<GridView>::template Model>
  {
    using Super = Functions::TypeErasureBase<VTKLocalFunctionInterface<GridView>,
                                             VTKLocalFunctionImpl<GridView>::template Model>;

    using Entity = typename GridView::template Codim<0>::Entity;
    using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

  public:
    template <class F, disableCopyMove<VTKLocalFunction, F> = 0>
    VTKLocalFunction (F&& f)
      : Super(std::forward<F>(f))
    {}

    VTKLocalFunction () = default;

    /// Bind the function to the grid entity
    void bind (Entity const& entity)
    {
      this->asInterface().bind(entity);
    }

    /// Unbind from the currently bound entity
    void unbind ()
    {
      this->asInterface().unbind();
    }

    /// Evaluate the `comp` component of the Range value at local coordinate `xi`
    double evaluate (int comp, LocalCoordinate const& xi) const
    {
      return this->asInterface().evaluate(comp, xi);
    }
  };

  // ---------------------------------------------------------------------------

  /// An abstract base class for GlobalFunctions
  template <class GridView>
  class VTKFunctionInterface
  {
  public:
    /// Create a local function
    virtual VTKLocalFunction<GridView> makeLocalFunction () const = 0;

    /// Virtual destructor
    virtual ~VTKFunctionInterface () = default;
  };


  template <class GridView>
  struct VTKFunctionImpl
  {
    template <class Wrapper>
    class Model : public Wrapper
    {
    public:
      using Wrapper::Wrapper;
      virtual VTKLocalFunction<GridView> makeLocalFunction () const override
      {
        return VTKLocalFunction<GridView>{localFunction(this->get())};
      }
    };
  };


  template <class GridView>
  class VTKFunction
      : public Functions::TypeErasureBase<VTKFunctionInterface<GridView>,
                                          VTKFunctionImpl<GridView>::template Model>
  {
    using Super = Functions::TypeErasureBase<VTKFunctionInterface<GridView>,
                                             VTKFunctionImpl<GridView>::template Model>;

    template <class F>
    using HasLocalFunction = decltype(localFunction(std::declval<F>()));

  public:
    template <class F>
    VTKFunction (F&& f, std::string name, int ncomps = 1)
      : Super(std::forward<F>(f))
      , name_(std::move(name))
      , ncomps_(ncomps > 1 ? 3 : 1)
    {
      assert(1 <= ncomps && ncomps <= 3);
      static_assert(Std::is_detected<HasLocalFunction,F>::value,
        "Requires A GridFunction to be passed to the VTKFunction.");
    }

    VTKFunction () = default;

    /// Create a LocalFunction
    friend VTKLocalFunction<GridView> localFunction (VTKFunction const& self)
    {
      return self.asInterface().makeLocalFunction();
    }

    /// Return a name associated with the function
    std::string name () const
    {
      return name_;
    }

    /// Return the number of components of the Range
    int ncomps () const
    {
      return ncomps_;
    }

  private:
    std::string name_;

    int ncomps_ = 1;
  };

}} // end namespace Dune::experimental
