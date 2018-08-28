#pragma once

#include <type_traits>

#include <dune/common/std/type_traits.hh>
#include <dune/functions/common/signature.hh>
#include <dune/grid/io/file/vtk/function.hh>

namespace Dune { namespace experimental
{
  /// An abstract base class for LocalFunctions
  template <class GridView>
  class VtkLocalFunctionInterface
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
    virtual ~VtkLocalFunctionInterface () = default;
  };


  /// Type erasure for dune-functions LocalFunction interface
  template <class GridView, class LocalFunction>
  class LocalFunctionWrapper
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


  /// Type erasure for Legacy VTKFunction
  template <class GridView>
  class VTKLocalFunctionWrapper
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


  template <class GridView>
  class VtkLocalFunction
  {
    using Self = VtkLocalFunction;
    using Entity = typename GridView::template Codim<0>::Entity;
    using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

    template <class LF, class E>
    using HasBind = decltype(std::declval<LF>().bind(std::declval<E>()));

  public:
    // Store wrapper around dune-function LocalFunction
    template <class LocalFct, disableCopyMove<Self, LocalFct> = 0,
      std::enable_if_t<Std::is_detected<HasBind, LocalFct, Entity>::value,int> = 0>
    VtkLocalFunction (LocalFct&& lf)
      : localFct_(std::make_unique<LocalFunctionWrapper<GridView, std::decay_t<LocalFct>>>(std::forward<LocalFct>(lf)))
    {}

    // store wrapper around legacy VTKFunction
    VtkLocalFunction (std::shared_ptr<VTKFunction<GridView> const> const& lf)
      : localFct_(std::make_unique<VTKLocalFunctionWrapper<GridView>>(lf))
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

  // ---------------------------------------------------------------------------

  /// An abstract base class for GlobalFunctions
  template <class GridView>
  class VtkFunctionInterface
  {
  public:
    /// Create a local function
    virtual VtkLocalFunction<GridView> makeLocalFunction () const = 0;

    /// Virtual destructor
    virtual ~VtkFunctionInterface () = default;
  };


  template <class GridView, class GridViewFunction>
  class GridViewFunctionWrapper
      : public VtkFunctionInterface<GridView>
  {
    using Self = GridViewFunctionWrapper;

  public:
    template <class GVFct, disableCopyMove<Self, GVFct> = 0>
    GridViewFunctionWrapper (GVFct&& gvFct)
      : gvFct_(std::forward<GVFct>(gvFct))
    {}

    virtual VtkLocalFunction<GridView> makeLocalFunction () const override
    {
      return VtkLocalFunction<GridView>{localFunction(gvFct_)};
    }

  private:
    GridViewFunction gvFct_;
  };


  template <class GridView>
  class VTKFunctionWrapper
      : public VtkFunctionInterface<GridView>
  {
  public:
    VTKFunctionWrapper (std::shared_ptr<VTKFunction<GridView> const> const& fct)
      : fct_(fct)
    {}

    virtual VtkLocalFunction<GridView> makeLocalFunction () const override
    {
      return VtkLocalFunction<GridView>{fct_};
    }

  private:
    std::shared_ptr<VTKFunction<GridView> const> fct_;
  };


  template <class GridView>
  class VtkFunction
  {
    template <class F>
    using HasLocalFunction = decltype(localFunction(std::declval<F>()));

    template <class F>
    using Signature = typename std::decay_t<F>::Signature;

  public:
    template <class F,
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0,
      class Range = typename Functions::SignatureTraits<Signature<F>>::Range>
    VtkFunction (F&& f, std::string name, int ncomps = 1,
                 Vtk::DataTypes type = Vtk::Map::type<Range>)
      : fct_(std::make_unique<GridViewFunctionWrapper<GridView,std::decay_t<F>>>(std::forward<F>(f)))
      , name_(std::move(name))
      , ncomps_(ncomps > 3 ? 9 : ncomps > 1 ? 3 : 1) // tensor, vector, or scalar
      , type_(type)
    {}

    template <class F,
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0,
      std::enable_if_t<not Std::is_detected<Signature,F>::value,int> = 0>
    VtkFunction (F&& f, std::string name, int ncomps = 1,
                 Vtk::DataTypes type = Vtk::FLOAT32)
      : fct_(std::make_unique<GridViewFunctionWrapper<GridView,std::decay_t<F>>>(std::forward<F>(f)))
      , name_(std::move(name))
      , ncomps_(ncomps > 3 ? 9 : ncomps > 1 ? 3 : 1) // tensor, vector, or scalar
      , type_(type)
    {}

    VtkFunction (std::shared_ptr<VTKFunction<GridView> const> const& fct,
                 std::string name, int ncomps = 1,
                 Vtk::DataTypes type = Vtk::FLOAT32)
      : fct_(std::make_unique<VTKFunctionWrapper<GridView>>(fct))
      , name_(std::move(name))
      , ncomps_(ncomps > 3 ? 9 : ncomps > 1 ? 3 : 1) // tensor, vector, or scalar
      , type_(type)
    {}

    VtkFunction () = default;

    /// Create a LocalFunction
    friend VtkLocalFunction<GridView> localFunction (VtkFunction const& self)
    {
      return self.fct_->makeLocalFunction();
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

    /// Return the VTK Datatype associated with the functions range type
    Vtk::DataTypes type () const
    {
      return type_;
    }

  private:
    std::shared_ptr<VtkFunctionInterface<GridView>> fct_;
    std::string name_;
    int ncomps_ = 1;
    Vtk::DataTypes type_;
  };

}} // end namespace Dune::experimental
