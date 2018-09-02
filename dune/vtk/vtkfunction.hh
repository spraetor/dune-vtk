#pragma once

#include <type_traits>

#include <dune/common/std/optional.hh>
#include <dune/common/std/type_traits.hh>

#include "vtklocalfunction.hh"
#include "vtktypes.hh"

namespace Dune
{
  template <class T, int N>
  class FieldVector;

  template <class T, int N, int M>
  class FieldMatrix;

  namespace Impl
  {
    template <class T, class = void>
    struct SizeImpl
        : std::integral_constant<int, 1> {};

    template <class T, int N>
    struct SizeImpl<FieldVector<T,N>>
        : std::integral_constant<int, N> {};

    template <class T, int N, int M>
    struct SizeImpl<FieldMatrix<T,N,M>>
        : std::integral_constant<int, N*M> {};
  }

  template <class T>
  constexpr int Size = Impl::SizeImpl<std::decay_t<T>>::value;


  template <class GridView>
  class VtkFunction
  {
    template <class F>
    using HasLocalFunction = decltype(localFunction(std::declval<F>()));

    using Domain = typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;

    template <class F>
    using Range = std::decay_t<decltype(std::declval<F>()(std::declval<Domain>()))>;

  public:
    /// Constructor VtkFunction from legacy VTKFunction
    VtkFunction (std::shared_ptr<VTKFunction<GridView> const> const& fct)
      : localFct_(fct)
      , name_(fct->name())
      , ncomps_(fct->ncomps())
      , type_(Vtk::FLOAT64)
    {}

    /// Construct VtkFunction from dune-functions GridFunction with Signature
    template <class F,
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0>
    VtkFunction (F&& fct, std::string name,
                 Std::optional<int> ncomps = {},
                 Std::optional<Vtk::DataTypes> type = {})
      : localFct_(localFunction(std::forward<F>(fct)))
      , name_(std::move(name))
    {
      using R = Range<decltype(localFunction(std::forward<F>(fct)))>;

      ncomps_ = ncomps ? *ncomps : Size<R>;
      type_ = type ? *type : Vtk::Map::type<R>;
    }

    VtkFunction () = default;

    /// Create a LocalFunction
    friend VtkLocalFunction<GridView> localFunction (VtkFunction const& self)
    {
      return self.localFct_;
    }

    /// Return a name associated with the function
    std::string name () const
    {
      return name_;
    }

    /// Return the number of components of the Range
    int ncomps () const
    {
      return ncomps_ > 3 ? 9 : ncomps_ > 1 ? 3 : 1; // tensor, vector, scalar
    }

    /// Return the VTK Datatype associated with the functions range type
    Vtk::DataTypes type () const
    {
      return type_;
    }

  private:
    VtkLocalFunction<GridView> localFct_;
    std::string name_;
    int ncomps_ = 1;
    Vtk::DataTypes type_ = Vtk::FLOAT32;
  };

} // end namespace Dune
