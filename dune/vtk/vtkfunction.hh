#pragma once

#include <type_traits>

#include <dune/common/std/type_traits.hh>
#include <dune/functions/common/signature.hh>

#include "vtklocalfunction.hh"

namespace Dune { namespace experimental
{
  template <class GridView>
  class VtkFunction
  {
    template <class F>
    using HasLocalFunction = decltype(localFunction(std::declval<F>()));

    template <class F>
    using Domain = typename std::decay_t<F>::EntitySet::GlobalCoordinate;

    template <class F>
    using Range = std::decay_t<decltype(std::declval<F>()(std::declval<Domain<F>>()))>;

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
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0,
      std::enable_if_t<Std::is_detected<Range,F>::value,int> = 0>
    VtkFunction (F&& fct, std::string name, int ncomps = 1, Vtk::DataTypes type = Vtk::Map::type<Range<F>>)
      : localFct_(localFunction(std::forward<F>(fct)))
      , name_(std::move(name))
      , ncomps_(ncomps)
      , type_(type)
    {}

    /// Construct VtkFunction from dune-functions GridFunction without Signature
    template <class F,
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0,
      std::enable_if_t<not Std::is_detected<Range,F>::value,int> = 0>
    VtkFunction (F const& fct, std::string name, int ncomps = 1, Vtk::DataTypes type = Vtk::FLOAT32)
      : localFct_(localFunction(std::forward<F>(fct)))
      , name_(std::move(name))
      , ncomps_(ncomps)
      , type_(type)
    {}

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
      return ncomps_;
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

}} // end namespace Dune::experimental
