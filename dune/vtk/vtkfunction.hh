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


  /// Wrapper class for functions allowing local evaluations.
  template <class GridView>
  class VtkFunction
  {
    template <class F>
    using HasLocalFunction = decltype(localFunction(std::declval<F>()));

    using Domain = typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;

    template <class F>
    using Range = std::decay_t<decltype(std::declval<F>()(std::declval<Domain>()))>;

  private:

    template <class T, int N>
    static auto sizeOfImpl (FieldVector<T,N> const&)
      -> std::integral_constant<int, N> { return {}; }

    template <class T, int N, int M>
    static auto sizeOfImpl (FieldMatrix<T,N,M> const&)
      -> std::integral_constant<int, N*M> { return {}; };

    static auto sizeOfImpl (...)
      -> std::integral_constant<int, 1> { return {}; }

    template <class T>
    static constexpr int sizeOf () { return decltype(sizeOfImpl(std::declval<T>()))::value; }

  public:
    /// Constructor VtkFunction from legacy VTKFunction
    /**
     * \param fct   The VTKFunction to wrap
     * \param type  The VTK datatype how to write the function values to the output [Vtk::FLOAT64]
     **/
    VtkFunction (std::shared_ptr<VTKFunction<GridView> const> const& fct,
                 Std::optional<Vtk::DataTypes> type = {})
      : localFct_(fct)
      , name_(fct->name())
      , ncomps_(fct->ncomps())
      , type_(type ? *type : Vtk::FLOAT64)
    {}

    /// Construct VtkFunction from dune-functions GridFunction with Signature
    // NOTE: Stores the localFunction(fct) by value.
    /**
     * \param fct     A Grid(View)-function, providing a `localFunction(fct)`
     * \param name    The name to use component identification in the VTK file
     * \param ncomps  Number of components of the pointwise data. Is extracted
     *                from the range type of the GridFunction if not given.
     * \param type    The \ref Vtk::DataTypes used in the output. E.g. FLOAT32,
     *                or FLOAT64. Is extracted from the range type of the
     *                GridFunction if not given.
     **/
    template <class F,
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0>
    VtkFunction (F&& fct, std::string name,
                 Std::optional<int> ncomps = {},
                 Std::optional<Vtk::DataTypes> type = {})
      : localFct_(localFunction(std::forward<F>(fct)))
      , name_(std::move(name))
    {
      using R = Range<decltype(localFunction(std::forward<F>(fct)))>;

      ncomps_ = ncomps ? *ncomps : sizeOf<R>();
      type_ = type ? *type : Vtk::Map::type<R>();
    }

    template <class F,
      std::enable_if_t<Std::is_detected<HasLocalFunction,F>::value, int> = 0>
    VtkFunction (F&& fct, Vtk::FieldInfo fieldInfo,
                 Std::optional<Vtk::DataTypes> type = {})
      : VtkFunction(std::forward<F>(fct), fieldInfo.name(), fieldInfo.ncomps(), type)
    {}

    VtkFunction () = default;

    /// Create a LocalFunction
    friend VtkLocalFunction<GridView> localFunction (VtkFunction const& self)
    {
      return self.localFct_;
    }

    /// Return a name associated with the function
    std::string const& name () const
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
