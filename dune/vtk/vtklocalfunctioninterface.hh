#pragma once

namespace Dune
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

} // end namespace Dune
