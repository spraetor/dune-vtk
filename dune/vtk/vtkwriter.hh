#pragma once

#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/writers/vtkrectilineargridwriter.hh>
#include <dune/vtk/writers/vtkstructuredgridwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#include <dune/vtk/datacollectors/spdatacollector.hh>
#endif

#include <dune/grid/geometrygrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>

namespace Dune { namespace experimental
{
  namespace Impl
  {
    // The default writer assumes an unstructured grid
    template <class GridView, class Grid>
    struct VtkWriterImpl
    {
      using type = VtkUnstructuredGridWriter<GridView>;
    };

#if HAVE_DUNE_SPGRID
    // A structured grid with constant spacing in x, y, and z direction.
    template <class GridView, class ct, int dim, template <int> class Ref, class Comm>
    struct VtkWriterImpl<GridView, SPGrid<ct,dim,Ref,Comm>>
    {
      using type = VtkImageDataWriter<GridView, SPDataCollector<GridView>>;
    };
#endif

    // A structured grid with constant spacing in x, y, and z direction.
    template <class GridView, int dim, class Coordinates>
    struct VtkWriterImpl<GridView, YaspGrid<dim,Coordinates>>
    {
      using type = VtkImageDataWriter<GridView, YaspDataCollector<GridView>>;
    };

    // A structured grid with coordinates in x, y, and z direction with arbitrary spacing
    template <class GridView, int dim, class ct>
    struct VtkWriterImpl<GridView, YaspGrid<dim,TensorProductCoordinates<ct,dim>>>
    {
      using type = VtkRectilinearGridWriter<GridView, YaspDataCollector<GridView>>;
    };

    // A transformed structured grid has structured connectivity but unstructured point
    // coordinates.
    template <class GridView, int dim, class Coordinates, class CoordFunction, class Allocator>
    struct VtkWriterImpl<GridView, GeometryGrid<YaspGrid<dim,Coordinates>, CoordFunction, Allocator>>
    {
      using type = VtkStructuredGridWriter<GridView, YaspDataCollector<GridView>>;
    };

  } // end namespace Impl


  /// Default choice for several grid types, uses the default data-collector.
  template <class GridView>
  using VtkWriter = Impl::VtkWriterImpl<GridView, typename GridView::Grid>;

}} // end namespace Dune::experimental
