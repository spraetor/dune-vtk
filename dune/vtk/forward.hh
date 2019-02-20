#pragma once

namespace Dune
{
  // forward declaration of all classes in dune-vtk

  template <class GridView, class Derived, class Partition = Partitions::InteriorBorder>
  class DataCollectorInterface;

  // @{ datacollectors
  template <class GridView, class Derived, class Partition = Partitions::InteriorBorder>
  class UnstructuredDataCollectorInterface;

  // @{ unstructured-datacollectors
  template <class GridView, class Partition = Partitions::InteriorBorder>
  class ContinuousDataCollector;

  template <class GridView, class Partition = Partitions::InteriorBorder>
  class DiscontinuousDataCollector;

  template <class GridView>
  class QuadraticDataCollector;
  // @} unstructured-datacollectors

  template <class GridView, class Derived>
  class StructuredDataCollectorInterface;

  namespace Impl
  {
    // Should be specialized for concrete structured grid
    template <class GridView, class Grid>
    struct StructuredDataCollectorImpl;
  }

  template <class GridView>
  using StructuredDataCollector = typename Impl::StructuredDataCollectorImpl<GridView, typename GridView::Grid>::type;

  // @{ structured-datacollectors
  template <class GridView>
  class SPDataCollector;

  template <class GridView>
  class YaspDataCollector;
  // @} structured-datacollectors

  // @} datacollectors

  template <class Grid, class Derived>
  class GridCreatorInterface;

  template <class GridCreator, class Derived>
  struct DerivedGridCreator;

  // @{ gridcreators
  template <class Grid>
  struct ContinuousGridCreator;

  template <class Grid>
  struct DiscontinuousGridCreator;

  template <class Grid>
  struct ParallelGridCreator;

  template <class Grid>
  struct SerialGridCreator;
  // @} gridcreators


  template <class Grid, class FilerReaderImp>
  class FileReader;

  template <class Grid, class GridCreator = ContinuousGridCreator<Grid>>
  class VtkReader;


  class FileWriter;

  // @{ filewriters
  template <class VtkWriter>
  class PvdWriter;

  template <class VtkWriter>
  class VtkTimeseriesWriter;

  template <class GridView, class DataCollector>
  class VtkWriterInterface;

  // @{ vtkwriters
  template <class GridView, class DataCollector = StructuredDataCollector<GridView>>
  class VtkImageDataWriter;

  template <class GridView, class DataCollector = StructuredDataCollector<GridView>>
  class VtkRectilinearGridWriter;

  template <class GridView, class DataCollector = StructuredDataCollector<GridView>>
  class VtkStructuredGridWriter;

  template <class GridView, class DataCollector = ContinuousDataCollector<GridView>>
  class VtkUnstructuredGridWriter;
  // @} vtkwriters

  // @} filewriters

} // end namespace Dune
