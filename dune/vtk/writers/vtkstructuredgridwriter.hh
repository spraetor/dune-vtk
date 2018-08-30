#pragma once

#include <array>
#include <iosfwd>
#include <map>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/vtkfunction.hh>
#include <dune/vtk/vtktypes.hh>
#include <dune/vtk/datacollectors/structureddatacollector.hh>

#include <dune/vtk/vtkwriterinterface.hh>

namespace Dune
{
  /// File-Writer for StructuredGrid VTK .vts files
  /**
   * Requirement:
   * - DataCollector must be a model of \ref StructuredDataCollector
   **/
  template <class GridView, class DataCollector = StructuredDataCollector<GridView>>
  class VtkStructuredGridWriter
      : public VtkWriterInterface<GridView, DataCollector>
  {
    static constexpr int dimension = GridView::dimension;

    using Super = VtkWriterInterface<GridView, DataCollector>;
    using pos_type = typename Super::pos_type;

  public:
    /// Constructor, stores the gridView
    VtkStructuredGridWriter (GridView const& gridView,
                             Vtk::FormatTypes format = Vtk::BINARY,
                             Vtk::DataTypes datatype = Vtk::FLOAT32)
      : Super(gridView, format, datatype)
    {}

  private:
    /// Write a serial VTK file in Unstructured format
    virtual void writeSerialFile (std::string const& filename) const override;

    /// Write a parallel VTK file `pfilename.pvtu` in Unstructured format,
    /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
    /// for [i] in [0,...,size).
    virtual void writeParallelFile (std::string const& pfilename, int size) const override;

    virtual std::string fileExtension () const override
    {
      return "vts";
    }

  private:
    using Super::dataCollector_;
    using Super::format_;
    using Super::datatype_;

    // attached data
    using Super::pointData_;
    using Super::cellData_;
  };

} // end namespace Dune

#include "vtkstructuredgridwriter.impl.hh"
