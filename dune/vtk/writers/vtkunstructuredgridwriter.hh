#pragma once

#include <array>
#include <iosfwd>
#include <map>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/vtkfunction.hh>
#include <dune/vtk/vtktypes.hh>
#include <dune/vtk/datacollectors/continuousdatacollector.hh>

#include <dune/vtk/vtkwriterinterface.hh>

namespace Dune
{
  /// File-Writer for VTK .vtu files
  /**
   * Requirement:
   * - DataCollector must be a model of \ref DataCollector
   **/
  template <class GridView, class DataCollector = ContinuousDataCollector<GridView>>
  class VtkUnstructuredGridWriter
      : public VtkWriterInterface<GridView, DataCollector>
  {
    template <class> friend class VtkTimeseriesWriter;

    static constexpr int dimension = GridView::dimension;

    using Super = VtkWriterInterface<GridView, DataCollector>;
    using pos_type = typename Super::pos_type;

  public:
    /// Constructor, stores the gridView
    VtkUnstructuredGridWriter (GridView const& gridView,
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

    /// Write a series of timesteps in one file
    void writeTimeseriesFile (std::string const& filename,
                              std::string const& filenameMesh,
                              std::vector<std::pair<double, std::string>> const& timesteps,
                              std::vector<std::uint64_t> const& blocksize) const;

    virtual std::string fileExtension () const override
    {
      return "vtu";
    }

    virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const override;

    // Write the element connectivity to the output stream `out`. In case
    // of binary format, stores the streampos of XML attributes "offset" in the
    // vector `offsets`.
    void writeCells (std::ofstream& oust,
                     std::vector<pos_type>& offsets,
                     Std::optional<std::size_t> timestep = {}) const;

  private:
    using Super::dataCollector_;
    using Super::format_;
    using Super::datatype_;

    // attached data
    using Super::pointData_;
    using Super::cellData_;
  };

} // end namespace Dune

#include "vtkunstructuredgridwriter.impl.hh"
