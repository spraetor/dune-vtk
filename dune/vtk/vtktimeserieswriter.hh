#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/vtktypes.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/uid.hh>

namespace Dune
{
  /// File-Writer for Vtk timeseries .vtu files
  /**
   * \tparam VtkWriter  Type of a FileWriter derived from \ref VtkWriterInterface that
   *                    additionally supports writeTimeseriesSerialFile() and writeTimeseriesParallelFile(),
   *                    e.g. \ref VtkUnstructuredGridWriter.
   **/
  template <class VtkWriter>
  class VtkTimeseriesWriter
      : public FileWriter
  {
  protected:
    using Self = VtkTimeseriesWriter;
    using pos_type = typename std::ostream::pos_type;

    template <class W>
    using HasWriteTimeseriesSerialFile = decltype(&W::writeTimeseriesSerialFile);
    static_assert(Std::is_detected<HasWriteTimeseriesSerialFile, VtkWriter>::value, "");

    template <class W>
    using HasWriteTimeseriesParallelFile = decltype(&W::writeTimeseriesParallelFile);
    static_assert(Std::is_detected<HasWriteTimeseriesParallelFile, VtkWriter>::value, "");

  public:
    /// Constructor, stores the gridView
    template <class... Args, disableCopyMove<Self, Args...> = 0>
    VtkTimeseriesWriter (Args&&... args)
      : vtkWriter_{std::forward<Args>(args)...}
    {
      assert(vtkWriter_.format_ != Vtk::ASCII && "Timeseries writer requires APPENDED mode");
      std::srand(std::time(nullptr));
      // put temporary file to /tmp directory
      tmpDir_ = filesystem::path("/tmp/vtktimeserieswriter_" + uid(10) + "/");
      assert( filesystem::exists("/tmp") );
      filesystem::create_directories(tmpDir_);
    }

    ~VtkTimeseriesWriter ()
    {
      std::remove(tmpDir_.string().c_str());
    }

    /// Write the attached data to the file with \ref Vtk::FormatTypes and \ref Vtk::DataTypes
    void writeTimestep (double time, std::string const& fn);

    /// Writes all timesteps to single timeseries file.
    // NOTE: requires an aforging call to \ref writeTimestep
    virtual void write (std::string const& fn) override;

    /// Attach point data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    VtkTimeseriesWriter& addPointData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addPointData(fct, std::forward<Args>(args)...);
      return *this;
    }

    /// Attach cell data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    VtkTimeseriesWriter& addCellData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addCellData(fct, std::forward<Args>(args)...);
      return *this;
    }

  protected:
    VtkWriter vtkWriter_;
    filesystem::path tmpDir_;

    bool initialized_ = false;

    // block size of attached data
    std::vector<std::uint64_t> blocks_;

    std::string filenameMesh_;
    std::vector<std::pair<double, std::string>> timesteps_;
  };

} // end namespace Dune

#include "vtktimeserieswriter.impl.hh"