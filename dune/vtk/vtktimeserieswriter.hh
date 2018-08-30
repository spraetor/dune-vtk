#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/vtktypes.hh>

namespace Dune
{
  /// File-Writer for Vtk .vtu files
  template <class VtkWriter>
  class VtkTimeseriesWriter
      : public FileWriter
  {
  protected:
    using Self = VtkTimeseriesWriter;
    using pos_type = typename std::ostream::pos_type;

  public:
    /// Constructor, stores the gridView
    template <class... Args, disableCopyMove<Self, Args...> = 0>
    VtkTimeseriesWriter (Args&&... args)
      : vtkWriter_{std::forward<Args>(args)...}
    {
      assert(vtkWriter_.format_ != Vtk::ASCII && "Timeseries writer requires APPENDED mode");
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

    bool initialized_ = false;

    // block size of attached data
    std::vector<std::uint64_t> blocksize_;

    std::string filenameMesh_;
    std::vector<std::pair<double, std::string>> timesteps_;
  };

} // end namespace Dune

#include "vtktimeserieswriter.impl.hh"
