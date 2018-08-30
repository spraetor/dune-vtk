#pragma once

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <dune/common/std/optional.hh>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/vtkfunction.hh>
#include <dune/vtk/vtktypes.hh>

namespace Dune
{
  /// File-Writer for Vtk .vtu files
  template <class VtkWriter>
  class VtkTimeseriesWriter
  {
  protected:
    using Self = VtkTimeseriesWriter;
    using pos_type = typename std::ostream::pos_type;

  public:
    /// Constructor, stores the gridView
    template <class... Args, disableCopyMove<Self, Args...> = 0>
    VtkTimeseriesWriter (Args&&... args)
      : vtkWriter_{std::forward<Args>(args)...}
    {}

    /// Write the attached data to the file with \ref Vtk::FormatTypes and \ref Vtk::DataTypes
    void writeTimestep (double time, std::string const& fn);

    // NOTE: requires a aforging call to \ref write
    void write (std::string const& fn);

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
