#pragma once

#include <iosfwd>
#include <string>
#include <vector>
#include <tuple>

#include <dune/vtk/vtktypes.hh>
// #include <dune/vtk/vtkwriterinterface.hh>

namespace Dune
{
  /// File-Writer for ParaView .pvd files
  template <class VtkWriter>
  class PvdWriter
  {
    using Self = PvdWriter;

    // static_assert(IsVtkWriter<VtkWriter>::value, "Writer must implement the VtkWriterInterface.");

  public:
    /// Constructor, creates a VtkWriter with constructor arguments forwarded
    template <class... Args, disableCopyMove<Self,Args...> = 0>
    explicit PvdWriter (Args&&... args)
      : vtkWriter_{std::forward<Args>(args)...}
    {
      format_ = vtkWriter_.getFormat();
      datatype_ = vtkWriter_.getDatatype();
    }

    /// Write the attached data to the file
    void write (double time, std::string const& fn) const;

    /// Attach point data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    PvdWriter& addPointData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addPointData(fct, std::forward<Args>(args)...);
      return *this;
    }

    /// Attach cell data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    PvdWriter& addCellData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addCellData(fct, std::forward<Args>(args)...);
      return *this;
    }

  protected:
    /// Write a series of vtk files in a .pvd ParaView Data file
    void writeFile (double time, std::string const& filename) const;

  protected:
    VtkWriter vtkWriter_;
    Vtk::FormatTypes format_;
    Vtk::DataTypes datatype_;

    mutable std::vector<std::pair<double, std::string>> timesteps_;
  };

} // end namespace Dune

#include "pvdwriter.impl.hh"
