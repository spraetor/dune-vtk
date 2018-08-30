#pragma once

#include <iomanip>

#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
void PvdWriter<W>
  ::write (double time, std::string const& fn, Vtk::FormatTypes format, Vtk::DataTypes datatype)
{
  format_ = format;
  datatype_ = datatype;

#ifndef HAVE_ZLIB
  if (format_ == Vtk::COMPRESSED)
    format_ = Vtk::BINARY;
#endif

  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  std::string ext = "." + vtkWriter_.getFileExtension();
  std::string filename = p.string() + "_t" + std::to_string(timeSeries_.size());

  int rank = 0;
  int num_ranks = 1;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  if (num_ranks > 1)
    ext = ".p" + vtkWriter_.getFileExtension();
#endif

  timeSeries_.emplace_back(time, filename + ext);
  vtkWriter_.write(filename + ext, format_, datatype_);

  if (rank == 0)
    writeFile(time, p.string() + ".pvd");
}


template <class W>
void PvdWriter<W>
  ::writeFile (double time, std::string const& filename) const
{
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  if (datatype_ == Vtk::FLOAT32)
    out << std::setprecision(std::numeric_limits<float>::digits10+2);
  else
    out << std::setprecision(std::numeric_limits<double>::digits10+2);

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile"
      << " type=\"Collection\""
      << " version=\"0.1\""
      << (format_ != Vtk::ASCII ? " byte_order=\"" << vtkWriter_.getEndian() << "\"" : "")
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  out << "<Collection>\n";

  // Write all timesteps
  for (auto const& timestep : timeSeries_) {
    out << "<DataSet"
        << " timestep=\"" << timestep.first << "\""
        << " part=\"0\""
        << " file=\"" << timestep.second << "\""
        << " />\n";
  }

  out << "</Collection>\n";
  out << "</VTKFile>";
}

} // end namespace Dune
