#pragma once

#include <iomanip>

#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
void PvdWriter<W>
  ::writeTimestep (double time, std::string const& fn, bool writeCollection) const
{
  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  std::string ext = "." + vtkWriter_.getFileExtension();
  std::string filename = p.string() + "_t" + std::to_string(timesteps_.size());

  int rank = 0;
  int num_ranks = 1;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  if (num_ranks > 1)
    ext = ".p" + vtkWriter_.getFileExtension();
#endif

  timesteps_.emplace_back(time, filename + ext);
  vtkWriter_.write(filename + ext);

  if (rank == 0 && writeCollection) {
    std::ofstream out(p.string() + ".pvd", std::ios_base::ate | std::ios::binary);
    assert(out.is_open());

    out.imbue(std::locale::classic());
    out << std::setprecision(datatype_ == Vtk::FLOAT32
      ? std::numeric_limits<float>::digits10+2
      : std::numeric_limits<double>::digits10+2);

    writeFile(out);
  }
}


template <class W>
void PvdWriter<W>
  ::write (std::string const& fn) const
{
  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  int rank = 0;
  int num_ranks = 1;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

  if (rank == 0) {
    std::ofstream out(p.string() + ".pvd", std::ios_base::ate | std::ios::binary);
    assert(out.is_open());

    out.imbue(std::locale::classic());
    out << std::setprecision(datatype_ == Vtk::FLOAT32
      ? std::numeric_limits<float>::digits10+2
      : std::numeric_limits<double>::digits10+2);

    writeFile(out);
  }
}


template <class W>
void PvdWriter<W>
  ::writeFile (std::ofstream& out) const
{
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile"
      << " type=\"Collection\""
      << " version=\"0.1\""
      << (format_ != Vtk::ASCII ? " byte_order=\"" + vtkWriter_.getEndian() + "\"" : "")
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  out << "<Collection>\n";

  // Write all timesteps
  for (auto const& timestep : timesteps_) {
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
