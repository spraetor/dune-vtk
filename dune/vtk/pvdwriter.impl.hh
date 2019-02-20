#pragma once

#include <iomanip>

#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
void PvdWriter<W>
  ::writeTimestep (double time, std::string const& fn, Std::optional<std::string> dir, bool writeCollection) const
{
  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();

  filesystem::path fn_dir = p;
  filesystem::path data_dir = dir ? filesystem::path(*dir) : fn_dir;
  filesystem::path rel_dir = filesystem::relative(data_dir, fn_dir);

  std::string pvd_fn = fn_dir.string() + '/' + name.string();
  std::string seq_fn = data_dir.string() + '/' + name.string() + "_t" + std::to_string(timesteps_.size());
  std::string rel_fn = rel_dir.string() + '/' + name.string() + "_t" + std::to_string(timesteps_.size());

  std::string ext = "." + vtkWriter_.getFileExtension();

  int commRank = vtkWriter_.comm().rank();
  int commSize = vtkWriter_.comm().size();
  if (commSize > 1)
    ext = ".p" + vtkWriter_.getFileExtension();

  timesteps_.emplace_back(time, rel_fn + ext);
  vtkWriter_.write(seq_fn + ext);

  if (commRank == 0 && writeCollection) {
    std::ofstream out(pvd_fn + ".pvd", std::ios_base::ate | std::ios::binary);
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
  ::write (std::string const& fn, Std::optional<std::string> /*dir*/) const
{
  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  int commRank = vtkWriter_.comm().rank();
  if (commRank == 0) {
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
