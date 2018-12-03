#pragma once

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#if HAVE_VTK_ZLIB
#include <zlib.h>
#endif

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
VtkTimeseriesWriter<W>::~VtkTimeseriesWriter ()
{
  if (initialized_) {
    int ec = std::remove(filenameMesh_.c_str());
    assert(ec == 0);
    for (auto const& timestep : timesteps_) {
      ec = std::remove(timestep.second.c_str());
      assert(ec == 0);
    }
  }
  std::remove(tmpDir_.string().c_str());
}


template <class W>
void VtkTimeseriesWriter<W>
  ::writeTimestep (double time, std::string const& fn, Std::optional<std::string> tmpDir, bool writeCollection) const
{
  auto name = filesystem::path(fn).stem();
  auto tmp = tmpDir ? filesystem::path(*tmpDir) : tmpDir_;
  tmp /= name.string();

  vtkWriter_.dataCollector_.update();

  std::string filenameBase = tmp.string();

  int rank = vtkWriter_.rank_;
  int numRanks = vtkWriter_.numRanks_;
  if (numRanks > 1)
    filenameBase = tmp.string() + "_p" + std::to_string(rank);

  if (!initialized_) {
    // write points and cells only once
    filenameMesh_ = filenameBase + ".mesh.vtkdata";
    std::ofstream out(filenameMesh_, std::ios_base::ate | std::ios::binary);
    vtkWriter_.writeGridAppended(out, blocks_);
    initialized_ = true;
  }

  std::string filenameData = filenameBase + "_t" + std::to_string(timesteps_.size()) + ".vtkdata";
  std::ofstream out(filenameData, std::ios_base::ate | std::ios::binary);
  vtkWriter_.writeDataAppended(out, blocks_);
  timesteps_.emplace_back(time, filenameData);

  if (writeCollection)
    write(fn);
}


template <class W>
void VtkTimeseriesWriter<W>
  ::write (std::string const& fn, Std::optional<std::string> dir) const
{
  assert( initialized_ );

  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();

  filesystem::path fn_dir = p;
  filesystem::path data_dir = dir ? filesystem::path(*dir) : fn_dir;
  filesystem::path rel_dir = filesystem::relative(data_dir, fn_dir);

  std::string serial_fn = fn_dir.string() + '/' + name.string() + "_ts";
  std::string parallel_fn = data_dir.string() + '/' + name.string() + "_ts";
  std::string rel_fn = rel_dir.string() + '/' + name.string() + "_ts";

  int rank = vtkWriter_.rank_;
  int numRanks = vtkWriter_.numRanks_;
  if (numRanks > 1)
    serial_fn += "_p" + std::to_string(rank);

  { // write serial file
    std::ofstream serial_out(serial_fn + "." + vtkWriter_.getFileExtension(),
                             std::ios_base::ate | std::ios::binary);
    assert(serial_out.is_open());

    serial_out.imbue(std::locale::classic());
    serial_out << std::setprecision(vtkWriter_.getDatatype() == Vtk::FLOAT32
      ? std::numeric_limits<float>::digits10+2
      : std::numeric_limits<double>::digits10+2);

    vtkWriter_.writeTimeseriesSerialFile(serial_out, filenameMesh_, timesteps_, blocks_);
  }

  if (numRanks > 1 && rank == 0) {
    // write parallel file
    std::ofstream parallel_out(parallel_fn + ".p" + vtkWriter_.getFileExtension(),
                               std::ios_base::ate | std::ios::binary);
    assert(parallel_out.is_open());

    parallel_out.imbue(std::locale::classic());
    parallel_out << std::setprecision(vtkWriter_.getDatatype() == Vtk::FLOAT32
      ? std::numeric_limits<float>::digits10+2
      : std::numeric_limits<double>::digits10+2);

    vtkWriter_.writeTimeseriesParallelFile(parallel_out, rel_fn, numRanks, timesteps_);
  }
}

} // end namespace Dune
