#pragma once

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
void VtkTimeseriesWriter<W>
  ::writeTimestep (double time, std::string const& fn)
{
  auto name = filesystem::path(fn).stem();
  auto tmp = tmpDir_;
  tmp /= name.string();

  vtkWriter_.dataCollector_.update();

  std::string filenameBase = tmp.string();

  int rank = 0;
  int num_ranks = 1;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    if (num_ranks > 1)
      filenameBase = tmp.string() + "_p" + std::to_string(rank);
  #endif

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
}


template <class W>
void VtkTimeseriesWriter<W>
  ::write (std::string const& fn)
{
  assert( initialized_ );
  assert( timesteps_.size() < 1000 ); // maximal number of allowed timesteps in timeseries file

  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  std::string filenameBase = p.string() + "_ts";

  int rank = 0;
  int num_ranks = 1;
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    if (num_ranks > 1)
      filenameBase = p.string() + "_ts_p" + std::to_string(rank);
#endif

  std::string filename = filenameBase + "." + vtkWriter_.getFileExtension();
  vtkWriter_.writeTimeseriesSerialFile(filename, filenameMesh_, timesteps_, blocks_);

#ifdef HAVE_MPI
  if (num_ranks > 1 && rank == 0)
    vtkWriter_.writeTimeseriesParallelFile(p.string() + "_ts", num_ranks, timesteps_);
#endif

  // remove all temporary data files
  int ec = std::remove(filenameMesh_.c_str());
  assert(ec == 0);
  for (auto const& timestep : timesteps_) {
    ec = std::remove(timestep.second.c_str());
    assert(ec == 0);
  }
}

} // end namespace Dune
