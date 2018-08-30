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
  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  vtkWriter_.dataCollector_.update();

  if (!initialized_) {
    // write points and cells only once
    filenameMesh_ = p.string() + ".mesh.vtkdata";
    std::ofstream file_mesh(filenameMesh_, std::ios_base::ate | std::ios::binary);

    if (vtkWriter_.datatype_ == Vtk::FLOAT32)
      blocksize_.push_back( vtkWriter_.template writePointsAppended<float>(file_mesh) );
    else
      blocksize_.push_back( vtkWriter_.template writePointsAppended<double>(file_mesh) );

    auto bs = vtkWriter_.writeCellsAppended(file_mesh);
    blocksize_.insert(blocksize_.end(), bs.begin(), bs.end());
    initialized_ = true;
  }

  std::string filenameData = p.string() + "_t" + std::to_string(timesteps_.size()) + ".vtkdata";
  std::ofstream file_data(filenameData, std::ios_base::ate | std::ios::binary);

  for (auto const& v : vtkWriter_.pointData_) {
    if (v.type() == Vtk::FLOAT32)
      blocksize_.push_back( vtkWriter_.template writeDataAppended<float>(file_data, v, W::POINT_DATA) );
    else
      blocksize_.push_back( vtkWriter_.template writeDataAppended<double>(file_data, v, W::POINT_DATA) );
  }
  for (auto const& v : vtkWriter_.cellData_) {
    if (v.type() == Vtk::FLOAT32)
      blocksize_.push_back( vtkWriter_.template writeDataAppended<float>(file_data, v, W::CELL_DATA) );
    else
      blocksize_.push_back( vtkWriter_.template writeDataAppended<double>(file_data, v, W::CELL_DATA) );
  }

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
  std::string filename = p.string() + ".ts." + vtkWriter_.getFileExtension();
  vtkWriter_.writeTimeseriesFile(filename, filenameMesh_, timesteps_, blocksize_);

  // remove all temporary data files
  int ec = std::remove(filenameMesh_.c_str());
  assert(ec == 0);
  for (auto const& timestep : timesteps_) {
    ec = std::remove(timestep.second.c_str());
    assert(ec == 0);
  }
}

} // end namespace Dune
