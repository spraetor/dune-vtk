#pragma once

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class GV, class DC>
void VtkUnstructuredGridWriter<GV,DC>
  ::writeSerialFile (std::string const& filename) const
{
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  if (format_ == Vtk::ASCII) {
    if (datatype_ == Vtk::FLOAT32)
      out << std::setprecision(std::numeric_limits<float>::digits10+2);
    else
      out << std::setprecision(std::numeric_limits<double>::digits10+2);
  }

  std::vector<pos_type> offsets; // pos => offset
  out << "<VTKFile"
      << " type=\"UnstructuredGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  out << "<UnstructuredGrid>\n";
  out << "<Piece"
      << " NumberOfPoints=\"" << dataCollector_.numPoints() << "\""
      << " NumberOfCells=\"" << dataCollector_.numCells() << "\""
      << ">\n";

  // Write data associated with grid points
  out << "<PointData" << this->getNames(pointData_) << ">\n";
  for (auto const& v : pointData_)
    this->writeData(out, offsets, v, Super::POINT_DATA);
  out << "</PointData>\n";

  // Write data associated with grid cells
  out << "<CellData" << this->getNames(cellData_) << ">\n";
  for (auto const& v : cellData_)
    this->writeData(out, offsets, v, Super::CELL_DATA);
  out << "</CellData>\n";

  // Write point coordinates
  out << "<Points>\n";
  this->writePoints(out, offsets);
  out << "</Points>\n";

  // Write element connectivity, types and offsets
  out << "<Cells>\n";
  writeCells(out, offsets);
  out << "</Cells>\n";

  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";

  std::vector<std::uint64_t> blocks; // size of i'th appended block
  pos_type appended_pos = 0;
  if (is_a(format_, Vtk::APPENDED)) {
    out << "<AppendedData encoding=\"raw\">\n_";
    appended_pos = out.tellp();
    for (auto const& v : pointData_) {
      if (v.type() == Vtk::FLOAT32)
        blocks.push_back( this->template writeDataAppended<float>(out, v, Super::POINT_DATA) );
      else
        blocks.push_back( this->template writeDataAppended<double>(out, v, Super::POINT_DATA) );
    }
    for (auto const& v : cellData_) {
      if (v.type() == Vtk::FLOAT32)
        blocks.push_back( this->template writeDataAppended<float>(out, v, Super::CELL_DATA) );
      else
        blocks.push_back( this->template writeDataAppended<double>(out, v, Super::CELL_DATA) );
    }

    if (datatype_ == Vtk::FLOAT32)
      blocks.push_back( this->template writePointsAppended<float>(out) );
    else
      blocks.push_back( this->template writePointsAppended<double>(out) );

    auto bs = writeCellsAppended(out);
    blocks.insert(blocks.end(), bs.begin(), bs.end());
    out << "</AppendedData>\n";
  }

  out << "</VTKFile>";

  // fillin offset values and block sizes
  if (is_a(format_, Vtk::APPENDED)) {
    pos_type offset = 0;
    for (std::size_t i = 0; i < offsets.size(); ++i) {
      out.seekp(offsets[i]);
      out << '"' << offset << '"';
      offset += pos_type(blocks[i]);
    }
  }
}


template <class GV, class DC>
void VtkUnstructuredGridWriter<GV,DC>
  ::writeParallelFile (std::string const& pfilename, int size) const
{
  std::string filename = pfilename + ".pvtu";
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  out << "<VTKFile"
      << " type=\"PUnstructuredGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  out << "<PUnstructuredGrid GhostLevel=\"0\">\n";

  // Write data associated with grid points
  out << "<PPointData" << this->getNames(pointData_) << ">\n";
  for (auto const& v : pointData_) {
    out << "<PDataArray"
        << " Name=\"" << v.name() << "\""
        << " type=\"" << to_string(v.type()) << "\""
        << " NumberOfComponents=\"" << v.ncomps() << "\""
        << " />\n";
  }
  out << "</PPointData>\n";

  // Write data associated with grid cells
  out << "<PCellData" << this->getNames(cellData_) << ">\n";
  for (auto const& v : cellData_) {
    out << "<PDataArray"
        << " Name=\"" << v.name() << "\""
        << " type=\"" << to_string(v.type()) << "\""
        << " NumberOfComponents=\"" << v.ncomps() << "\""
        << " />\n";
  }
  out << "</PCellData>\n";

  // Write points
  out << "<PPoints>\n";
  out << "<PDataArray"
      << " type=\"" << to_string(datatype_) << "\""
      << " NumberOfComponents=\"3\""
      << " />\n";
  out << "</PPoints>\n";

  // Write piece file references
  for (int p = 0; p < size; ++p) {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + this->fileExtension();
    out << "<Piece Source=\"" << piece_source << "\" />\n";
  }

  out << "</PUnstructuredGrid>\n";
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkUnstructuredGridWriter<GV,DC>
  ::writeTimeseriesFile (std::string const& filename,
                         std::string const& filenameMesh,
                         std::vector<std::pair<double, std::string>> const& timesteps,
                         std::vector<std::uint64_t> const& blocksize) const
{
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  assert(is_a(format_, Vtk::APPENDED));
  if (datatype_ == Vtk::FLOAT32)
    out << std::setprecision(std::numeric_limits<float>::digits10+2);
  else
    out << std::setprecision(std::numeric_limits<double>::digits10+2);

  std::vector<std::vector<pos_type>> offsets(timesteps.size()); // pos => offset
  out << "<VTKFile"
      << " type=\"UnstructuredGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  out << "<UnstructuredGrid"
      << " TimeValues=\"";
  {
    std::size_t i = 0;
    for (auto const& timestep : timesteps)
      out << timestep.first << (++i % 6 != 0 ? ' ' : '\n');
  }
  out << "\">\n";

  out << "<Piece"
      << " NumberOfPoints=\"" << dataCollector_.numPoints() << "\""
      << " NumberOfCells=\"" << dataCollector_.numCells() << "\""
      << ">\n";

  // Write point coordinates
  out << "<Points>\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    out << "<DataArray"
        << " type=\"" << to_string(datatype_) << "\""
        << " NumberOfComponents=\"3\""
        << " TimeStep=\"" << i << "\""
        << " format=\"appended\""
        << " offset=";
    offsets[i].push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
  out << "</Points>\n";

  // Write element connectivity, types and offsets
  out << "<Cells>\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\""
        << " TimeStep=\"" << i << "\""
        << " offset=";
    offsets[i].push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\""
        << " TimeStep=\"" << i << "\""
        << " offset=";
    offsets[i].push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\""
        << " TimeStep=\"" << i << "\""
        << " offset=";
    offsets[i].push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
  out << "</Cells>\n";

  const std::size_t shift = offsets[0].size(); // number of blocks to write the grid

  // Write data associated with grid points
  out << "<PointData" << this->getNames(pointData_) << ">\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    for (auto const& v : pointData_) {
      out << "<DataArray"
          << " Name=\"" << v.name() << "\""
          << " TimeStep=\"" << i << "\""
          << " type=\"" << to_string(datatype_) << "\""
          << " NumberOfComponents=\"" << v.ncomps() << "\""
          << " format=\"appended\""
          << " offset=";
      offsets[i].push_back(out.tellp());
      out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
      out << "/>\n";
    }
  }
  out << "</PointData>\n";

  // Write data associated with grid cells
  out << "<CellData" << this->getNames(cellData_) << ">\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    for (auto const& v : cellData_) {
      out << "<DataArray"
          << " Name=\"" << v.name() << "\""
          << " TimeStep=\"" << i << "\""
          << " type=\"" << to_string(datatype_) << "\""
          << " NumberOfComponents=\"" << v.ncomps() << "\""
          << " format=\"appended\""
          << " offset=";
      offsets[i].push_back(out.tellp());
      out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
      out << "/>\n";
    }
  }
  out << "</CellData>\n";

  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";

  std::vector<std::uint64_t> blocks; // size of i'th appended block
  pos_type appended_pos = 0;
  out << "<AppendedData encoding=\"raw\">\n_";
  appended_pos = out.tellp();

  std::ifstream file_mesh(filenameMesh, std::ios_base::in | std::ios_base::binary);
  out << file_mesh.rdbuf();
  file_mesh.close();
  assert( std::uint64_t(out.tellp()) == std::accumulate(blocksize.begin(), std::next(blocksize.begin(),shift), std::uint64_t(appended_pos)) );

  for (auto const& timestep : timesteps) {
    std::ifstream file(timestep.second, std::ios_base::in | std::ios_base::binary);
    out << file.rdbuf();
  }
  out << "</AppendedData>\n";

  out << "</VTKFile>";

  // write correct offsets in file.
  pos_type offset = 0;
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    offset = 0;
    auto const& off = offsets[i];

    // write mesh data offsets
    for (std::size_t j = 0; j < shift; ++j) {
      out.seekp(off[j]);
      out << '"' << offset << '"';
      offset += pos_type(blocksize[j]);
    }
  }

  std::size_t j = shift;
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    auto const& off = offsets[i];

    for (std::size_t k = shift; k < off.size(); ++k) {
      out.seekp(off[k]);
      out << '"' << offset << '"';
      offset += pos_type(blocksize[j++]);
    }
  }
}


template <class GV, class DC>
void VtkUnstructuredGridWriter<GV,DC>
  ::writeCells (std::ofstream& out, std::vector<pos_type>& offsets) const
{
  if (format_ == Vtk::ASCII) {
    auto cells = dataCollector_.cells();
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    std::size_t i = 0;
    for (auto const& c : cells.connectivity)
      out << c << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    i = 0;
    for (auto const& o : cells.offsets)
      out << o << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    i = 0;
    for (auto const& t : cells.types)
      out << int(t) << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
  }
  else { // Vtk::APPENDED format
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
std::array<std::uint64_t,3> VtkUnstructuredGridWriter<GV,DC>
  ::writeCellsAppended (std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto cells = dataCollector_.cells();

  // write conncetivity, offsets, and types
  std::uint64_t bs0 = this->writeAppended(out, cells.connectivity);
  std::uint64_t bs1 = this->writeAppended(out, cells.offsets);
  std::uint64_t bs2 = this->writeAppended(out, cells.types);

  return {bs0, bs1, bs2};
}


} // end namespace Dune
