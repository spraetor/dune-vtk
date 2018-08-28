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
void VtkRectilinearGridWriter<GV,DC>
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
      << " type=\"RectilinearGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  auto const& wholeExtent = dataCollector_.wholeExtent();
  out << "<RectilinearGrid"
      << " WholeExtent=\"" << join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << ">\n";

  dataCollector_.writeLocalPiece([&out](auto const& extent) {
    out << "<Piece Extent=\"" << join(extent.begin(), extent.end()) << "\">\n";
  });

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

  // Write point coordinates for x, y, and z ordinate
  out << "<Coordinates>\n";
  writeCoordinates(out, offsets);
  out << "</Coordinates>\n";

  out << "</Piece>\n";
  out << "</RectilinearGrid>\n";

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

    if (datatype_ == Vtk::FLOAT32) {
      auto bs = writeCoordinatesAppended<float>(out);
      blocks.insert(blocks.end(), bs.begin(), bs.end());
    } else {
      auto bs = writeCoordinatesAppended<double>(out);
      blocks.insert(blocks.end(), bs.begin(), bs.end());
    }
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
void VtkRectilinearGridWriter<GV,DC>
  ::writeParallelFile (std::string const& pfilename, int /*size*/) const
{
  std::string filename = pfilename + ".p" + this->fileExtension();
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  out << "<VTKFile"
      << " type=\"PRectilinearGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  auto const& wholeExtent = dataCollector_.wholeExtent();
  out << "<PRectilinearGrid"
      << " GhostLevel=\"" << dataCollector_.ghostLevel() << "\""
      << " WholeExtent=\"" << join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << ">\n";

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
        << " type=\"" <<  to_string(v.type()) << "\""
        << " NumberOfComponents=\"" << v.ncomps() << "\""
        << " />\n";
  }
  out << "</PCellData>\n";

  // Write point coordinates for x, y, and z ordinate
  out << "<PCoordinates>\n";
  out << "<PDataArray Name=\"x\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "<PDataArray Name=\"y\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "<PDataArray Name=\"z\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "</PCoordinates>\n";

  // Write piece file references
  dataCollector_.writePieces([&out,pfilename,ext=this->fileExtension()](int p, auto const& extent, bool write_extent)
  {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + ext;
    out << "<Piece Source=\"" << piece_source << "\"";
    if (write_extent)
      out << " Extent=\"" << join(extent.begin(), extent.end()) << "\"";
     out << " />\n";
  });

  out << "</PRectilinearGrid>\n";
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkRectilinearGridWriter<GV,DC>
  ::writeCoordinates (std::ofstream& out, std::vector<pos_type>& offsets) const
{
  std::string names = "xyz";
  if (format_ == Vtk::ASCII) {
    auto coordinates = dataCollector_.template coordinates<double>();
    for (std::size_t d = 0; d < 3; ++d) {
      out << "<DataArray type=\"" << to_string(datatype_) << "\" Name=\"" << names[d] << "\" format=\"ascii\">\n";
      std::size_t i = 0;
      for (auto const& c : coordinates[d])
        out << c << (++i % 6 != 0 ? ' ' : '\n');
      out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
    }
  }
  else { // Vtk::APPENDED format
    for (std::size_t j = 0; j < 3; ++j) {
      out << "<DataArray type=\"" << to_string(datatype_) << "\" Name=\"" << names[j] << "\" format=\"appended\"";
      out << " offset=";
      offsets.push_back(out.tellp());
      out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
      out << "/>\n";
    }
  }
}


template <class GV, class DC>
  template <class T>
std::array<std::uint64_t,3> VtkRectilinearGridWriter<GV,DC>
  ::writeCoordinatesAppended (std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto coordinates = dataCollector_.template coordinates<T>();

  // write conncetivity, offsets, and types
  std::uint64_t bs0 = this->writeAppended(out, coordinates[0]);
  std::uint64_t bs1 = this->writeAppended(out, coordinates[1]);
  std::uint64_t bs2 = this->writeAppended(out, coordinates[2]);

  return {bs0, bs1, bs2};
}

} // end namespace Dune
