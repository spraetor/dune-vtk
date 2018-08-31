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

  // Write point coordinates for x, y, and z ordinate
  out << "<Coordinates>\n";
  writeCoordinates(out, offsets);
  out << "</Coordinates>\n";

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

  out << "</Piece>\n";
  out << "</RectilinearGrid>\n";

  this->writeAppended(out, offsets);
  out << "</VTKFile>";
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

  // Write point coordinates for x, y, and z ordinate
  out << "<PCoordinates>\n";
  out << "<PDataArray Name=\"x\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "<PDataArray Name=\"y\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "<PDataArray Name=\"z\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "</PCoordinates>\n";

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
  ::writeCoordinates (std::ofstream& out, std::vector<pos_type>& offsets,
                      Std::optional<std::size_t> timestep) const
{
  std::string names = "xyz";
  if (format_ == Vtk::ASCII) {
    auto coordinates = dataCollector_.template coordinates<double>();
    for (std::size_t d = 0; d < 3; ++d) {
      out << "<DataArray type=\"" << to_string(datatype_) << "\" Name=\"" << names[d] << "\" format=\"ascii\"";
      if (timestep)
        out << " TimeStep=\"" << *timestep << "\"";
      out << ">\n";
      std::size_t i = 0;
      for (auto const& c : coordinates[d])
        out << c << (++i % 6 != 0 ? ' ' : '\n');
      out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
    }
  }
  else { // Vtk::APPENDED format
    for (std::size_t j = 0; j < 3; ++j) {
      out << "<DataArray type=\"" << to_string(datatype_) << "\" Name=\"" << names[j] << "\" format=\"appended\"";
      if (timestep)
        out << " TimeStep=\"" << *timestep << "\"";
      out << " offset=";
      offsets.push_back(out.tellp());
      out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
      out << "/>\n";
    }
  }
}


template <class GV, class DC>
void VtkRectilinearGridWriter<GV,DC>
  ::writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  // write coordinates along axis
  if (datatype_ == Vtk::FLOAT32) {
    auto coordinates = dataCollector_.template coordinates<float>();
    blocks.push_back(this->writeValuesAppended(out, coordinates[0]));
    blocks.push_back(this->writeValuesAppended(out, coordinates[1]));
    blocks.push_back(this->writeValuesAppended(out, coordinates[2]));
  } else {
    auto coordinates = dataCollector_.template coordinates<double>();
    blocks.push_back(this->writeValuesAppended(out, coordinates[0]));
    blocks.push_back(this->writeValuesAppended(out, coordinates[1]));
    blocks.push_back(this->writeValuesAppended(out, coordinates[2]));
  }
}

} // end namespace Dune
