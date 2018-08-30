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
void VtkStructuredGridWriter<GV,DC>
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
      << " type=\"StructuredGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  auto const& wholeExtent = dataCollector_.wholeExtent();
  out << "<StructuredGrid WholeExtent=\"" << join(wholeExtent.begin(), wholeExtent.end()) << "\">\n";

  dataCollector_.writeLocalPiece([&out](auto const& extent) {
    out << "<Piece Extent=\"" << join(extent.begin(), extent.end()) << "\">\n";
  });

  // Write point coordinates
  out << "<Points>\n";
  this->writePoints(out, offsets);
  out << "</Points>\n";

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
  out << "</StructuredGrid>\n";

  this->writeAppended(out, offsets);
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkStructuredGridWriter<GV,DC>
  ::writeParallelFile (std::string const& pfilename, int /*size*/) const
{
  std::string filename = pfilename + ".p" + this->fileExtension();
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  out << "<VTKFile"
      << " type=\"PStructuredGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  auto const& wholeExtent = dataCollector_.wholeExtent();
  out << "<PStructuredGrid"
      << " GhostLevel=\"" << dataCollector_.ghostLevel() << "\""
      << " WholeExtent=\"" << join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << ">\n";

  // Write points
  out << "<PPoints>\n";
  out << "<PDataArray"
      << " type=\"" << to_string(datatype_) << "\""
      << " NumberOfComponents=\"3\""
      << " />\n";
  out << "</PPoints>\n";

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

  // Write piece file references
  dataCollector_.writePieces([&out,pfilename,ext=this->fileExtension()](int p, auto const& extent, bool write_extent)
  {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + ext;
    out << "<Piece Source=\"" << piece_source << "\"";
    if (write_extent)
      out << " Extent=\"" << join(extent.begin(), extent.end()) << "\"";
    out << " />\n";
  });

  out << "</PStructuredGrid>\n";
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkStructuredGridWriter<GV,DC>
  ::writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  // write points
  if (datatype_ == Vtk::FLOAT32) {
    auto points = dataCollector_.template points<float>();
    blocks.push_back(this->writeValuesAppended(out, points));
  } else {
    auto points = dataCollector_.template points<double>();
    blocks.push_back(this->writeValuesAppended(out, points));
  }
}

} // end namespace Dune
