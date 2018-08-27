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

namespace Dune { namespace experimental {

template <class GV, class DC>
void VtkImageDataWriter<GV,DC>
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
  auto const& origin = dataCollector_.origin();
  auto const& spacing = dataCollector_.spacing();
  out << "<ImageData"
      << " WholeExtent=\"" << join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << " Origin=\"" << join(origin.begin(), origin.end()) << "\""
      << " Spacing=\"" << join(spacing.begin(), spacing.end()) << "\""
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

  out << "</Piece>\n";
  out << "</ImageData>\n";

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
void VtkImageDataWriter<GV,DC>
  ::writeParallelFile (std::string const& pfilename, int /*size*/) const
{
  std::string filename = pfilename + ".p" + this->fileExtension();
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  assert(out.is_open());

  out << "<VTKFile"
      << " type=\"StructuredGrid\""
      << " version=\"1.0\""
      << " byte_order=\"" << this->getEndian() << "\""
      << " header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  auto const& wholeExtent = dataCollector_.wholeExtent();
  auto const& origin = dataCollector_.origin();
  auto const& spacing = dataCollector_.spacing();
  out << "<PImageData"
      << " GhostLevel=\"" << dataCollector_.ghostLevel() << "\""
      << " WholeExtent=\"" << join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << " Origin=\"" << join(origin.begin(), origin.end()) << "\""
      << " Spacing=\"" << join(spacing.begin(), spacing.end()) << "\""
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

  // Write piece file references
  dataCollector_.writePieces([&out,pfilename,ext=this->fileExtension()](int p, auto const& extent, bool write_extent)
  {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + ext;
    out << "<Piece Source=\"" << piece_source << "\"";
    if (write_extent)
      out << " Extent=\"" << join(extent.begin(), extent.end()) << "\"";
     out << " />\n";
  });

  out << "</PImageData>\n";
  out << "</VTKFile>";
}

}} // end namespace Dune::experimental
