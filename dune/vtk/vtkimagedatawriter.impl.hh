#pragma once

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include "utility/enum.hh"
#include "utility/filesystem.hh"
#include "utility/string.hh"

namespace Dune { namespace experimental {

template <class GV, class DC>
void VtkImageDataWriter<GV,DC>
  ::writeSerialFile (std::string const& filename) const
{
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  if (format_ == Vtk::ASCII) {
    if (datatype_ == Vtk::FLOAT32)
      out << std::setprecision(std::numeric_limits<float>::digits10+2);
    else
      out << std::setprecision(std::numeric_limits<double>::digits10+2);
  }

  std::vector<pos_type> offsets; // pos => offset
  out << "<VTKFile type=\"ImageData\" version=\"1.0\" "
      << "byte_order=\"" << this->getEndian() << "\" header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<ImageData WholeExtent=\"";
  auto const& wholeExtent = dataCollector_.wholeExtent();
  for (int i = 0; i < 3; ++i) {
    out << (i == 0 ? "" : " ") << wholeExtent[2*i] << " " << wholeExtent[2*i+1];
  }
  out << "\" Origin=\"";
  auto const& origin = dataCollector_.origin();
  for (int i = 0; i < 3; ++i) {
    out << (i == 0 ? "" : " ") << origin[i];
  }
  out << "\" Spacing=\"";
  auto const& spacing = dataCollector_.spacing();
  for (int i = 0; i < 3; ++i) {
    out << (i == 0 ? "" : " ") << spacing[i];
  }
  out << "\">\n";

  dataCollector_.writeLocalPiece([&out](std::array<int,6> const& extent)
  {
    out << "<Piece Extent=\"";
    for (int i = 0; i < 3; ++i) {
      out << (i == 0 ? "" : " ") << extent[2*i] << " " << extent[2*i+1];
    }
    out << "\">\n";
  });

  { // Write data associated with grid points
    auto scalar = std::find_if(pointData_.begin(), pointData_.end(), [](auto const& v) { return v.ncomps() == 1; });
    auto vector = std::find_if(pointData_.begin(), pointData_.end(), [](auto const& v) { return v.ncomps() > 1; });
    out << "<PointData" << (scalar != pointData_.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
                        << (vector != pointData_.end() ? " Vectors=\"" + vector->name() + "\"" : "")
                        << ">\n";
    for (auto const& v : pointData_)
      this->writeData(out, offsets, v, Super::POINT_DATA);
    out << "</PointData>\n";
  }

  { // Write data associated with grid cells
    auto scalar = std::find_if(cellData_.begin(), cellData_.end(), [](auto const& v) { return v.ncomps() == 1; });
    auto vector = std::find_if(cellData_.begin(), cellData_.end(), [](auto const& v) { return v.ncomps() > 1; });
    out << "<CellData" << (scalar != cellData_.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
                       << (vector != cellData_.end() ? " Vectors=\"" + vector->name() + "\"" : "")
                       << ">\n";
    for (auto const& v : cellData_)
      this->writeData(out, offsets, v, Super::CELL_DATA);
    out << "</CellData>\n";
  }
  out << "</Piece>\n";
  out << "</ImageData>\n";

  std::vector<std::uint64_t> blocks; // size of i'th appended block
  pos_type appended_pos = 0;
  if (is_a(format_, Vtk::APPENDED)) {
    out << "<AppendedData encoding=\"raw\">\n_";
    appended_pos = out.tellp();
    for (auto const& v : pointData_) {
      if (datatype_ == Vtk::FLOAT32)
        blocks.push_back( this->template writeDataAppended<float>(out, v, Super::POINT_DATA) );
      else
        blocks.push_back( this->template writeDataAppended<double>(out, v, Super::POINT_DATA) );
    }
    for (auto const& v : cellData_) {
      if (datatype_ == Vtk::FLOAT32)
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
  ::writeParallelFile (std::string const& pfilename, int size) const
{
  std::string filename = pfilename + ".p" + this->fileExtension();
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);

  out << "<VTKFile type=\"PImageData\" version=\"1.0\" "
      << "byte_order=\"" << this->getEndian() << "\" header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<PImageData GhostLevel=\"0\" WholeExtent=\"";
  auto const& wholeExtent = dataCollector_.wholeExtent();
  for (int i = 0; i < 3; ++i) {
    out << (i == 0 ? "" : " ") << wholeExtent[2*i] << " " << wholeExtent[2*i+1];
  }
  out << "\" Origin=\"";
  auto const& origin = dataCollector_.origin();
  for (int i = 0; i < 3; ++i) {
    out << (i == 0 ? "" : " ") << origin[i];
  }
  out << "\" Spacing=\"";
  auto const& spacing = dataCollector_.spacing();
  for (int i = 0; i < 3; ++i) {
    out << (i == 0 ? "" : " ") << spacing[i];
  }
  out << "\">\n";

  { // Write data associated with grid points
    auto scalar = std::find_if(pointData_.begin(), pointData_.end(), [](auto const& v) { return v.ncomps() == 1; });
    auto vector = std::find_if(pointData_.begin(), pointData_.end(), [](auto const& v) { return v.ncomps() > 1; });
    out << "<PPointData" << (scalar != pointData_.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
                         << (vector != pointData_.end() ? " Vectors=\"" + vector->name() + "\"" : "")
                         << ">\n";
    for (auto const& v : pointData_) {
      out << "<PDataArray Name=\"" << v.name() << "\" type=\"" << Vtk::Map::from_datatype[datatype_] << "\""
          << " NumberOfComponents=\"" << v.ncomps() << "\" />\n";
    }
    out << "</PPointData>\n";
  }

  { // Write data associated with grid cells
    auto scalar = std::find_if(cellData_.begin(), cellData_.end(), [](auto const& v) { return v.ncomps() == 1; });
    auto vector = std::find_if(cellData_.begin(), cellData_.end(), [](auto const& v) { return v.ncomps() > 1; });
    out << "<PCellData" << (scalar != cellData_.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
                        << (vector != cellData_.end() ? " Vectors=\"" + vector->name() + "\"" : "")
                        << ">\n";
    for (auto const& v : cellData_) {
      out << "<PDataArray Name=\"" << v.name() << "\" type=\"" << Vtk::Map::from_datatype[datatype_] << "\""
          << " NumberOfComponents=\"" << v.ncomps() << "\" />\n";
    }
    out << "</PCellData>\n";
  }

  // Write piece file references
  dataCollector_.writePieces([&out,pfilename,ext=this->fileExtension()](int p, std::array<int,6> const& extent, bool write_extent)
  {
    out << "<Piece Source=\"" << pfilename << "_p" << std::to_string(p) << "." << ext << "\"";
    if (write_extent) {
      out << " Extent=\"";
      for (int i = 0; i < 3; ++i) {
        out << (i == 0 ? "" : " ") << extent[2*i] << " " << extent[2*i+1];
      }
      out << "\"";
     }
     out << " />\n";
  });

  out << "</PImageData>\n";
  out << "</VTKFile>";
}

}} // end namespace Dune::experimental
