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
void VtkUnstructuredGridWriter<GV,DC>
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
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
      << "byte_order=\"" << this->getEndian() << "\" header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<UnstructuredGrid>\n";
  out << "<Piece NumberOfPoints=\"" << dataCollector_.numPoints() << "\" "
      << "NumberOfCells=\"" << dataCollector_.numCells() << "\">\n";

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

  // Write point coordinates
  out << "<Points>\n";
  this->writePoints(out, offsets);
  out << "</Points>\n";

  // Write element connectivity, types and offsets
  out << "<Cells>\n";
  this->writeCells(out, offsets);
  out << "</Cells>\n";

  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";

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
    if (datatype_ == Vtk::FLOAT32)
      blocks.push_back( this->template writePointsAppended<float>(out) );
    else
      blocks.push_back( this->template writePointsAppended<double>(out) );
    auto bs = this->writeCellsAppended(out);
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

  out << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" "
      << "byte_order=\"" << this->getEndian() << "\" header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<PUnstructuredGrid GhostLevel=\"0\">\n";

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

  // Write points
  out << "<PPoints>\n";
  out << "<PDataArray type=\"" << Vtk::Map::from_datatype[datatype_] << "\" NumberOfComponents=\"3\" />\n";
  out << "</PPoints>\n";

  // Write piece file references
  for (int i = 0; i < size; ++i) {
    out << "<Piece Source=\"" << pfilename << "_p" << std::to_string(i) << "." << this->fileExtension() << "\" />\n";
  }

  out << "</PUnstructuredGrid>\n";
  out << "</VTKFile>";
}

}} // end namespace Dune::experimental
