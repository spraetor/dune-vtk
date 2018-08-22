#pragma once

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

#include "utility/enum.hh"
#include "utility/filesystem.hh"
#include "utility/string.hh"

namespace Dune { namespace experimental {

template <class GV, class DC>
void VtkWriter<GV,DC>::write (std::string const& fn, Vtk::FormatTypes format, Vtk::DataTypes datatype)
{
  format_ = format;
  datatype_ = datatype;

#ifndef HAVE_ZLIB
  if (format_ == Vtk::COMPRESSED) {
    std::cout << "Dune is compiled without compression. Falling back to BINARY VTK output!\n";
    format_ = Vtk::BINARY;
  }
#endif

  std::string filename = fn;

#ifdef DUNE_HAS_MPI
  int rank = -1;
  int num_ranks = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  if (num_ranks > 1) {
    auto p = filesystem::path(fn);
    auto name = p.stem();
    p.remove_filename();
    p /= name.string() + "_p" + std::to_string(rank) + ".vtu";
    filename = p.string();
  }
#endif
  writeImpl(filename);
}


template <class GV, class DC>
void VtkWriter<GV,DC>::writeImpl (std::string const& filename) const
{
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  if (format_ == Vtk::ASCII) {
    if (datatype_ == Vtk::FLOAT32)
      out << std::setprecision(std::numeric_limits<float>::digits10+2);
    else
      out << std::setprecision(std::numeric_limits<double>::digits10+2);
  }

  dataCollector_.update();

  std::vector<pos_type> offsets; // pos => offset
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
      << "byte_order=\"" << getEndian() << "\" header_type=\"UInt64\""
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
      writeData(out, offsets, v, POINT_DATA);
    out << "</PointData>\n";
  }

  { // Write data associated with grid cells
    auto scalar = std::find_if(cellData_.begin(), cellData_.end(), [](auto const& v) { return v.ncomps() == 1; });
    auto vector = std::find_if(cellData_.begin(), cellData_.end(), [](auto const& v) { return v.ncomps() > 1; });
    out << "<CellData" << (scalar != cellData_.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
                       << (vector != cellData_.end() ? " Vectors=\"" + vector->name() + "\"" : "")
                       << ">\n";
    for (auto const& v : cellData_)
      writeData(out, offsets, v, CELL_DATA);
    out << "</CellData>\n";
  }

  // Write point coordinates
  out << "<Points>\n";
  writePoints(out, offsets);
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
      if (datatype_ == Vtk::FLOAT32)
        blocks.push_back( writeDataAppended<float>(out, v, POINT_DATA) );
      else
        blocks.push_back( writeDataAppended<double>(out, v, POINT_DATA) );
    }
    for (auto const& v : cellData_) {
      if (datatype_ == Vtk::FLOAT32)
        blocks.push_back( writeDataAppended<float>(out, v, CELL_DATA) );
      else
        blocks.push_back( writeDataAppended<double>(out, v, CELL_DATA) );
    }
    if (datatype_ == Vtk::FLOAT32)
      blocks.push_back( writePointsAppended<float>(out) );
    else
      blocks.push_back( writePointsAppended<double>(out) );
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
void VtkWriter<GV,DC>::writeData (std::ofstream& out, std::vector<pos_type>& offsets,
                                     GlobalFunction const& fct, PositionTypes type) const
{
  out << "<DataArray Name=\"" << fct.name() << "\" type=\"" << Vtk::Map::from_datatype[datatype_] << "\""
      << " NumberOfComponents=\"" << fct.ncomps() << "\" format=\"" << (format_ == Vtk::ASCII ? "ascii\">\n" : "appended\"");

  if (format_ == Vtk::ASCII) {
    std::size_t i = 0;
    if (type == POINT_DATA) {
      auto data = dataCollector_.template pointData<double>(fct);
      for (auto const& v : data)
        out << v << (++i % 6 != 0 ? ' ' : '\n');
    } else {
      auto data = dataCollector_.template cellData<double>(fct);
      for (auto const& v : data)
        out << v << (++i % 6 != 0 ? ' ' : '\n');
    }
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkWriter<GV,DC>::writePoints (std::ofstream& out, std::vector<pos_type>& offsets) const
{
  out << "<DataArray type=\"" << Vtk::Map::from_datatype[datatype_] << "\""
      << " NumberOfComponents=\"3\" format=\"" << (format_ == Vtk::ASCII ? "ascii\">\n" : "appended\"");

  if (format_ == Vtk::ASCII) {
    auto points = dataCollector_.template points<double>();
    std::size_t i = 0;
    for (auto const& v : points)
      out << v << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkWriter<GV,DC>::writeCells (std::ofstream& out, std::vector<pos_type>& offsets) const
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


namespace Impl {

template <class T>
std::uint64_t writeValuesToBuffer (std::size_t max_num_values, unsigned char* buffer,
                                   std::vector<T> const& vec, std::size_t shift)
{
  std::size_t num_values = std::min(max_num_values, vec.size()-shift);
  std::uint64_t bs = num_values*sizeof(T);
  std::memcpy(buffer, (unsigned char*)(vec.data()+shift), std::size_t(bs));
  return bs;
}


template <class OStream>
std::uint64_t writeCompressed (unsigned char const* buffer, unsigned char* buffer_out,
                               std::uint64_t bs, std::uint64_t cbs, int level, OStream& outb)
{
#ifdef HAVE_ZLIB
  uLongf uncompressed_space = uLongf(bs);
  uLongf compressed_space = uLongf(cbs);

  Bytef* out = reinterpret_cast<Bytef*>(buffer_out);
  Bytef const* in = reinterpret_cast<Bytef const*>(buffer);

  if (compress2(out, &compressed_space, in, uncompressed_space, level) != Z_OK) {
    std::cerr << "Zlib error while compressing data.\n";
    std::abort();
  } else {
    outb.write((char*)out, compressed_space);
  }

  return compressed_space;
#else
  std::cerr << "Can not call writeCompressed without compression enabled!\n";
  std::abort();
  return 0;
#endif
}

} // end namespace Impl


template <class GV, class DC>
  template <class T>
std::uint64_t VtkWriter<GV,DC>::writeAppended (std::ofstream& out, std::vector<T> const& values) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");
  pos_type begin_pos = out.tellp();

  std::uint64_t size = values.size() * sizeof(T);

  std::uint64_t num_full_blocks = size / block_size;
  std::uint64_t last_block_size = size % block_size;
  std::uint64_t num_blocks = num_full_blocks + (last_block_size > 0 ? 1 : 0);

  // write block-size(s)
  std::uint64_t zero = 0;
  if (format_ == Vtk::COMPRESSED) {
    out.write((char*)&num_blocks, sizeof(std::uint64_t));
    out.write((char*)&block_size, sizeof(std::uint64_t));
    out.write((char*)&last_block_size, sizeof(std::uint64_t));
    for (std::uint64_t i = 0; i < num_blocks; ++i)
      out.write((char*)&zero, sizeof(std::uint64_t));
  } else {
    out.write((char*)&size, sizeof(std::uint64_t));
  }

  std::uint64_t compressed_block_size = block_size + (block_size + 999)/1000 + 12;
  std::vector<unsigned char> buffer(block_size);
  std::vector<unsigned char> buffer_out;
  std::size_t num_values = block_size / sizeof(T);

  std::vector<std::uint64_t> cbs(std::size_t(num_blocks), 0); // compressed block sizes
  for (std::size_t i = 0; i < std::size_t(num_blocks); ++i) {
    std::uint64_t bs = Impl::writeValuesToBuffer<T>(num_values, buffer.data(), values, i*num_values);

    if (format_ == Vtk::COMPRESSED) {
      buffer_out.resize(std::size_t(compressed_block_size));
      cbs[i] = Impl::writeCompressed(buffer.data(), buffer_out.data(), bs,
                                     compressed_block_size, compression_level, out);
    } else
      out.write((char*)buffer.data(), bs);
  }

  pos_type end_pos = out.tellp();
  if (format_ == Vtk::COMPRESSED) {
    out.seekp(begin_pos + std::int64_t(3*sizeof(std::uint64_t)));
    out.write((char*)cbs.data(), num_blocks*sizeof(std::uint64_t));
    out.seekp(end_pos);
  }

  return std::uint64_t(end_pos - begin_pos);
}


template <class GV, class DC>
  template <class T>
std::uint64_t VtkWriter<GV,DC>::writeDataAppended (std::ofstream& out, GlobalFunction const& fct, PositionTypes type) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  if (type == POINT_DATA) {
    auto data = dataCollector_.template pointData<T>(fct);
    return writeAppended(out, data);
  } else {
    auto data = dataCollector_.template cellData<T>(fct);
    return writeAppended(out, data);
  }
}


template <class GV, class DC>
  template <class T>
std::uint64_t VtkWriter<GV,DC>::writePointsAppended (std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto points = dataCollector_.template points<T>();
  return writeAppended(out, points);
}


template <class GV, class DC>
std::array<std::uint64_t,3> VtkWriter<GV,DC>::writeCellsAppended (std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto cells = dataCollector_.cells();

  // write conncetivity, offsets, and types
  std::uint64_t bs0 = writeAppended(out, cells.connectivity);
  std::uint64_t bs1 = writeAppended(out, cells.offsets);
  std::uint64_t bs2 = writeAppended(out, cells.types);

  return {bs0, bs1, bs2};
}

}} // end namespace Dune::experimental
