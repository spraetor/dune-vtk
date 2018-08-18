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

namespace Dune {

template <class GridView>
void VtkWriter<GridView>::write(std::string const& fn, Vtk::FormatTypes format, Vtk::DataTypes datatype)
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


template <class GridView>
void VtkWriter<GridView>::writeImpl(std::string const& filename) const
{
  std::ofstream out(filename, std::ios_base::ate | std::ios::binary);
  if (format_ == Vtk::ASCII) {
    if (datatype_ == Vtk::FLOAT32)
      out << std::setprecision(std::numeric_limits<float>::digits10+1);
    else
      out << std::setprecision(std::numeric_limits<double>::digits10+1);
  }

  std::vector<pos_type> offsets; // pos => offset
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
      << "byte_order=\"" << getEndian() << "\" header_type=\"UInt64\""
      << (format_ == Vtk::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<UnstructuredGrid>\n";
  out << "<Piece NumberOfPoints=\"" << gridView_.size(dimension) << "\" "
      << "NumberOfCells=\"" << gridView_.size(0) << "\">\n";

  // TODO: detect vector and scalars by ncomps() property of data
  out << "<PointData" << (!vertexData_.empty() ? " Scalars=\"" + vertexData_.front().name() + "\"" : "") << ">\n";
  for (auto const& v : vertexData_)
    writeData(out, offsets, v, Vtk::VERTEX_DATA);
  out << "</PointData>\n";

  // TODO: detect vector and scalars by ncomps() property of data
  out << "<CellData" << (!cellData_.empty() ? " Scalars=\"" + cellData_.front().name() + "\"" : "") << ">\n";
  for (auto const& v : cellData_)
    writeData(out, offsets, v, Vtk::CELL_DATA);
  out << "</CellData>\n";

  out << "<Points>\n";
  writePoints(out, offsets);
  out << "</Points>\n";

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
    for (auto const& v : vertexData_) {
      if (datatype_ == Vtk::FLOAT32)
        blocks.push_back( writeDataAppended<float>(out, v, Vtk::VERTEX_DATA) );
      else
        blocks.push_back( writeDataAppended<double>(out, v, Vtk::VERTEX_DATA) );
    }
    for (auto const& v : cellData_) {
      if (datatype_ == Vtk::FLOAT32)
        blocks.push_back( writeDataAppended<float>(out, v, Vtk::CELL_DATA) );
      else
        blocks.push_back( writeDataAppended<double>(out, v, Vtk::CELL_DATA) );
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


// @{ implementation details

template <class T, class GridView>
std::vector<T> getPoints(GridView const& gridView)
{
  const int dim = GridView::dimension;
  std::vector<T> data(gridView.size(dim) * 3);
  auto const& indexSet = gridView.indexSet();
  for (auto const& vertex : vertices(gridView)) {
    std::size_t idx = 3 * indexSet.index(vertex);
    auto v = vertex.geometry().center();
    for (std::size_t j = 0; j < v.size(); ++j)
      data[idx + j] = T(v[j]);
    for (std::size_t j = v.size(); j < 3u; ++j)
      data[idx + j] = T(0);
  }
  return data;
}

template <class T, class GridView, class GlobalFunction>
std::vector<T> getVertexData(GridView const& gridView, GlobalFunction const& fct)
{
  const int dim = GridView::dimension;
  std::vector<T> data(gridView.size(dim) * fct.ncomps());
  auto const& indexSet = gridView.indexSet();
  auto localFct = localFunction(fct);
  for (auto const& e : elements(gridView)) {
    localFct.bind(e);
    Vtk::CellType cellType{e.type()};
    auto refElem = referenceElement(e.geometry());
    for (int j = 0; j < e.subEntities(dim); ++j) {
      std::size_t idx = fct.ncomps() * indexSet.subIndex(e,cellType.localIndex(j),dim);
      for (int comp = 0; comp < fct.ncomps(); ++comp)
        data[idx + comp] = T(localFct.evaluate(comp, refElem.position(cellType.localIndex(j),dim)));
    }
    localFct.unbind();
  }
  return data;
}

template <class T, class GridView, class GlobalFunction>
std::vector<T> getCellData(GridView const& gridView, GlobalFunction const& fct)
{
  const int dim = GridView::dimension;
  std::vector<T> data(gridView.size(0) * fct.ncomps());
  auto const& indexSet = gridView.indexSet();
  auto localFct = localFunction(fct);
  for (auto const& e : elements(gridView)) {
    localFct.bind(e);
    auto refElem = referenceElement(e.geometry());
    std::size_t idx = fct.ncomps() * indexSet.index(e);
    for (int comp = 0; comp < fct.ncomps(); ++comp)
      data[idx + comp] = T(localFct.evaluate(comp, refElem.position(0,0)));
    localFct.unbind();
  }
  return data;
}

// @}


template <class GridView>
void VtkWriter<GridView>::writeData(std::ofstream& out, std::vector<pos_type>& offsets,
                          GlobalFunction const& fct, Vtk::PositionTypes type) const
{
  out << "<DataArray Name=\"" << fct.name() << "\" type=\"" << (datatype_ == Vtk::FLOAT32 ? "Float32" : "Float64") << "\""
      << " NumberOfComponents=\"" << fct.ncomps() << "\" format=\"" << (format_ == Vtk::ASCII ? "ascii\">\n" : "appended\"");

  if (format_ == Vtk::ASCII) {
    std::size_t i = 0;
    if (type == Vtk::VERTEX_DATA) {
      auto data = getVertexData<double>(gridView_, fct);
      for (auto const& v : data)
        out << v << (++i % 6 != 0 ? ' ' : '\n');
    } else {
      auto data = getCellData<double>(gridView_, fct);
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


template <class GridView>
void VtkWriter<GridView>::writePoints(std::ofstream& out, std::vector<pos_type>& offsets) const
{
  out << "<DataArray type=\"" << (datatype_ == Vtk::FLOAT32 ? "Float32" : "Float64") << "\""
      << " NumberOfComponents=\"3\" format=\"" << (format_ == Vtk::ASCII ? "ascii\">\n" : "appended\"");

  if (format_ == Vtk::ASCII) {
    auto points = getPoints<double>(gridView_);
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


template <class GridView>
void VtkWriter<GridView>::writeCells(std::ofstream& out, std::vector<pos_type>& offsets) const
{
  auto const& indexSet = gridView_.indexSet();
  if (format_ == Vtk::ASCII) {
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    std::size_t i = 0;
    for (auto const& c : elements(gridView_)) {
      auto cellType = getType(c.type());
      for (int j = 0; j < c.subEntities(dimension); ++j)
        out << indexSet.subIndex(c,cellType.localIndex(j),dimension) << (++i % 6 != 0 ? ' ' : '\n');
    }
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    i = 0;
    std::size_t old_o = 0;
    for (auto const& c : elements(gridView_))
      out << (old_o += c.subEntities(dimension)) << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    i = 0;
    for (auto const& c : elements(gridView_))
      out << int(getType(c.type()).type()) << (++i % 6 != 0 ? ' ' : '\n');
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


// @{ implementation details

template <class T>
std::uint64_t write_values_to_buffer(std::size_t max_num_values, unsigned char* buffer,
                                     std::vector<T> const& vec, std::size_t shift)
{
  std::size_t num_values = std::min(max_num_values, vec.size()-shift);
  std::uint64_t bs = num_values*sizeof(T);
  std::memcpy(buffer, (unsigned char*)(vec.data()+shift), std::size_t(bs));
  return bs;
}


template <class OStream>
std::uint64_t write_compressed(unsigned char const* buffer, unsigned char* buffer_out,
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
  std::cerr << "Can not call write_compressed without compression enabled!\n";
  std::abort();
  return 0;
#endif
}

// @}

template <class GridView>
  template <class T>
std::uint64_t VtkWriter<GridView>::writeAppended(std::ofstream& out, std::vector<T> const& values) const
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
    std::uint64_t bs = write_values_to_buffer<T>(num_values, buffer.data(), values, i*num_values);

    if (format_ == Vtk::COMPRESSED) {
      buffer_out.resize(std::size_t(compressed_block_size));
      cbs[i] = write_compressed(buffer.data(), buffer_out.data(), bs,
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


template <class GridView>
  template <class T>
std::uint64_t VtkWriter<GridView>::writeDataAppended(std::ofstream& out, GlobalFunction const& fct, Vtk::PositionTypes type) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  if (type == Vtk::VERTEX_DATA) {
    auto data = getVertexData<T>(gridView_, fct);
    return writeAppended(out, data);
  } else {
    auto data = getCellData<T>(gridView_, fct);
    return writeAppended(out, data);
  }
}


template <class GridView>
  template <class T>
std::uint64_t VtkWriter<GridView>::writePointsAppended(std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto points = getPoints<T>(gridView_);
  return writeAppended(out, points);
}


template <class GridView>
std::array<std::uint64_t,3> VtkWriter<GridView>::writeCellsAppended(std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto const& indexSet = gridView_.indexSet();
  auto types = indexSet.types(0);
  int maxVertices = std::accumulate(types.begin(), types.end(), 1, [](int m, GeometryType t) {
    auto refElem = referenceElement<double,dimension>(t);
    return std::max(m, refElem.size(dimension));
  });

  std::vector<std::int64_t> connectivity;
  std::vector<std::int64_t> cell_offsets;
  std::vector<std::uint8_t> cell_types;

  connectivity.reserve(gridView_.size(0) * maxVertices);
  cell_offsets.reserve(gridView_.size(0));
  cell_types.reserve(gridView_.size(0));

  std::int64_t old_o = 0;
  for (auto const& c : elements(gridView_)) {
    auto cellType = getType(c.type());
    for (int j = 0; j < c.subEntities(dimension); ++j)
      connectivity.push_back( std::int64_t(indexSet.subIndex(c,cellType.localIndex(j),dimension)) );
    cell_offsets.push_back(old_o += c.subEntities(dimension));
    cell_types.push_back(cellType.type());
  }

  // write conncetivity
  std::uint64_t bs0 = writeAppended(out, connectivity);

  // write offsets
  std::uint64_t bs1 = writeAppended(out, cell_offsets);

  // write cell types
  std::uint64_t bs2 = writeAppended(out, cell_types);

  return {bs0, bs1, bs2};
}

} // end namespace Dec
