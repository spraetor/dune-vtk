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

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::write (std::string const& fn, Vtk::FormatTypes format, Vtk::DataTypes datatype)
{
  format_ = format;
  datatype_ = datatype;

#ifndef HAVE_ZLIB
  if (format_ == Vtk::COMPRESSED) {
    std::cout << "Dune is compiled without compression. Falling back to BINARY VTK output!\n";
    format_ = Vtk::BINARY;
  }
#endif

  dataCollector_.update();

  auto p = filesystem::path(fn);
  auto name = p.stem();
  p.remove_filename();
  p /= name.string();

  std::string filename = p.string() + "." + fileExtension();

#ifdef HAVE_MPI
  int rank = -1;
  int num_ranks = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  if (num_ranks > 1) {
    filename = p.string() + "_p" + std::to_string(rank) + "." + fileExtension();

    writeSerialFile(filename);
    if (rank == 0) {
      writeParallelFile(p.string(), num_ranks);
    }
  } else {
    writeSerialFile(filename);
  }
#endif
}


template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::writeData (std::ofstream& out, std::vector<pos_type>& offsets,
               VtkFunction const& fct, PositionTypes type) const
{
  out << "<DataArray Name=\"" << fct.name() << "\" type=\"" << to_string(fct.type()) << "\""
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
void VtkWriterInterface<GV,DC>
  ::writePoints (std::ofstream& out, std::vector<pos_type>& offsets) const
{
  out << "<DataArray type=\"" << to_string(datatype_) << "\""
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
  template <class T>
std::uint64_t VtkWriterInterface<GV,DC>
  ::writeDataAppended (std::ofstream& out, VtkFunction const& fct, PositionTypes type) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  if (type == POINT_DATA) {
    auto data = dataCollector_.template pointData<T>(fct);
    return this->writeAppended(out, data);
  } else {
    auto data = dataCollector_.template cellData<T>(fct);
    return this->writeAppended(out, data);
  }
}


template <class GV, class DC>
  template <class T>
std::uint64_t VtkWriterInterface<GV,DC>
  ::writePointsAppended (std::ofstream& out) const
{
  assert(is_a(format_, Vtk::APPENDED) && "Function should by called only in appended mode!\n");

  auto points = dataCollector_.template points<T>();
  return this->writeAppended(out, points);
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
std::uint64_t VtkWriterInterface<GV,DC>
  ::writeAppended (std::ofstream& out, std::vector<T> const& values) const
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

} // end namespace Dune
