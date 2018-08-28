#include <sstream>
#include <fstream>
#include <iterator>
#include <string>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#include <dune/common/classname.hh>

#include "utility/filesystem.hh"
#include "utility/string.hh"

namespace Dune {

template <class Grid, class Creator>
void VtkReader<Grid,Creator>::readFromFile (std::string const& filename)
{
  // check whether file exists!
  if (!filesystem::exists(filename))
    DUNE_THROW(IOError, "File " << filename << " does not exist!");

  std::ifstream input(filename, std::ios_base::in | std::ios_base::binary);

  std::string compressor = "";
  std::string data_name = "", data_format = "";
  Vtk::DataTypes data_type = Vtk::UNKNOWN;
  std::size_t data_components = 0;
  std::uint64_t data_offset = 0;

  Sections section = NO_SECTION;
  for (std::string line; std::getline(input, line); ) {
    ltrim(line);

    if (isSection(line, "VTKFile", section)) {
      bool closed = false;
      auto attr = parseXml(line, closed);

      if (!attr["type"].empty())
        assert(attr["type"] == "UnstructuredGrid");
      if (!attr["version"].empty())
        assert(std::stod(attr["version"]) == 1.0);
      if (!attr["byte_order"].empty())
        assert(attr["byte_order"] == "LittleEndian");
      if (!attr["header_type"].empty())
        assert(attr["header_type"] == "UInt64");
      if (!attr["compressor"].empty()) {
        compressor = attr["compressor"];
        assert(compressor == "vtkZLibDataCompressor"); // only ZLib compression supported
      }

      // std::cout << "<VTKFile type='" << attr["type"] << "' version='" << attr["version"] << "' header_type='" << attr["header_type"] << "' byte_order='" << attr["byte_order"] << "' compressor='" << attr["compressor"] << "'>\n";

      section = VTK_FILE;
    }
    else if (isSection(line, "/VTKFile", section, VTK_FILE))
      section = NO_SECTION;
    else if (isSection(line, "UnstructuredGrid", section, VTK_FILE))
      section = UNSTRUCTURED_GRID;
    else if (isSection(line, "/UnstructuredGrid", section, UNSTRUCTURED_GRID))
      section = VTK_FILE;
    else if (isSection(line, "Piece", section, UNSTRUCTURED_GRID)) {
      bool closed = false;
      auto attr = parseXml(line, closed);

      assert(attr.count("NumberOfPoints") > 0 && attr.count("NumberOfCells") > 0);
      numberOfPoints_ = std::stoul(attr["NumberOfPoints"]);
      numberOfCells_ = std::stoul(attr["NumberOfCells"]);
      section = PIECE;
    }
    else if (isSection(line, "/Piece", section, PIECE))
      section = UNSTRUCTURED_GRID;
    else if (isSection(line, "PointData", section, PIECE))
      section = POINT_DATA;
    else if (isSection(line, "/PointData", section, POINT_DATA))
      section = PIECE;
    else if (isSection(line, "CellData", section, PIECE))
      section = CELL_DATA;
    else if (isSection(line, "/CellData", section, CELL_DATA))
      section = PIECE;
    else if (isSection(line, "Points", section, PIECE))
      section = POINTS;
    else if (isSection(line, "/Points", section, POINTS))
      section = PIECE;
    else if (isSection(line, "Cells", section, PIECE))
      section = CELLS;
    else if (isSection(line, "/Cells", section, CELLS))
      section = PIECE;
    else if (line.substr(1,9) == "DataArray") {
      bool closed = false;
      auto attr = parseXml(line, closed);

      data_type = Vtk::Map::to_datatype[attr["type"]];

      if (!attr["Name"].empty())
        data_name = to_lower(attr["Name"]);
      else if (section == POINTS)
        data_name = "points";

      data_components = 1;
      if (!attr["NumberOfComponents"].empty())
        data_components = std::stoul(attr["NumberOfComponents"]);

      // determine FormatType
      data_format = to_lower(attr["format"]);
      if (data_format == "appended") {
        format_ = !compressor.empty() ? Vtk::COMPRESSED : Vtk::BINARY;
      } else {
        format_ = Vtk::ASCII;
      }

      // Offset makes sense in appended mode only
      data_offset = 0;
      if (!attr["offset"].empty()) {
        data_offset = std::stoul(attr["offset"]);
        assert(data_format == "appended");
      }

      // Store attributes of DataArray
      dataArray_[data_name] = {data_type, data_components, data_offset};

      // std::cout << "<DataArray type='" << attr["type"] << "' Name='" << data_name << "' NumberOfComponents='" << attr["NumberOfComponents"] << "' format='" << data_format << "' offset='" << attr["offset"] << "' " << (closed ? "/" : "") << ">\n";

      // Skip section in appended mode
      if (data_format == "appended") {
        if (!closed) {
          while (std::getline(input, line)) {
            ltrim(line);
            if (line.substr(1,10) == "/DataArray")
              break;
          }
        }
        continue;
      }

      if (section == POINT_DATA)
        section = PD_DATA_ARRAY;
      else if (section == POINTS)
        section = POINTS_DATA_ARRAY;
      else if (section == CELL_DATA)
        section = CD_DATA_ARRAY;
      else if (section == CELLS)
        section = CELLS_DATA_ARRAY;
      else
        DUNE_THROW(Exception, "Wrong section for <DataArray>");
    }
    else if (line.substr(1,10) == "/DataArray") {
      if (section == PD_DATA_ARRAY)
        section = POINT_DATA;
      else if (section == POINTS_DATA_ARRAY)
        section = POINTS;
      else if (section == CD_DATA_ARRAY)
        section = CELL_DATA;
      else if (section == CELLS_DATA_ARRAY)
        section = CELLS;
      else
        DUNE_THROW(Exception, "Wrong section for </DataArray>");
    }
    else if (isSection(line, "AppendedData", section, VTK_FILE)) {
      bool closed = false;
      auto attr = parseXml(line, closed);
      if (!attr["encoding"].empty())
        assert(attr["encoding"] == "raw"); // base64 encoding not supported

      offset0_ = findAppendedDataPosition(input);
      if (dataArray_["points"].type == Vtk::FLOAT32)
        readPointsAppended<float>(input);
      else
        readPointsAppended<double>(input);

      readCellsAppended(input);
      section = NO_SECTION; // finish reading after appended section
    }
    else if (isSection(line, "/AppendedData", section, APPENDED_DATA))
      section = VTK_FILE;

    switch (section) {
      case PD_DATA_ARRAY:
        if (data_type == Vtk::FLOAT32) {
          std::vector<float> values;
          section = readPointData(input, values, data_name);
        } else if (data_type == Vtk::FLOAT64) {
          std::vector<double> values;
          section = readPointData(input, values, data_name);
        }
        break;
      case POINTS_DATA_ARRAY:
        section = readPoints(input);
        break;
      case CD_DATA_ARRAY:
        if (data_type == Vtk::FLOAT32) {
          std::vector<float> values;
          section = readCellData(input, values, data_name);
        } else if (data_type == Vtk::FLOAT64) {
          std::vector<double> values;
          section = readCellData(input, values, data_name);
        }
        break;
      case CELLS_DATA_ARRAY:
        section = readCells(input, data_name);
        break;
      default:
        // do nothing
        break;
    }

    if (section == NO_SECTION)
      break;
  }

  if (section != NO_SECTION)
    DUNE_THROW(IOError, "VTK-File is incomplete. It must end with </VTKFile>!");

  createGrid();
}


// @{ implementation detail
/**
 * Read ASCII data from `input` stream into vector `values`
 * \param max_size  Upper bound for the number of values
 * \param section   Current XML section you are reading in
 * \param parent_section   XML Section to return when current `section` is finished.
 **/
template <class IStream, class T, class Sections>
Sections readDataArray (IStream& input, std::vector<T>& values, std::size_t max_size,
                        Sections section, Sections parent_section)
{
  values.reserve(max_size);
  using S = std::conditional_t<(sizeof(T) <= 1), std::uint16_t, T>; // problem when reading chars as ints

  std::size_t idx = 0;
  for (std::string line; std::getline(input, line);) {
    trim(line);
    if (line.substr(1,10) == "/DataArray")
      return parent_section;

    std::istringstream stream(line);
    S value;
    for (; stream >> value; idx++)
      values.push_back(T(value));
  }

  return section;
}
// @}


template <class Grid, class Creator>
typename VtkReader<Grid,Creator>::Sections
VtkReader<Grid,Creator>::readPoints (std::ifstream& input)
{
  using T = typename GlobalCoordinate::value_type;
  assert(numberOfPoints_ > 0);
  assert(dataArray_["points"].components == 3u);

  std::vector<T> point_values;
  auto sec = readDataArray(input, point_values, 3*numberOfPoints_, POINTS_DATA_ARRAY, POINTS);
  assert(sec == POINTS);
  assert(point_values.size() == 3*numberOfPoints_);

  // extract points from continuous values
  GlobalCoordinate p;
  vec_points.reserve(numberOfPoints_);
  std::size_t idx = 0;
  for (std::size_t i = 0; i < numberOfPoints_; ++i) {
    for (std::size_t j = 0; j < p.size(); ++j)
      p[j] = point_values[idx++];
    idx += (3u - p.size());
    vec_points.push_back(p);
  }

  return sec;
}


template <class Grid, class Creator>
  template <class T>
void VtkReader<Grid,Creator>::readPointsAppended (std::ifstream& input)
{
  assert(numberOfPoints_ > 0);
  assert(dataArray_["points"].components == 3u);

  std::vector<T> point_values;
  readAppended(input, point_values, dataArray_["points"].offset);
  assert(point_values.size() == 3*numberOfPoints_);

  // extract points from continuous values
  GlobalCoordinate p;
  vec_points.reserve(numberOfPoints_);
  std::size_t idx = 0;
  for (std::size_t i = 0; i < numberOfPoints_; ++i) {
    for (std::size_t j = 0; j < p.size(); ++j)
      p[j] = T(point_values[idx++]);
    idx += (3u - p.size());
    vec_points.push_back(p);
  }
}


template <class Grid, class Creator>
typename VtkReader<Grid,Creator>::Sections
VtkReader<Grid,Creator>::readCells (std::ifstream& input, std::string name)
{
  Sections sec = CELLS_DATA_ARRAY;

  assert(numberOfCells_ > 0);
  if (name == "types") {
    sec = readDataArray(input, vec_types, numberOfCells_, CELLS_DATA_ARRAY, CELLS);
    assert(vec_types.size() == numberOfCells_);
  } else if (name == "offsets") {
    sec = readDataArray(input, vec_offsets, numberOfCells_, CELLS_DATA_ARRAY, CELLS);
    assert(vec_offsets.size() == numberOfCells_);
  } else if (name == "connectivity") {
    std::size_t max_size = 0;
    int max_vertices = (Grid::dimension == 1 ? 2 : Grid::dimension == 2 ? 4 : 8);
    if (!vec_offsets.empty())
      max_size = vec_offsets.back();
    else
      max_size = numberOfCells_ * max_vertices;
    sec = readDataArray(input, vec_connectivity, max_size, CELLS_DATA_ARRAY, CELLS);
  }

  return sec;
}


template <class Grid, class Creator>
void VtkReader<Grid,Creator>::readCellsAppended (std::ifstream& input)
{
  assert(numberOfCells_ > 0);
  auto types_data = dataArray_["types"];
  auto dataArray_data = dataArray_["offsets"];
  auto connectivity_data = dataArray_["connectivity"];

  assert(types_data.type == Vtk::UINT8);
  readAppended(input, vec_types, types_data.offset);
  assert(vec_types.size() == numberOfCells_);

  assert(dataArray_data.type == Vtk::INT64);
  readAppended(input, vec_offsets, dataArray_data.offset);
  assert(vec_offsets.size() == numberOfCells_);

  assert(connectivity_data.type == Vtk::INT64);
  readAppended(input, vec_connectivity, connectivity_data.offset);
  assert(vec_connectivity.size() == std::size_t(vec_offsets.back()));
}


// @{ implementation detail
/**
 * Read compressed data into `buffer_in`, uncompress it and store the result in
 * the concrete-data-type `buffer`
 * \param bs     Size of the uncompressed data
 * \param cbs    Size of the compressed data
 * \param input  Stream to read from.
 **/
template <class T, class IStream>
void read_compressed (T* buffer, unsigned char* buffer_in,
                      std::uint64_t bs, std::uint64_t cbs, IStream& input)
{
#ifdef HAVE_ZLIB
  uLongf uncompressed_space = uLongf(bs);
  uLongf compressed_space = uLongf(cbs);

  Bytef* compressed_buffer = reinterpret_cast<Bytef*>(buffer_in);
  Bytef* uncompressed_buffer = reinterpret_cast<Bytef*>(buffer);

  input.read((char*)(compressed_buffer), compressed_space);
  assert(uLongf(input.gcount()) == compressed_space);

  if (uncompress(uncompressed_buffer, &uncompressed_space, compressed_buffer, compressed_space) != Z_OK) {
    std::cerr << "Zlib error while uncompressing data.\n";
    std::abort();
  }
  assert(uLongf(bs) == uncompressed_space);
#else
  std::cerr << "Can not call read_compressed without compression enabled!\n";
  std::abort();
#endif
}
// @}


template <class Grid, class Creator>
  template <class T>
void VtkReader<Grid,Creator>::readAppended (std::ifstream& input, std::vector<T>& values, std::uint64_t offset)
{
  input.seekg(offset0_ + offset);

  std::uint64_t size = 0;

  std::uint64_t num_blocks = 0;
  std::uint64_t block_size = 0;
  std::uint64_t last_block_size = 0;
  std::vector<std::uint64_t> cbs; // compressed block sizes

  // read total size / block-size(s)
  if (format_ == Vtk::COMPRESSED) {
    input.read((char*)&num_blocks, sizeof(std::uint64_t));
    input.read((char*)&block_size, sizeof(std::uint64_t));
    input.read((char*)&last_block_size, sizeof(std::uint64_t));

    assert(block_size % sizeof(T) == 0);

    // total size of the uncompressed data
    size = block_size * (num_blocks-1) + last_block_size;

    // size of the compressed blocks
    cbs.resize(num_blocks);
    input.read((char*)cbs.data(), num_blocks*sizeof(std::uint64_t));
  } else {
    input.read((char*)&size, sizeof(std::uint64_t));
  }
  assert(size > 0 && (size % sizeof(T)) == 0);
  values.resize(size / sizeof(T));

  if (format_ == Vtk::COMPRESSED) {
    // upper bound for compressed block-size
    std::uint64_t compressed_block_size = block_size + (block_size + 999)/1000 + 12;
    // number of values in the full blocks
    std::size_t num_values = block_size / sizeof(T);

    std::vector<unsigned char> buffer_in(compressed_block_size);
    for (std::size_t i = 0; i < std::size_t(num_blocks); ++i) {
      std::uint64_t bs = i < std::size_t(num_blocks-1) ? block_size : last_block_size;
      read_compressed(values.data() + i*num_values, buffer_in.data(), bs, cbs[i], input);
    }
  } else {
    input.read((char*)(values.data()), size);
    assert(input.gcount() == std::streamsize(size));
  }
}


template <class Grid, class Creator>
void VtkReader<Grid,Creator>::createGrid () const
{
  assert(vec_points.size() == numberOfPoints_);
  assert(vec_types.size() == numberOfCells_);
  assert(vec_offsets.size() == numberOfCells_);

  Creator::create(*factory_, vec_points, vec_types, vec_offsets, vec_connectivity);
}

// Assume input already read the line <AppendedData ...>
template <class Grid, class Creator>
std::uint64_t VtkReader<Grid,Creator>::findAppendedDataPosition (std::ifstream& input) const
{
  char c;
  while (input.get(c) && std::isblank(c)) { /*do nothing*/ }

  std::uint64_t offset = input.tellg();
  if (c != '_')
    --offset; // if char is not '_', assume it is part of the data.

  return offset;
}


template <class Grid, class Creator>
std::map<std::string, std::string> VtkReader<Grid,Creator>::parseXml (std::string const& line, bool& closed)
{
  closed = false;
  std::map<std::string, std::string> attr;

  Sections sec = NO_SECTION;
  bool escape = false;

  std::string name = "";
  std::string value = "";
  for (auto c : line) {
    switch (sec) {
    case NO_SECTION:
      if (std::isalpha(c) || c == '_') {
        name.clear();
        sec = XML_NAME;
        name.push_back(c);
      } else if (c == '/') {
        closed = true;
      }
      break;
    case XML_NAME:
      if (std::isalpha(c) || c == '_')
        name.push_back(c);
      else
        sec = (c == '=' ? XML_NAME_ASSIGN : NO_SECTION);
      break;
    case XML_NAME_ASSIGN:
      value.clear();
      escape = false;
      assert( c == '"' && "Format error!" );
      sec = XML_VALUE;
      break;
    case XML_VALUE:
      if (c == '"' && !escape) {
        attr[name] = value;
        sec = NO_SECTION;
      } else if (c == '\\' && !escape) {
        escape = true;
      }  else {
        value.push_back(c);
        escape = false;
      }
      break;
    default:
      assert(false && "Format error!");
    }
  }

  return attr;
}

} // end namespace Dune
