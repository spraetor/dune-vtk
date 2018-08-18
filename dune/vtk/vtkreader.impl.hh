#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
//#include <regex>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#include "utility/filesystem.hh"
#include "utility/string.hh"

namespace Dune {

template <class Grid>
void VtkReader<Grid>::read(GridFactory<Grid>& factory, std::string const& filename)
{
  // check whether file exists!
  if (!filesystem::exists(filename))
    DUNE_THROW(IOError, "File " << filename << " does not exist!");

  std::ifstream input(filename, std::ios_base::in | std::ios_base::binary);

  std::string data_name = "",
              data_format = "";
  Vtk::DataTypes data_type = Vtk::UNKNOWN;
  std::size_t data_components = 0;
  std::size_t data_offset = 0;

  std::string file_type = "",
              byte_order = "",
              header_type = "",
              compressor = "",
              encoding = "";
  double version = 0.0;

  Sections section = NO_SECTION;
  for (std::string line; std::getline(input, line); ) {
    ltrim(line);

    if (isSection(line, "VTKFile", section)) {

      auto attr = parseXml(line);
      file_type = attr["type"];
      if (file_type != "UnstructuredGrid")
        DUNE_THROW(NotImplemented, "Only UnstructuredGrid format implemented. Found: " << file_type);
      if (!attr["version"].empty())
        version = std::stod(attr["version"]);
      byte_order = attr["byte_order"];
      header_type = attr["header_type"];
      if (!attr["compressor"].empty())
        compressor = attr["compressor"];
      section = VTK_FILE;
    }
    else if (isSection(line, "/VTKFile", section, VTK_FILE))
      section = NO_SECTION;
    else if (isSection(line, "UnstructuredGrid", section, VTK_FILE))
      section = UNSTRUCTURED_GRID;
    else if (isSection(line, "/UnstructuredGrid", section, UNSTRUCTURED_GRID))
      section = VTK_FILE;
    else if (isSection(line, "Piece", section, UNSTRUCTURED_GRID)) {
      auto attr = parseXml(line);
      numVertices_ = std::stol(attr["NumberOfPoints"]);
      numCells_ = std::stol(attr["NumberOfCells"]);
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
      auto attr = parseXml(line);

      data_type = Vtk::Map::datatype[attr["type"]];
      data_name = attr["Name"];
      if (!attr["NumberOfComponents"].empty())
        data_components = std::stoul(attr["NumberOfComponents"]);

      // determine FormatType
      data_format = attr["format"];
      if (data_format == "appended") {
        if (!compressor.empty())
          format_ = Vtk::COMPRESSED;
        else
          format_ = Vtk::BINARY;
      } else {
        format_ = Vtk::ASCII;
      }

      if (!attr["offset"].empty()) {
        data_offset = std::stoul(attr["offset"]);
        assert(data_format == "appended");
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
      auto attr = parseXml(line);
      encoding = attr["encoding"];
      if (encoding != "raw")
        DUNE_THROW(NotImplemented, "Binary encoding != raw not implemented.");
      offset0_ = input.tellg() + 1;
      section = APPENDED_DATA;
    }
    else if (isSection(line, "/AppendedData", section, APPENDED_DATA))
      section = VTK_FILE;

    switch (section) {
      case PD_DATA_ARRAY:
        section = readPointData(input, factory, data_name, data_type, data_components, data_format, data_offset);
        break;
      case POINTS_DATA_ARRAY:
        section = readPoints(input, factory, data_name, data_type, data_components, data_format, data_offset);
        break;
      case CD_DATA_ARRAY:
        section = readCellData(input, factory, data_name, data_type, data_components, data_format, data_offset);
        break;
      case CELLS_DATA_ARRAY:
        section = readCells(input, factory, data_name, data_type, data_format, data_offset);
        break;
      case APPENDED_DATA:
        if (offsets_["points"].first == Vtk::FLOAT32)
          readPointsAppended<float>(input, factory);
        else
          readPointsAppended<double>(input, factory);

        readCellsAppended(input);
        section = NO_SECTION; // finish reading after appended section
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

  createGrid(factory);
}


template <class IStream, class T, class Sections>
Sections read_data_array(IStream& input, std::vector<T>& values, std::size_t max_size,
                         Sections section, Sections parent_section)
{
  values.reserve(max_size);
  using S = std::conditional_t<sizeof(T) <= 8, std::uint16_t, T>;

  std::size_t idx = 0;
  std::string line;
  while (std::getline(input, line)) {
    trim(line);
    if (line.substr(0,12) == std::string("</DataArray>"))
      return parent_section;

    std::istringstream stream(line);
    S value;
    for (; stream >> value; idx++)
      values.push_back(T(value));
  }

  return section;
}


template <class Grid>
typename VtkReader<Grid>::Sections
VtkReader<Grid>::readPoints(std::ifstream& input, GridFactory<Grid>& factory,
                            std::string name, Vtk::DataTypes type,
                            std::size_t nComponents, std::string format, std::uint64_t offset)
{
  if (format == "appended") {
    offsets_["points"] = {type, offset};
    return POINTS;
  }

  assert(numVertices_ > 0);
  using T = typename GlobalCoordinate::value_type;
  std::vector<T> point_values;
  auto sec = read_data_array(input, point_values, 3*numVertices_, POINTS_DATA_ARRAY, POINTS);
  assert(sec == POINTS);
  assert(point_values.size() == 3*numVertices_);

  // extract points from continuous values
  GlobalCoordinate p;
  std::size_t idx = 0;
  for (std::size_t i = 0; i < numVertices_; ++i) {
    for (std::size_t j = 0; j < p.size(); ++j)
      p[j] = point_values[idx++];
    idx += (3u - p.size());
    factory.insertVertex(p);
  }

  return sec;
}


template <class Grid>
  template <class T>
void VtkReader<Grid>::readPointsAppended(std::ifstream& input, GridFactory<Grid>& factory)
{
  assert(numVertices_ > 0);
  auto offset_data = offsets_["points"];
  std::vector<T> point_values;
  readAppended(input, point_values, offset_data.second);
  assert(point_values.size() == 3*numVertices_);

  // extract points from continuous values
  GlobalCoordinate p;
  std::size_t idx = 0;
  for (std::size_t i = 0; i < numVertices_; ++i) {
    for (std::size_t j = 0; j < p.size(); ++j)
      p[j] = T(point_values[idx++]);
    idx += (3u - p.size());

    factory.insertVertex(p);
  }
}


template <class Grid>
typename VtkReader<Grid>::Sections
VtkReader<Grid>::readCells(std::ifstream& input, GridFactory<Grid>& factory,
                           std::string name, Vtk::DataTypes type, std::string format, std::uint64_t offset)
{
  if (format == "appended") {
    offsets_[name] = {type, offset};
    return CELLS;
  }

  Sections sec = CELLS_DATA_ARRAY;

  assert(numCells_ > 0);
  if (name == "types") {
    sec = read_data_array(input, vec_types, numCells_, CELLS_DATA_ARRAY, CELLS);
    assert(vec_types.size() == numCells_);
  } else if (name == "offsets") {
    sec = read_data_array(input, vec_offsets, numCells_, CELLS_DATA_ARRAY, CELLS);
    assert(vec_offsets.size() == numCells_);
  } else if (name == "connectivity") {
    std::size_t max_size = 0;
    int max_vertices = (Grid::dimension == 1 ? 2 : Grid::dimension == 2 ? 4 : 8);
    if (!vec_offsets.empty())
      max_size = vec_offsets.back();
    else
      max_size = numCells_ * max_vertices;
    sec = read_data_array(input, vec_connectivity, max_size, CELLS_DATA_ARRAY, CELLS);
  }

  return sec;
}


template <class Grid>
void VtkReader<Grid>::readCellsAppended(std::ifstream& input)
{
  assert(numCells_ > 0);
  auto types_data = offsets_["types"];
  auto offsets_data = offsets_["offsets"];
  auto connectivity_data = offsets_["connectivity"];

  assert(types_data.first == Vtk::UINT8);
  readAppended(input, vec_types, types_data.second);
  assert(vec_types.size() == numCells_);

  assert(offsets_data.first == Vtk::INT64);
  readAppended(input, vec_offsets, offsets_data.second);
  assert(vec_offsets.size() == numCells_);

  assert(connectivity_data.first == Vtk::INT64);
  readAppended(input, vec_connectivity, connectivity_data.second);
  assert(vec_connectivity.size() == vec_offsets.back());
}


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

  if (uncompress(uncompressed_buffer, &uncompressed_space, compressed_buffer, compressed_space) != Z_OK) {
    std::cerr << "Zlib error while uncompressing data.\n";
    std::abort();
  }
#else
  std::cerr << "Can not call read_compressed without compression enabled!\n";
  std::abort();
#endif
}


template <class Grid>
  template <class T>
void VtkReader<Grid>::readAppended (std::ifstream& input, std::vector<T>& values, std::uint64_t offset)
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
  }
}


template <class Grid>
void VtkReader<Grid>::createGrid(GridFactory<Grid>& factory) const
{
  assert(vec_types.size() == vec_offsets.size());
  std::size_t idx = 0;
  for (std::size_t i = 0; i < vec_types.size(); ++i) {
    if (Vtk::Map::type.count(vec_types[i]) == 0)
      DUNE_THROW(Exception, "Unknown ElementType: " << vec_types[i]);
    auto type = Vtk::Map::type[vec_types[i]];
    Vtk::CellType cellType{type};

    std::size_t nNodes = vec_offsets[i] - (i == 0 ? 0 : vec_offsets[i-1]);
    assert(nNodes > 0);
    std::vector<unsigned int> vtk_cell; vtk_cell.reserve(nNodes);
    for (std::size_t j = 0; j < nNodes; ++j)
      vtk_cell.push_back( vec_connectivity[idx++] );

    // apply index permutation
    std::vector<unsigned int> cell(nNodes);
    for (std::size_t j = 0; j < nNodes; ++j)
      cell[j] = vtk_cell[cellType.localIndex(j)];

    factory.insertElement(type,cell);
  }
}


template <class Grid>
std::map<std::string, std::string> VtkReader<Grid>::parseXml(std::string const& line)
{
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
