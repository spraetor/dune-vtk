#pragma once

#include <iosfwd>
#include <map>

#include "filereader.hh"
#include "vtktypes.hh"

namespace Dune
{
  /// File-Reader for Vtk .vtu files
  /**
   * Reads .vtu files and constructs a grid from the cells stored in the file
   * Additionally, stored data can be read.
   *
   * Assumption on the file structure: Each XML tag must be on a separate line.
   **/
  template <class Grid>
  class VtkReader
      : public FileReader<Grid, VtkReader<Grid>>
  {
    // Sections visited during the xml parsing
    enum Sections {
      NO_SECTION = 0, VTK_FILE, UNSTRUCTURED_GRID, PIECE, POINT_DATA, PD_DATA_ARRAY, CELL_DATA, CD_DATA_ARRAY,
      POINTS, POINTS_DATA_ARRAY, CELLS, CELLS_DATA_ARRAY, APPENDED_DATA, XML_NAME, XML_NAME_ASSIGN, XML_VALUE
    };

    using Entity = typename Grid::template Codim<0>::Entity;
    using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;

  public:
    /// Constructor. Stores a pointer to the GridFactory.
    VtkReader (GridFactory<Grid>& factory)
      : factory_(&factory)
    {}

    /// Read the grid from file with `filename` into the GridFactory `factory`
    void readFromFile (std::string const& filename);

    /// Implementation of \ref FileReader interface
    static void readFactoryImpl (GridFactory<Grid>& factory, std::string const& filename)
    {
      VtkReader reader{factory};
      reader.readFromFile(filename);
    }

  private:
    // Read values stored on the cells with name `name`
    template <class T>
    Sections readCellData (std::ifstream& /*input*/,
                           std::vector<T>& /*values*/,
                           std::string /*name*/,
                           Vtk::DataTypes /*type*/,
                           std::size_t /*nComponents*/,
                           std::string /*format*/,
                           std::uint64_t /*offset*/)
    {
      /* does not read anything */
      return CD_DATA_ARRAY;
    }

    template <class T>
    Sections readPointData (std::ifstream& /*input*/,
                            std::vector<T>& /*values*/,
                            std::string /*name*/,
                            Vtk::DataTypes /*type*/,
                            std::size_t /*nComponents*/,
                            std::string /*format*/,
                            std::uint64_t /*offset*/)
    {
      /* does not read anything */
      return PD_DATA_ARRAY;
    }

    // Read vertex coordinates from `input` stream and store in into `factory`
    Sections readPoints (std::ifstream& input,
                         std::string name,
                         Vtk::DataTypes type,
                         std::size_t nComponents,
                         std::string format,
                         std::uint64_t offset);

    template <class T>
    void readPointsAppended (std::ifstream& input);

    // Read cell type, cell offsets and connectivity from `input` stream
    Sections readCells (std::ifstream& input,
                        std::string name,
                        Vtk::DataTypes type,
                        std::string format,
                        std::uint64_t offset);

    void readCellsAppended (std::ifstream& input);

    // Read data from appended section in vtk file, starting from `offset`
    template <class T>
    void readAppended (std::ifstream& input, std::vector<T>& values, std::uint64_t offset);

    // Test whether line belongs to section
    bool isSection (std::string line,
                    std::string key,
                    Sections current,
                    Sections parent = NO_SECTION)
    {
      bool result = line.substr(1, key.length()) == key;
      if (result && current != parent)
        DUNE_THROW(Exception , "<" << key << "> in wrong section." );
      return result;
    }

    // Read attributes from current xml tag
    std::map<std::string, std::string> parseXml(std::string const& line, bool& closed);

    // Construct a grid using the GridFactory `factory` and the read vectors
    // \ref vec_types, \ref vec_offsets, and \ref vec_connectivity
    void createGrid() const;

  private:
    GridFactory<Grid>* factory_;

    /// Data format, i.e. ASCII, BINARY or COMPRESSED. Read from xml attributes.
    Vtk::FormatTypes format_;

    // Temporary data to construct the grid elements
    std::vector<std::uint8_t> vec_types; //< VTK cell type ID
    std::vector<std::int64_t> vec_offsets; //< offset of vertices of cell
    std::vector<std::int64_t> vec_connectivity; //< vertex indices of cell

    std::size_t numCells_; //< Number of cells in the grid
    std::size_t numVertices_; // Number of vertices in the grid

    // offset information for appended data
    // map Name -> {DataType,Offset}
    std::map<std::string, std::pair<Vtk::DataTypes,std::uint64_t>> offsets_;

    /// Offset of beginning of appended data
    std::uint64_t offset0_;
  };

} // end namespace Dune

#include "vtkreader.impl.hh"
