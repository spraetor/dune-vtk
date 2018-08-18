#pragma once

#include <iosfwd>
#include <map>

#include "filereader.hh"
#include "vtktypes.hh"

namespace Dune
{
  /// File-Reader for Vtk .vtu files
  template <class Grid>
  class VtkReader
      : public FileReader<Grid>
  {
    enum Sections {
      NO_SECTION = 0, VTK_FILE, UNSTRUCTURED_GRID, PIECE, POINT_DATA, PD_DATA_ARRAY, CELL_DATA, CD_DATA_ARRAY,
      POINTS, POINTS_DATA_ARRAY, CELLS, CELLS_DATA_ARRAY, APPENDED_DATA, XML_NAME, XML_NAME_ASSIGN, XML_VALUE
    };

    using Entity = typename Grid::template Codim<0>::Entity;
    using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;

  public:
    /// Constructor
    VtkReader ()
      : buffer_(block_size)
    {}

    /// read the file
    virtual void read (GridFactory<Grid>& factory, std::string const& filename) override;
    using FileReader<Grid>::read;

  private:
    Sections readCellData (std::ifstream&,
                           GridFactory<Grid>& /*factory*/,
                           std::string /*name*/,
                           Vtk::DataTypes /*type*/,
                           std::size_t /*nComponents*/,
                           std::string /*format*/,
                           std::size_t /*offset*/)
    {
      /* does not read anything */
      return CD_DATA_ARRAY;
    }

    Sections readPointData (std::ifstream&,
                            GridFactory<Grid>& /*factory*/,
                            std::string /*name*/,
                            Vtk::DataTypes /*type*/,
                            std::size_t /*nComponents*/,
                            std::string /*format*/,
                            std::size_t /*offset*/)
    {
      /* does not read anything */
      return PD_DATA_ARRAY;
    }

    Sections readPoints (std::ifstream& input,
                         GridFactory<Grid>& factory,
                         std::string name,
                         Vtk::DataTypes type,
                         std::size_t nComponents,
                         std::string format,
                         std::uint64_t offset);

    template <class T>
    void readPointsAppended (std::ifstream& input,
                             GridFactory<Grid>& factory);

    Sections readCells (std::ifstream& input,
                        GridFactory<Grid>& factory,
                        std::string name,
                        Vtk::DataTypes type,
                        std::string format,
                        std::uint64_t offset);

    void readCellsAppended (std::ifstream& input);

    template <class T>
    void readAppended (std::ifstream& input, std::vector<T>& values, std::uint64_t offset);

    inline bool isSection (std::string line,
                           std::string key,
                           Sections current,
                           Sections parent = NO_SECTION)
    {
      bool result = line.substr(1, key.length()) == key;
      if (result && current != parent)
        DUNE_THROW(Exception , "<" << key << "> in wrong section." );
      return result;
    }

    std::map<std::string, std::string> parseXml(std::string const& line);

    void createGrid(GridFactory<Grid>& factory) const;

  private:
    Vtk::FormatTypes format_;

    std::vector<std::uint8_t> vec_types;
    std::vector<std::int64_t> vec_offsets;
    std::vector<std::int64_t> vec_connectivity;

    std::size_t numCells_;
    std::size_t numVertices_;
    std::size_t numData_;

    // map Name -> {DataType,Offset}
    std::map<std::string, std::pair<Vtk::DataTypes,std::uint64_t>> offsets_;
    std::uint64_t offset0_;

    std::size_t const block_size = 1024*32;
    std::vector<unsigned char> buffer_;
  };

} // end namespace Dune

#include "vtkreader.impl.hh"
