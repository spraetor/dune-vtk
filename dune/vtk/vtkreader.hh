#pragma once

#include <iosfwd>
#include <map>

#include <dune/vtk/filereader.hh>
#include <dune/vtk/forward.hh>
#include <dune/vtk/vtktypes.hh>

// default GridCreator
#include <dune/vtk/gridcreators/continuousgridcreator.hh>

namespace Dune
{
  /// File-Reader for Vtk unstructured .vtu files
  /**
   * Reads .vtu files and constructs a grid from the cells stored in the file
   * Additionally, stored data can be read.
   *
   * Assumption on the file structure: Each XML tag must be on a separate line.
   **/
  template <class Grid, class GridCreator>
  class VtkReader
      : public FileReader<Grid, VtkReader<Grid, GridCreator>>
  {
    // Sections visited during the xml parsing
    enum Sections {
      NO_SECTION = 0, VTK_FILE, UNSTRUCTURED_GRID, PIECE, POINT_DATA, PD_DATA_ARRAY, CELL_DATA, CD_DATA_ARRAY,
      POINTS, POINTS_DATA_ARRAY, CELLS, CELLS_DATA_ARRAY, APPENDED_DATA, XML_NAME, XML_NAME_ASSIGN, XML_VALUE
    };

    struct DataArrayAttributes
    {
      Vtk::DataTypes type;
      std::size_t components = 1;
      std::uint64_t offset = 0;
    };

    using Entity = typename Grid::template Codim<0>::Entity;
    using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;

  public:
    /// Constructor. Creates a new GridCreator with the passed factory
    template <class... Args>
    explicit VtkReader (Args&&... args)
      : creatorStorage_(std::make_unique<GridCreator>(std::forward<Args>(args)...))
      , creator_(*creatorStorage_)
    {}

    /// Constructor. Stores a references to the passed creator
    VtkReader (GridCreator& creator)
      : creator_(creator)
    {}

    // disable copy and move operations
    VtkReader(VtkReader const&) = delete;
    VtkReader(VtkReader&&) = delete;
    VtkReader& operator=(VtkReader const&) = delete;
    VtkReader& operator=(VtkReader&&) = delete;

    /// Read the grid from file with `filename` into the GridFactory \ref factory_
    /**
     * \param filename  The name of the input file
     * \param create    If `false`, only fill internal data structures, if `true`, also create the grid. [true]
     **/
    void readFromFile (std::string const& filename, bool create = true);

    /// Implementation of \ref FileReader interface
    static void readFactoryImpl (GridFactory<Grid>& factory, std::string const& filename)
    {
      VtkReader reader{factory};
      reader.readFromFile(filename);
    }

    /// Advanced read methods
    /// @{

    /// Read the grid from an input stream, referring to a .vtu file, into the GridFactory \ref factory_
    /**
     * \param input   A STL input stream to read the VTK file from.
     * \param create  If `false`, only fill internal data structures, if `true`, also create the grid. [true]
     **/
    void readSerialFileFromStream (std::ifstream& input, bool create = true);

    /// Read the grid from and input stream, referring to a .pvtu file, into the GridFactory \ref factory_
    /**
     * \param input   A STL input stream to read the VTK file from.
     * \param create  If `false`, only fill internal data structures, if `true`, also create the grid. [true]
     **/
    void readParallelFileFromStream (std::ifstream& input, int rank, int size, bool create = true);

    /// Construct a grid using the GridCreator
    // NOTE: requires the internal data structures to be filled by an aforgoing call to readFromFile
    void createGrid(bool insertPieces = true);

    /// @}

    /// Return the filenames of parallel pieces
    std::vector<std::string> const& pieces () const
    {
      return pieces_;
    }

  private:
    // Read values stored on the cells with name `name`
    template <class T>
    Sections readCellData (std::ifstream& /*input*/, std::vector<T>& /*values*/, std::string /*name*/)
    {
      /* does not read anything */
      return CD_DATA_ARRAY;
    }

    template <class T>
    Sections readPointData (std::ifstream& /*input*/, std::vector<T>& /*values*/, std::string /*name*/)
    {
      /* does not read anything */
      return PD_DATA_ARRAY;
    }

    // Read vertex coordinates from `input` stream and store in into `factory`
    Sections readPoints (std::ifstream& input, std::string name);

    template <class T>
    void readPointsAppended (std::ifstream& input);

    // Read cell type, cell offsets and connectivity from `input` stream
    Sections readCells (std::ifstream& input, std::string name);

    void readCellsAppended (std::ifstream& input);

    // Read data from appended section in vtk file, starting from `offset`
    template <class T>
    void readAppended (std::ifstream& input, std::vector<T>& values, std::uint64_t offset);

    // Test whether line belongs to section
    bool isSection (std::string line,
                    std::string key,
                    Sections current,
                    Sections parent = NO_SECTION) const
    {
      bool result = line.substr(1, key.length()) == key;
      if (result && current != parent)
        DUNE_THROW(Exception , "<" << key << "> in wrong section." );
      return result;
    }

    // Find beginning of appended binary data
    std::uint64_t findAppendedDataPosition (std::ifstream& input) const;

    // Read attributes from current xml tag
    std::map<std::string, std::string> parseXml(std::string const& line, bool& closed);

    // clear all vectors
    void clear ();

    auto comm () const
    {
      return MPIHelper::getCollectiveCommunication();
    }

  private:
    std::unique_ptr<GridCreator> creatorStorage_ = nullptr;
    GridCreator& creator_;

    /// Data format, i.e. ASCII, BINARY or COMPRESSED. Read from xml attributes.
    Vtk::FormatTypes format_;

    // Temporary data to construct the grid elements
    std::vector<GlobalCoordinate> vec_points;
    std::vector<std::uint64_t> vec_point_ids;   //< Global unique vertex ID
    std::vector<std::uint8_t> vec_types;        //< VTK cell type ID
    std::vector<std::int64_t> vec_offsets;      //< offset of vertices of cell
    std::vector<std::int64_t> vec_connectivity; //< vertex indices of cell

    std::size_t numberOfCells_ = 0;   //< Number of cells in the grid
    std::size_t numberOfPoints_ = 0;  //< Number of vertices in the grid

    // offset information for appended data
    // map Name -> {DataType,NumberOfComponents,Offset}
    std::map<std::string, DataArrayAttributes> dataArray_;

    // vector of filenames of parallel pieces
    std::vector<std::string> pieces_;

    /// Offset of beginning of appended data
    std::uint64_t offset0_ = 0;
  };

} // end namespace Dune

#include "vtkreader.impl.hh"
