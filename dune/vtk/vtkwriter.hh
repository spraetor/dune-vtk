#pragma once

#include <array>
#include <iosfwd>
#include <map>

#include "filewriter.hh"
#include "vtkfunction.hh"
#include "vtktypes.hh"

namespace Dune
{
  /// File-Writer for Vtk .vtu files
  template <class GridView>
  class VtkWriter
      : public FileWriter
  {
    static constexpr int dimension = GridView::dimension;

    using GlobalFunction = VTKFunction<GridView>;
    using LocalFunction = VTKLocalFunction<GridView>;
    using pos_type = typename std::ostream::pos_type;

  public:
    /// Constructor, stores the gridView
    VtkWriter (GridView const& gridView)
      : gridView_(gridView)
    {}

    /// Write the attached data to the file
    virtual void write (std::string const& fn) override
    {
      write(fn, Vtk::BINARY);
    }

    /// Write the attached data to the file with \ref Vtk::FormatTypes and \ref Vtk::DataTypes
    virtual void write (std::string const& fn,
                        Vtk::FormatTypes format,
                        Vtk::DataTypes datatype = Vtk::FLOAT32);

    /// Attach vertex data to the writer
    template <class GridViewFunction>
    VtkWriter& addVertexData (GridViewFunction const& gridViewFct,
                              std::string const& name = {},
                              int ncomps = 1)
    {
      vertexData_.emplace_back(gridViewFct, name, ncomps);
      return *this;
    }

    /// Attach cell data to the writer
    template <class GridViewFunction>
    VtkWriter& addCellData (GridViewFunction const& gridViewFct,
                            std::string const& name = {},
                            int ncomps = 1)
    {
      cellData_.emplace_back(gridViewFct, name, ncomps);
      return *this;
    }

  private:
    // Write \ref vertexData_ and \ref cellData_ with set \ref format_ and
    // \ref datatype_ to file given by filename
    void writeImpl (std::string const& filename) const;

    // Write the vertex or cell values given by the grid function `fct` to the
    // output stream `out`. In case of binary format, stores the streampos of XML
    // attributes "offset" in the vector `offsets`.
    void writeData (std::ofstream& out,
                    std::vector<pos_type>& offsets,
                    GlobalFunction const& fct,
                    Vtk::PositionTypes type) const;

    // Write the coordinates of the vertices to the output stream `out`. In case
    // of binary format, stores the streampos of XML attributes "offset" in the
    // vector `offsets`.
    void writePoints (std::ofstream& out,
                      std::vector<pos_type>& offsets) const;

    // Write the element connectivity to the output stream `out`. In case
    // of binary format, stores the streampos of XML attributes "offset" in the
    // vector `offsets`.
    void writeCells (std::ofstream& oust,
                     std::vector<pos_type>& offsets) const;

    // Write the `values` in blocks (possibly compressed) to the output
    // stream `out`. Return the written block size.
    template <class T>
    std::uint64_t writeAppended (std::ofstream& out,
                                 std::vector<T> const& values) const;

    // Collect vertex or cell data (depending on \ref PositionTypes) and pass
    // the resulting vector to \ref writeAppended.
    template <class T>
    std::uint64_t writeDataAppended (std::ofstream& out,
                                     GlobalFunction const& localFct,
                                     Vtk::PositionTypes type) const;

    // Collect vertex positions and pass the resulting vector to \ref writeAppended.
    template <class T>
    std::uint64_t writePointsAppended (std::ofstream& out) const;

    // Collect element connectivity, offsets and element types, and pass the
    // resulting vectors to \ref writeAppended.
    std::array<std::uint64_t,3> writeCellsAppended (std::ofstream& out) const;

    // Returns endianness
    std::string getEndian () const
    {
      short i = 1;
      return (reinterpret_cast<char*>(&i)[1] == 1 ? "BigEndian" : "LittleEndian");
    }

    Vtk::CellType getType (GeometryType const& t) const
    {
      return {t};
    }

  private:
    GridView gridView_;

    std::string filename_;
    Vtk::FormatTypes format_;
    Vtk::DataTypes datatype_;

    // attached data
    std::vector<GlobalFunction> vertexData_;
    std::vector<GlobalFunction> cellData_;

    std::size_t const block_size = 1024*32;
    int compression_level = -1; // in [0,9], -1 ... use default value
  };

} // end namespace Dune

#include "vtkwriter.impl.hh"
