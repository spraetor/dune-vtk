#pragma once

#include <array>
#include <iosfwd>
#include <map>

#include "datacollector.hh"
#include "filewriter.hh"
#include "vtkfunction.hh"
#include "vtktypes.hh"

namespace Dune { namespace experimental
{
  /// File-Writer for Vtk .vtu files
  template <class GridView, class DataCollector>
  class VtkWriter
      : public FileWriter
  {
  protected:
    static constexpr int dimension = GridView::dimension;

    using GlobalFunction = VTKFunction<GridView>;
    using LocalFunction = VTKLocalFunction<GridView>;
    using pos_type = typename std::ostream::pos_type;

    enum PositionTypes {
      POINT_DATA,
      CELL_DATA
    };

  public:
    /// Constructor, stores the gridView
    VtkWriter (GridView const& gridView)
      : dataCollector_(gridView)
    {}

    /// Write the attached data to the file
    virtual void write (std::string const& fn) override
    {
      write(fn, Vtk::BINARY);
    }

    /// Write the attached data to the file with \ref Vtk::FormatTypes and \ref Vtk::DataTypes
    void write (std::string const& fn,
                Vtk::FormatTypes format,
                Vtk::DataTypes datatype = Vtk::FLOAT32);

    /// Attach point data to the writer
    template <class GridViewFunction>
    VtkWriter& addPointData (GridViewFunction const& gridViewFct,
                              std::string const& name = {},
                              int ncomps = 1)
    {
      pointData_.emplace_back(gridViewFct, name, ncomps);
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

  protected:
    /// Write a serial VTK file in Unstructured format
    virtual void writeSerialFile (std::string const& filename) const = 0;

    /// Write a parallel VTK file `pfilename.pvtu` in Unstructured format,
    /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
    /// for [i] in [0,...,size).
    virtual void writeParallelFile (std::string const& pfilename, int size) const = 0;

    /// Return the file extension of the serial file (not including the dot)
    virtual std::string fileExtension () const = 0;

    // Write the point or cell values given by the grid function `fct` to the
    // output stream `out`. In case of binary format, stores the streampos of XML
    // attributes "offset" in the vector `offsets`.
    void writeData (std::ofstream& out,
                    std::vector<pos_type>& offsets,
                    GlobalFunction const& fct,
                    PositionTypes type) const;

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

    // Collect point or cell data (depending on \ref PositionTypes) and pass
    // the resulting vector to \ref writeAppended.
    template <class T>
    std::uint64_t writeDataAppended (std::ofstream& out,
                                     GlobalFunction const& localFct,
                                     PositionTypes type) const;

    // Collect point positions and pass the resulting vector to \ref writeAppended.
    template <class T>
    std::uint64_t writePointsAppended (std::ofstream& out) const;

    // Collect element connectivity, offsets and element types, and pass the
    // resulting vectors to \ref writeAppended.
    std::array<std::uint64_t,3> writeCellsAppended (std::ofstream& out) const;

    // Write the `values` in blocks (possibly compressed) to the output
    // stream `out`. Return the written block size.
    template <class T>
    std::uint64_t writeAppended (std::ofstream& out, std::vector<T> const& values) const;

    // Returns endianness
    std::string getEndian () const
    {
      short i = 1;
      return (reinterpret_cast<char*>(&i)[1] == 1 ? "BigEndian" : "LittleEndian");
    }

  protected:
    mutable DataCollector dataCollector_;

    std::string filename_;
    Vtk::FormatTypes format_;
    Vtk::DataTypes datatype_;

    // attached data
    std::vector<GlobalFunction> pointData_;
    std::vector<GlobalFunction> cellData_;

    std::size_t const block_size = 1024*32;
    int compression_level = -1; // in [0,9], -1 ... use default value
  };

}} // end namespace Dune::experimental

#include "vtkwriter.impl.hh"
