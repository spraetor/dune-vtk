#pragma once

#include <array>
#include <iosfwd>
#include <map>

#include <dune/common/std/optional.hh>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/vtkfunction.hh>
#include <dune/vtk/vtktypes.hh>

namespace Dune { namespace experimental
{
  /// File-Writer for Vtk .vtu files
  template <class GridView, class DataCollector>
  class VtkWriterInterface
      : public FileWriter
  {
  protected:
    static constexpr int dimension = GridView::dimension;

    using VtkFunction = Dune::experimental::VtkFunction<GridView>;
    using pos_type = typename std::ostream::pos_type;

    enum PositionTypes {
      POINT_DATA,
      CELL_DATA
    };

  public:
    /// Constructor, stores the gridView
    VtkWriterInterface (GridView const& gridView)
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

    /// Attach point data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    VtkWriterInterface& addPointData (Function const& fct, Args&&... args)
    {
      pointData_.emplace_back(fct, std::forward<Args>(args)...);
      return *this;
    }

    /// Attach cell data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    VtkWriterInterface& addCellData (Function const& fct, Args&&... args)
    {
      cellData_.emplace_back(fct, std::forward<Args>(args)...);
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
                    VtkFunction const& fct,
                    PositionTypes type) const;

    // Collect point or cell data (depending on \ref PositionTypes) and pass
    // the resulting vector to \ref writeAppended.
    template <class T>
    std::uint64_t writeDataAppended (std::ofstream& out,
                                     VtkFunction const& fct,
                                     PositionTypes type) const;

    // Write the coordinates of the vertices to the output stream `out`. In case
    // of binary format, stores the streampos of XML attributes "offset" in the
    // vector `offsets`.
    void writePoints (std::ofstream& out,
                      std::vector<pos_type>& offsets) const;

    // Collect point positions and pass the resulting vector to \ref writeAppended.
    template <class T>
    std::uint64_t writePointsAppended (std::ofstream& out) const;

    // Write the `values` in blocks (possibly compressed) to the output
    // stream `out`. Return the written block size.
    template <class T>
    std::uint64_t writeAppended (std::ofstream& out, std::vector<T> const& values) const;

    /// Return PointData/CellData attributes for the name of the first scalar/vector/tensor DataArray
    std::string getNames (std::vector<VtkFunction> const& data) const
    {
      auto scalar = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.ncomps() == 1; });
      auto vector = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.ncomps() == 3; });
      auto tensor = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.ncomps() == 9; });
      return (scalar != data.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
           + (vector != data.end() ? " Vectors=\"" + vector->name() + "\"" : "")
           + (tensor != data.end() ? " Tensors=\"" + tensor->name() + "\"" : "");
    }

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
    std::vector<VtkFunction> pointData_;
    std::vector<VtkFunction> cellData_;

    std::size_t const block_size = 1024*32;
    int compression_level = -1; // in [0,9], -1 ... use default value
  };

}} // end namespace Dune::experimental

#include "vtkwriterinterface.impl.hh"
