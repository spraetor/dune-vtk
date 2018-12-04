#pragma once

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <dune/common/std/optional.hh>
#include <dune/vtk/filewriter.hh>
#include <dune/vtk/forward.hh>
#include <dune/vtk/vtkfunction.hh>
#include <dune/vtk/vtktypes.hh>

namespace Dune
{
  /// Interface for file writers for the Vtk XML file formats
  template <class GridView, class DataCollector>
  class VtkWriterInterface
      : public FileWriter
  {
    template <class> friend class VtkTimeseriesWriter;
    template <class> friend class PvdWriter;

  protected:
    static constexpr int dimension = GridView::dimension;

    using VtkFunction = Dune::VtkFunction<GridView>;
    using pos_type = typename std::ostream::pos_type;

    enum PositionTypes {
      POINT_DATA,
      CELL_DATA
    };

  public:
    /// Constructor, stores the gridView
    VtkWriterInterface (GridView const& gridView,
                        Vtk::FormatTypes format = Vtk::BINARY,
                        Vtk::DataTypes datatype = Vtk::FLOAT32)
      : dataCollector_(gridView)
      , rank_(gridView.comm().rank())
      , numRanks_(gridView.comm().size())
      , format_(format)
      , datatype_(datatype)
    {
#if !HAVE_VTK_ZLIB
      if (format_ == Vtk::COMPRESSED) {
        std::cout << "Dune is compiled without compression. Falling back to BINARY VTK output!\n";
        format_ = Vtk::BINARY;
      }
#endif
    }

    /// Write the attached data to the file
    /**
     * \param fn   Filename of the VTK file. May contain a directory and any file extension.
     * \param dir  The optional parameter specifies the directory of the partition files.
     **/
    virtual void write (std::string const& fn, Std::optional<std::string> dir = {}) const override;

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

  private:
    /// Write a serial VTK file in Unstructured format
    virtual void writeSerialFile (std::ofstream& out) const = 0;

    /// Write a parallel VTK file `pfilename.pvtu` in Unstructured format,
    /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
    /// for [i] in [0,...,size).
    virtual void writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const = 0;

    /// Return the file extension of the serial file (not including the dot)
    virtual std::string fileExtension () const = 0;

    /// Write points and cells in raw/compressed format to output stream
    virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const = 0;

  protected:
    // Write the point or cell values given by the grid function `fct` to the
    // output stream `out`. In case of binary format, stores the streampos of XML
    // attributes "offset" in the vector `offsets`.
    void writeData (std::ofstream& out,
                    std::vector<pos_type>& offsets,
                    VtkFunction const& fct,
                    PositionTypes type,
                    Std::optional<std::size_t> timestep = {}) const;

    // Write points-data and cell-data in raw/compressed format to output stream
    void writeDataAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const;

    // Write the coordinates of the vertices to the output stream `out`. In case
    // of binary format, stores the streampos of XML attributes "offset" in the
    // vector `offsets`.
    void writePoints (std::ofstream& out,
                      std::vector<pos_type>& offsets,
                      Std::optional<std::size_t> timestep = {}) const;

    // Write Appended section and fillin offset values to XML attributes
    void writeAppended (std::ofstream& out, std::vector<pos_type> const& offsets) const;

    // Write the `values` in blocks (possibly compressed) to the output
    // stream `out`. Return the written block size.
    template <class T>
    std::uint64_t writeValuesAppended (std::ofstream& out, std::vector<T> const& values) const;

    template <class T>
    void writeValuesAscii (std::ofstream& out, std::vector<T> const& values) const;

    void writeHeader (std::ofstream& out, std::string const& type) const;

    /// Return PointData/CellData attributes for the name of the first scalar/vector/tensor DataArray
    std::string getNames (std::vector<VtkFunction> const& data) const;

    // Returns endianness
    std::string getEndian () const
    {
      short i = 1;
      return (reinterpret_cast<char*>(&i)[1] == 1 ? "BigEndian" : "LittleEndian");
    }

    // provide accessor to \ref fileExtension virtual method
    std::string getFileExtension () const
    {
      return fileExtension();
    }

    Vtk::FormatTypes getFormat () const
    {
      return format_;
    }

    Vtk::DataTypes getDatatype () const
    {
      return datatype_;
    }

  protected:
    mutable DataCollector dataCollector_;

    // the rank and size of the gridView collective communication
    int rank_ = 0;
    int numRanks_ = 1;

    Vtk::FormatTypes format_;
    Vtk::DataTypes datatype_;

    // attached data
    std::vector<VtkFunction> pointData_;
    std::vector<VtkFunction> cellData_;

    std::size_t const block_size = 1024*32;
    int compression_level = -1; // in [0,9], -1 ... use default value
  };


  template <class Writer>
  struct IsVtkWriter
  {
    template <class GV, class DC>
    static std::uint16_t test(VtkWriterInterface<GV,DC> const&);
    static std::uint8_t  test(...); // fall-back overload

    static constexpr bool value = sizeof(test(std::declval<Writer>())) > sizeof(std::uint8_t);
  };

} // end namespace Dune

#include "vtkwriterinterface.impl.hh"
