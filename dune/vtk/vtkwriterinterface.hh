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
  /**
   * \tparam GridView       Model of Dune::GridView
   * \tparam DataCollector  Model of \ref DataCollectorInterface
   **/
  template <class GridView, class DataCollector>
  class VtkWriterInterface
      : public FileWriter
  {
    template <class> friend class VtkTimeseriesWriter;
    template <class> friend class PvdWriter;

  protected:
    static constexpr int dimension = GridView::dimension;

    using VtkFunction = Dune::VtkFunction<GridView>;
    using Communicator = CollectiveCommunication<typename MPIHelper::MPICommunicator>;
    using pos_type = typename std::ostream::pos_type;

    enum PositionTypes {
      POINT_DATA,
      CELL_DATA
    };

  public:
    /// \brief Constructor, passes the gridView to the DataCollector
    /**
     * Creates a new VtkWriterInterface for the provided GridView. Initializes a
     * DataCollector that is used to collect point coordinates, cell connectivity and
     * data values.
     *
     * \param gridView  Implementation of Dune::GridView
     * \param format    Format of the VTK file, either Vtk::BINARY, Vtk::ASCII, or Vtk::COMPRESSED
     * \param datatype  Output datatype used for coordinates and other global floating point values
     **/
    VtkWriterInterface (GridView const& gridView,
                        Vtk::FormatTypes format = Vtk::BINARY,
                        Vtk::DataTypes datatype = Vtk::FLOAT32)
      : dataCollector_(gridView)
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

    /// \brief Write the attached data to the file
    /**
     * \param fn   Filename of the VTK file. May contain a directory and any file extension.
     * \param dir  The optional parameter specifies the directory of the partition files for parallel writes.
     **/
    virtual void write (std::string const& fn, Std::optional<std::string> dir = {}) const override;

    /// \brief Attach point data to the writer
    /**
     * Attach a global function to the writer that will be evaluated at grid points
     * (vertices and higher order points). The global function must be
     * assignable to the function wrapper \ref VtkFunction. Additional argument
     * for output datatype and number of components can be bassed. See \ref VtkFunction
     * Constructor for possible arguments.
     **/
    template <class Function, class... Args>
    VtkWriterInterface& addPointData (Function const& fct, Args&&... args)
    {
      pointData_.emplace_back(fct, std::forward<Args>(args)...);
      return *this;
    }

    /// \brief Attach cell data to the writer
    /**
     * Attach a global function to the writer that will be evaluated at cell centers.
     * The global function must be assignable to the function wrapper \ref VtkFunction.
     * Additional argument for output datatype and number of components can be bassed.
     * See \ref VtkFunction Constructor for possible arguments.
     **/
    template <class Function, class... Args>
    VtkWriterInterface& addCellData (Function const& fct, Args&&... args)
    {
      cellData_.emplace_back(fct, std::forward<Args>(args)...);
      return *this;
    }

  private:
    /// Write a serial VTK file in Unstructured format
    virtual void writeSerialFile (std::ofstream& out) const = 0;

    /// Write a parallel VTK file `pfilename.pvtx` in XML format,
    /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
    /// for [i] in [0,...,size).
    virtual void writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const = 0;

    /// Return the file extension of the serial file (not including the dot)
    virtual std::string fileExtension () const = 0;

    /// Write points and cells in raw/compressed format to output stream
    virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const = 0;

  protected:
    // Write the point or cell values given by the grid function `fct` to the
    // output stream `out`. In case of binary format, append the streampos of XML
    // attributes "offset" to the vector `offsets`.
    void writeData (std::ofstream& out,
                    std::vector<pos_type>& offsets,
                    VtkFunction const& fct,
                    PositionTypes type,
                    Std::optional<std::size_t> timestep = {}) const;

    // Write point-data and cell-data in raw/compressed format to output stream
    void writeDataAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const;

    // Write the coordinates of the vertices to the output stream `out`. In case
    // of binary format, appends the streampos of XML attributes "offset" to the
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

    // Write the `values` in a space and newline separated list of ascii representations.
    // The precision is controlled by the datatype and numerical_limits::digits10.
    template <class T>
    void writeValuesAscii (std::ofstream& out, std::vector<T> const& values) const;

    // Write the XML file header of a VTK file `<VTKFile ...>`
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

    // Returns the VTK file format initialized in the constructor
    Vtk::FormatTypes getFormat () const
    {
      return format_;
    }

    // Returns the global datatype used for coordinates and other global float values
    Vtk::DataTypes getDatatype () const
    {
      return datatype_;
    }

    // Return the global MPI communicator.
    auto comm () const
    {
      return MPIHelper::getCollectiveCommunication();
    }

  protected:
    mutable DataCollector dataCollector_;

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
