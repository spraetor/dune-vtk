# Dune-Vtk
File reader and writer for the VTK Format

## Summary
Provides structured and unstructured file writers for the VTK XML File Formats 
that can be opened in the popular ParaView visualization application. Additionally
a file reader is provided to import VTK files into Dune grid and data objects.

## Usage
The VTK writer works similar to the dune-grid `VTKWriter`. It needs to be bound 
to a GridView and then data can be added to the points or cells in the grid.
Points are not necessarily grid vertices, but any coordinates placed inside the 
grid cells, so the data must be provided as GridViewFunction to allow the local
evaluation in arbitrary local coordinates.

General interface of a VtkWriter
```c++
template <class GridView, class DataCollector = DefaultDataCollector<GridView>>
class Vtk[Type]Writer
{
public:
  // Constructor
  Vtk[Type]Writer(GridView, Vtk::FormatTypes = Vtk::BINARY, Vtk::DataTypes = Vtk::FLOAT32);
  
  // Bind data to the writer
  Vtk[Type]Writer& addPointData(Function [, std::string name, int numComponents, Vtk::FormatTypes]);
  Vtk[Type]Writer& addCellData(Function [, std::string name, int numComponents, Vtk::FormatTypes]);
  
  // Write file with filename
  void write(std::string filename);
};
```
where `Function` is either a `GridViewFunction`, i.e. supports `bind()`, `unbind()`, and `localFunction(Function)`, or is a legacy `VTKFunction` from Dune-Grid. The optional parameters `name`, `numComponents` and `format` may be given for a `GridViewFunction`.

The parameter `Vtk::FormatTypes` is one of `Vtk::ASCII`, `Vtk::BINARY`, or `Vtk::COMPRESSED` and `Vtk::DataTypes` is one of `Vtk::FLOAT32`, or `Vtk::FLOAT64`. The `[Type]` of a VtkWriter is one of `UnstructuredGrid`, `StructuredGrid`, `RectilinearGrid`, `ImageData`, or `Timeseries`, see below for details. A `DataCollector` may be specified to control how point and cell values are extracted from the `GridView` and the bound data. See `dune/vtk/datacollectors/` of a list of poissible types. The default datacollector extracts a connected grid with continuous data, where points are grid vertices.

See also the `src/` directory for more examples.

## Comparison with Dune::VTKWriter
In Dune-Grid there is a VTK writer available, that is a bit different from the
proposed one. A comparions:

| **Property**       | **Dune-Grid** | **Dune-Vtk** |
| ------------------ | :-----------: | :----------: |
| VTK version        | 0.1           | 1.0          |
| UnstructuredGrid   | **x**         | **x**        |
| PolyData           | (1d)          | -            |
| StructuredGrid     | -             | **x**        |
| RectilinearGrid    | -             | **x**        |
| ImageData          | -             | **x**        |
| ASCII              | **x**         | **x**        |
| BASE64             | **x**         | -            |
| APPENDED_RAW       | **x**         | **x**        |
| APPENDED_BASE64    | **x**         | -            |
| BASE64_COMPRESSED  | -             | -            |
| APPENDED_COMPRESSED| -             | **x**        |
| Parallel files     | **x**         | **x**        |
| Conforming Data    | **x**         | **x**        |
| NonConforming Data | **x**         | **x**        |
| Quadratic Data     | -             | **x**        |
| Subdivided Data    | **x**         | -            |
| Sequence (PVD)     | **x**         | **x**        |
| Timeseries         | -             | **x**        |

## Writers and Readers
Dune-Vtk provides nearly all file formats specified in VTK + 2 time series formats: 
PVD and VTK-Timeseries.

### VtkUnstructuredGridWriter
Implements a VTK file format for unstructured grids with arbitrary element types 
in 1d, 2d, and 3d. Coordinates are specified explicitly and a connectivity table + 
element types are specified for all grid elements (of codim 0). Can be used with 
all Dune grid types.

### VtkStructuredGridWriter
Implements a writer for grid composed of cube elements (lines, pixels, voxels) with 
local numbering similar to Dunes `cube(d)` numbering. The coordinates of the vertices 
can be arbitrary but the connectivity is implicitly given and equals that of 
`Dune::YaspGrid` or `Dune::SPGrid`. Might be chosen as writer for a transformed 
structured grid, using, e.g., a `GeometryGrid` meta-grid. See `src/geometrygrid.cc` 
for an example.

### VtkRectilinearGridWriter
Rectilinear grids are tensor-product grids with given coordinates along the x, y, 
and z axes. Therefore, the grid must allow to extract these 1d coordinates and a 
specialization for a `StructuredDataCollector` must be provided, that implements 
the `ordinates()` function. By default, it assumes constant grid spacing starting 
from a lower left corner. For `YaspGrid` a specialization is implemented if the 
coordinates type is `TensorProductCoordinates`. See `src/structuredgridwriter.cc` 
for an example.

### VtkImageDataWriter
The *most structured* grid is composed of axis-parallel cube elements with constant 
size along each axis. The is implemented in the VtkImageDataWriter. A specialization 
of the `StructuredDataCollector` must implement `origin()` for the lower left corner, 
`wholeExtent()` for the range of cell numbers along each axis in the global grid, 
`extent()` for the range in the local grid, and `spacing()` for the constant grid 
spacing in each direction.

### PvdWriter
A sequence writer, i.e. a collection of timestep files, in the ParaView Data (PVD) 
format. Supports all VtkWriters for the timestep output. In each timestep a collection 
(.pvd) file is created.

### VtkTimseriesWriter
A timeseries is a collection of timesteps stored in one file, instead of separate 
files for each timestep value. Since in the `Vtk::APPENDED` mode, the data is written 
as binary blocks in the appended section of the file and references by an offset 
in the XML DataArray attributes, it allows to reuse written data. An example of 
usage is when the grid points and cells do not change over time, but just the 
point-/cell-data. Then, the grid is written only once and the data is just appended.

Timeseries file are create a bit differently from other Vtk file. There, in the 
first write the grid points and cells are stored in a separate file, and in each 
timestep just the data is written also to temporary files. When you need the timeseries 
file, these stored temporaries are collected and combined to one VTK file. Thus, 
only the minimum amount of data is written in each timestep. The intermediate files 
are stored, by default, in a `/tmp` folder, with (hopefully) fast write access.

### VtkReader
Read in unstructured grid files (.vtu files) and create a new grid, using a GridFactory.
The reader allows to create the grid in multiple ways, by providing a `GridCreator`
template parameter. The `ContinuousGridCreator` reads the connectivity of the grid
as it is and assumes that the elements are already connected correctly. On the other
hand, a `DiscontinuousGridCreator` reconnects separated elements, by identifying 
matching coordinates of the cell vertices.

The VtkReader supports grid creation in parallel. If a partition file .pvtu is 
provided, all partitions can be read by 1. one processor and distributed later on
(`SerialGridCreator`) or read directly in parallel (`ParallelGridCreator`). The later
is currently only available in dune-alugrid 2.6.

**TODO:**

- Provide an interface to read the points-/cell-data from the file
- Extent the implementation to other file formats
