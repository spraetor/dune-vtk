# Dune-Vtk
File reader and writer for the VTK Format

## Summary
Provides structured and unstructured file writers for the VTK XML File Formats 
that can be opened in the popular ParaView visualization application. Additionally
a file reader is provided to import VTK files into Dune grid and data objects.

## Usage
The VTK writer works similar to the dune-grid `VTKWriter`. It needs to be bound 
to a GridView and then data can be added to the points or cells in the grid.
Points are not necessarily grid vertices, but any coordinates places inside the 
grid cells, so the data must be provided as GridViewFunction to allow the local
evaluation in arbitrary local coordinates.

```c++
Grid grid = ...;
VtkWriter<typename Grid::LeafGridView> vtkWriter(grid.leafGridView());

auto fct = makeAnalyticGridViewFunction([](auto const& x) {
  return std::sin(x[0]);
});

vtkWriter.addPointData(fct, "u_points");
vtkWriter.addCellData(fct, "u_cells");

vtkWriter.write("filename.vtu", Vtk::ASCII);
```

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
| APPEDED_COMPRESSED | -             | **x**        |
| Parallel files     | **x**         | **x**        |
| Conforming Data    | **x**         | **x**        |
| NonConforming Data | **x**         | **x**        |
| Quadratic Data     | -             | **x**        |
| Subdivided Data    | **x**         | -            |

## Example
Create a fine structured YaspGrid and write it in parallel with 8 cores with
different write modi:

### 1. Create the grid

```c++
const int dim = 3;
FieldVector<double,dim> upperRight; upperRight = 1.0;
auto numElements = filledArray<dim,int>(16);
YaspGrid<dim> grid(upperRight,numElements,0,0); // no overlap
grid.globalRefine(3);
auto gridView = grid.leafGridView();
```

### 2. Create a function to write

```c++
auto fct = makeAnalyticGridViewFunction([](auto const& x) {
  return std::sin(10*x[0]) * std::cos(10*x[1]) + std::sin(10*x[2]);
}, gridView);
```

### 3. Write the grid and data to file

```c++
Vtk[FORMAT]Writer<decltype(gridView)> vtkWriter(gridView);
vtkWriter.addPointData(fct, "fct");
vtkWriter.write("filename.vtu", [FORMAT_TYPE], [DATA_TYPE]);
```

where `FORMAT` one of 

- `UnstructuredGrid`,
- `StructuredGrid` (structured connectivity, arbitrary coordinates),
- `RectilinearGrid` (structured connectivity, tensor-product coordinates), or 
- `ImageData` (structured connectivity, axis-parallel coordinates with constant grid-spacing),

`FORMAT_TYPE` one of

- `Vtk::ASCII` (inline ascii format),
- `Vtk::BINARY` (appended raw format), or 
- `Vtk::COMPRESSED` (appended compressed raw format), 

and `DATA_TYPE` one of 

- `Vtk::FLOAT32` (single precision), or 
- `Vtk::FLOAT64` (double precision).

We measure the file size (per processor) and the memory requirement of ParaView
to visualize the data.

| **Setup**                           | **Filesize** | **Memory** |
| ----------------------------------- | ------------ | ---------- |
| UnstructuredGrid, ASCII, FLOAT32    | 26M | 330M |
| *UnstructuredGrid, ASCII, FLOAT64*  | *29M* | *360M* |
| UnstructuredGrid, BINAR, FLOAT32    | 23M | 330M |
| UnstructuredGrid, BINAR, FLOAT64    | 27M | 360M |
| UnstructuredGrid, COMPR, FLOAT32    | 4.5M | 330M |
| UnstructuredGrid, COMPR, FLOAT64    | 5.7M | 360M |
| StructuredGrid, ASCII, FLOAT32      | 10M | 34M |
| StructuredGrid, ASCII, FLOAT64      | 13M | 67M |
| StructuredGrid, BINAR, FLOAT32      | 4.2M | 34M |
| StructuredGrid, BINAR, FLOAT64      | 8.4M | 67M | 
| StructuredGrid, COMPR, FLOAT32      | 1.4M | 34M |
| StructuredGrid, COMPR, FLOAT64      | 2.6M | 67M |
| RectilinearGrid, ASCII, FLOAT32     | 3.0M | 8.4M |
| RectilinearGrid, ASCII, FLOAT64     | 5.3M | 17M |
| RectilinearGrid, BINAR, FLOAT32     | 1.1M | 8.4M |
| RectilinearGrid, BINAR, FLOAT64     | 2.1M | 17M |
| RectilinearGrid, COMPR, FLOAT32     | 970K | 8.4M |
| RectilinearGrid, COMPR, FLOAT64     | 2.0M | 17M |
| ImageData, ASCII, FLOAT32           | 3.0M | 8.4M |
| ImageData, ASCII, FLOAT64           | 5.3M | 17M |
| ImageData, BINAR, FLOAT32           | 1.1M | 8.4M |
| ImageData, BINAR, FLOAT64           | 2.1M | 17M |
| **ImageData, COMPR, FLOAT32**       | **970K** | **8.4M** |
| ImageData, COMPR, FLOAT64           | 2.0M | 17M |
