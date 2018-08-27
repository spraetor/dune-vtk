# Dune-Vtk
File reader and writer for the VTK Format

## Summary
Provides structured and unstructured file writers for the VTK XML File Formats 
that can be opened in the popular ParaView visualization application. Additionally
a file reader is provided to import VTK files into Dune grid and data objects.

## Requirements
For the management of the grids the `dune-grid` module is required. Additionally
for the data evaluation, the `dune-functions` module is used. 

### Optional modules
For tests and examples `dune-spgrid` and `dune-polygongrid` are suggested, that 
support structured grid data or special element types.

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

| Property          | Dune-Grid | Dune-Vtk |
+-------------------+-----------+----------+
| VTK version       | 0.1       | 1.0      |
| UnstructuredGrid  | x         | x        |
| PolyData          | x         | -        |
| StructuredGrid    | -         | x        |
| RectilinearGrid   | -         | x        |
| ImageData         | -         | x        |
| ASCII             | x         | x        |
| BASE64            | x         | -        |
| APPENDED_RAW      | x         | x        |
| APPENDED_BASE64   | x         | -        |
| BASE64_COMPRESSED | -         | -        |
| APPEDED_COMPRESSED| -         | x        |
| Parallel files    | x         | x        |
