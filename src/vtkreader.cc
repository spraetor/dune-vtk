// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/filledarray.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::experimental;

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 3;

  using GridType = UGGrid<dim>;
  using GridView = typename GridType::LeafGridView;
  {
    FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim,unsigned int>(4);
    auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
    auto& grid = *gridPtr;

    GridView gridView = grid.leafGridView();

    VtkWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("test_ascii_float32.vtu", Vtk::ASCII);
    vtkWriter.write("test_binary_float32.vtu", Vtk::BINARY);
    vtkWriter.write("test_compressed_float64.vtu", Vtk::COMPRESSED, Vtk::FLOAT64);
  }

  {
    auto gridPtr = VtkReader<GridType>::read("test_ascii_float32.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float32_2.vtu", Vtk::ASCII);
  }

  {
    auto gridPtr = VtkReader<GridType>::read("test_binary_float32.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float32_3.vtu", Vtk::ASCII);
  }

  {
    auto gridPtr = VtkReader<GridType>::read("test_compressed_float64.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float64_3.vtu", Vtk::ASCII, Vtk::FLOAT64);
  }

  {
    auto gridPtr = VtkReader<GridType,ConnectedGridCreator>::read("test_ascii_float32.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float32_4.vtu", Vtk::ASCII);
  }

  if (filesystem::exists("paraview_3d.vtu")) {
    using GridType3d = UGGrid<3>;
    using GridView3d = typename GridType3d::LeafGridView;
    auto gridPtr = VtkReader<GridType3d>::read("paraview_3d.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView3d> vtkWriter(grid.leafGridView());
    vtkWriter.write("paraview_3d_ascii.vtu", Vtk::ASCII, Vtk::FLOAT64);
  }

  if (filesystem::exists("paraview_2d.vtu")) {
    std::cout << "paraview_2d_ascii...\n";
    using GridType2d = UGGrid<2>;
    using GridView2d = typename GridType2d::LeafGridView;
    auto gridPtr = VtkReader<GridType2d>::read("paraview_2d.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView2d> vtkWriter(grid.leafGridView());
    vtkWriter.write("paraview_2d_ascii.vtu", Vtk::ASCII, Vtk::FLOAT64);
  }
}
