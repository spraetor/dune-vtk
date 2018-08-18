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

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 3;

  using GridType = UGGrid<dim>;
  using GridView = typename GridType::LeafGridView;
  {
    FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim,unsigned int>(1);
    auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
    auto& grid = *gridPtr;

    GridView gridView = grid.leafGridView();

    VtkWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("test_ascii_float32.vtu", Vtk::ASCII);
    vtkWriter.write("test_binary_float32.vtu", Vtk::BINARY);
    vtkWriter.write("test_compressed_float64.vtu", Vtk::COMPRESSED, Vtk::FLOAT64);
  }

  {
    VtkReader<GridType> vtkReader{};
    auto gridPtr = vtkReader.read("test_ascii_float32.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float32_2.vtu", Vtk::ASCII);
  }

  {
    VtkReader<GridType> vtkReader{};
    auto gridPtr = vtkReader.read("test_binary_float32.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float32_3.vtu", Vtk::ASCII);
  }

  {
    VtkReader<GridType> vtkReader{};
    auto gridPtr = vtkReader.read("test_compressed_float64.vtu");
    auto& grid = *gridPtr;

    VtkWriter<GridView> vtkWriter(grid.leafGridView());
    vtkWriter.write("test_ascii_float64_3.vtu", Vtk::ASCII, Vtk::FLOAT64);
  }
}
