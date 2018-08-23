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

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/vtk/vtkstructuredgridwriter.hh>
#include <dune/vtk/vtkimagedatawriter.hh>

using namespace Dune;
using namespace Dune::experimental;
using namespace Dune::Functions;

#define GRID_TYPE 1

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 3;

  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(8);
  GridType grid(upperRight,numElements);

  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid.leafGridView();

  using namespace BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());

  std::vector<double> p1function(basis.dimension());

  interpolate(basis, p1function, [](auto const& x) {
    return 100*x[0] + 10*x[1] + 1*x[2];
  });

  // write discrete global-basis function
  auto p1FctWrapped = makeDiscreteGlobalBasisFunction<double>(basis, p1function);

  {
    using Writer = VtkStructuredGridWriter<GridView>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(p1FctWrapped, "p1");
    vtkWriter.addCellData(p1FctWrapped, "p0");

    // write analytic function
    auto p1Analytic = makeAnalyticGridViewFunction([](auto const& x) {
      return std::sin(10*x[0])*std::cos(10*x[1])+std::sin(10*x[2]);
    }, gridView);

    vtkWriter.addPointData(p1Analytic, "analytic");

    vtkWriter.write("sg_ascii_float32.vts", Vtk::ASCII);
    vtkWriter.write("sg_binary_float32.vts", Vtk::BINARY);
    vtkWriter.write("sg_compressed_float32.vts", Vtk::COMPRESSED);
    vtkWriter.write("sg_ascii_float64.vts", Vtk::ASCII, Vtk::FLOAT64);
    vtkWriter.write("sg_binary_float64.vts", Vtk::BINARY, Vtk::FLOAT64);
    vtkWriter.write("sg_compressed_float64.vts", Vtk::COMPRESSED, Vtk::FLOAT64);
  }

  {
    using Writer = VtkImageDataWriter<GridView>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(p1FctWrapped, "p1");
    vtkWriter.addCellData(p1FctWrapped, "p0");

    // write analytic function
    auto p1Analytic = makeAnalyticGridViewFunction([](auto const& x) {
      return std::sin(10*x[0])*std::cos(10*x[1])+std::sin(10*x[2]);
    }, gridView);

    vtkWriter.addPointData(p1Analytic, "analytic");

    vtkWriter.write("id_ascii_float32.vti", Vtk::ASCII);
    vtkWriter.write("id_binary_float32.vti", Vtk::BINARY);
    vtkWriter.write("id_compressed_float32.vti", Vtk::COMPRESSED);
    vtkWriter.write("id_ascii_float64.vti", Vtk::ASCII, Vtk::FLOAT64);
    vtkWriter.write("id_binary_float64.vti", Vtk::BINARY, Vtk::FLOAT64);
    vtkWriter.write("id_compressed_float64.vti", Vtk::COMPRESSED, Vtk::FLOAT64);
  }
}
