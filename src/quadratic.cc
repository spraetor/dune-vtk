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

#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/datacollectors/quadraticdatacollector.hh>

using namespace Dune;
using namespace Dune::experimental;
using namespace Dune::Functions;

#define GRID_TYPE 2

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 3;

#if GRID_TYPE == 1
  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(4);
  GridType grid(upperRight,numElements);
#elif GRID_TYPE == 2
  using GridType = UGGrid<dim>;
  FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,unsigned int>(4);
  auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);
  auto& grid = *gridPtr;
#endif

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

  using Writer = VtkUnstructuredGridWriter<GridView, QuadraticDataCollector<GridView>>;
  Writer vtkWriter(gridView);
  vtkWriter.addPointData(p1FctWrapped, "p1");
  vtkWriter.addCellData(p1FctWrapped, "p0");

  // write analytic function
  auto p1Analytic = makeAnalyticGridViewFunction([](auto const& x) {
    return std::sin(10*x[0])*std::cos(10*x[1])+std::sin(10*x[2]);
  }, gridView);

  vtkWriter.addPointData(p1Analytic, "analytic");

  vtkWriter.write("p2_ascii_float32.vtu", Vtk::ASCII);
  vtkWriter.write("p2_binary_float32.vtu", Vtk::BINARY);
  vtkWriter.write("p2_compressed_float32.vtu", Vtk::COMPRESSED);
  vtkWriter.write("p2_ascii_float64.vtu", Vtk::ASCII, Vtk::FLOAT64);
  vtkWriter.write("p2_binary_float64.vtu", Vtk::BINARY, Vtk::FLOAT64);
  vtkWriter.write("p2_compressed_float64.vtu", Vtk::COMPRESSED, Vtk::FLOAT64);
}
