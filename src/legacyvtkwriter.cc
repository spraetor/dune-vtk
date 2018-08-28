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
#include <dune/vtk/legacyvtkfunction.hh>

using namespace Dune;
using namespace Dune::Functions;

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

  using P1Function = P1VTKFunction<GridView,std::vector<double>>;
  std::shared_ptr<VTKFunction<GridView> const> p1FctWrapped(new P1Function(gridView, p1function, "p1"));

  using Writer = VtkUnstructuredGridWriter<GridView>;
  Writer vtkWriter(gridView);
  vtkWriter.addPointData(p1FctWrapped);
  vtkWriter.write("test_ascii_float32.vtu", Vtk::ASCII);
}
