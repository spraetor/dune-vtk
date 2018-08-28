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

#include <dune/grid/yaspgrid.hh>

#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

#include <dune/vtk/datacollectors/continuousdatacollector.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/datacollectors/quadraticdatacollector.hh>

using namespace Dune;
using namespace Dune::Functions;

template <class DataCollector, class GridView, class Fct1, class Fct2>
void write_dc (std::string prefix, GridView const& gridView, Fct1 const& fct1, Fct2 const& fct2)
{
  VtkUnstructuredGridWriter<GridView, DataCollector> vtkWriter(gridView);
  vtkWriter.addPointData(fct1, "p1");
  vtkWriter.addCellData(fct1, "p0");
  vtkWriter.addPointData(fct2, "q1");
  vtkWriter.addCellData(fct2, "q0");

  vtkWriter.write(prefix + "_" + std::to_string(GridView::dimension) + "d_ascii.vtu",
    Vtk::ASCII, Vtk::FLOAT32);
}

template <class GridView>
void write (std::string prefix, GridView const& gridView)
{
  using namespace BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());

  FieldVector<double,GridView::dimension> c;
  if (GridView::dimension > 0) c[0] = 11.0;
  if (GridView::dimension > 1) c[1] = 7.0;
  if (GridView::dimension > 2) c[2] = 3.0;

  std::vector<double> vec(basis.dimension());
  interpolate(basis, vec, [&c](auto const& x) { return c.dot(x); });

  // write discrete global-basis function
  auto p1Interpol = makeDiscreteGlobalBasisFunction<double>(basis, vec);

  // write analytic function
  auto p1Analytic = makeAnalyticGridViewFunction([&c](auto const& x) { return c.dot(x); }, gridView);

  write_dc<ContinuousDataCollector<GridView>>(prefix + "_continuous", gridView, p1Interpol, p1Analytic);
  write_dc<DiscontinuousDataCollector<GridView>>(prefix + "_discontinuous", gridView, p1Interpol, p1Analytic);
  write_dc<QuadraticDataCollector<GridView>>(prefix + "_quadratic", gridView, p1Interpol, p1Analytic);
}

template <int I>
using int_ = std::integral_constant<int,I>;

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  Hybrid::forEach(std::make_tuple(int_<1>{}, int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(4);
    GridType grid(upperRight, numElements);
    write("yasp", grid.leafGridView());
  });
}
