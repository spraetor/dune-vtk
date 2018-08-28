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

#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::experimental;
using namespace Dune::Functions;

using TestCases = std::set<std::tuple<std::string,Vtk::FormatTypes,Vtk::DataTypes>>;
static TestCases test_cases = {
  {"ascii32", Vtk::ASCII, Vtk::FLOAT32},
  {"bin32", Vtk::BINARY, Vtk::FLOAT32},
  {"zlib32", Vtk::COMPRESSED, Vtk::FLOAT32},
  {"ascii64", Vtk::ASCII, Vtk::FLOAT64},
  {"bin64", Vtk::BINARY, Vtk::FLOAT64},
  {"zlib64", Vtk::COMPRESSED, Vtk::FLOAT64},
};

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

  VtkWriter<GridView> vtkWriter(gridView);
  vtkWriter.addPointData(p1Interpol, "p1");
  vtkWriter.addCellData(p1Interpol, "p0");
  vtkWriter.addPointData(p1Analytic, "q1");
  vtkWriter.addCellData(p1Analytic, "q0");
  for (auto const& test_case : test_cases) {
    vtkWriter.write(prefix + "_" + std::to_string(GridView::dimension) + "d_" + std::get<0>(test_case) + ".vtu",
      std::get<1>(test_case), std::get<2>(test_case));
  }
}

template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

#ifdef HAVE_UG
  // Test VtkWriter for UGGrid
  Hybrid::forEach(std::make_tuple(int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = UGGrid<dim.value>;
    {
      FieldVector<double,dim.value> lowerLeft; lowerLeft = 0.0;
      FieldVector<double,dim.value> upperRight; upperRight = 1.0;
      auto numElements = filledArray<dim.value,unsigned int>(4);
      auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);

      write("ug", gridPtr->leafGridView());
    }
  });
#endif

  // Test VtkWriter for YaspGrid
  Hybrid::forEach(std::make_tuple(int_<1>{}, int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(4);
    GridType grid(upperRight, numElements);
    write("yasp", grid.leafGridView());
  });
}