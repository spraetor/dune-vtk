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
#include <dune/vtk/datacollectors/yaspdatacollector.hh>
#include <dune/vtk/datacollectors/spdatacollector.hh>

using namespace Dune;
using namespace Dune::experimental;
using namespace Dune::Functions;

namespace Impl_
{
  template <class GridView, class Grid>
  struct StructuredDataCollector;

#if HAVE_DUNE_SPGRID
  template<class GridView, class ct, int dim, template< int > class Ref, class Comm>
  struct StructuredDataCollector<GridView, SPGrid<ct,dim,Ref,Comm>>
  {
    using type = SPDataCollector<GridView>;
  };
#endif

  template<class GridView, int dim, class Coordinates>
  struct StructuredDataCollector<GridView, YaspGrid<dim,Coordinates>>
  {
    using type = YaspDataCollector<GridView>;
  };
}

template <class GridView>
using StructuredDataCollector = typename Impl_::StructuredDataCollector<GridView, typename GridView::Grid>::type;


template <int dim>
using int_ = std::integral_constant<int,dim>;

template <class GridView>
void write(std::string prefix, GridView const& gridView)
{
  using namespace BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());

  std::vector<double> p1function(basis.dimension());
  interpolate(basis, p1function, [](auto const& x) {
    return 100*x[0] + 10*x[1] + 1*x[2];
  });

  auto fct1 = makeDiscreteGlobalBasisFunction<double>(basis, p1function);
  auto fct2 = makeAnalyticGridViewFunction([](auto const& x) {
    return std::sin(10*x[0])*std::cos(10*x[1])+std::sin(10*x[2]);
  }, gridView);

  {
    using Writer = VtkStructuredGridWriter<GridView, StructuredDataCollector<GridView>>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(fct1, "p1");
    vtkWriter.addCellData(fct1, "p0");
    vtkWriter.addPointData(fct2, "analytic");
    vtkWriter.write(prefix + "sg_ascii_float32.vts", Vtk::ASCII);
  }

  {
    using Writer = VtkImageDataWriter<GridView, StructuredDataCollector<GridView>>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(fct1, "p1");
    vtkWriter.addCellData(fct1, "p0");
    vtkWriter.addPointData(fct2, "analytic");
    vtkWriter.write(prefix + "id_ascii_float32.vti", Vtk::ASCII);
  }
}

template <int dim>
void write_yaspgrid(std::integral_constant<int,dim>)
{
  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(8);
  GridType grid(upperRight,numElements,0,0);

  write("yasp_" + std::to_string(dim) + "d_", grid.leafGridView());
}

template <int dim>
void write_spgrid(std::integral_constant<int,dim>)
{
#if HAVE_DUNE_SPGRID
  using GridType = SPGrid<double,dim, SPIsotropicRefinement>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(8);
  GridType grid(SPDomain<double,dim>::unitCube(),numElements);

  write("sp_" + std::to_string(dim) + "d_", grid.leafGridView());
#endif
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  Hybrid::forEach(std::make_tuple(int_<1>{}, int_<2>{}, int_<3>{}), [](auto const dim)
  {
    write_yaspgrid(dim);
    write_spgrid(dim);
  });
}
