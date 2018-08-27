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

#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/writers/vtkrectilineargridwriter.hh>
#include <dune/vtk/writers/vtkstructuredgridwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>
#include <dune/vtk/datacollectors/spdatacollector.hh>

using namespace Dune;
using namespace Dune::experimental;
using namespace Dune::Functions;

template <int dim>
using int_ = std::integral_constant<int,dim>;

template <class GridView>
void write(std::string prefix, GridView const& gridView)
{
  auto fct2 = makeAnalyticGridViewFunction([](auto const& x) -> float {
    return std::sin(10*x[0]) * (x.size() > 1 ? std::cos(10*x[1]) : 1) + (x.size() > 2 ? std::sin(10*x[2]) : 0);
  }, gridView);

  {
    using Writer = VtkImageDataWriter<GridView>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(fct2, "analytic");
    vtkWriter.write(prefix + "id_ascii_float32.vti", Vtk::ASCII, Vtk::FLOAT32);
  }

  {
    using Writer = VtkRectilinearGridWriter<GridView>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(fct2, "analytic");
    vtkWriter.write(prefix + "rg_ascii_float32.vtr", Vtk::ASCII, Vtk::FLOAT32);
  }

  {
    using Writer = VtkStructuredGridWriter<GridView>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(fct2, "analytic");
    vtkWriter.write(prefix + "sg_ascii_float32.vts", Vtk::ASCII, Vtk::FLOAT32);
  }

  {
    using Writer = VtkUnstructuredGridWriter<GridView>;
    Writer vtkWriter(gridView);
    vtkWriter.addPointData(fct2, "analytic");
    vtkWriter.write(prefix + "ug_ascii_float32.vts", Vtk::ASCII, Vtk::FLOAT32);
  }
}

template <int dim>
void write_yaspgrid(std::integral_constant<int,dim>)
{
  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(16);
  GridType grid(upperRight,numElements,0,0);
  grid.globalRefine(3);

  write("yasp_" + std::to_string(dim) + "d_", grid.leafGridView());
}

template <int dim>
void write_spgrid(std::integral_constant<int,dim>)
{
#if HAVE_DUNE_SPGRID
  using GridType = SPGrid<double,dim, SPIsotropicRefinement>;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,int>(12);
  GridType grid(SPDomain<double,dim>::unitCube(),numElements);
  // grid.globalRefine(1);

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
