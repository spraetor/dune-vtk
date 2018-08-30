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

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/vtk/vtktimeserieswriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;
using namespace Dune::Functions;

template <class GridView>
void write (std::string prefix, GridView const& gridView)
{
  FieldVector<double,GridView::dimension> c{11.0, 7.0, 3.0};
  auto p1Analytic = makeAnalyticGridViewFunction([&c](auto const& x) -> float { return c.dot(x); }, gridView);

  using Writer = VtkUnstructuredGridWriter<GridView>;
  VtkTimeseriesWriter<Writer> seriesWriter(gridView, Vtk::BINARY, Vtk::FLOAT32);
  seriesWriter.addPointData(p1Analytic, "q1");
  seriesWriter.addCellData(p1Analytic, "q0");
  std::string filename = prefix + "_" + std::to_string(GridView::dimension) + "d_binary32.vtu";
  for (double t = 0.0; t < 5; t += 0.5) {
    seriesWriter.writeTimestep(t, filename);
  }
  seriesWriter.write(filename);

  Writer vtkWriter(gridView, Vtk::BINARY, Vtk::FLOAT32);
  vtkWriter.addPointData(p1Analytic, "q1");
  vtkWriter.addCellData(p1Analytic, "q0");
  vtkWriter.write(filename);
}

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using GridType = YaspGrid<3>;
  FieldVector<double,3> upperRight; upperRight = 1.0;
  auto numElements = filledArray<3,int>(4);
  GridType grid(upperRight, numElements);
  write("yasp", grid.leafGridView());
}