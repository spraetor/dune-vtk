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

#include <dune/polygongrid/grid.hh>
#include <dune/polygongrid/gridfactory.hh>

#include <dune/vtk/vtkunstructuredgridwriter.hh>

using namespace Dune;
using namespace Dune::Functions;

using GridType = Dune::PolygonGrid<double>;

std::unique_ptr<GridType> createArbitraryGrid ()
{
  const std::vector< Dune::FieldVector< double, 2 > > vertices
    = { { 0.0, 0.0 }, { 0.5, 0.0 }, { 1.0, 0.0 },
        { 0.0, 0.4 }, { 0.5, 0.2 }, { 0.7, 0.4 }, { 1.0, 0.4 },
        { 0.0, 0.7 }, { 0.7, 0.6 }, { 1.0, 0.6 },
        { 0.0, 1.0 }, { 0.3, 1.0 }, { 1.0, 1.0 } };
  const std::vector< std::vector< unsigned int > > polys
    = { { 0, 1, 4, 3 }, { 1, 2, 6, 5, 4 }, { 3, 4, 5, 8, 11, 7 }, { 5, 6, 9, 8 }, { 7, 11, 10 }, { 8, 9, 12, 11 } };

  Dune::GridFactory<GridType> factory;
  for( const auto &vertex : vertices )
    factory.insertVertex( vertex );
  for( const auto &poly : polys )
    factory.insertElement( Dune::GeometryTypes::none( 2 ), poly );
  return std::unique_ptr<GridType>( factory.createGrid() );
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  auto gridPtr = createArbitraryGrid();

  using GridView = typename GridType::LeafGridView;
  GridView gridView = gridPtr->leafGridView();

  using Writer = VtkUnstructuredGridWriter<GridView>;
  Writer vtkWriter(gridView);
  auto p1Analytic = makeAnalyticGridViewFunction([](auto const& x) {
    return std::sin(10*x[0])*std::cos(10*x[1]);
  }, gridView);

  vtkWriter.addPointData(p1Analytic, "p1");
  vtkWriter.addCellData(p1Analytic, "p0");

  vtkWriter.write("poly_ascii_float32.vtu", Vtk::ASCII);
  vtkWriter.write("poly_binary_float32.vtu", Vtk::BINARY);
  vtkWriter.write("poly_compressed_float32.vtu", Vtk::COMPRESSED);
  vtkWriter.write("poly_ascii_float64.vtu", Vtk::ASCII, Vtk::FLOAT64);
  vtkWriter.write("poly_binary_float64.vtu", Vtk::BINARY, Vtk::FLOAT64);
  vtkWriter.write("poly_compressed_float64.vtu", Vtk::COMPRESSED, Vtk::FLOAT64);
}
