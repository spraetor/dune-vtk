// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstring>
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/filledarray.hh>
#include <dune/common/timer.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/vtkwriter.hh>

using namespace Dune;

using TestCasesOld = std::set<std::tuple<std::string,VTK::OutputType,VTK::DataMode>>;
static TestCasesOld test_cases_old = {
  {"asciiC", VTK::ascii, VTK::conforming},
//  {"asciiNC", VTK::ascii, VTK::nonconforming},
//  {"base64C", VTK::appendedbase64, VTK::conforming},
//  {"base64NC", VTK::appendedbase64, VTK::nonconforming},
  {"binC", VTK::appendedraw, VTK::conforming},
//  {"binaryNC", VTK::appendedraw, VTK::nonconforming}
};


using TestCasesNew = std::set<std::tuple<std::string,experimental::Vtk::FormatTypes,experimental::Vtk::DataTypes>>;
static TestCasesNew test_cases_new = {
  {"ascii32", experimental::Vtk::ASCII, experimental::Vtk::FLOAT32},
  {"bin32", experimental::Vtk::BINARY, experimental::Vtk::FLOAT32},
  // {"zlib32", experimental::Vtk::COMPRESSED, experimental::Vtk::FLOAT32},
  // {"ascii64", experimental::Vtk::ASCII, experimental::Vtk::FLOAT64},
  // {"bin64", experimental::Vtk::BINARY, experimental::Vtk::FLOAT64},
  // {"zlib64", experimental::Vtk::COMPRESSED, experimental::Vtk::FLOAT64}
};

template <class GridView>
void writer_old (GridView const& gridView)
{
  Timer t;
  for (auto const& test_case : test_cases_old) {
    t.reset();
    VTKWriter<GridView> vtkWriter(gridView, std::get<2>(test_case));
    vtkWriter.write("/tmp/writer_old_" + std::get<0>(test_case) + ".vtu",
      std::get<1>(test_case));
    std::cout << "  time (writer_old_" + std::get<0>(test_case) + ") = " << t.elapsed() << "\n";
  }
}

template <class GridView>
void writer_new (GridView const& gridView)
{
  Timer t;
  experimental::VtkWriter<GridView> vtkWriter(gridView);
  for (auto const& test_case : test_cases_new) {
    t.reset();
    vtkWriter.write("/tmp/writer_new_" + std::get<0>(test_case) + ".vtu",
      std::get<1>(test_case), std::get<2>(test_case));
    std::cout << "  time (writer_new_" + std::get<0>(test_case) + ") = " << t.elapsed() << "\n";
  }
}

template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  TestSuite test{};

#ifdef HAVE_UG
  // Test VtkWriter for UGGrid
  std::cout << "UGGrid\n";
  Hybrid::forEach(std::make_tuple(int_<2>{}, int_<3>{}), [&test](auto dim)
  {
    using GridType = UGGrid<dim.value>;
    FieldVector<double,dim.value> lowerLeft; lowerLeft = 0.0;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,unsigned int>(10);
    auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);

    std::cout << "DIMENSION " << dim.value << "\n";
    writer_old(gridPtr->leafGridView());
    writer_new(gridPtr->leafGridView());
  });
#endif

  // Test VtkWriter for YaspGrid
  std::cout << "YaspGrid\n";
  Hybrid::forEach(std::make_tuple(int_<1>{}, int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(10);
    GridType grid(upperRight, numElements);

    std::cout << "DIMENSION " << dim.value << "\n";
    writer_old(grid.leafGridView());
    writer_new(grid.leafGridView());
  });

  return test.exit();
}
