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
#include <dune/common/test/testsuite.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>

using namespace Dune;

// see https://stackoverflow.com/questions/6163611/compare-two-files
bool compare_files (std::string const& fn1, std::string const& fn2)
{
  std::ifstream in1(fn1, std::ios::binary);
  std::ifstream in2(fn2, std::ios::binary);

  std::ifstream::pos_type size1 = in1.seekg(0, std::ifstream::end).tellg();
  in1.seekg(0, std::ifstream::beg);

  std::ifstream::pos_type size2 = in2.seekg(0, std::ifstream::end).tellg();
  in2.seekg(0, std::ifstream::beg);

  if (size1 != size2)
    return false;

  static const std::size_t BLOCKSIZE = 4096;
  std::size_t remaining = size1;

  while (remaining) {
    char buffer1[BLOCKSIZE], buffer2[BLOCKSIZE];
    std::size_t size = std::min(BLOCKSIZE, remaining);

    in1.read(buffer1, size);
    in2.read(buffer2, size);

    if (0 != std::memcmp(buffer1, buffer2, size))
      return false;

    remaining -= size;
  }

  return true;
}

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
void writer_test (GridView const& gridView)
{
  VtkUnstructuredGridWriter<GridView> vtkWriter(gridView);
  for (auto const& test_case : test_cases) {
    vtkWriter.write("/tmp/reader_writer_test_" + std::get<0>(test_case) + ".vtu",
      std::get<1>(test_case), std::get<2>(test_case));
  }
}

template <class Grid, class Test>
void reader_test (Test& test)
{
  for (auto const& test_case : test_cases) {
    auto grid = VtkReader<Grid>::read("/tmp/reader_writer_test_" + std::get<0>(test_case) + ".vtu");
    VtkUnstructuredGridWriter<typename Grid::LeafGridView> vtkWriter(grid->leafGridView());
    vtkWriter.write("/tmp/reader_writer_test_" + std::get<0>(test_case) + "_2.vtu",
      std::get<1>(test_case), std::get<2>(test_case));
    test.check(compare_files("/tmp/reader_writer_test_" + std::get<0>(test_case) + ".vtu",
                             "/tmp/reader_writer_test_" + std::get<0>(test_case) + "_2.vtu"));
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
  Hybrid::forEach(std::make_tuple(int_<2>{}, int_<3>{}), [&test](auto dim)
  {
    using GridType = UGGrid<dim.value>;
    {
      FieldVector<double,dim.value> lowerLeft; lowerLeft = 0.0;
      FieldVector<double,dim.value> upperRight; upperRight = 1.0;
      auto numElements = filledArray<dim.value,unsigned int>(4);
      auto gridPtr = StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, numElements);

      writer_test(gridPtr->leafGridView());
    }

    reader_test<GridType>(test);
  });
#endif

  // Test VtkWriter for YaspGrid
  Hybrid::forEach(std::make_tuple(int_<1>{}, int_<2>{}, int_<3>{}), [](auto dim)
  {
    using GridType = YaspGrid<dim.value>;
    FieldVector<double,dim.value> upperRight; upperRight = 1.0;
    auto numElements = filledArray<dim.value,int>(4);
    GridType grid(upperRight, numElements);
    writer_test(grid.leafGridView());
  });

  return test.exit();
}
