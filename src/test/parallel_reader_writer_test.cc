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

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/gridcreators/parallelgridcreator.hh>
#include <dune/vtk/gridcreators/serialgridcreator.hh>

using namespace Dune;

// see https://stackoverflow.com/questions/6163611/compare-two-files
bool compare_files (std::string const& fn1, std::string const& fn2)
{
  std::ifstream in1(fn1, std::ios::binary);
  std::ifstream in2(fn2, std::ios::binary);
  if (!in1 || !in2) {
    std::cout << "can not find file " << fn1 << " or file " << fn2 << "\n";
    return false;
  }

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

using TestCase = std::tuple<std::string,Vtk::FormatTypes,Vtk::DataTypes>;
using TestCases = std::set<TestCase>;
static TestCases test_cases = {
  {"ascii32", Vtk::ASCII, Vtk::FLOAT32},
  {"bin32", Vtk::BINARY, Vtk::FLOAT32},
  {"zlib32", Vtk::COMPRESSED, Vtk::FLOAT32},
  {"ascii64", Vtk::ASCII, Vtk::FLOAT64},
  {"bin64", Vtk::BINARY, Vtk::FLOAT64},
  {"zlib64", Vtk::COMPRESSED, Vtk::FLOAT64},
};

template <class Test>
void compare (Test& test, filesystem::path const& dir, filesystem::path const& name)
{
  test.check(compare_files(dir.string() + '/' + name.string() + ".vtu",
                           dir.string() + '/' + name.string() + "_2.vtu"));
}

template <class GridView>
void writer_test (GridView const& gridView)
{
  for (auto const& test_case : test_cases) {
    VtkUnstructuredGridWriter<GridView> vtkWriter(gridView, std::get<1>(test_case), std::get<2>(test_case));
    vtkWriter.write("parallel_reader_writer_test_" + std::get<0>(test_case) + ".vtu");
  }
}

template <class G> struct IsALUGrid : std::false_type {};
#if DUNE_VERSION_GT(DUNE_GRID,2,6) && HAVE_DUNE_ALUGRID
template<int dim, int dimworld, Dune::ALUGridElementType elType, Dune::ALUGridRefinementType refineType, class Comm>
struct IsALUGrid<Dune::ALUGrid<dim,dimworld,elType,refineType,Comm>> : std::true_type {};
#endif

template <class Grid, class Creator>
void reader_test (MPIHelper& mpi, TestSuite& test)
{
  std::string ext = ".vtu";
  if (mpi.size() > 1)
    ext = ".pvtu";

  TestCase test_case = {"ascii32", Vtk::ASCII, Vtk::FLOAT32};

  GridFactory<Grid> factory;
  VtkReader<Grid, Creator> reader{factory};
  reader.readFromFile("parallel_reader_writer_test_" + std::get<0>(test_case) + ext);

  std::unique_ptr<Grid> grid = Hybrid::ifElse(IsALUGrid<Grid>{},
    [&](auto id) { return id(factory).createGrid(std::true_type{}); },
    [&](auto id) { return id(factory).createGrid(); });
  std::vector<std::string> pieces1 = grid->comm().size() > 1 ?
    reader.pieces() :
    std::vector<std::string>{"parallel_reader_writer_test_" + std::get<0>(test_case) + ".vtu"};

  VtkUnstructuredGridWriter<typename Grid::LeafGridView> vtkWriter(grid->leafGridView(),
    std::get<1>(test_case), std::get<2>(test_case));
  vtkWriter.write("parallel_reader_writer_test_" + std::get<0>(test_case) + "_2.vtu");

  GridFactory<Grid> factory2;
  VtkReader<Grid, Creator> reader2{factory2};
  reader2.readFromFile("parallel_reader_writer_test_" + std::get<0>(test_case) + "_2" + ext, false);
  std::vector<std::string> pieces2 = grid->comm().size() > 1 ?
    reader.pieces() :
    std::vector<std::string>{"parallel_reader_writer_test_" + std::get<0>(test_case) + "_2.vtu"};

  test.check(pieces1.size() == pieces2.size(), "pieces1.size == pieces2.size");
  for (std::size_t i = 0; i < pieces1.size(); ++i)
    test.check(compare_files(pieces1[i], pieces2[i]));
}


template <class Grid, class Creator>
void reader_writer_test(MPIHelper& mpi, TestSuite& test, std::string const& testName)
{
  std::cout << "== " << testName << "\n";
  const int dim = Grid::dimension;
  FieldVector<double,dim> lowerLeft; lowerLeft = 0.0;
  FieldVector<double,dim> upperRight; upperRight = 1.0;
  auto numElements = filledArray<dim,unsigned int>(4);
  auto gridPtr = StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, numElements);
  gridPtr->loadBalance();

  writer_test(gridPtr->leafGridView());
  MPI_Barrier(MPI_COMM_WORLD);
  reader_test<Grid, Creator>(mpi,test);
}

#if HAVE_DUNE_ALUGRID
  template <int dim>
  using ALUGridType = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
#endif

template <int I>
using int_ = std::integral_constant<int,I>;

int main (int argc, char** argv)
{
  auto& mpi = Dune::MPIHelper::instance(argc, argv);

  TestSuite test{};

#if HAVE_UG
  reader_writer_test<UGGrid<2>, SerialGridCreator<UGGrid<2>>>(mpi, test, "UGGrid<2>");
  reader_writer_test<UGGrid<3>, SerialGridCreator<UGGrid<3>>>(mpi, test, "UGGrid<3>");
#endif

#if HAVE_DUNE_ALUGRID
  reader_writer_test<ALUGridType<2>, SerialGridCreator<ALUGridType<2>>>(mpi, test, "ALUGridType<2>");
  reader_writer_test<ALUGridType<3>, SerialGridCreator<ALUGridType<3>>>(mpi, test, "ALUGridType<3>");

// #if DUNE_VERSION_LT(DUNE_GRID,2,7)
  reader_writer_test<ALUGridType<2>, ParallelGridCreator<ALUGridType<2>>>(mpi, test, "ALUGridType<2, Parallel>");
  reader_writer_test<ALUGridType<3>, ParallelGridCreator<ALUGridType<3>>>(mpi, test, "ALUGridType<3, Parallel>");
// #endif
#endif

  return test.exit();
}
