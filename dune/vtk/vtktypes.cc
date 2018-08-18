#include "vtktypes.hh"

#include <iostream>

namespace Dune {
namespace Vtk {

std::map<std::uint8_t, GeometryType> Map::type = {
  { 1, GeometryTypes::vertex },
  { 3, GeometryTypes::line },
  { 5, GeometryTypes::triangle },
  { 9, GeometryTypes::quadrilateral },
  {10, GeometryTypes::tetrahedron },
  {12, GeometryTypes::hexahedron },
  {13, GeometryTypes::prism },
  {14, GeometryTypes::pyramid },
};

std::map<std::string, DataTypes> Map::datatype = {
  {"Int8",    INT8},
  {"UInt8",   UINT8},
  {"Int16",   INT16},
  {"UInt16",  UINT16},
  {"Int32",   INT32},
  {"UInt32",  UINT32},
  {"Int64",   INT64},
  {"UInt64",  UINT64},
  {"Float32", FLOAT32},
  {"Float64", FLOAT64}
};


CellType::CellType (GeometryType const& t)
{
  if (t.isVertex()) {
    type_ = 1;
    permutation_ = {0};
  }
  else if (t.isLine()) {
    type_ = 3;
    permutation_ = {0,1};
  }
  else if (t.isTriangle()) {
    type_ = 5;
    permutation_ = {0,1,2};
  }
  else if (t.isQuadrilateral()) {
    type_ = 9;
    permutation_ = {0,1,2,3};
  }
  else if (t.isTetrahedron()) {
    type_ = 10;
    permutation_ = {0,1,2,3};
  }
  else if (t.isHexahedron()) {
    type_ = 12;
    permutation_ = {0,1,3,2,4,5,7,6};
  }
  else if (t.isPrism()) {
    type_ = 13;
    permutation_ = {0,2,1,3,5,4};
  }
  else if (t.isPyramid()) {
    type_ = 14;
    permutation_ = {0,1,3,2,4};
  }
  else {
    std::cerr << "Geometry Type not supported by VTK!\n";
    std::abort();
  }
}

}} // end namespace Vtk::Dune
