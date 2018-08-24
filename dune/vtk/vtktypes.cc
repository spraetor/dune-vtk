#include "vtktypes.hh"

#include <iostream>

namespace Dune { namespace experimental {
namespace Vtk {

std::map<std::uint8_t, GeometryType> Map::from_type = {
  {VERTEX,     GeometryTypes::vertex },
  {LINE,       GeometryTypes::line },
  {TRIANGLE,   GeometryTypes::triangle },
  {QUAD,       GeometryTypes::quadrilateral },
  {TETRA,      GeometryTypes::tetrahedron },
  {HEXAHEDRON, GeometryTypes::hexahedron },
  {WEDGE,      GeometryTypes::prism },
  {PYRAMID,    GeometryTypes::pyramid },
};

std::map<std::string, DataTypes> Map::to_datatype = {
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

std::map<DataTypes, std::string> Map::from_datatype = {
  {INT8,    "Int8"},
  {UINT8,   "UInt8"},
  {INT16,   "Int16"},
  {UINT16,  "UInt16"},
  {INT32,   "Int32"},
  {UINT32,  "UInt32"},
  {INT64,   "Int64"},
  {UINT64,  "UInt64"},
  {FLOAT32, "Float32"},
  {FLOAT64, "Float64"}
};

CellType::CellType (GeometryType const& t, CellParametrization parametrization)
  : noPermutation_(true)
{
  if (parametrization == LINEAR) {
    if (t.isVertex()) {
      type_ = VERTEX;
      permutation_ = {0};
    }
    else if (t.isLine()) {
      type_ = LINE;
      permutation_ = {0,1};
    }
    else if (t.isTriangle()) {
      type_ = TRIANGLE;
      permutation_ = {0,1,2};
    }
    else if (t.isQuadrilateral()) {
      type_ = QUAD;
      permutation_ = {0,1,3,2};
      noPermutation_ = false;
    }
    else if (t.isTetrahedron()) {
      type_ = TETRA;
      permutation_ = {0,1,2,3};
    }
    else if (t.isHexahedron()) {
      type_ = HEXAHEDRON;
      permutation_ = {0,1,3,2,4,5,7,6};
      noPermutation_ = false;
    }
    else if (t.isPrism()) {
      type_ = WEDGE;
      permutation_ = {0,2,1,3,5,4};
      noPermutation_ = false;
    }
    else if (t.isPyramid()) {
      type_ = PYRAMID;
      permutation_ = {0,1,3,2,4};
      noPermutation_ = false;
    }
    else if (t.isNone() && t.dim() == 1) {
      type_ = LINE;
      permutation_ = {0,1};
    }
    else if (t.isNone() && t.dim() == 2) {
      type_ = POLYGON;
      permutation_ = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    }
    else {
      std::cerr << "Geometry Type not supported by VTK!\n";
      std::abort();
    }
  } else if (parametrization == QUADRATIC) {
    if (t.isLine()) {
      type_ = QUADRATIC_EDGE;
      permutation_ = {0,1, 0};
    }
    else if (t.isTriangle()) {
      type_ = QUADRATIC_TRIANGLE;
      permutation_ = {0,1,2, 0,2,1};
      noPermutation_ = false;
    }
    else if (t.isQuadrilateral()) {
      type_ = QUADRATIC_QUAD;
      permutation_ = {0,1,3,2, 3,1,0,2}; // maybe need inverse permutation as well
      noPermutation_ = false;
    }
    else if (t.isTetrahedron()) {
      type_ = QUADRATIC_TETRA;
      permutation_ = {0,1,2,3, 0,2,1,3,4,5};
      noPermutation_ = false;
    }
    else if (t.isHexahedron()) {
      type_ = QUADRATIC_HEXAHEDRON;
      permutation_ = {0,1,3,2,4,5,7,6, 6,5,7,4,10,9,11,8,0,1,3,2}; // maybe need inverse permutation as well
      noPermutation_ = false;
    }
    else {
      std::cerr << "Geometry Type not supported by VTK!\n";
      std::abort();
    }
  }
}

}}} // end namespace Dune::experimental::Vtk
