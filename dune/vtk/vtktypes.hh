#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <dune/geometry/type.hh>

namespace Dune { namespace experimental
{
  namespace Vtk
  {
    enum FormatTypes {
      ASCII      = 1<<0,
      BINARY     = 1<<1,
      COMPRESSED = 1<<2,
      APPENDED = BINARY | COMPRESSED
    };

    enum DataTypes {
      UNKNOWN = 0,
      INT8, UINT8,
      INT16, UINT16,
      INT32, UINT32,
      INT64, UINT64,
      FLOAT32 = 32,
      FLOAT64 = 64
    };

    enum CellParametrization {
      LINEAR,
      QUADRATIC
    };

    enum CellTypes : std::uint8_t {
      // Linear VTK cell types
      VERTEX         = 1,
      /* POLY_VERTEX    = 2, // not supported */
      LINE           = 3,
      /* POLY_LINE      = 4, // not supported */
      TRIANGLE       = 5,
      /* TRIANGLE_STRIP = 6, // not supported */
      POLYGON        = 7,
      /* PIXEL          = 8, // not supported */
      QUAD           = 9,
      TETRA          = 10,
      /* VOXEL          = 11, // not supported */
      HEXAHEDRON     = 12,
      WEDGE          = 13,
      PYRAMID        = 14,
      // Quadratic VTK cell types
      QUADRATIC_EDGE       = 21,
      QUADRATIC_TRIANGLE   = 22,
      QUADRATIC_QUAD       = 23,
      QUADRATIC_TETRA      = 24,
      QUADRATIC_HEXAHEDRON = 25
    };

    struct Map
    {
      static std::map<std::uint8_t, GeometryType> from_type; // VTK Cell type -> Dune::GeometryType
      static std::map<std::string, DataTypes> to_datatype; // String -> DataTypes
      static std::map<DataTypes, std::string> from_datatype; // DataTypes -> String
    };


    /// Mapping of Dune geometry types to VTK cell types
    class CellType
    {
    public:
      CellType (GeometryType const& t, CellParametrization = LINEAR);

      /// Return VTK Cell type
      std::uint8_t type () const
      {
        return type_;
      }

      /// Return a permutation of Dune elemenr vertices to conform to VTK element numbering
      int permutation (int idx) const
      {
        return permutation_[idx];
      }

      bool noPermutation () const
      {
        return noPermutation_;
      }

    private:
      std::uint8_t type_;
      std::vector<int> permutation_;
      bool noPermutation_;
    };

  } // end namespace Vtk
}} // end namespace Dune::experimental
