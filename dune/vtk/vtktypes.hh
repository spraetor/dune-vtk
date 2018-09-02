#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <dune/common/ftraits.hh>
#include <dune/geometry/type.hh>

namespace Dune
{
  namespace Vtk
  {
    enum FormatTypes {
      ASCII      = 1<<0,
      BINARY     = 1<<1,
      COMPRESSED = 1<<2,
      APPENDED = BINARY | COMPRESSED
    };
    std::string to_string (FormatTypes);

    enum DataTypes {
      UNKNOWN = 0,
      INT8, UINT8,
      INT16, UINT16,
      INT32, UINT32,
      INT64, UINT64,
      FLOAT32 = 32,
      FLOAT64 = 64
    };
    std::string to_string (DataTypes);

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
    GeometryType to_geometry (std::uint8_t);

    struct Map
    {
      static std::map<std::string, DataTypes> to_datatype; // String -> DataTypes

      template <class T> struct Type {};

      static constexpr DataTypes typeImpl (Type<std::int8_t>)   { return INT8; }
      static constexpr DataTypes typeImpl (Type<std::uint8_t>)  { return UINT8; }
      static constexpr DataTypes typeImpl (Type<std::int16_t>)  { return INT16; }
      static constexpr DataTypes typeImpl (Type<std::uint16_t>) { return UINT16; }
      static constexpr DataTypes typeImpl (Type<std::int32_t>)  { return INT32; }
      static constexpr DataTypes typeImpl (Type<std::uint32_t>) { return UINT32; }
      static constexpr DataTypes typeImpl (Type<std::int64_t>)  { return INT64; }
      static constexpr DataTypes typeImpl (Type<std::uint64_t>) { return UINT64; }
      static constexpr DataTypes typeImpl (Type<float>)       { return FLOAT32; }
      static constexpr DataTypes typeImpl (Type<double>)      { return FLOAT64; }
      static constexpr DataTypes typeImpl (Type<long double>) { return FLOAT64; }

      template <class T>
      static constexpr DataTypes type = typeImpl(Type<typename FieldTraits<T>::field_type>{});
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
} // end namespace Dune
