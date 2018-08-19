#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

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

    enum DataTypes {
      UNKNOWN = 0,
      INT8, UINT8,
      INT16, UINT16,
      INT32, UINT32,
      INT64, UINT64,
      FLOAT32 = 32,
      FLOAT64 = 64
    };

    enum PositionTypes {
      VERTEX_DATA,
      CELL_DATA
    };

    enum ContinuityTypes {
      CONFORMING,
      NONCONFORMING
    };

    struct Map
    {
      static std::map<std::uint8_t, GeometryType> type; // VTK Cell type -> Dune::GeometryType
      static std::map<std::string, DataTypes> datatype; // String -> DataTypes
    };


    /// Mapping of Dune geometry types to VTK cell types
    class CellType
    {
    public:
      CellType (GeometryType const& t);

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
