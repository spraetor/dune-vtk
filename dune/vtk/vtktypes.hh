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

    inline std::size_t size(DataTypes type)
    {
      switch (type) {
        case INT8:    return 8;
        case UINT8:   return 8;
        case INT16:   return 16;
        case UINT16:  return 16;
        case INT32:   return 32;
        case UINT32:  return 32;
        case INT64:   return 64;
        case UINT64:  return 64;
        case FLOAT32: return 32;
        case FLOAT64: return 64;
        default:
          std::abort();
          return 0;
      }
    }

    template <class T>
    void cast_value(char const* data, DataTypes type, T& value)
    {
      switch (type) {
        case INT8:
          value = *reinterpret_cast<std::int8_t const*>(data);
          break;
        case UINT8:
          value = *reinterpret_cast<std::uint8_t const*>(data);
          break;
        case INT16:
          value = *reinterpret_cast<std::int16_t const*>(data);
          break;
        case UINT16:
          value = *reinterpret_cast<std::uint16_t const*>(data);
          break;
        case INT32:
          value = *reinterpret_cast<std::int32_t const*>(data);
          break;
        case UINT32:
          value = *reinterpret_cast<std::uint32_t const*>(data);
          break;
        case INT64:
          value = *reinterpret_cast<std::int64_t const*>(data);
          break;
        case UINT64:
          value = *reinterpret_cast<std::uint64_t const*>(data);
          break;
        case FLOAT32:
          value = *reinterpret_cast<float const*>(data);
          break;
        case FLOAT64:
          value = *reinterpret_cast<double const*>(data);
          break;
        default:
          std::abort();
      }
    }

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
      static std::map<std::size_t, GeometryType> type;
      static std::map<std::string, DataTypes> datatype;
    };


    /// Mapping of Dune geometry types to VTK cell types
    class CellType
    {
    public:
      CellType (GeometryType const& t);

      std::uint8_t type () const
      {
        return type_;
      }

      int localIndex (int idx) const
      {
        return permutation_[idx];
      }

    private:
      std::uint8_t type_;
      std::vector<int> permutation_;
    };

  } // end namespace Vtk
} // end namespace Dune
