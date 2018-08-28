#pragma once

#include <memory>
#include <string>
#include <utility>

namespace Dune
{
  template <class Grid, class FilerReaderImp>
  class FileReader
  {
  private:
    // type of underlying implementation, for internal use only
    using Implementation = FilerReaderImp;

    /// \brief An accessor class to call protected members of reader implementations.
    struct Accessor : public Implementation
    {
      template <class... Args>
      static std::unique_ptr<Grid> readImp (Args&&... args)
      {
        return Implementation::readImpl(std::forward<Args>(args)...);
      }

      template <class... Args>
      static void readFactoryImpl (Args&&... args)
      {
        return Implementation::readFactoryImpl(std::forward<Args>(args)...);
      }
    };

  public:
    /// Reads the grid from a file with filename and returns a unique_ptr to the created grid.
    /// Redirects to concrete implementation of derivated class.
    template <class... Args>
    static std::unique_ptr<Grid> read (const std::string &filename, Args&&... args)
    {
      return Accessor::readImpl(filename, std::forward<Args>(args)...);
    }

    /// Reads the grid from a file with filename into a grid-factory.
    /// Redirects to concrete implementation of derivated class.
    template <class... Args>
    static void read (GridFactory<Grid> &factory, const std::string &filename, Args&&... args)
    {
      Accessor::readFactoryImpl(factory, filename, std::forward<Args>(args)...);
    }

  protected: // default implementations

    // Default implementation, redirects to factory read implementation.
    template <class... Args>
    static std::unique_ptr<Grid> readImpl (const std::string &filename, Args&&... args)
    {
      GridFactory<Grid> factory;
      read(factory, filename, std::forward<Args>(args)...);

      return std::unique_ptr<Grid>{ factory.createGrid() };
    }

    // Default implementation for reading into grid-factory: produces a runtime-error.
    template <class... Args>
    static void readFactoryImpl (GridFactory<Grid> &/*factory*/, const std::string &/*filename*/,
                                 Args&&... /*args*/)
    {
      DUNE_THROW(NotImplemented,
        "GridReader using a factory argument not implemented for concrete reader implementation.");
    }
  };

} // end namespace Dune
