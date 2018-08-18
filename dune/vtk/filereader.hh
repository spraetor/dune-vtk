#pragma once

#include <string>

namespace Dune
{
  template <class Grid>
  class FileReader
  {
  public:
    /// Virtual destructor
    virtual ~FileReader() = default;

    /// Reads the grid from a file with `filename` and returns a unique_ptr to the created grid.
    virtual std::unique_ptr<Grid> read(std::string const& filename)
    {
      GridFactory<Grid> factory;
      read(factory, filename);
      return std::unique_ptr<Grid>{ factory.createGrid() };
    }

    /// Reads the grid from a file with `filename` into a grid-factory.
    virtual void read(GridFactory<Grid>& factory, std::string const& filename) = 0;
  };

} // end namespace Dune
