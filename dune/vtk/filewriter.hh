#pragma once

#include <string>

namespace Dune
{
  class FileWriter
  {
  public:
    /// Virtual destructor
    virtual ~FileWriter () = default;

    /// Write to file given by `filename`
    virtual void write (std::string const& filename) = 0;
  };

} // end namespace Dune
