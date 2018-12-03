#pragma once

#include <string>
#include <dune/common/std/optional.hh>

namespace Dune
{
  class FileWriter
  {
  public:
    /// Virtual destructor
    virtual ~FileWriter () = default;

    /// Write to file given by `filename` and (optionally) store additional data in `dataDir`
    virtual void write (std::string const& filename, Std::optional<std::string> dataDir = {}) const = 0;
  };

} // end namespace Dune
