#pragma once

#include <string>
#include <vector>

#include "string.hh"

namespace Dune
{
  namespace filesystem
  {
    // A minimalistic filesystem class
    class path
        : public std::vector<std::string>
    {
      using Super = std::vector<std::string>;
      using iterator = Super::iterator;
      using const_iterator = Super::const_iterator;

    public:
#ifdef _WIN32
      static constexpr char preferred_separator = '\\';
#else
      static constexpr char preferred_separator = '/';
#endif

    public:
      path() = default;

      // NOTE: implicit conversion is allowed here
      template <class String>
      path(String const& p)
        : original(p)
      {
        split(p);
      }

      template <class InputIt>
      path(InputIt it, InputIt end_it)
        : Super(it, end_it)
      {
        original = this->string();
      }

      template <class String>
      path(std::initializer_list<String> const& list)
        : path(list.begin(), list.end())
      {}

      /// Removes filename path component
      path& remove_filename()
      {
        this->pop_back();
        return *this;
      }

      /// Returns the path of the parent path
      path parent_path() const
      {
        return empty() ? path() : path(begin(), --end());
      }

      /// Returns filename path component
      path filename() const
      {
        return empty() ? path() : path(back());
      }

      /// Returns the stem path component
      path stem() const;

      /// Returns the file extension path component
      path extension() const;

      /// Return the path as string
      std::string string() const;

      /// \brief Return whether a path is an absolute path.
      /** In Linux, test whether the path starts with `/`, in Windows whether it starts
        * with `[a-z]:\\`.
        **/
      static bool is_absolute(std::string p);

      bool is_absolute() const { return is_absolute(original); }

      bool is_relative() const { return !is_absolute(); }

      /// Check whether path is a regular file
      bool is_file() const;

      /// Check whether path is a regular file
      bool is_directory() const;

      /// Lexicographically compares two paths
      bool operator==(path const& p)
      {
        return this->string() == p.string();
      }

      /// Appends elements to the path
      path& operator/=(path const& p);

      /// output of the path
      template <class CharT, class Traits>
      friend std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& out, path const& p)
      {
        out << '"' << p.string() << '"';
        return out;
      }

    protected:

      // split the path string into names separated by a `/`, remove relative directories,
      // like `.` or `..`, if possible.
      void split(std::string p);

    private:
      std::string original = "";
    };

    /// Test whether the path is a valid (existing and accessible) file / directory
    bool exists(path const&);

    /// Create directory and non existing parent directories.
    bool create_directories(path const&);

    /// Returns the current path
    path current_path();

    /// Find the path of `a` relative to directory of `b`
    path relative(path const& a, path const& b);

  } // end namespace filesystem
} // end namespace Dune
