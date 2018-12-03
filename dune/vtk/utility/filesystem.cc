#include "filesystem.hh"

#ifdef _WIN32
  #include <direct.h>   // _mkdir
  #define GET_CURRENT_DIR _getcwd
#else
  #include <unistd.h>
  #define GET_CURRENT_DIR getcwd
#endif

#include <errno.h>    // errno, ENOENT, EEXIST
#include <stdio.h>  // defines FILENAME_MAX
#include <sys/stat.h>
#include <sys/types.h>

#include <cctype>
#include <iostream>
#include <fstream>
#include <string>

template <class... Args>
void inline _ignore_(Args&&...) {}

namespace Dune { namespace filesystem {

std::string path::string() const
{
  if (empty())
    return ".";

  auto it = begin();
  auto result = *it;
  for (++it; it != end(); ++it)
    result += preferred_separator + *it;
  return result;
}


void path::split(std::string p)
{
  std::string separators = "/\\";
  bool relative = true;

  trim(p);
  Dune::split(p.begin(), p.end(), separators.begin(), separators.end(),
    [this,&relative](auto first, auto end)
    {
      auto token = std::string(first, end);

      if ((!token.empty() && token != "." && token != "..") || (token.empty() && this->empty())) {
        this->push_back(token);
        relative = false;
      } else if (token == "..") {
        if (relative || this->empty()) {
          this->push_back(token);
        }
        else {
          this->pop_back();
        }
      }
    });
}


path path::stem() const
{
  auto f = filename().string();
  auto pos = f.find_last_of('.');
  if (f == "." || f == ".." || pos == std::string::npos)
    return {f};
  else
    return {f.substr(0,pos)};
}


path path::extension() const
{
  auto f = filename().string();
  auto pos = f.find_last_of('.');
  if (f == "." || f == ".." || pos == std::string::npos)
    return {};
  else
    return {f.substr(pos)};
}


bool path::is_absolute(std::string p)
{
  if (p[0] == '/')
    return true;

  // c:\ or z:/
  if (std::isalpha(p[0]) && p[1] == ':' && (p[2] == '/' || p[2] == '\\'))
    return true;

  return false;
}


path& path::operator/=(path const& p)
{
  insert(end(), p.begin(), p.end());
  original += preferred_separator + p.original;
  return *this;
}


bool path::is_file() const
{
  std::string p = this->string();
  struct stat info;
  return stat(p.c_str(), &info) == 0 && (info.st_mode & S_IFREG) != 0;
}


bool path::is_directory() const
{
  std::string p = this->string();
  struct stat info;
  return stat(p.c_str(), &info) == 0 && (info.st_mode & S_IFDIR) != 0;
}


path current_path()
{
  char cwd_[FILENAME_MAX];
  _ignore_(GET_CURRENT_DIR(cwd_, sizeof(cwd_)));
  std::string cwd(cwd_);
  return { trim(cwd) };
}


bool exists(path const& p)
{
  return p.is_file() || p.is_directory();
}


bool create_directories(path const& p)
{
  if (p.is_directory())
    return true;

  auto parent = p.parent_path();
  if (!parent.empty() && !parent.is_directory())
    create_directories(parent);

#ifdef _WIN32
  int ret = _mkdir(p.string().c_str());
#else
  mode_t mode = 0755;
  int ret = mkdir(p.string().c_str(), mode);
#endif
  if (ret == 0)
    return true;

  switch (errno)
  {
    case ENOENT:
      std::cerr << "parent didn't exist. Should not happen, since parent directory created before!\n";
      std::abort();
      return false;
      break;
    case EEXIST:
      return true;
      break;
    default:
      return false;
  }
}

path relative(path const& a, path const& b)
{
  // find common base path
  auto a_it = a.begin();
  auto b_it = b.begin();
  for (; a_it != a.end() && b_it != b.end(); ++a_it, ++b_it) {
    if (*a_it != *b_it)
      break;
  }

  // combine remaining parts of a to result path
  path rel(".");
  for (; a_it != a.end(); ++a_it)
    rel /= *a_it;

  return rel;
}

} } // end namespace Dune::filesystem
