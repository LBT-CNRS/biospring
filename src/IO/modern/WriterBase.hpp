#ifndef __WRITER_BASE_HPP__
#define __WRITER_BASE_HPP__

#include <fstream>
#include <string>

#include "utils.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

class WriterBase
{
  protected:
    std::string _path;
    std::ofstream _ostream;

  public:
    WriterBase() {}
    WriterBase(const std::string & path) : _path(path) {}

    void set_path(const std::string & path) { _path = path; }

    // Open file for writing.
    // Dies if something goes wrong.
    void safe_open()
    {
        if (!_ostream.is_open())
            biospring::utils::file::openwrite(_path, _ostream);
    }

    bool is_open() { return _ostream.is_open(); }
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __WRITER_BASE_HPP__