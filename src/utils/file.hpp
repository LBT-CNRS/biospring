// File manipulation utilities.

#ifndef __UTILS_FILE_HPP__
#define __UTILS_FILE_HPP__

#include "logging.h"

#include <fstream>
#include <sys/stat.h>

namespace biospring
{
namespace utils
{
namespace file
{

// Opens a file in read mode.
// Dies if:
//     - cannot open it,
//     - file is empty (optionnal).
inline void openread(const std::string & path, std::ifstream & infile, bool diesIfEmpty = false)
{
    if (path.empty())
    {
        logging::die("openread: empty file name");
    }
    infile.open(path);
    if (not infile)
    {
        logging::die("can't open file: '%s'", path.c_str());
    }
    if (diesIfEmpty and infile.peek() == std::ifstream::traits_type::eof())
    {
        logging::die("%s: empty file", path.c_str());
    }
}

// Opens a file in write mode.
// Dies if:
//     - cannot open it,
//     - file is empty (optionnal).
inline void openwrite(const std::string & path, std::ofstream & outfile)
{
    outfile.open(path);
    if (not outfile)
    {
        logging::die("can't open file: '%s'", path.c_str());
    }
}

// Returns true if a file exists.
inline bool exists(const std::string & path)
{
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}


} // namespace file
} // namespace utils
} // namespace biospring

#endif // __UTILS_FILE_HPP__