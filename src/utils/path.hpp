// Path manipulation utilities.

#ifndef __UTILS_PATH_HPP__
#define __UTILS_PATH_HPP__

#include <string>

namespace biospring
{
namespace utils
{
namespace path
{

// Returns true if a path has an extension.
inline bool hasExtension(const std::string & s)
{
    return s.find(".") != std::string::npos;
}

// Returns true if a path does not have an extension.
inline bool hasNoExtension(const std::string & s)
{
    return not hasExtension(s);
}


// Splits a path into a pair (root, extension).
inline std::pair<std::string, std::string> splitExtension(const std::string & s)
{
    std::string prefix, ext;
    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos)
    {
        prefix = s.substr(0, i);
        ext = s.substr(i + 1, s.length() - i);
    }
    else
    {
        prefix = s;
    }
    std::pair<std::string, std::string> out;
    out.first = prefix;
    out.second = ext;
    return out;
}

// Returns a path's extension.
inline std::string getExtension(std::string path)
{
    size_t n = path.find_last_of('.');
    if (n == std::string::npos)
        return "";
    return path.substr(path.find_last_of('.') + 1, path.size());
}

} // namespace path
} // namespace utils
} // namespace biospring



#endif // __UTILS_PATH_HPP__