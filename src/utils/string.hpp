#ifndef __UTILS_STRING_HPP__
#define __UTILS_STRING_HPP__

#include <algorithm>
#include <cstdarg>
#include <iterator>
#include <regex>
#include <sstream>
#include <string>
#include <vector>


namespace biospring
{
namespace utils
{
namespace string
{

// default implementation
template <typename T> struct TypeName
{
    static const char* Get() { return typeid(T).name(); }
};

template <> struct TypeName<int>
{
    static const char* Get() { return "int"; }
};

template <> struct TypeName<size_t>
{
    static const char* Get() { return "size_t"; }
};

template <> struct TypeName<float>
{
    static const char* Get() { return "float"; }
};

template <> struct TypeName<std::string>
{
    static const char* Get() { return "std::string"; }
};


// Returns the concatenation of strings present in input container and separated by delimiter.
inline std::string join(const std::vector<std::string> & tokens, const std::string & delimiter = "")
{
    std::string s = "";
    for (std::vector<std::string>::const_iterator p = tokens.begin(); p != tokens.end(); ++p)
    {
        s += *p;
        if (p != tokens.end() - 1)
            s += delimiter;
    }
    return s;
}

// Returns a vector of the words in the string s using a delimiter.
inline std::vector<std::string> split(const std::string & s, const std::string & delimiter)
{
    std::string buffer(s);
    std::vector<std::string> tokens;
    size_t pos = 0;
    while ((pos = buffer.find(delimiter)) != std::string::npos)
    {
        tokens.push_back(buffer.substr(0, pos));
        buffer.erase(0, pos + delimiter.length());
    }
    tokens.push_back(buffer);
    return tokens;
}

// Returns a vector of the words in the string s using the
// space as delimiter.
inline std::vector<std::string> split(const std::string & s)
{
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
         std::back_inserter<std::vector<std::string>>(tokens));
    return tokens;
}

inline std::string toupper(const std::string & s)
{
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::toupper);
    return out;
}

inline std::string tolower(const std::string & s)
{
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

// String formatting.
// From http://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
inline std::string format(const std::string fmt, ...)
{
    int size = 100;
    std::string str;
    va_list ap;
    while (1) {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
        va_end(ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
    return str;
}

//
// String trim functions.
// Source: https://www.techiedelight.com/trim-string-cpp-remove-leading-trailing-spaces/
//

// Trim from left.
inline std::string ltrim(const std::string & s) { return std::regex_replace(s, std::regex("^\\s+"), std::string("")); }

// Trim from right.
inline std::string rtrim(const std::string & s) { return std::regex_replace(s, std::regex("\\s+$"), std::string("")); }

// Trim from both sides.
inline std::string trim(const std::string & s) { return ltrim(rtrim(s)); }

// Converts a string to any type.
// Adapted from http://forums.codeguru.com/showthread.php?231054-C-String-How-to-convert-a-string-into-a-numeric-type
template <class T> inline bool from_string(T & t, const std::string & s)
{
    std::istringstream iss(s);
    return !(iss >> t).fail();
}

template <> inline bool from_string(bool & t, const std::string & s)
{
    bool success = true;
    const std::string str = tolower(s);
    if (str == "true" or str == "1" or str == "on" or str == "yes")
        t = true;
    else if (str == "false" or str == "0" or str == "off" or str == "no")
        t = false;
    else
        success = false;
    return success;
}

template <class T> inline T parse_string(const std::string & s)
{
    T value;
    if (not from_string(value, s))
        throw format("cannot convert string \"%s\" to %s", s.c_str(), TypeName<T>::Get()).c_str();
    return value;
}


} // namespace string
} // namespace utils
} // namespace biospring

#endif // __UTILS_STRING_HPP__