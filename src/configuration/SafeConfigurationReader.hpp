#ifndef __SAFE_CONFIGURATION_READER_H__
#define __SAFE_CONFIGURATION_READER_H__

#include "Configuration.hpp"
#include "IO/ReaderBase.h"

#include <set>

namespace biospring
{
namespace configuration
{

// A SafeConfigurationReader reads a configuration file and dies if it reads a parameter
// that has not been defined in its internal configuration.
// Therefore, it should be initialized using a default configuration.
class SafeConfigurationReader : public ReaderBase
{
  public:
    SafeConfigurationReader(const Configuration & config) : ReaderBase(), _config(config) {}
    SafeConfigurationReader(const std::string & path, const Configuration & config) : ReaderBase(path), _config(config)
    {
    }
    SafeConfigurationReader(const char * const path, const Configuration & config) : ReaderBase(path), _config(config)
    {
    }

    const Configuration & getConfiguration() const { return _config; }

    void read();

  protected:
    Configuration _config;
    std::set<std::string> _initializedParameters;

    struct _Token
    {
        std::string name;
        std::string value;
    };

    static bool _isCommentLine(const std::string & line) { return line.substr(0, 1) == "#"; }
    static bool _isEmptyLine(const std::string & line) { return line.size() == 0; }
    static bool _isDataLine(const std::string & line) { return not(_isCommentLine(line) or _isEmptyLine(line)); }

    bool _isDuplicateParameter(const std::string & name) { return _initializedParameters.count(name); }
    void _registerParameter(const _Token & param);

    static _Token _splitLine(const std::string & line, const size_t lineid);
};

} // namespace configuration
} // namespace biospring

#endif // __SAFE_CONFIGURATION_READER_H__