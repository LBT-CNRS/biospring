
#include "SafeConfigurationReader.hpp"
#include "../logging.h"
#include "../utils/string.hpp"

namespace biospring
{
namespace configuration
{

void SafeConfigurationReader::read()
{
    // Opens file for reading, dies if error occurs.
    safeOpen();

    // Actual reading.
    std::string buffer;
    for (size_t lineid = 1; std::getline(_instream, buffer); lineid++)
    {
        utils::string::trim(buffer);

        if (_isDataLine(buffer))
        {
            auto param = _splitLine(buffer, lineid);

            if (not _config.exists(param.name))
            {
                logging::die("SafeConfigurationReader: line %d: invalid parameter \"%s\"", lineid, param.name.c_str());
            }

            if (_isDuplicateParameter(param.name))
            {
                logging::die("SafeConfigurationReader: line %d: duplicate parameter error: \"%s\" has been found "
                             "multiple times in the parameter file",
                             lineid, param.name.c_str());
            }

            _registerParameter(param);
        }
    }
    close();
}

void SafeConfigurationReader::_registerParameter(const _Token & param)
{
    // Stores in internal list of initialized parameters for duplicate parameter detection.
    _initializedParameters.insert(param.name);

    // Sets configuration.
    _config.setFromString(param.name, param.value);
}

SafeConfigurationReader::_Token SafeConfigurationReader::_splitLine(const std::string & line, const size_t lineid)
{
    size_t loc = line.find("=");
    if (loc == string::npos)
    {
        logging::die("SafeConfigurationReader: line %d: syntax error: expects param = value (found \"%s\")", lineid,
                     line.c_str());
    }

    _Token tok;
    tok.name = utils::string::trim(line.substr(0, loc));
    tok.value = utils::string::trim(line.substr(loc + 1, line.size()));

    if (tok.name.empty() or tok.value.empty())
    {
        logging::die("SafeConfigurationReader: line %d: syntax error: expects param = value (found \"%s\")", lineid,
                     line.c_str());
    }

    return tok;
}

} // namespace configuration
} // namespace biospring
