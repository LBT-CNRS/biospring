#include "IO/PeptideSiteRuleReader.h"
#include "IO/LobeSpec.h"

#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace io
{

void PeptideSiteRuleReader::_parse_line(const std::string & line, size_t line_id)
{
    const auto tokens = utils::string::split(line);
    if (tokens.size() < 3 || tokens.size() > 4)
        logging::die("PeptideSiteRuleReader: line %d: invalid number of tokens (expected 3 or 4, found %d) -- "
                     "'<resname> <atomname> <antecedent> [[~|^]<plane>[:angle]]'",
                     static_cast<int>(line_id), static_cast<int>(tokens.size()));

    PeptideSite site;
    site.antecedent = tokens[2];
    if (tokens.size() == 4)
        site.antecedent2 =
            parse_lobe_spec(tokens[3], site.lobeMode, site.lobeAngle, "PeptideSiteRuleReader", line_id);

    const auto key = std::make_pair(tokens[0], tokens[1]);
    if (tokens[0] == "*" && site.antecedent == "*")
        logging::die("PeptideSiteRuleReader: line %d: '*' is a residue wildcard, not an atom or an antecedent",
                     static_cast<int>(line_id));
    if (_sites.count(key))
        logging::die("PeptideSiteRuleReader: line %d: %s:%s is declared twice", static_cast<int>(line_id),
                     tokens[0].c_str(), tokens[1].c_str());
    _sites[key] = site;
}

void PeptideSiteRuleReader::read()
{
    safeOpen();

    std::string buffer;
    size_t line_id = 0;
    while (_instream)
    {
        line_id++;
        std::getline(_instream, buffer);
        buffer = biospring::utils::string::trim(buffer);
        if (!buffer.empty() && buffer[0] != '#')
            _parse_line(buffer, line_id);
    }
    close();

    logging::info("PeptideSiteRuleReader: read %d reactive site(s) from %s.", static_cast<int>(_sites.size()),
                  _filename.c_str());
}

} // namespace io
} // namespace biospring
