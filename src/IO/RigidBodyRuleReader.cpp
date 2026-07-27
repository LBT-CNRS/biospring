
#include "IO/RigidBodyRuleReader.h"
#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace rigidbodygroup
{

RigidBodyRule RigidBodyRuleReader::parse_rule(const std::string & line, size_t line_id) const
{
    const auto & tokens = utils::string::split(line);

    if (tokens.size() < 3)
    {
        logging::die("RigidBodyRuleReader: line %d: invalid number of tokens (expected at least 3, found %d)",
                     line_id, tokens.size());
    }

    std::string rule_name = tokens[0];
    std::string residue_name = tokens[1];

    RigidBodyRule rule(rule_name, residue_name);

    for (size_t i = 2; i < tokens.size(); i++)
        rule.addAtom(tokens[i]);

    return rule;
}

void RigidBodyRuleReader::read()
{
    safeOpen();

    std::string buffer;
    size_t line_id = 0;
    while (_instream)
    {
        line_id++;
        std::getline(_instream, buffer);
        buffer = utils::string::trim(buffer);
        if (!buffer.empty() && buffer[0] != '#')
        {
            RigidBodyRule rule = parse_rule(buffer, line_id);
            _rules.append(rule);
        }
    }
    close();
}

} // namespace rigidbodygroup
} // namespace biospring
