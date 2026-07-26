
#include "IO/ReduceRuleReader.h"
#include "logging.h"
#include "utils/string.hpp"

namespace logging = biospring::logging;

namespace biospring
{
namespace reduce
{

ReduceRule ReduceRuleReader::parse_rule(const std::string & line, size_t line_id) const
{
    const auto & tokens = utils::string::split(line);

    // Dies if the line is misformatted.
    if (tokens.size() < 3)
    {
        logging::die("ReduceRuleReader: line %d: invalid number of tokens (expected at least 3, found %d)", line_id,
                     tokens.size());
    }

    std::string rule_name = tokens[0];
    std::string residue_name = tokens[1];

    // Rule names are not required to be globally unique: the same grain name
    // (e.g. the nucleic acid backbone beads in CGZaccharias.grp) can
    // legitimately be reused across different residues that reduce to the
    // same physical bead type. Rules are always looked up by residue name
    // (see ReduceRuleContainer::get_rules_for_residue), never by name alone.

    // Initializes the rule using its name and residue name.
    ReduceRule rule(rule_name, residue_name);

    // Iterates over the atom names composing the grain in the line.
    for (size_t i = 2; i < tokens.size(); i++)
    {
        if (tokens[i] == "*")
            rule.setAllAtom(true);
        else
        {
            std::string atom_name = tokens[i];
            rule.addAtom(atom_name);
        }
    }
    return rule;
}

void ReduceRuleReader::read()
{
    safeOpen();

    std::string buffer;
    size_t line_id = 0;
    while (_instream)
    {
        line_id++;
        std::getline(_instream, buffer);
        biospring::utils::string::trim(buffer);
        if (buffer[0] != '#' && !buffer.empty())
        {
            ReduceRule rule = parse_rule(buffer, line_id);
            _rules.append(rule);
        }
    }
    close();
}

namespace legacy
{

void ReduceRuleReader::read()
{
    safeOpen();

    if (_reduce == nullptr)
        _reduce = new biospring::reduce::legacy::Reduce();

    std::string buffer;
    size_t lineid = 0;
    while (_instream)
    {
        lineid++;
        std::getline(_instream, buffer);
        biospring::utils::string::trim(buffer);
        if (buffer[0] != '#' and not buffer.empty())
        {
            vector<std::string> tokens = biospring::utils::string::split(buffer);

            if (tokens.size() < 3)
            {
                biospring::logging::die(
                    "ReduceRuleReader: line %d: invalid number of tokens (expected at least 3, found %d)", lineid,
                    tokens.size());
            }
            else
            {
                std::string rulename = tokens[0];
                std::string resname = tokens[1];

                _reduce->addReduceRuleName(rulename);
                _reduce->addResidueName(resname);

                biospring::reduce::ReduceRule rule(rulename, resname);

                for (size_t i = 2; i < tokens.size(); i++)
                {
                    if (tokens[i] == "*")
                    {
                        rule.setAllAtom(true);
                    }
                    else
                    {
                        std::string atomname = tokens[i];
                        _reduce->addAtomName(atomname);
                        rule.addAtom(atomname);
                    }
                }
                _reduce->addReduceRule(rule);
                tokens.clear();
            }
        }
    }
    _reduce->checkOneParticlePerRule();
    _instream.close();
}

} // namespace legacy
} // namespace reduce
} // namespace biospring