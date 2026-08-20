#include "IO/DonorAcceptorRuleReader.h"

#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace io
{

void DonorAcceptorRuleReader::_parse_line(const std::string & line, size_t line_id)
{
    const auto tokens = utils::string::split(line);

    if (tokens.size() != 4 && tokens.size() != 5)
    {
        logging::die("DonorAcceptorRuleReader: line %d: invalid number of tokens (expected 4 or 5, found %d)",
                     static_cast<int>(line_id), static_cast<int>(tokens.size()));
    }

    const std::string & resname = tokens[0];
    const std::string & atomname = tokens[1];

    DonorAcceptorRole entry;
    if (!utils::string::from_string(entry.donorCapacity, tokens[2]))
        logging::die("DonorAcceptorRuleReader: line %d: invalid donor capacity '%s' (expected a count)",
                     static_cast<int>(line_id), tokens[2].c_str());
    if (!utils::string::from_string(entry.acceptorCapacity, tokens[3]))
        logging::die("DonorAcceptorRuleReader: line %d: invalid acceptor capacity '%s' (expected a count)",
                     static_cast<int>(line_id), tokens[3].c_str());
    if (tokens.size() == 5)
        entry.antecedent = tokens[4];

    _roles[{resname, atomname}] = entry;
}

void DonorAcceptorRuleReader::read()
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
}

void DonorAcceptorRuleReader::tagParticles(spn::SpringNetwork & spn) const
{
    unsigned nb_donors = 0;
    unsigned nb_acceptors = 0;
    unsigned nb_directed = 0;
    unsigned nb_missing_antecedent = 0;

    // (chain, residue id, atom name) -> particle index, so a rule's
    // antecedent column can be resolved inside its own residue. Built once:
    // the alternative is a rescan of every particle per antecedent.
    std::map<std::tuple<std::string, int, std::string>, unsigned> by_atom;
    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        const spn::Particle & p = spn.getParticle(i);
        by_atom.emplace(std::make_tuple(p.getChainName(), p.getResId(), p.getName()), i);
    }

    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        spn::Particle & p = spn.getParticle(i);
        const auto it = _roles.find({p.getResName(), p.getName()});
        if (it == _roles.end())
            continue;

        p.setDonorCapacity(it->second.donorCapacity);
        p.setAcceptorCapacity(it->second.acceptorCapacity);
        if (it->second.donorCapacity > 0)
            ++nb_donors;
        if (it->second.acceptorCapacity > 0)
            ++nb_acceptors;

        if (it->second.antecedent.empty())
            continue;
        const auto anc = by_atom.find(std::make_tuple(p.getChainName(), p.getResId(), it->second.antecedent));
        if (anc == by_atom.end())
        {
            ++nb_missing_antecedent;
            continue;
        }
        p.setAntecedentIndex(static_cast<int>(anc->second));
        ++nb_directed;
    }

    logging::info("Hydrogen bond tagging: %d donor(s), %d acceptor(s) among %d particles, %d with a resolved "
                  "antecedent (directional).",
                  nb_donors, nb_acceptors, spn.getNumberOfParticles(), nb_directed);
    if (nb_missing_antecedent > 0)
        logging::warning("DonorAcceptorRuleReader: %d particle(s) name an antecedent absent from their residue -- "
                         "left undirected (distance-only), not an error.",
                         nb_missing_antecedent);
}

} // namespace io
} // namespace biospring
