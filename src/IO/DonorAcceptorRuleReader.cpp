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

    if (tokens.size() != 4)
    {
        logging::die("DonorAcceptorRuleReader: line %d: invalid number of tokens (expected 4, found %d)",
                     static_cast<int>(line_id), static_cast<int>(tokens.size()));
    }

    const std::string & resname = tokens[0];
    const std::string & atomname = tokens[1];

    DonorAcceptorRole entry;
    if (!utils::string::from_string(entry.isDonor, tokens[2]))
        logging::die("DonorAcceptorRuleReader: line %d: invalid donor flag '%s' (expected 0 or 1)",
                     static_cast<int>(line_id), tokens[2].c_str());
    if (!utils::string::from_string(entry.isAcceptor, tokens[3]))
        logging::die("DonorAcceptorRuleReader: line %d: invalid acceptor flag '%s' (expected 0 or 1)",
                     static_cast<int>(line_id), tokens[3].c_str());

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

    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        spn::Particle & p = spn.getParticle(i);
        const auto it = _roles.find({p.getResName(), p.getName()});
        if (it == _roles.end())
            continue;

        if (it->second.isDonor)
        {
            p.setDonor(true);
            ++nb_donors;
        }
        if (it->second.isAcceptor)
        {
            p.setAcceptor(true);
            ++nb_acceptors;
        }
    }

    logging::info("Hydrogen bond donor/acceptor tagging: %d donor(s), %d acceptor(s) found among %d particles.",
                  nb_donors, nb_acceptors, spn.getNumberOfParticles());
}

} // namespace io
} // namespace biospring
