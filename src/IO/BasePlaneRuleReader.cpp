#include "IO/BasePlaneRuleReader.h"

#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace io
{

void BasePlaneRuleReader::_parse_line(const std::string & line, size_t line_id)
{
    const auto tokens = utils::string::split(line);
    if (tokens.size() != 4)
        logging::die("BasePlaneRuleReader: line %d: invalid number of tokens (expected 4, found %d) -- "
                     "'<resname> <atom0> <atom1> <atom2>'",
                     static_cast<int>(line_id), static_cast<int>(tokens.size()));

    // Three DISTINCT atoms, or the cross product is zero and the base has no
    // normal at all. Worth refusing here rather than at the first force
    // evaluation, where it would surface as a silently skipped residue.
    for (unsigned a = 1; a < 4; ++a)
        for (unsigned b = a + 1; b < 4; ++b)
            if (tokens[a] == tokens[b])
                logging::die("BasePlaneRuleReader: line %d: %s names '%s' twice; three atoms of a base plane "
                             "must be distinct",
                             static_cast<int>(line_id), tokens[0].c_str(), tokens[a].c_str());

    if (_planes.count(tokens[0]))
        logging::die("BasePlaneRuleReader: line %d: %s is declared twice", static_cast<int>(line_id),
                     tokens[0].c_str());
    _planes[tokens[0]] = {tokens[1], tokens[2], tokens[3]};
}

void BasePlaneRuleReader::read()
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

    logging::info("BasePlaneRuleReader: read %d base plane(s) from %s.", static_cast<int>(_planes.size()),
                  _filename.c_str());
}

void BasePlaneRuleReader::buildFrames(spn::SpringNetwork & spn) const
{
    // (chain, residue id, atom name) -> particle, so a rule's three atoms are
    // resolved inside their own residue. Built once, like
    // DonorAcceptorRuleReader::tagParticles does it, rather than rescanning
    // every particle per atom named.
    std::map<std::tuple<std::string, unsigned, std::string>, unsigned> by_atom;
    // The residues present, in particle order, so the frames come out in a
    // deterministic order whatever the map's iteration order is: the stacking
    // energy is a sum of floats and its value must not depend on that.
    std::vector<std::tuple<std::string, unsigned, std::string>> residues;
    std::set<std::tuple<std::string, unsigned, std::string>> seen;
    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        const spn::Particle & p = spn.getParticle(i);
        by_atom.emplace(std::make_tuple(p.getChainName(), p.getResId(), p.getName()), i);
        const auto residue = std::make_tuple(p.getChainName(), p.getResId(), p.getResName());
        if (seen.insert(residue).second)
            residues.push_back(residue);
    }

    std::vector<spn::SpringNetwork::BaseFrame> frames;
    // Which frame each PARTICLE's residue has, for a caller holding an atom
    // rather than a residue -- the planarity factor is asked about a donor or
    // an acceptor, which is not generally one of the three plane atoms.
    std::vector<int> frame_of_particle(spn.getNumberOfParticles(), -1);
    std::map<std::pair<std::string, unsigned>, int> frame_of_residue;
    unsigned nb_incomplete = 0;
    for (const auto & residue : residues)
    {
        const auto rule = _planes.find(std::get<2>(residue));
        if (rule == _planes.end())
            continue;

        spn::SpringNetwork::BaseFrame frame;
        bool complete = true;
        for (unsigned k = 0; k < 3; ++k)
        {
            const auto found =
                by_atom.find(std::make_tuple(std::get<0>(residue), std::get<1>(residue), rule->second[k]));
            if (found == by_atom.end())
            {
                complete = false;
                break;
            }
            frame.atoms[k] = found->second;
        }
        if (!complete)
        {
            nb_incomplete++;
            continue;
        }
        frame_of_residue[std::make_pair(std::get<0>(residue), std::get<1>(residue))] =
            static_cast<int>(frames.size());
        frames.push_back(frame);
    }

    // One pass, after the frames are known: scanning every particle per framed
    // residue instead would be 4757 x 100 000 comparisons on a ribosome.
    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        const spn::Particle & p = spn.getParticle(i);
        const auto found = frame_of_residue.find(std::make_pair(p.getChainName(), p.getResId()));
        if (found != frame_of_residue.end())
            frame_of_particle[i] = found->second;
    }

    if (nb_incomplete)
        logging::warning("BasePlaneRuleReader: %u residue(s) the table names are missing one of their three "
                         "plane atoms and take no part in the stacking term.",
                         nb_incomplete);
    logging::info("BasePlaneRuleReader: %d base frame(s) installed.", static_cast<int>(frames.size()));
    spn.setBaseFrames(frames);
    spn.setBaseFrameOfParticle(frame_of_particle);
}

} // namespace io
} // namespace biospring
