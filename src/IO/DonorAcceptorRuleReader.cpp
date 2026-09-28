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

    if (tokens.size() < 4 || tokens.size() > 6)
    {
        logging::die("DonorAcceptorRuleReader: line %d: invalid number of tokens (expected 4 to 6, found %d)",
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
    if (tokens.size() >= 5)
        entry.antecedent = tokens[4];
    // A second antecedent turns the direction into the bisector's opposite,
    // which is exact for a planar sp2 site with two heavy neighbours -- see
    // DonorAcceptorRole. A '~' prefix says the atom is a PLANE reference
    // instead, for a site with one heavy neighbour and two hydrogens or two
    // lone pairs, whose directions are the two lobes at +/- 62 degrees.
    if (tokens.size() == 6)
    {
        std::string a2 = tokens[5];
        const char marker = a2[0];
        if (marker == '~' || marker == '^')
        {
            entry.lobeMode = marker == '~' ? spn::Particle::HBOND_LOBES_IN_PLANE
                                           : spn::Particle::HBOND_LOBES_OUT_OF_PLANE;
            a2 = a2.substr(1);
            // An angle after a colon overrides the default, which is 62 degrees for
            // the in-plane form and nothing for the out-of-plane one -- a face
            // approach has no canonical angle, so the file has to say.
            const size_t colon = a2.find(':');
            if (colon != std::string::npos)
            {
                if (!utils::string::from_string(entry.lobeAngle, a2.substr(colon + 1)))
                    logging::die("DonorAcceptorRuleReader: line %d: invalid lobe angle '%s'",
                                 static_cast<int>(line_id), a2.substr(colon + 1).c_str());
                a2 = a2.substr(0, colon);
            }
            else if (marker == '^')
                logging::die("DonorAcceptorRuleReader: line %d: '^%s' needs an explicit angle, "
                             "as '^%s:105' -- an out-of-plane approach has no default",
                             static_cast<int>(line_id), a2.c_str(), a2.c_str());
            if (a2.empty())
                logging::die("DonorAcceptorRuleReader: line %d: '%c' names no plane atom",
                             static_cast<int>(line_id), marker);
            if (entry.lobeAngle <= 0.0f || entry.lobeAngle >= 180.0f)
                logging::die("DonorAcceptorRuleReader: line %d: lobe angle %g is outside (0, 180)",
                             static_cast<int>(line_id), static_cast<double>(entry.lobeAngle));
        }
        entry.antecedent2 = a2;
    }

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
    unsigned nb_bisector = 0;
    unsigned nb_lobes = 0;

    // (chain, residue id, atom name) -> particle index, so a rule's
    // antecedent column can be resolved inside its own residue. Built once:
    // the alternative is a rescan of every particle per antecedent.
    std::map<std::tuple<std::string, int, std::string>, unsigned> by_atom;
    // Second index, keyed on the atom name with its residue-code prefix
    // stripped: amber.grp names ALA's carbonyl AC and Arg's RC, so "the
    // previous residue's C" cannot be written as one fixed type. Dropping the
    // first character recovers the plain PDB name (AC -> C, ACA -> CA), which
    // is what a cross-residue antecedent names.
    std::map<std::tuple<std::string, int, std::string>, unsigned> by_plain;
    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        const spn::Particle & p = spn.getParticle(i);
        by_atom.emplace(std::make_tuple(p.getChainName(), p.getResId(), p.getName()), i);
        const std::string & n = p.getName();
        if (n.size() > 1)
            by_plain.emplace(std::make_tuple(p.getChainName(), p.getResId(), n.substr(1)), i);
    }

    // Resolves one antecedent name, which may carry a CHARMM-style "-"/"+"
    // prefix for the previous/next residue -- the same convention .rbody and
    // .bi.ff use. A prefixed name is matched on the plain atom name, since
    // the type of an atom in another residue depends on that residue.
    auto resolve = [&](const spn::Particle & self, const std::string & name) -> int {
        if (name.empty())
            return -1;
        const char prefix = name[0];
        if (prefix == '-' || prefix == '+')
        {
            const int rid = self.getResId() + (prefix == '+' ? 1 : -1);
            const auto f = by_plain.find(std::make_tuple(self.getChainName(), rid, name.substr(1)));
            return f == by_plain.end() ? -1 : static_cast<int>(f->second);
        }
        const auto f = by_atom.find(std::make_tuple(self.getChainName(), self.getResId(), name));
        return f == by_atom.end() ? -1 : static_cast<int>(f->second);
    };

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
        const int anc = resolve(p, it->second.antecedent);
        if (anc < 0)
        {
            ++nb_missing_antecedent;
            continue;
        }
        p.setAntecedentIndex(anc);
        ++nb_directed;

        // The second antecedent is optional and its absence is not an error:
        // a chain's first residue has no previous one, and the site simply
        // falls back to the single-antecedent direction.
        const int anc2 = resolve(p, it->second.antecedent2);
        if (anc2 >= 0)
        {
            p.setAntecedentIndex2(anc2);
            p.setLobes(it->second.lobeMode, it->second.lobeAngle);
            if (it->second.lobeMode != spn::Particle::HBOND_LOBES_NONE)
                ++nb_lobes;
            else
                ++nb_bisector;
        }
    }

    logging::info("Hydrogen bond tagging: %d donor(s), %d acceptor(s) among %d particles, %d with a resolved "
                  "antecedent (directional), %d of them with a bisector (two heavy neighbours, one direction), "
                  "%d with two lobes (one heavy neighbour, two directions).",
                  nb_donors, nb_acceptors, spn.getNumberOfParticles(), nb_directed, nb_bisector, nb_lobes);
    if (nb_missing_antecedent > 0)
        logging::warning("DonorAcceptorRuleReader: %d particle(s) name an antecedent absent from their residue -- "
                         "left undirected (distance-only), not an error.",
                         nb_missing_antecedent);
}

} // namespace io
} // namespace biospring
