#include "IO/DonorAcceptorRuleReader.h"
#include "IO/LobeSpec.h"
#include "IO/ReduceRuleReader.h"

#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace io
{

void DonorAcceptorRuleReader::_parse_line(const std::string & line, size_t line_id)
{
    auto tokens = utils::string::split(line);

    // GROUP <name> <welldepth> <equilibrium> <width>: a named set of Morse
    // parameters. See HydrogenBondGroup for why a PAIR has to opt into it at
    // both ends.
    if (!tokens.empty() && tokens[0] == "GROUP")
    {
        if (tokens.size() != 5)
            logging::die("DonorAcceptorRuleReader: line %d: GROUP takes a name and three numbers "
                         "(welldepth kJ.mol-1, equilibrium A, width A-1), found %d tokens",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));
        HydrogenBondGroup group;
        group.name = tokens[1];
        if (!utils::string::from_string(group.welldepth, tokens[2]) ||
            !utils::string::from_string(group.equilibrium, tokens[3]) ||
            !utils::string::from_string(group.width, tokens[4]))
            logging::die("DonorAcceptorRuleReader: line %d: GROUP '%s' has a parameter that is not a number",
                         static_cast<int>(line_id), group.name.c_str());
        if (group.welldepth <= 0.0f || group.equilibrium <= 0.0f || group.width <= 0.0f)
            logging::die("DonorAcceptorRuleReader: line %d: GROUP '%s' needs all three parameters > 0",
                         static_cast<int>(line_id), group.name.c_str());
        for (const auto & known : _groups)
            if (known.name == group.name)
                logging::die("DonorAcceptorRuleReader: line %d: GROUP '%s' is declared twice",
                             static_cast<int>(line_id), group.name.c_str());
        _groups.push_back(group);
        return;
    }

    // A '@name' token says which GROUP this site belongs to. Taken out before
    // anything else, so the columns that follow keep their positions and a
    // table written before groups existed reads exactly the same.
    std::string group;
    for (size_t k = tokens.size(); k-- > 0;)
        if (!tokens[k].empty() && tokens[k][0] == '@')
        {
            if (!group.empty())
                logging::die("DonorAcceptorRuleReader: line %d: a site may name only one group",
                             static_cast<int>(line_id));
            group = tokens[k].substr(1);
            if (group.empty())
                logging::die("DonorAcceptorRuleReader: line %d: '@' names no group", static_cast<int>(line_id));
            tokens.erase(tokens.begin() + static_cast<long>(k));
        }

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
        entry.antecedent2 =
            parse_lobe_spec(tokens[5], entry.lobeMode, entry.lobeAngle, "DonorAcceptorRuleReader", line_id);

    entry.group = group;
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

void DonorAcceptorRuleReader::setNaming(const std::string & path)
{
    if (path.empty())
        return;
    reduce::ReduceRuleReader naming(path);
    naming.read();
    _translation = naming.rules();
    _hastranslation = true;
}

std::string DonorAcceptorRuleReader::_plainName(const spn::Particle & p) const
{
    if (_hastranslation)
    {
        // A one-atom rule IS a renaming: its name is the type, its single atom
        // the original. A rule with several atoms is a real coarse-grain grain
        // and has no plain equivalent, so it is left alone.
        for (const auto & rule : _translation.get_rules_for_residue(p.getResName()))
            if (rule.getName() == p.getName() && rule.getNumberOfAtoms() == 1)
                return *rule.getAtomNames().begin();
        return p.getName();
    }
    const std::string & n = p.getName();
    return n.size() > 1 ? n.substr(1) : n;
}

const DonorAcceptorRole * DonorAcceptorRuleReader::_roleFor(const spn::Particle & p) const
{
    const std::string plain = _plainName(p);
    // Four lookups, most specific first: the residue's own entry under the name
    // the topology carries, then under the plain name, then the wildcard's two.
    for (const std::string & res : {p.getResName(), std::string("*")})
        for (const std::string & atom : {p.getName(), plain})
        {
            const auto it = _roles.find({res, atom});
            if (it != _roles.end())
                return &it->second;
        }
    return nullptr;
}

void DonorAcceptorRuleReader::tagParticles(spn::SpringNetwork & spn) const
{
    // The named Morse parameters, in the order they were declared: a site's
    // group index is 1-based into this, and 0 means the force field's own.
    std::vector<spn::SpringNetwork::HydrogenBondParameters> groups;
    for (const HydrogenBondGroup & g : _groups)
        groups.push_back({g.welldepth, g.equilibrium, g.width});
    spn.setHydrogenBondGroups(groups);

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
        by_plain.emplace(std::make_tuple(p.getChainName(), p.getResId(), _plainName(p)), i);
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
        // The carried name first, then the plain one. BOTH are needed: a table
        // keyed by amber.grp's types writes this antecedent 'ACA', a table in
        // plain names writes it 'CA', and the same table has to work on a
        // reduced topology and an unreduced one. Trying only the carried name
        // silently left every unprefixed antecedent of a plain table
        // unresolved -- measured on the three-resolution ribosome: 1459
        // directional sites instead of 1548, with no error of any kind.
        const auto key = std::make_tuple(self.getChainName(), self.getResId(), name);
        const auto own = by_atom.find(key);
        if (own != by_atom.end())
            return static_cast<int>(own->second);
        const auto plain = by_plain.find(key);
        return plain == by_plain.end() ? -1 : static_cast<int>(plain->second);
    };

    for (unsigned i = 0; i < spn.getNumberOfParticles(); ++i)
    {
        spn::Particle & p = spn.getParticle(i);
        const DonorAcceptorRole * role = _roleFor(p);
        if (role == nullptr)
            continue;
        const DonorAcceptorRole & entry = *role;

        p.setDonorCapacity(entry.donorCapacity);
        p.setAcceptorCapacity(entry.acceptorCapacity);
        if (!entry.group.empty())
        {
            int index = 0;
            for (size_t g = 0; g < _groups.size(); ++g)
                if (_groups[g].name == entry.group)
                    index = static_cast<int>(g) + 1;
            if (index == 0)
                logging::die("DonorAcceptorRuleReader: %s:%s names group '%s', which no GROUP line declares",
                             p.getResName().c_str(), p.getName().c_str(), entry.group.c_str());
            p.setHydrogenBondGroup(index);
        }
        if (entry.donorCapacity > 0)
            ++nb_donors;
        if (entry.acceptorCapacity > 0)
            ++nb_acceptors;

        if (entry.antecedent.empty())
            continue;
        const int anc = resolve(p, entry.antecedent);
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
        const int anc2 = resolve(p, entry.antecedent2);
        if (anc2 >= 0)
        {
            p.setAntecedentIndex2(anc2);
            p.setLobes(entry.lobeMode, entry.lobeAngle);
            if (entry.lobeMode != spn::Particle::HBOND_LOBES_NONE)
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
