#include "IO/BondedForceFieldReader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>

#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace rigidbodygroup
{

// The .bi.ff keyword naming each dihedral family, in DihedralFamily (i.e.
// spn::SpringNetwork::DihedralFamilyIndex) order. This is the file format's
// own spelling of the family list, so it lives with the reader rather than
// with the enum -- SpringNetwork's DIHEDRAL_FAMILY_NAMES are a third
// spelling of the same list, the .nc variable prefixes.
static constexpr const char * DIHEDRAL_FAMILY_KEYWORDS[spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT] = {
    "PHI",              // proper dihedral, backbone phi
    "PSI",              // proper dihedral, backbone psi
    "OMEGA",            // proper dihedral, backbone omega (peptide bond)
    "SIDECHAIN",        // proper dihedral, chi1-4
    "PLANARITY",        // improper dihedral, ring/guanidinium planarity
    "NUCLEIC_BACKBONE", // proper dihedral, alpha..zeta
    "NUCLEIC_CHI",      // proper dihedral, glycosidic torsion + RNA's 2'-OH rotor
    "NUCLEIC_SUGAR"     // proper dihedral, the four furanose ring bonds (the pucker)
};

// The build-time flag each family answers to (see buildSprings's
// enableDihedral* parameters and -dihedralbackbone/-dihedralsidechain/
// -dihedralplanarity in pdb2spn-cli.cpp). Not one flag per family: PHI/PSI/
// OMEGA are built together under -dihedralbackbone and only become
// independent at runtime, via their own .msp toggles.
enum DihedralBuildGate
{
    GATE_BACKBONE = 0,
    GATE_SIDECHAIN,
    GATE_PLANARITY,
    GATE_COUNT
};
static constexpr DihedralBuildGate DIHEDRAL_FAMILY_GATES[spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT] = {
    GATE_BACKBONE,  // PHI
    GATE_BACKBONE,  // PSI
    GATE_BACKBONE,  // OMEGA
    GATE_SIDECHAIN, // SIDECHAIN
    GATE_PLANARITY, // PLANARITY
    // The nucleic families split along the same line: what lies on the
    // chain itself is backbone, what hangs off it is sidechain. The
    // furanose counts as backbone -- the chain runs THROUGH the ring
    // (C4'-C3'), so building the phosphate torsions without the ring bonds
    // that carry delta would leave the backbone half-restrained. They stay
    // a separate FAMILY (own energy channel, own runtime toggle, since the
    // pucker has to be measurable on its own) without needing a separate
    // build flag.
    GATE_BACKBONE,  // NUCLEIC_BACKBONE
    GATE_SIDECHAIN, // NUCLEIC_CHI -- the base's own orientation, plus the 2'-OH rotor
    GATE_BACKBONE   // NUCLEIC_SUGAR
};

void BondedForceFieldReader::_parse_line(const std::string & line, size_t line_id)
{
    const auto tokens = utils::string::split(line);

    if (tokens.empty())
        return;

    const std::string & type = tokens[0];

    if (type == "TORSIONTABLE")
    {
        // <id> <bins> then 2*(bins+1) numbers, energy and torque interleaved.
        if (tokens.size() < 3)
            logging::die("BondedForceFieldReader: line %d: TORSIONTABLE expects an id and a bin count",
                         static_cast<int>(line_id));
        unsigned id = 0, bins = 0;
        if (!utils::string::from_string(id, tokens[1]) || !utils::string::from_string(bins, tokens[2]) || bins == 0)
            logging::die("BondedForceFieldReader: line %d: TORSIONTABLE has an invalid id or bin count",
                         static_cast<int>(line_id));
        if (tokens.size() != 3 + 2 * (bins + 1))
            logging::die("BondedForceFieldReader: line %d: TORSIONTABLE %u declares %u bins, so %d values were "
                         "expected and %d found",
                         static_cast<int>(line_id), id, bins, static_cast<int>(2 * (bins + 1)),
                         static_cast<int>(tokens.size()) - 3);
        if (_torsiontables.size() <= id)
            _torsiontables.resize(id + 1);
        TorsionTableEntry & tab = _torsiontables[id];
        tab.bins = bins;
        tab.energy.resize(bins + 1);
        tab.torque.resize(bins + 1);
        for (unsigned b = 0; b <= bins; ++b)
            if (!utils::string::from_string(tab.energy[b], tokens[3 + 2 * b]) ||
                !utils::string::from_string(tab.torque[b], tokens[4 + 2 * b]))
                logging::die("BondedForceFieldReader: line %d: TORSIONTABLE %u has a non-numeric sample at %u",
                             static_cast<int>(line_id), id, b);
    }
    else if (type == "TORSION")
    {
        if (tokens.size() != 8)
            logging::die("BondedForceFieldReader: line %d: TORSION expects 8 tokens (type resname family "
                         "atom1 atom2 atom3 atom4 table), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));
        TorsionEntry entry;
        entry.resname = tokens[1];
        unsigned family_id = spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT;
        for (unsigned f = 0; f < spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT; ++f)
            if (tokens[2] == DIHEDRAL_FAMILY_KEYWORDS[f])
                family_id = f;
        if (family_id == spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT)
            logging::die("BondedForceFieldReader: line %d: invalid TORSION family '%s'", static_cast<int>(line_id),
                         tokens[2].c_str());
        entry.family = static_cast<DihedralFamily>(family_id);
        for (unsigned k = 0; k < 4; ++k)
            entry.atoms[k] = tokens[3 + k];
        if (!utils::string::from_string(entry.table, tokens[7]))
            logging::die("BondedForceFieldReader: line %d: invalid TORSION table index '%s'",
                         static_cast<int>(line_id), tokens[7].c_str());
        _torsion.push_back(entry);
    }
    else
    {
        logging::die("BondedForceFieldReader: line %d: unknown entry type '%s' (expected "
                     "TORSION or TORSIONTABLE)",
                     static_cast<int>(line_id), type.c_str());
    }
}

void BondedForceFieldReader::read()
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
            _parse_line(buffer, line_id);
    }
    close();

    logging::info("BondedForceFieldReader: read %zu torsion(s) over %zu table(s).", _torsion.size(),
                  _torsiontables.size());
}

std::vector<BondedForceFieldReader::ResidueParticleIndices>
BondedForceFieldReader::_group_particles_by_residue(const topology::Topology & topology) const
{
    std::vector<ResidueParticleIndices> residues;
    if (topology.number_of_particles() == 0)
        return residues;

    ResidueParticleIndices current = {0};

    for (size_t i = 1; i < topology.number_of_particles(); i++)
    {
        const auto & previous_properties = topology.get_particle(current.back()).properties();
        const auto & properties = topology.get_particle(i).properties();

        const bool same_residue =
            properties.residue_id() == previous_properties.residue_id() &&
            properties.chain_name() == previous_properties.chain_name();

        if (same_residue)
            current.push_back(i);
        else
        {
            residues.push_back(current);
            current = {i};
        }
    }
    residues.push_back(current);

    return residues;
}

topology::Particle * BondedForceFieldReader::_resolve_atom(const std::string & atomname,
                                                           const std::vector<ResidueParticleIndices> & residues,
                                                           size_t index, topology::Topology & topology,
                                                           const reduce::ReduceRuleContainer * translation) const
{
    std::string name = atomname;
    size_t target_index = index;

    if (!name.empty() && (name[0] == '+' || name[0] == '-'))
    {
        const bool is_next = name[0] == '+';
        name = name.substr(1);

        if (is_next && index + 1 >= residues.size())
            return nullptr;
        if (!is_next && index == 0)
            return nullptr;

        target_index = is_next ? index + 1 : index - 1;

        const std::string & this_chain = topology.get_particle(residues[index][0]).properties().chain_name();
        const std::string & other_chain = topology.get_particle(residues[target_index][0]).properties().chain_name();
        if (this_chain != other_chain)
            return nullptr;
    }

    for (size_t particle_index : residues[target_index])
    {
        topology::Particle & p = topology.get_particle(particle_index);
        if (p.properties().name() == name)
            return &p;
    }

    // Not found under its original name: if a reduction renamed the
    // particles (e.g. amber.grp: CA -> ACA for ALA), try the translated name.
    if (translation != nullptr)
    {
        const std::string & target_resname =
            topology.get_particle(residues[target_index][0]).properties().residue_name();
        const std::string translated = _translate(*translation, target_resname, name);

        if (!translated.empty())
            for (size_t particle_index : residues[target_index])
            {
                topology::Particle & p = topology.get_particle(particle_index);
                if (p.properties().name() == translated)
                    return &p;
            }
    }

    return nullptr;
}

std::string BondedForceFieldReader::_translate(const reduce::ReduceRuleContainer & translation,
                                               const std::string & resname, const std::string & atomname) const
{
    reduce::ReduceRuleContainer rules = translation.get_rules_for_residue(resname);

    for (const auto & rule : rules)
        if (rule.hasAtomNamed(atomname))
            return rule.getName();

    return "";
}

void BondedForceFieldReader::_check_translation_is_one_atom_per_rule(
    const reduce::ReduceRuleContainer & translation) const
{
    for (const auto & rule : translation)
        if (rule.number_of_atoms() != 1)
            logging::die("BondedForceFieldReader: naming-translation file must have exactly one atom per rule (rule "
                         "'%s' for residue '%s' has %zu) -- it must be an all-atom identity mapping (like "
                         "amber.grp), not a coarse-grain reduction.",
                         rule.getName().c_str(), rule.getResidueName().c_str(), rule.number_of_atoms());
}

void BondedForceFieldReader::buildSprings(topology::Topology & topology,
                                          const reduce::ReduceRuleContainer * translation,
                                          bool enableDihedralBackbone, bool enableDihedralSidechain,
                                          bool enableDihedralPlanarity) const
{
    if (translation != nullptr)
        _check_translation_is_one_atom_per_rule(*translation);

    std::vector<ResidueParticleIndices> residues = _group_particles_by_residue(topology);

    unsigned nb_torsions_applied = 0;

    // The tables travel with the torsions and are given once, before any of
    // them: a torsion is useless without the curve it indexes.
    if (!_torsiontables.empty())
    {
        std::vector<topology::Topology::TorsionTable> tabs;
        tabs.reserve(_torsiontables.size());
        for (const TorsionTableEntry & e : _torsiontables)
        {
            topology::Topology::TorsionTable t;
            t.bins = e.bins;
            t.energy = e.energy;
            t.torque = e.torque;
            tabs.push_back(std::move(t));
        }
        topology.set_torsion_tables(std::move(tabs));
    }

    // A torsion spanning two residues is written under each of them, so the
    // same four atoms arrive twice; applying it twice would double its torque.
    // Keyed on the quadruplet either way round, a torsion and its reverse
    // being one torsion.
    std::set<std::array<size_t, 4>> seen_torsions;

    for (size_t index = 0; index < residues.size(); index++)
    {
        const std::string resname = topology.get_particle(residues[index][0]).properties().residue_name();

        for (const TorsionEntry & entry : _torsion)
        {
            if (entry.resname != resname)
                continue;
            const bool gates[GATE_COUNT] = {enableDihedralBackbone, enableDihedralSidechain, enableDihedralPlanarity};
            if (!gates[DIHEDRAL_FAMILY_GATES[entry.family]])
                continue;

            topology::Particle * p[4] = {nullptr, nullptr, nullptr, nullptr};
            bool resolved = true;
            for (unsigned k = 0; k < 4 && resolved; ++k)
            {
                p[k] = _resolve_atom(entry.atoms[k], residues, index, topology, translation);
                resolved = p[k] != nullptr;
            }
            if (!resolved)
                continue; // a chain terminus: the torsion simply does not exist there

            std::array<size_t, 4> q{p[0]->unique_id(), p[1]->unique_id(), p[2]->unique_id(), p[3]->unique_id()};
            const std::array<size_t, 4> rev{q[3], q[2], q[1], q[0]};
            if (rev < q)
                q = rev;
            if (!seen_torsions.insert(q).second)
                continue;

            topology.add_torsion(entry.family, *p[0], *p[1], *p[2], *p[3], entry.table);
            nb_torsions_applied++;
        }
    }

    logging::info("BondedForceFieldReader: applied %u torsion(s) over %zu table(s).", nb_torsions_applied,
                  _torsiontables.size());
}

} // namespace rigidbodygroup
} // namespace biospring
