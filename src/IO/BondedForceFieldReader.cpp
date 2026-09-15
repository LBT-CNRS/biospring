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

        if (type == "DIHEDRAL")
    {
        // 11 tokens when the entry names its own axis, 9 when it does not:
        // every file written before the axis existed stays valid, and so does
        // every entry whose spring has a ghost at one end to ask instead.
        if (tokens.size() != 9 && tokens.size() != 11)
            logging::die("BondedForceFieldReader: line %d: DIHEDRAL expects 9 tokens (type name resname family "
                         "atom_ref atom_rotant d0 k dc_offset), or 11 with its axis (... axis_b axis_c), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));

        DihedralEntry entry;
        entry.resname = tokens[2];
        const std::string & family = tokens[3];
        unsigned family_id = spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT;
        for (unsigned f = 0; f < spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT; ++f)
            if (family == DIHEDRAL_FAMILY_KEYWORDS[f])
            {
                family_id = f;
                break;
            }
        if (family_id == spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT)
        {
            std::string expected;
            for (unsigned f = 0; f < spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT; ++f)
                expected += (f == 0 ? "" : ", ") + std::string(DIHEDRAL_FAMILY_KEYWORDS[f]);
            logging::die("BondedForceFieldReader: line %d: invalid DIHEDRAL family '%s' (expected one of %s)",
                         static_cast<int>(line_id), family.c_str(), expected.c_str());
        }
        entry.family = static_cast<DihedralFamily>(family_id);
        entry.atom_ref = tokens[4];
        entry.atom_rotant = tokens[5];
        if (!utils::string::from_string(entry.d0, tokens[6]))
            logging::die("BondedForceFieldReader: line %d: invalid d0 '%s'", static_cast<int>(line_id),
                         tokens[6].c_str());
        if (!utils::string::from_string(entry.k, tokens[7]))
            logging::die("BondedForceFieldReader: line %d: invalid k '%s'", static_cast<int>(line_id),
                         tokens[7].c_str());
        if (!utils::string::from_string(entry.dc_offset, tokens[8]))
            logging::die("BondedForceFieldReader: line %d: invalid dc_offset '%s'", static_cast<int>(line_id),
                         tokens[8].c_str());
        if (tokens.size() == 11)
        {
            entry.axis_b = tokens[9];
            entry.axis_c = tokens[10];
        }
        _dihedral.push_back(entry);
    }
    else if (type == "TORSIONTABLE")
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
    else if (type == "GHOSTPARTICLE")
    {
        if (tokens.size() != 9)
            logging::die("BondedForceFieldReader: line %d: GHOSTPARTICLE expects 9 tokens (type name resname atom_B "
                         "atom_C atom_ref r theta delta), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));

        GhostParticleEntry entry;
        entry.name = tokens[1];
        entry.resname = tokens[2];
        entry.atom_B = tokens[3];
        entry.atom_C = tokens[4];
        entry.atom_ref = tokens[5];
        if (!utils::string::from_string(entry.r, tokens[6]))
            logging::die("BondedForceFieldReader: line %d: invalid r '%s'", static_cast<int>(line_id),
                         tokens[6].c_str());
        if (!utils::string::from_string(entry.theta_deg, tokens[7]))
            logging::die("BondedForceFieldReader: line %d: invalid theta '%s'", static_cast<int>(line_id),
                         tokens[7].c_str());
        if (!utils::string::from_string(entry.delta_deg, tokens[8]))
            logging::die("BondedForceFieldReader: line %d: invalid delta '%s'", static_cast<int>(line_id),
                         tokens[8].c_str());
        _ghostparticles.push_back(entry);
    }
    else
    {
        logging::die("BondedForceFieldReader: line %d: unknown entry type '%s' (expected "
                     "GHOSTPARTICLE, DIHEDRAL, TORSION or TORSIONTABLE)",
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

    logging::info("BondedForceFieldReader: read %zu ghost particle, %zu dihedral rule(s) and %zu torsion(s) "
                 "over %zu table(s).",
                 _ghostparticles.size(), _dihedral.size(), _torsion.size(), _torsiontables.size());
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

topology::Spring * BondedForceFieldReader::_find_spring(topology::SpringCollection & collection,
                                                         const topology::Particle & p1,
                                                         const topology::Particle & p2) const
{
    if (!collection.exists(p1, p2))
        return nullptr;

    try
    {
        return &collection.at_uid(topology::Spring::generate_uid(p1, p2));
    }
    catch (const std::out_of_range &)
    {
        return &collection.at_uid(topology::Spring::generate_uid(p2, p1));
    }
}

bool BondedForceFieldReader::_existing_equilibrium(topology::Topology & topology, const topology::Particle & p1,
                                                    const topology::Particle & p2, double & equilibrium) const
{
    topology::Spring * existing = _find_spring(topology.springs(), p1, p2);
    if (existing == nullptr)
        return false;
    equilibrium = existing->equilibrium();
    return true;
}

unsigned BondedForceFieldReader::_create_ghost_particles(topology::Topology & topology,
                                                          std::vector<ResidueParticleIndices> & residues,
                                                          size_t index, const std::string & resname,
                                                          const reduce::ReduceRuleContainer * translation) const
{
    unsigned nb_created = 0;
    for (const GhostParticleEntry & entry : _ghostparticles)
    {
        if (entry.resname != resname)
            continue;

        topology::Particle * anchor_b = _resolve_atom(entry.atom_B, residues, index, topology, translation);
        topology::Particle * anchor_c = _resolve_atom(entry.atom_C, residues, index, topology, translation);
        topology::Particle * anchor_ref = _resolve_atom(entry.atom_ref, residues, index, topology, translation);
        if (anchor_b == nullptr || anchor_c == nullptr || anchor_ref == nullptr)
        {
            logging::warning("BondedForceFieldReader: GHOSTPARTICLE %s %s has an unresolved anchor (%s/%s/%s), "
                             "skipped.",
                             resname.c_str(), entry.name.c_str(), entry.atom_B.c_str(), entry.atom_C.c_str(),
                             entry.atom_ref.c_str());
            continue;
        }

        // Belongs to the same residue the GHOSTPARTICLE rule matched
        // (residues[index][0]'s own identity), not necessarily anchor_B's:
        // a ghost's anchors may include a +/- cross-residue atom, but the
        // ghost itself is a property of the residue the .bi.ff rule was
        // written for.
        const topology::ParticleProperties & residue_properties = topology.get_particle(residues[index][0]).properties();

        topology::ParticleProperties properties;
        properties.set_name(entry.name);
        properties.set_residue_name(residue_properties.residue_name());
        properties.set_chain_name(residue_properties.chain_name());
        properties.set_residue_id(residue_properties.residue_id());
        properties.set_static(true);
        properties.set_mass(0.0f);

        topology::Particle ghost_template(properties);
        // Every GHOSTPARTICLE line in a .bi.ff is a dihedral-ring ghost:
        // the rotated image of its reference atom.
        topology.add_ghost_particle(ghost_template, *anchor_b, *anchor_c, *anchor_ref, entry.r, entry.theta_deg,
                                    entry.delta_deg,
                                    static_cast<unsigned>(spn::GhostPlacement::AxisRotation));

        residues[index].push_back(topology.number_of_particles() - 1);
        nb_created++;
    }
    return nb_created;
}

void BondedForceFieldReader::_add_or_combine_dihedral_spring(topology::SpringCollection & collection,
                                                              topology::Particle & p1, topology::Particle & p2,
                                                              double equilibrium, double stiffness,
                                                              double dc_offset, topology::Particle * axis_b,
                                                              topology::Particle * axis_c) const
{
    // Two different Fourier-term ghost-spring groups for the same axis may
    // legitimately pick the same real substituent pair (see
    // doc/BondedForceFieldSprings.md, "Choosing substituents when fewer are
    // needed than exist" -- the choice of which real atoms a secondary term
    // uses is explicit per case, not guaranteed disjoint from the dominant
    // term's grid). Since both are quadratic in the same pair distance,
    // their combined contribution is *exactly* one equivalent spring:
    // 0.5*k1*(x-d1)^2 + 0.5*k2*(x-d2)^2 = 0.5*(k1+k2)*(x-d_combined)^2 + const,
    // with d_combined = (k1*d1 + k2*d2) / (k1+k2) -- so a collision is
    // combined here rather than rejected or silently overwritten.
    //
    // NOTE: that "+const" (= k1*k2/(2*(k1+k2)) * (d1-d2)^2, from completing
    // the square) is a genuine part of the two original springs' combined
    // energy that the merged single-spring representation cannot reproduce
    // on its own -- currently NOT folded into dc_offset below (unlike the
    // ring-construction artifact, which is). Harmless today: no current
    // DIHEDRAL axis ever produces two groups targeting the same real pair
    // (verified empirically, zero collisions across all generated data), so
    // this branch is never actually exercised -- but if a future axis ever
    // does hit it, the combined spring's reported energy would be short by
    // exactly that dropped constant.
    topology::Spring * existing = _find_spring(collection, p1, p2);

    if (existing == nullptr)
    {
        topology::Spring & added = collection.add_spring(p1, p2, equilibrium, stiffness);
        added.set_dc_offset(dc_offset);
        if (axis_b != nullptr && axis_c != nullptr)
            added.set_axis(axis_b->unique_id(), axis_c->unique_id());
    }
    else
    {
        const double k1 = existing->stiffness();
        const double k2 = stiffness;
        const double combined_equilibrium = (k1 * existing->equilibrium() + k2 * equilibrium) / (k1 + k2);
        existing->set_equilibrium(combined_equilibrium);
        existing->set_stiffness(k1 + k2);
        existing->set_dc_offset(existing->dc_offset() + dc_offset);
        // Both entries describe the same real pair, so they turn the same
        // bond: the first axis seen is the axis, and a later one only fills a
        // gap it left.
        if (!existing->has_axis() && axis_b != nullptr && axis_c != nullptr)
            existing->set_axis(axis_b->unique_id(), axis_c->unique_id());
    }
}

void BondedForceFieldReader::_retune_or_add_spring(topology::Topology & topology, topology::Particle & p1,
                                                    topology::Particle & p2, double equilibrium, double stiffness,
                                                    const char * kind) const
{
    topology::Spring * existing = _find_spring(topology.springs(), p1, p2);

    if (existing != nullptr)
    {
        existing->set_equilibrium(equilibrium);
        existing->set_stiffness(stiffness);
    }
    else
    {
        logging::info("BondedForceFieldReader: no pre-existing --rigidbody spring for %s pair %s-%s, adding a new "
                     "one.",
                     kind, p1.properties().name().c_str(), p2.properties().name().c_str());
        topology.add_spring(p1, p2, equilibrium, stiffness);
    }
}

size_t BondedForceFieldReader::countExpectedGhostParticles(const topology::Topology & topology) const
{
    const std::vector<ResidueParticleIndices> residues = _group_particles_by_residue(topology);
    size_t count = 0;
    for (const ResidueParticleIndices & residue : residues)
    {
        if (residue.empty())
            continue;
        const std::string resname = topology.get_particle(residue[0]).properties().residue_name();
        for (const GhostParticleEntry & entry : _ghostparticles)
            if (entry.resname == resname)
                count++;
    }
    return count;
}

void BondedForceFieldReader::buildSprings(topology::Topology & topology,
                                          const reduce::ReduceRuleContainer * translation,
                                          bool enableDihedralBackbone, bool enableDihedralSidechain,
                                          bool enableDihedralPlanarity) const
{
    if (translation != nullptr)
        _check_translation_is_one_atom_per_rule(*translation);

    // Non-const: _create_ghost_particles appends each newly-created ghost's
    // index to residues[index] so later DIHEDRAL entries in the same
    // residue can resolve it exactly like a real atom (see its own
    // comment in the header).
    std::vector<ResidueParticleIndices> residues = _group_particles_by_residue(topology);

    unsigned nb_dihedral_applied = 0;
    unsigned nb_ghostparticles_created = 0;
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

    // Only worth creating ghost particles at all if some dihedral family is
    // actually enabled -- otherwise every DIHEDRAL entry that would use them
    // is skipped anyway (see the family_enabled check below), and they'd
    // just be dead, unused particles sitting in the topology.
    const bool enableGhostParticles = enableDihedralBackbone || enableDihedralSidechain;

    for (size_t index = 0; index < residues.size(); index++)
    {
        const std::string resname = topology.get_particle(residues[index][0]).properties().residue_name();

        if (enableGhostParticles)
            nb_ghostparticles_created += _create_ghost_particles(topology, residues, index, resname, translation);

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

        for (const DihedralEntry & entry : _dihedral)
        {
            if (entry.resname != resname)
                continue;

            // PHI/PSI/OMEGA are all still gated by the single
            // enableDihedralBackbone flag at build time (see buildSprings's
            // own comment) -- they only gain independent control at runtime,
            // via SpringNetwork's dihedral.phi/psi/omega .msp settings.
            // Hence the indirection: several families share one build flag,
            // so the mapping is a table (DIHEDRAL_FAMILY_GATES) rather than
            // one flag per family.
            const bool build_gates[GATE_COUNT] = {enableDihedralBackbone, enableDihedralSidechain,
                                                  enableDihedralPlanarity};
            if (!build_gates[DIHEDRAL_FAMILY_GATES[entry.family]])
                continue;

            topology::Particle * p_ref = _resolve_atom(entry.atom_ref, residues, index, topology, translation);
            topology::Particle * p_rot = _resolve_atom(entry.atom_rotant, residues, index, topology, translation);
            if (p_ref == nullptr || p_rot == nullptr)
                continue;

            // A ghost spring is always a new addition (see DihedralEntry's
            // comment): never a retune of a real --rigidbody spring, so
            // this goes straight to the family's own collection instead of
            // _retune_or_add_spring -- combining in place if another
            // Fourier-term group already added a spring for this exact pair
            // (see _add_or_combine_dihedral_spring).
            // The axis, when the entry names one, resolved exactly like any
            // other atom so "+N"/"-C" work here too. A spring whose axis atom
            // is missing (a chain terminus) keeps no axis rather than a wrong
            // one: unfiltered is a degraded model, a wrong axis is a broken one.
            topology::Particle * p_axis_b = nullptr;
            topology::Particle * p_axis_c = nullptr;
            if (!entry.axis_b.empty() && !entry.axis_c.empty())
            {
                p_axis_b = _resolve_atom(entry.axis_b, residues, index, topology, translation);
                p_axis_c = _resolve_atom(entry.axis_c, residues, index, topology, translation);
            }
            _add_or_combine_dihedral_spring(topology.dihedral_springs(entry.family), *p_ref, *p_rot, entry.d0, entry.k,
                                            entry.dc_offset, p_axis_b, p_axis_c);
            nb_dihedral_applied++;
        }
    }

    logging::info("BondedForceFieldReader: applied %u dihedral spring(s), created %u ghost particle(s), "
                 "applied %u torsion(s) over %zu table(s).",
                 nb_dihedral_applied, nb_ghostparticles_created, nb_torsions_applied, _torsiontables.size());
}

} // namespace rigidbodygroup
} // namespace biospring
