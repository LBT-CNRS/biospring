#include "IO/BondedForceFieldReader.h"

#include <cmath>
#include <stdexcept>

#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace rigidbodygroup
{

void BondedForceFieldReader::_parse_line(const std::string & line, size_t line_id)
{
    const auto tokens = utils::string::split(line);

    if (tokens.empty())
        return;

    const std::string & type = tokens[0];

    if (type == "STRETCH")
    {
        if (tokens.size() != 7)
            logging::die("BondedForceFieldReader: line %d: STRETCH expects 7 tokens (type name resname atom1 atom2 "
                         "r0 k), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));

        StretchEntry entry;
        entry.resname = tokens[2];
        entry.atom1 = tokens[3];
        entry.atom2 = tokens[4];
        if (!utils::string::from_string(entry.r0, tokens[5]))
            logging::die("BondedForceFieldReader: line %d: invalid r0 '%s'", static_cast<int>(line_id),
                         tokens[5].c_str());
        if (!utils::string::from_string(entry.k, tokens[6]))
            logging::die("BondedForceFieldReader: line %d: invalid k '%s'", static_cast<int>(line_id),
                         tokens[6].c_str());
        _stretch.push_back(entry);
    }
    else if (type == "BEND")
    {
        if (tokens.size() != 8)
            logging::die("BondedForceFieldReader: line %d: BEND expects 8 tokens (type name resname atom1 atom2 "
                         "atom3 theta0 k), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));

        BendEntry entry;
        entry.resname = tokens[2];
        entry.atom1 = tokens[3];
        entry.atom2 = tokens[4];
        entry.atom3 = tokens[5];
        if (!utils::string::from_string(entry.theta0_deg, tokens[6]))
            logging::die("BondedForceFieldReader: line %d: invalid theta0 '%s'", static_cast<int>(line_id),
                         tokens[6].c_str());
        if (!utils::string::from_string(entry.k, tokens[7]))
            logging::die("BondedForceFieldReader: line %d: invalid k '%s'", static_cast<int>(line_id),
                         tokens[7].c_str());
        _bend.push_back(entry);
    }
    else if (type == "DIHEDRAL")
    {
        if (tokens.size() != 9)
            logging::die("BondedForceFieldReader: line %d: DIHEDRAL expects 9 tokens (type name resname family "
                         "atom_ref atom_rotant d0 k dc_offset), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));

        DihedralEntry entry;
        entry.resname = tokens[2];
        const std::string & family = tokens[3];
        if (family == "PHI")
            entry.family = DihedralFamily::PHI;
        else if (family == "PSI")
            entry.family = DihedralFamily::PSI;
        else if (family == "OMEGA")
            entry.family = DihedralFamily::OMEGA;
        else if (family == "SIDECHAIN")
            entry.family = DihedralFamily::SIDECHAIN;
        else if (family == "PLANARITY")
            entry.family = DihedralFamily::PLANARITY;
        else
            logging::die("BondedForceFieldReader: line %d: invalid DIHEDRAL family '%s' (expected PHI, PSI, "
                         "OMEGA, SIDECHAIN or PLANARITY)",
                         static_cast<int>(line_id), family.c_str());
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
        _dihedral.push_back(entry);
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
        logging::die("BondedForceFieldReader: line %d: unknown entry type '%s' (expected STRETCH, BEND, "
                     "GHOSTPARTICLE or DIHEDRAL)",
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

    logging::info("BondedForceFieldReader: read %zu stretching, %zu bending, %zu ghost particle and %zu dihedral "
                 "rule(s).",
                 _stretch.size(), _bend.size(), _ghostparticles.size(), _dihedral.size());
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
        topology.add_ghost_particle(ghost_template, *anchor_b, *anchor_c, *anchor_ref, entry.r, entry.theta_deg,
                                    entry.delta_deg);

        residues[index].push_back(topology.number_of_particles() - 1);
        nb_created++;
    }
    return nb_created;
}

void BondedForceFieldReader::_add_or_combine_dihedral_spring(topology::SpringCollection & collection,
                                                              topology::Particle & p1, topology::Particle & p2,
                                                              double equilibrium, double stiffness,
                                                              double dc_offset) const
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
        collection.add_spring(p1, p2, equilibrium, stiffness).set_dc_offset(dc_offset);
    }
    else
    {
        const double k1 = existing->stiffness();
        const double k2 = stiffness;
        const double combined_equilibrium = (k1 * existing->equilibrium() + k2 * equilibrium) / (k1 + k2);
        existing->set_equilibrium(combined_equilibrium);
        existing->set_stiffness(k1 + k2);
        existing->set_dc_offset(existing->dc_offset() + dc_offset);
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
        // Every applied BEND entry creates exactly 2 vertex-anchored ghosts
        // (see buildSprings) -- counted here too (a safe overcount, same as
        // the dihedral entries above: this doesn't know or care whether
        // -bending is actually enabled for this run).
        for (const BendEntry & entry : _bend)
            if (entry.resname == resname)
                count += 2;
    }
    return count;
}

void BondedForceFieldReader::buildSprings(topology::Topology & topology,
                                          const reduce::ReduceRuleContainer * translation, bool enableStretch,
                                          bool enableBend, bool enableDihedralBackbone,
                                          bool enableDihedralSidechain) const
{
    if (translation != nullptr)
        _check_translation_is_one_atom_per_rule(*translation);

    if (enableBend && !enableStretch)
        logging::die("BondedForceFieldReader: BEND requires STRETCH to be enabled too -- bending's 1-3 conversion "
                     "needs the real (AMBER) 1-2 bond lengths STRETCH provides; without them it would silently use "
                     "whichever arbitrary equilibrium --rigidbody happened to set instead, an inconsistent result.");

    // Non-const: _create_ghost_particles appends each newly-created ghost's
    // index to residues[index] so later DIHEDRAL entries in the same
    // residue can resolve it exactly like a real atom (see its own
    // comment in the header).
    std::vector<ResidueParticleIndices> residues = _group_particles_by_residue(topology);

    unsigned nb_stretch_applied = 0;
    unsigned nb_bend_applied = 0;
    unsigned nb_bend_skipped = 0;
    unsigned nb_dihedral_applied = 0;
    unsigned nb_ghostparticles_created = 0;

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

        if (enableStretch)
            for (const StretchEntry & entry : _stretch)
            {
                if (entry.resname != resname)
                    continue;

                topology::Particle * p1 = _resolve_atom(entry.atom1, residues, index, topology, translation);
                topology::Particle * p2 = _resolve_atom(entry.atom2, residues, index, topology, translation);
                if (p1 == nullptr || p2 == nullptr)
                    continue;

                // Zeroes (never removes) the existing --rigidbody spring's
                // stiffness for this real pair, but keeps its equilibrium
                // at the real AMBER r0: BEND's own _existing_equilibrium
                // lookup below still finds it there (only topology.springs()
                // is ever searched), and the real 1-2 pair stays a
                // spring-neighbour (nonbonded exclusion unaffected) even
                // though it no longer contributes any energy on its own.
                _retune_or_add_spring(topology, *p1, *p2, entry.r0, 0.0, "stretching");
                // The actual physics: a brand new spring, in its own
                // collection, between the same two real atoms -- see
                // _stretch_springs' own comment for why this is no longer a
                // plain in-place retune.
                // A real 1-2 bond at a chain boundary can legitimately be
                // described twice -- once by each flanking residue's own
                // STRETCH entries (e.g. a residue's own "C +N" rule and the
                // next residue's/cap's own "N -C" rule for that exact same
                // peptide bond, see ProteinAtomBonded.bi.ff's NME_N_C entry).
                // Both describe the same real AMBER bond with the same real
                // parameters, so the second occurrence is a genuine, expected
                // duplicate (unlike a plain naming collision) -- skipped here
                // rather than added again, which would otherwise either
                // throw (SpringCollection::add_spring forbids a duplicate
                // pair) or, if that check were bypassed, silently double the
                // effective stiffness for that one bond.
                if (_find_spring(topology.stretch_springs(), *p1, *p2) == nullptr)
                {
                    topology.add_stretch_spring(*p1, *p2, entry.r0, entry.k);
                    nb_stretch_applied++;
                }
            }

        if (enableBend)
            for (const BendEntry & entry : _bend)
            {
                if (entry.resname != resname)
                    continue;

                topology::Particle * p1 = _resolve_atom(entry.atom1, residues, index, topology, translation);
                topology::Particle * p2 = _resolve_atom(entry.atom2, residues, index, topology, translation);
                topology::Particle * p3 = _resolve_atom(entry.atom3, residues, index, topology, translation);
                if (p1 == nullptr || p2 == nullptr || p3 == nullptr)
                    continue;

                double r12 = 0.0;
                double r23 = 0.0;
                const bool has12 = _existing_equilibrium(topology, *p1, *p2, r12);
                const bool has23 = _existing_equilibrium(topology, *p2, *p3, r23);
                if (!has12 || !has23)
                {
                    logging::warning("BondedForceFieldReader: BEND %s %s-%s-%s has no existing 1-2 spring for "
                                     "r12/r23 (STRETCH must cover both real bonds first), skipped.",
                                     resname.c_str(), entry.atom1.c_str(), entry.atom2.c_str(), entry.atom3.c_str());
                    nb_bend_skipped++;
                    continue;
                }

                const double theta0 = entry.theta0_deg * M_PI / 180.0;

                const double r13sq = r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * std::cos(theta0);
                const double r13 = std::sqrt(r13sq);

                // k13 = K_theta * r13^2 / (r12 * r23 * sin(theta0))^2 --
                // matching curvatures (d^2E/dtheta^2) at equilibrium between
                // topology/Spring.hpp's 0.5*k13*(r13-r13_0)^2 and
                // 0.5*K_theta*(theta-theta0)^2, both already in the same
                // 0.5-prefactor convention (entry.k is meant to be sourced
                // directly from OpenMM's HarmonicAngleForce, which already
                // uses this convention -- same as STRETCH's k. No extra
                // factor of 2 here: that would only apply if entry.k were
                // instead sourced from a raw AMBER parm*.dat file, whose K
                // has no 1/2 baked in).
                //
                // This formula assumes r12/r23 stay fixed at their AMBER r0
                // -- true only at the calibration point, not during real
                // dynamics, where the real 1-2 bonds independently flex via
                // their own STRETCH springs. Connecting this spring directly
                // between the two real 1-3 atoms would let that real
                // bond-stretch leak into the measured r13 (found and
                // quantified via a real MD-trajectory validation: ~60%
                // systematic energy excess, ~95% of it from exactly this
                // leakage, only ~5% genuine linearization error) -- see
                // _bend_springs' own comment. Anchoring two ghosts at the
                // FIXED r12/r23 instead, each along the REAL, current bond
                // *direction*, removes it: this spring's own d0/k13 stay
                // exactly this formula, only what it connects to changes.
                const double denom = r12 * r23 * std::sin(theta0);
                const double k13 = entry.k * r13sq / (denom * denom);

                const topology::ParticleProperties & residue_properties =
                    topology.get_particle(residues[index][0]).properties();
                topology::ParticleProperties ghost_properties;
                ghost_properties.set_residue_name(residue_properties.residue_name());
                ghost_properties.set_chain_name(residue_properties.chain_name());
                ghost_properties.set_residue_id(residue_properties.residue_id());
                ghost_properties.set_static(true);
                ghost_properties.set_mass(0.0f);

                ghost_properties.set_name("GB" + entry.atom1);
                topology::Particle ghost_a_template(ghost_properties);
                // theta_deg=0: places the ghost exactly at B + r*normalize(C-B)
                // (the Ref-dependent azimuthal term is multiplied by
                // sin(0)=0, so p3 here only needs to be non-degenerate, not
                // meaningful on its own -- see spn::GhostParticle::computePosition).
                topology::Particle & ghost_a =
                    topology.add_ghost_particle(ghost_a_template, *p2, *p1, *p3, r12, 0.0, 0.0);

                ghost_properties.set_name("GB" + entry.atom3);
                topology::Particle ghost_c_template(ghost_properties);
                topology::Particle & ghost_c =
                    topology.add_ghost_particle(ghost_c_template, *p2, *p3, *p1, r23, 0.0, 0.0);

                // Same zero-out-but-keep treatment as STRETCH above, for the
                // same reason (nonbonded exclusion for the real 1-3 pair
                // must survive even though the real spring no longer carries
                // any of the actual restraint).
                _retune_or_add_spring(topology, *p1, *p3, r13, 0.0, "bending");
                topology.add_bend_spring(ghost_a, ghost_c, r13, k13);
                nb_bend_applied++;
                nb_ghostparticles_created += 2;
            }

        for (const DihedralEntry & entry : _dihedral)
        {
            if (entry.resname != resname)
                continue;

            // PLANARITY has no enabling flag yet (see header comment) --
            // never applied until that future work adds one. PHI/PSI/OMEGA
            // are all still gated by the single enableDihedralBackbone flag
            // at build time (see buildSprings's own comment) -- they only
            // gain independent control at runtime, via SpringNetwork's
            // dihedral.phi/psi/omega .msp settings.
            const bool is_backbone = entry.family == DihedralFamily::PHI || entry.family == DihedralFamily::PSI ||
                                     entry.family == DihedralFamily::OMEGA;
            const bool family_enabled = (is_backbone && enableDihedralBackbone) ||
                                        (entry.family == DihedralFamily::SIDECHAIN && enableDihedralSidechain);
            if (!family_enabled)
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
            switch (entry.family)
            {
            case DihedralFamily::PHI:
                _add_or_combine_dihedral_spring(topology.dihedral_phi_springs(), *p_ref, *p_rot, entry.d0, entry.k,
                                                entry.dc_offset);
                break;
            case DihedralFamily::PSI:
                _add_or_combine_dihedral_spring(topology.dihedral_psi_springs(), *p_ref, *p_rot, entry.d0, entry.k,
                                                entry.dc_offset);
                break;
            case DihedralFamily::OMEGA:
                _add_or_combine_dihedral_spring(topology.dihedral_omega_springs(), *p_ref, *p_rot, entry.d0, entry.k,
                                                entry.dc_offset);
                break;
            case DihedralFamily::SIDECHAIN:
                _add_or_combine_dihedral_spring(topology.dihedral_sidechain_springs(), *p_ref, *p_rot, entry.d0,
                                                entry.k, entry.dc_offset);
                break;
            case DihedralFamily::PLANARITY:
                _add_or_combine_dihedral_spring(topology.dihedral_planarity_springs(), *p_ref, *p_rot, entry.d0,
                                                entry.k, entry.dc_offset);
                break;
            }
            nb_dihedral_applied++;
        }
    }

    logging::info("BondedForceFieldReader: applied %u stretching, %u bending and %u dihedral spring(s), created %u "
                 "ghost particle(s) (%u bending rule(s) skipped for missing STRETCH cross-reference).",
                 nb_stretch_applied, nb_bend_applied, nb_dihedral_applied, nb_ghostparticles_created, nb_bend_skipped);
}

} // namespace rigidbodygroup
} // namespace biospring
