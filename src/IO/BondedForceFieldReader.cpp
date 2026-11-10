#include "IO/BondedForceFieldReader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

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
        const size_t shifted = _tablebase + id;
        if (_torsiontables.size() <= shifted)
            _torsiontables.resize(shifted + 1);
        TorsionTableEntry & tab = _torsiontables[shifted];
        // A SECOND TORSIONTABLE WITH THE SAME ID IS NOT AN OVERRIDE, IT IS A
        // CONCATENATION ACCIDENT, AND IT USED TO BE SILENT. Each generated
        // file numbers its tables from 0 -- ProteinAtomTorsion.bi.ff has 32
        // (0..31), RNAAtomTorsion.bi.ff 39, DNAAtomTorsion.bi.ff 44 -- so
        // `cat` of any two of them collides on every id the smaller one uses,
        // and the entry that happens to come last wins. Measured on
        // 078.Ribosome70S's r70.bi.ff, which is protein + RNA concatenated:
        // 71 TORSIONTABLE lines, 39 distinct ids, 32 collisions, and so ALL
        // 894 protein records read a nucleic table instead of their own.
        //
        //   arginine's guanidinium chi, NE-CZ-NH1-HH11
        //     amber99sb X-CA-N2-X   n=2, k = 10.0416 kJ.mol-1, phase 180
        //     protein table 11      n=2, k = 10.0416, phase 180   correct
        //     RNA table 11          n=2, k = 41.8399, phase 180   4.17x
        //
        // That axis carries 4 such paths, so it ended up 348.7 kJ.mol-1 deep
        // where AMBER asks 168.2, and being both the stiffest torsion of the
        // network and one acting on 1.008 Da hydrogens, it is what threw them
        // out: 079's production run at 0.25 fs lost its arginine HH atoms
        // between step 400 and 450 and a spring reached 5e8 A. The mis-tabled
        // records are not even detectable from the energies -- the sum looks
        // plausible, 23897.7 kJ.mol-1 against the correct 13613.1 on one
        // 1476-torsion fragment.
        //
        // Nothing legitimate redeclares a table: the generators emit each id
        // once. A mixed system passes its files to -bondedinteraction
        // SEPARATELY instead, and read() shifts each one's ids past the
        // previous file's -- so this only ever fires within a single file,
        // which is where it means a real mistake.
        if (tab.bins != 0)
            logging::die("BondedForceFieldReader: line %d: TORSIONTABLE %u is declared twice. Each generated "
                         ".bi.ff numbers its tables from 0, so concatenating two of them collides -- pass them "
                         "to -bondedinteraction as separate files instead",
                         static_cast<int>(line_id), id);
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
        // "<resname>" or "<resname>:NTER|CTER". The residue name stays AMBER's
        // own in both cases (see TerminalPosition).
        if (!topology::parse_residue_field(tokens[1], entry.resname, entry.position))
            logging::die("BondedForceFieldReader: line %d: '%s' is not a residue name with at most one terminal "
                         "qualifier (expected '<resname>', '<resname>:NTER' or '<resname>:CTER')",
                         static_cast<int>(line_id), tokens[1].c_str());
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
        // Same shift as the TORSIONTABLE above, so a record keeps naming its
        // own file's table whichever position that file was read in.
        entry.table += static_cast<unsigned>(_tablebase);
        _torsion.push_back(entry);
    }
    else if (type == "STIFFNESS")
    {
        if (tokens.size() != 5)
            logging::die("BondedForceFieldReader: line %d: STIFFNESS expects 5 tokens (type resname "
                         "atom1 atom2 k), found %d",
                         static_cast<int>(line_id), static_cast<int>(tokens.size()));
        StiffnessEntry entry;
        if (!topology::parse_residue_field(tokens[1], entry.resname, entry.position))
            logging::die("BondedForceFieldReader: line %d: '%s' is not a residue name with at most one terminal "
                         "qualifier (expected '<resname>', '<resname>:NTER' or '<resname>:CTER')",
                         static_cast<int>(line_id), tokens[1].c_str());
        entry.atoms[0] = tokens[2];
        entry.atoms[1] = tokens[3];
        if (!utils::string::from_string(entry.stiffness, tokens[4]) || entry.stiffness <= 0.0f)
            logging::die("BondedForceFieldReader: line %d: invalid STIFFNESS constant '%s' (expected a positive "
                         "number of kJ.mol-1.A-2)",
                         static_cast<int>(line_id), tokens[4].c_str());
        _stiffness.push_back(entry);
    }
    else
    {
        logging::die("BondedForceFieldReader: line %d: unknown entry type '%s' (expected "
                     "TORSION, TORSIONTABLE or STIFFNESS)",
                     static_cast<int>(line_id), type.c_str());
    }
}

void BondedForceFieldReader::read()
{
    safeOpen();

    _tablebase = _torsiontables.size();
    const size_t firsttorsion = _torsion.size();

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

    // A record pointing at a table this file never declared used to land on a
    // default-constructed entry, which has no samples at all. With several
    // files it would be worse than that: the index is shifted by _tablebase,
    // so an id one past the end reaches into the NEXT file's tables and a
    // negative-looking one into the previous file's. Checking here is what
    // makes the shift safe.
    for (size_t i = firsttorsion; i < _torsion.size(); ++i)
        if (_torsion[i].table >= _torsiontables.size() || _torsiontables[_torsion[i].table].bins == 0)
            logging::die("BondedForceFieldReader: %s: TORSION %s %s names table %zu, which the file never "
                         "declares",
                         getFileName().c_str(), _torsion[i].resname.c_str(), _torsion[i].atoms[0].c_str(),
                         _torsion[i].table - _tablebase);

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
                                          bool enableDihedralPlanarity, BondStiffness bondStiffness) const
{
    if (translation != nullptr)
        _check_translation_is_one_atom_per_rule(*translation);

    std::vector<ResidueParticleIndices> residues = _group_particles_by_residue(topology);

    unsigned nb_torsions_applied = 0;
    // STIFFNESS retunes, it does not create: these count what it found and
    // what it did not, so a record set that does not match the mesh is
    // reported rather than silently half-applied.
    unsigned nb_stiffness_applied = 0;
    unsigned nb_stiffness_unsprung = 0;
    unsigned nb_stiffness_skipped_hydrogen = 0;
    std::set<size_t> retuned;
    double stiffness_min = 0.0;
    double stiffness_max = 0.0;

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

    // Which (residue, position) pairs the file actually carries, so a position
    // it says nothing about can fall back to the interior records.
    std::set<std::pair<std::string, TerminalPosition>> declared;
    for (const TorsionEntry & entry : _torsion)
        declared.insert({entry.resname, entry.position});
    for (const StiffnessEntry & entry : _stiffness)
        declared.insert({entry.resname, entry.position});

    // Chain ends that took the interior records, reported once rather than per
    // residue, and split by cause -- the two mean very different things. On a
    // .bi.ff with no qualifier at all, every chain end is undeclared, which is
    // the previous behaviour and not worth a line each. A chain end that IS
    // declared and still falls back is a carved-out fragment, which is a fact
    // about the structure the reader is the only place to notice.
    unsigned nb_undeclared = 0;
    unsigned nb_carved = 0;

    for (size_t index = 0; index < residues.size(); index++)
    {
        const std::string resname = topology.get_particle(residues[index][0]).properties().residue_name();
        const std::string & chain = topology.get_particle(residues[index][0]).properties().chain_name();
        const bool is_nterm =
            index == 0 || topology.get_particle(residues[index - 1][0]).properties().chain_name() != chain;
        const bool is_cterm = index + 1 == residues.size() ||
                              topology.get_particle(residues[index + 1][0]).properties().chain_name() != chain;

        // And being at the end of a chain is not the same as being a free
        // terminus: a carved-out fragment -- 056's lone methionine, 4 of 077's
        // chain ends -- has neither the ammonium nor the carboxylate, and its
        // backbone torsions are the interior ones. The residue's own atoms say
        // which it is; see topology::is_nterminal_only_atom. A reduction has
        // renamed them by now (H1 -> AH1), so the ORIGINAL name is the one to
        // ask, and it is set by exactly the reductions this reader accepts --
        // one atom per rule, hence one original per grain.
        bool has_nterminal_atom = false;
        bool has_cterminal_atom = false;
        for (size_t particle_index : residues[index])
        {
            const auto & properties = topology.get_particle(particle_index).properties();
            const std::string atom =
                properties.original_name().empty() ? properties.name() : properties.original_name();
            has_nterminal_atom |= topology::is_nterminal_only_atom(atom);
            has_cterminal_atom |= topology::is_cterminal_only_atom(atom);
        }

        const TerminalPosition position = topology::wanted_position(resname, is_nterm, is_cterm, has_nterminal_atom,
                                                                    has_cterminal_atom, declared);
        if ((is_nterm || is_cterm) && position == POSITION_INTERIOR)
        {
            // Reaching the interior set at a declared end can only mean the
            // atom test said fragment: wanted_position would have returned the
            // terminal one otherwise.
            const bool declared_here = (is_nterm && declared.count({resname, POSITION_NTER}) != 0) ||
                                       (is_cterm && declared.count({resname, POSITION_CTER}) != 0);
            if (declared_here)
                nb_carved++;
            else
                nb_undeclared++;
        }

        for (const TorsionEntry & entry : _torsion)
        {
            if (entry.resname != resname || entry.position != position)
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

        // The bond force constants. Nothing is created: the spring is already
        // there, and this replaces the uniform value with AMBER's own. A pair
        // with no spring means the mesh and these records disagree about the
        // topology, which is worth reporting and never worth fixing by adding
        // a bond RigidBodyBuilder chose not to make.
        if (bondStiffness == BondStiffness::NONE)
            continue;
        for (const StiffnessEntry & entry : _stiffness)
        {
            if (entry.resname != resname || entry.position != position)
                continue;
            topology::Particle * p1 = _resolve_atom(entry.atoms[0], residues, index, topology, translation);
            topology::Particle * p2 = _resolve_atom(entry.atoms[1], residues, index, topology, translation);
            if (p1 == nullptr || p2 == nullptr)
                continue; // a chain end: that bond does not exist here
            // A hydrogen's bond is the shortest period in the mesh, so
            // retuning it is what forces the timestep down; HEAVY leaves it
            // at --stiffness on purpose. See BondStiffness.
            if (bondStiffness == BondStiffness::HEAVY &&
                (p1->properties().mass() < 2.0f || p2->properties().mass() < 2.0f))
            {
                nb_stiffness_skipped_hydrogen++;
                continue;
            }
            // Either order: a spring is stored under whichever pair created
            // it, and the record names the bond the way the force field does.
            auto & by_uid = topology.springs().by_uid();
            auto it = by_uid.find(topology::Spring::generate_uid(*p1, *p2));
            if (it == by_uid.end())
                it = by_uid.find(topology::Spring::generate_uid(*p2, *p1));
            if (it == by_uid.end())
            {
                nb_stiffness_unsprung++;
                continue;
            }
            // A bond across the chain link is written under both residues, so
            // the same spring arrives twice with the same constant. Counting
            // SPRINGS rather than records is what makes the reported number
            // mean "this many springs now carry AMBER's value".
            if (!retuned.insert(it->second).second)
                continue;
            if (nb_stiffness_applied == 0 || entry.stiffness < stiffness_min)
                stiffness_min = entry.stiffness;
            if (nb_stiffness_applied == 0 || entry.stiffness > stiffness_max)
                stiffness_max = entry.stiffness;
            topology.springs()[it->second].set_stiffness(entry.stiffness);
            nb_stiffness_applied++;
        }
    }

    logging::info("BondedForceFieldReader: applied %u torsion(s) over %zu table(s).", nb_torsions_applied,
                  _torsiontables.size());
    if (nb_undeclared != 0)
        logging::info("BondedForceFieldReader: %u chain end(s) kept the interior torsion parameters -- the file "
                      "declares no :NTER/:CTER records for them.",
                      nb_undeclared);
    if (nb_carved != 0)
        logging::info("BondedForceFieldReader: %u chain end(s) kept the interior torsion parameters although the "
                      "file does declare them -- those residues carry neither the ammonium hydrogens nor OXT, so "
                      "they are cuts through a longer chain, not free ends.",
                      nb_carved);
    if (bondStiffness != BondStiffness::NONE)
    {
        logging::info("BondedForceFieldReader: retuned %u spring(s) to AMBER's own bond stiffness, %.1f to %.1f "
                      "kJ.mol-1.A-2. No spring was created: the mesh already carries every bond.",
                      nb_stiffness_applied, stiffness_min, stiffness_max);
        if (nb_stiffness_skipped_hydrogen != 0)
            logging::info("BondedForceFieldReader: %u bond(s) to a hydrogen left at the uniform stiffness "
                          "(-bondstiffness heavy). Those are the shortest periods in the mesh, so retuning them "
                          "is what would force the timestep down.",
                          nb_stiffness_skipped_hydrogen);
        if (nb_stiffness_unsprung != 0)
            logging::warning("BondedForceFieldReader: %u bond(s) named by a STIFFNESS record carry no spring, so "
                             "their constant was not applied. The mesh and the bonded force field disagree about "
                             "the topology -- check that -rigidbody and -bondedinteraction describe the same "
                             "chemistry.",
                             nb_stiffness_unsprung);
    }
}

} // namespace rigidbodygroup
} // namespace biospring
