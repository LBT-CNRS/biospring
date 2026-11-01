#ifndef __IO_BONDEDFORCEFIELDREADER_H__
#define __IO_BONDEDFORCEFIELDREADER_H__

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "IO/ReaderBase.h"
#include "reduce/ReduceRuleContainer.hpp"
#include "topology.hpp"

namespace biospring
{
namespace rigidbodygroup
{

// Which family a TORSION entry belongs to -- a grouping label, used to
// route it at build time (see buildSprings's enableDihedral* parameters)
// and to switch it at runtime (see spn::SpringNetwork's per-family .msp
// settings). Aliases spn::SpringNetwork::DihedralFamilyIndex.
using DihedralFamily = spn::SpringNetwork::DihedralFamilyIndex;

// Where in its chain a residue has to sit for a record to apply to it.
//
// AMBER types a free chain end from its own residue template (NALA/CALA, not
// ALA), and the retyping is not cosmetic: the N-terminal nitrogen is an
// ammonium N3 rather than an amide N, and the C-terminal carbonyl oxygen is a
// delocalised carboxylate O2 rather than O. amber99sb then assigns NO torsion
// at all to N3-CT-C-N, HP-CT-C-O or H1-CT-C-O2, while assigning one to the
// identical interior quadruplets -- so a record carrying the interior
// parameters invents a term where the force field has none (measured on
// 071: a psi well 17.8 kJ.mol-1 deep and 20.3 kJ.mol-1.rad-1 of torque, on a
// residue whose psi AMBER leaves free), and the terms AMBER does add there,
// the ammonium rotor and the carboxylate improper, are missing.
//
// So the residue column of a TORSION record may carry a terminal qualifier,
// "ALA:NTER" or "ALA:CTER", and a residue is parameterised from the records
// matching its own position: the interior set is NOT also applied at a
// terminus. Which is why the generator has to make every terminal-sensitive
// torsion belong to the residue whose end it is -- see the note on
// `_wanted_position`.
enum TerminalPosition
{
    POSITION_INTERIOR = 0, // no qualifier: a residue with neighbours on both sides
    POSITION_NTER,         // ":NTER", a free N-terminus (NH3+), AMBER's NXXX
    POSITION_CTER          // ":CTER", a free C-terminus (COO-), AMBER's CXXX
};

// One AMBER torsion ("TORSION" line): the four atoms it acts on and the table
// that holds its energy and torque. Nothing is summed over substituent pairs,
// so no phase origin is chosen and no sign convention is guessed -- the record
// is AMBER's own quadruplet, and the table AMBER's own V(phi) sampled once.
struct TorsionEntry
{
    std::string resname;
    TerminalPosition position = POSITION_INTERIOR;
    DihedralFamily family;
    std::string atoms[4];
    unsigned table;
};

// One tabulated parameter set ("TORSIONTABLE" line): energy in kJ.mol-1 and
// torque (-dV/dphi) in kJ.mol-1.rad-1, sampled at bins + 1 points over
// [-pi, pi]. Torsions sharing a parameter set share a table -- the protein
// force field has 891 torsions and 30 tables.
struct TorsionTableEntry
{
    unsigned bins = 0;
    std::vector<float> energy;
    std::vector<float> torque;
};

// Parses a .bi.ff file ("bonded interaction .force field" -- see the
// .nbi.ff/.bi.ff naming note in data/forcefield/*.nbi.ff): AMBER's torsion
// terms, tabulated.
//
//   TORSIONTABLE <id> <bins> <E_0> <T_0> ...   over phi in [-180, 180],
//       E in kJ.mol-1 and T = -dV/dphi in kJ.mol-1.rad-1
//   TORSION <resname>[:NTER|:CTER] <FAMILY> <a1> <a2> <a3> <a4> <table>
//
// Atom names may carry a "+"/"-" prefix (CHARMM-style: next/previous
// residue), same convention as .rbody. The residue name may carry a terminal
// qualifier (see TerminalPosition): the residue name itself stays AMBER's
// own, like every other data file here -- amber.grp and .rbody both spell a
// terminal alanine "ALA" too.
//
// RigidBodyBuilder must have already created a spring for every real bond and
// valence angle. Those keep --rigidbody's uniform --stiffness; this reader
// adds no spring at all, only the torsions, which act on the four atoms they
// name through the exact four-atom gradient.

class BondedForceFieldReader : public ReaderBase
{
  public:
    BondedForceFieldReader() : ReaderBase() {}
    BondedForceFieldReader(const std::string & path) : ReaderBase(path) {}
    BondedForceFieldReader(const char * const path) : ReaderBase(path) {}

    /// Reads the current file, APPENDING to whatever earlier calls read. Call
    /// it once per .bi.ff, with setFileName() in between: the files are
    /// per-polymer-type (protein, RNA, DNA) and a mixed system needs several,
    /// but each numbers its tables from 0, so they must not be concatenated by
    /// hand -- this shifts each file's ids instead. See `_tablebase`.
    void read();

    // Applies every dihedral parameter described by the file onto
    // `topology`.
    //
    // This is designed to always run on top of a --rigidbody-built topology:
    // RigidBodyBuilder already creates a spring for every real bond (pairs
    // involving a single-vertex group's own vertex) and every real valence
    // angle (the group's other pairs) -- that is exactly the topology this
    // reader needs, already validated to leave the right hinge axes free.
    // Bonds and angles stay at --rigidbody's uniform --stiffness: this
    // reader adds dihedral wells on top of that mesh and nothing else.
    //
    // `translation` mirrors RigidBodyBuilder's own naming-translation table:
    // when -grp/--grp has renamed particles (e.g. CA -> ACA for Ala), the
    // .bi.ff file still uses original atom names, so this lets them be
    // resolved anyway. Must be an all-atom identity mapping (one atom per
    // rule, like amber.grp), same requirement as RigidBodyBuilder.
    //
    // `enableDihedralBackbone`/`enableDihedralSidechain` independently
    // select which families of the .bi.ff file are actually applied (see
    // -dihedralbackbone/-dihedralsidechain in pdb2spn-cli.cpp, where
    // -dihedral is a convenience alias setting them together) -- none are
    // applied unless explicitly requested. This is a build-time decision: a
    // family not requested here never gets a spring created at all, so it
    // costs nothing (no particle, no NetCDF entry) -- unrelated to
    // SpringNetwork's own dihedral.* .msp settings, a runtime debug on/off
    // for whichever families WERE built (see Configuration.hpp).
    //
    // `enableDihedralPlanarity` gates the PLANARITY entries (improper
    // dihedrals -- aromatic-ring/His hub planarity, generated since the
    // aromatic rings' rigid cliques were split into per-vertex hinge
    // groups; see -dihedralplanarity in pdb2spn-cli.cpp, folded into the
    // -dihedral convenience alias like the other dihedral families).
    void buildSprings(topology::Topology & topology, const reduce::ReduceRuleContainer * translation,
                      bool enableDihedralBackbone, bool enableDihedralSidechain,
                      bool enableDihedralPlanarity) const;

  protected:
    std::vector<TorsionEntry> _torsion;
    std::vector<TorsionTableEntry> _torsiontables;

    // What this file's table id 0 becomes, so several files can be read into
    // one reader. Each generated .bi.ff numbers its own tables from 0 -- the
    // protein, RNA and DNA files all do -- so reading a second one has to
    // shift its ids past the first's, and shift its TORSION records' table
    // column by the same amount. Set at the top of every read(); zero for the
    // first file, so a single-file read is byte-identical to before.
    size_t _tablebase = 0;

    using ResidueParticleIndices = std::vector<size_t>;

    void _parse_line(const std::string & line, size_t line_id);

    std::vector<ResidueParticleIndices> _group_particles_by_residue(const topology::Topology & topology) const;

    // Resolves a single (possibly +/- prefixed) atom name relative to the
    // residue at `index`. Returns nullptr if the neighbour residue does not
    // exist, belongs to a different chain, or does not contain that atom
    // (directly, or once translated through `translation`).
    topology::Particle * _resolve_atom(const std::string & atomname,
                                       const std::vector<ResidueParticleIndices> & residues, size_t index,
                                       topology::Topology & topology,
                                       const reduce::ReduceRuleContainer * translation) const;

    // Translates `atomname` for a residue named `resname` through
    // `translation` (e.g. "CA" -> "ACA" for ALA). Returns an empty string if
    // no translation is found.
    std::string _translate(const reduce::ReduceRuleContainer & translation, const std::string & resname,
                           const std::string & atomname) const;

    // Dies if `translation` contains a rule with more than one atom -- must
    // be an all-atom identity mapping, never a real coarse-grain reduction.
    void _check_translation_is_one_atom_per_rule(const reduce::ReduceRuleContainer & translation) const;

    // Which records a residue is entitled to, from where it sits in its chain
    // and what the file actually declares.
    //
    // The position is read off the chain the same way `_resolve_atom` reads
    // "+"/"-": a residue is a terminus when the neighbouring group in particle
    // order belongs to another chain, or does not exist. So a chain BROKEN in
    // the middle -- a gap in the numbering, a missing loop -- is not seen as a
    // terminus, exactly as a "+N" there still resolves to the next residue
    // present. That is a property of the whole file format, not of this
    // function.
    //
    // `declared` is the set of (resname, position) pairs the file carries, and
    // a position the file says nothing about falls back to the interior
    // records. Without that fallback a .bi.ff generated before the qualifier
    // existed -- the DNA and RNA files still are -- would silently lose every
    // torsion at every chain end instead of keeping today's behaviour.
    //
    // A residue that is both the first and the last of its chain (a free amino
    // acid) is treated as an N-terminus: AMBER has no template for one, so
    // there is nothing to generate the real combination from. The ammonium
    // rotor is then right, and so are both suppressions, since HP-CT-C-O and
    // H1-CT-C-O2 are equally absent from the force field; only the
    // carboxylate's own improper is missed.
    TerminalPosition _wanted_position(const std::string & resname, bool is_nterm, bool is_cterm,
                                      const std::set<std::pair<std::string, TerminalPosition>> & declared) const;

};

} // namespace rigidbodygroup
} // namespace biospring

#endif // __IO_BONDEDFORCEFIELDREADER_H__
