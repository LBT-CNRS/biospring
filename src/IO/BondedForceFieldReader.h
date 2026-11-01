#ifndef __IO_BONDEDFORCEFIELDREADER_H__
#define __IO_BONDEDFORCEFIELDREADER_H__

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "IO/ReaderBase.h"
#include "reduce/ReduceRuleContainer.hpp"
#include "topology.hpp"
#include "topology/TerminalPosition.hpp"

namespace biospring
{
namespace rigidbodygroup
{

// Which family a TORSION entry belongs to -- a grouping label, used to
// route it at build time (see buildSprings's enableDihedral* parameters)
// and to switch it at runtime (see spn::SpringNetwork's per-family .msp
// settings). Aliases spn::SpringNetwork::DihedralFamilyIndex.
using DihedralFamily = spn::SpringNetwork::DihedralFamilyIndex;

// Where in its chain a residue has to sit for a record to apply to it, and the
// qualifier that spells it. Shared with the .grp reader, which carries the same
// qualifier with the same semantics -- see topology/TerminalPosition.hpp, where
// both the rule and the reason for it live.
using topology::TerminalPosition;
using topology::POSITION_CTER;
using topology::POSITION_INTERIOR;
using topology::POSITION_NTER;

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
};

} // namespace rigidbodygroup
} // namespace biospring

#endif // __IO_BONDEDFORCEFIELDREADER_H__
