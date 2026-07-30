#ifndef __IO_BONDEDFORCEFIELDREADER_H__
#define __IO_BONDEDFORCEFIELDREADER_H__

#include <string>
#include <vector>

#include "IO/ReaderBase.h"
#include "reduce/ReduceRuleContainer.hpp"
#include "topology.hpp"

namespace biospring
{
namespace rigidbodygroup
{

// One real 1-2 bond ("STRETCH" line): equilibrium distance + stiffness,
// translated from a real force field (e.g. AMBER's BOND parameters).
struct StretchEntry
{
    std::string resname, atom1, atom2;
    double r0;
    double k;
};

// One real valence angle ("BEND" line, vertex = atom2): converted at build
// time into a 1-3 distance spring (law of cosines + curvature-matching
// stiffness conversion), using the *already-retuned* 1-2 springs' current
// equilibrium for r12/r23 (see buildSprings: STRETCH runs first, so those
// springs already carry the real r0 by the time BEND needs them -- reading
// the live spring avoids having to guess which residue's STRETCH entry
// owns a cross-residue bond, which differs depending on which side of the
// peptide bond you look from).
struct BendEntry
{
    std::string resname, atom1, atom2, atom3;
    double theta0_deg;
    double k;
};

// Which physical axis a DIHEDRAL entry's ghost spring belongs to -- purely a
// grouping label (see DihedralEntry), used only so the springs it produces
// can be routed to the matching Topology collection and, from there,
// selectively applied at build time (see buildSprings's enableDihedral*
// parameters and -dihedralbackbone/-dihedralsidechain in pdb2spn-cli.cpp).
enum class DihedralFamily
{
    BACKBONE,  // proper dihedral, phi/psi
    SIDECHAIN, // proper dihedral, chi1-4
    PLANARITY  // improper dihedral, ring/guanidinium planarity
};

// One ghost spring ("DIHEDRAL" line) contributing to one Fourier term of one
// real AMBER proper or improper torsion. `atom_ref`/`atom_rotant` are a real
// substituent on each side of the rotation axis (never the axis atoms
// themselves) -- several DIHEDRAL lines make up one (M,N) ghost-spring
// group, several groups make up one fully-decomposed real dihedral (see
// doc/BondedForceFieldSprings.md for the full derivation). Unlike
// StretchEntry/BendEntry, a DihedralEntry never corresponds to a real 1-2
// bond: it is always a new spring between two atoms with no direct chemical
// bond between them.
struct DihedralEntry
{
    std::string resname, atom_ref, atom_rotant;
    DihedralFamily family;
    double d0;
    double k;
};

// Parses a .bi.ff file ("bonded interaction .force field" -- see the
// .nbi.ff/.bi.ff naming note in data/forcefield/*.nbi.ff): lines of
// "STRETCH <name> <resname> <atom1> <atom2> <r0_A> <k_kJ.mol-1.A-2>",
// "BEND <name> <resname> <atom1> <atom2> <atom3> <theta0_deg>
// <k_kJ.mol-1.rad-2>", or "DIHEDRAL <name> <resname>
// <BACKBONE|SIDECHAIN|PLANARITY> <atom_ref> <atom_rotant> <d0_A>
// <k_kJ.mol-1.A-2>". Atom names may carry a "+"/"-" prefix (CHARMM-style:
// next/previous residue), same convention as .rbody.
//
// TODO(ring planarity): ring groups in ProteinAtomRigidGroups.rbody
// (aromatic side chains, guanidinium, proline ring) also carry pairwise
// springs beyond the real 1-2/1-3 pairs -- the ones spanning further across
// the ring, which are what actually keep it planar (fixing every pairwise
// distance in a group removes all internal degrees of freedom, including
// out-of-plane pucker -- a distance-geometry fact, not a dedicated
// planarity term). Neither STRETCH nor BEND has a rule type for those, so
// they stay at the uniform --stiffness value until DIHEDRAL PLANARITY
// entries are generated for the residue in question (see
// doc/BondedForceFieldSprings.md, Section 3.2) -- Proline is a deliberate
// exception, its ring genuinely puckers in real AMBER (no improper term
// restrains it), so no PLANARITY entry should ever be generated for it.
//
// This is a pdb2spn-time mechanism (see -bondedinteraction/--bondedinteraction
// in pdb2spn-cli.cpp), always used together with -rigidbody/--rigidbody:
// RigidBodyBuilder must have already created a spring for every real bond
// and valence angle (see buildSprings), which this retunes in place with
// real parameters instead of the uniform "rigid" ones.
class BondedForceFieldReader : public ReaderBase
{
  public:
    BondedForceFieldReader() : ReaderBase() {}
    BondedForceFieldReader(const std::string & path) : ReaderBase(path) {}
    BondedForceFieldReader(const char * const path) : ReaderBase(path) {}

    void read();

    // Applies every stretching/bending parameter described by the file onto
    // `topology`.
    //
    // This is designed to always run on top of a --rigidbody-built topology:
    // RigidBodyBuilder already creates a spring for every real bond (pairs
    // involving a single-vertex group's own vertex) and every real valence
    // angle (the group's other pairs) -- that is exactly the topology this
    // reader needs, already validated to leave the right hinge axes free. So
    // the normal case is to *retune* the spring RigidBodyBuilder already
    // created for a given pair in place (topology::SpringCollection::exists
    // + at_uid), replacing its rigid/as-loaded parameters with the real
    // ones. Only if a pair named in the file has no existing spring yet
    // (fallback, not the expected case) is a new one added via
    // Topology::add_spring.
    //
    // A fallback-added spring is "virtual" in that it doesn't correspond to
    // any real chemical bond (the two atoms aren't directly bonded) -- but
    // it isn't arbitrary: it's the same 1-3 distance-spring proxy BEND
    // always uses for a real valence angle, just newly added instead of
    // retuning an existing --rigidbody one. The one case seen so far
    // (Proline's "-C"-N-CD angle: its ring closes onto N via CD instead of
    // an amide H, and P_PHI/P_RING in ProteinAtomRigidGroups.rbody don't
    // already connect "-C" to CD) carries real physical meaning once BEND
    // supplies the actual AMBER angle constant -- it's what realistically
    // restricts Proline's phi (ring pucker allows some variation, unlike a
    // free hinge). Deliberately NOT added directly to ProteinAtomRigidGroups.rbody
    // instead: with the uniform --stiffness there and no real angle
    // constant, it would freeze phi completely rather than restrict it
    // realistically -- worse than leaving phi free, since it would assert a
    // single arbitrary value with artificial confidence instead of just not
    // modelling the restriction.
    //
    // `translation` mirrors RigidBodyBuilder's own naming-translation table:
    // when -grp/--grp has renamed particles (e.g. CA -> ACA for Ala), the
    // .bi.ff file still uses original atom names, so this lets them be
    // resolved anyway. Must be an all-atom identity mapping (one atom per
    // rule, like amber.grp), same requirement as RigidBodyBuilder.
    //
    // DIHEDRAL entries are handled differently from STRETCH/BEND: a ghost
    // spring never corresponds to a real chemical bond, so it is always a
    // new addition (Topology::add_dihedral_backbone_spring /
    // add_dihedral_sidechain_spring / add_dihedral_planarity_spring,
    // depending on the entry's family), never a retune of an existing
    // --rigidbody spring.
    //
    // `enableStretch`/`enableBend`/`enableDihedralBackbone`/
    // `enableDihedralSidechain` independently select which categories of
    // the .bi.ff file are actually applied (see -stretching/-bending/
    // -dihedralbackbone/-dihedralsidechain in pdb2spn-cli.cpp, where
    // -dihedral is a pure convenience alias that sets both dihedral flags
    // at once -- this reader only ever sees the two resolved booleans) --
    // none are applied unless explicitly requested. This is a build-time
    // decision: once written to a .nc topology, springs are just springs,
    // with no separate per-category runtime switch in the .msp. Combining
    // either dihedral flag alone (without stretch/bend) gives the
    // "intermediate model" (see doc/BondedForceFieldSprings.md, Section 3.2
    // vs 3.1): bonds/angles stay at --rigidbody's uniform value, real
    // dihedral wells are layered on top anyway.
    //
    // PLANARITY entries (improper dihedrals) have no enabling flag yet --
    // not generated by the script at all currently, see doc/
    // BondedForceFieldSprings.md's Known Limitations -- so they are never
    // applied regardless of these two flags; a future -dihedralplanarity
    // flag/parameter should be added alongside that work, not before.
    void buildSprings(topology::Topology & topology, const reduce::ReduceRuleContainer * translation,
                      bool enableStretch, bool enableBend, bool enableDihedralBackbone,
                      bool enableDihedralSidechain) const;

  protected:
    std::vector<StretchEntry> _stretch;
    std::vector<BendEntry> _bend;
    std::vector<DihedralEntry> _dihedral;

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

    // Looks up the spring between `p1`/`p2` in `collection`, if any.
    // `SpringCollection::exists()` checks both orderings of the pair
    // (generate_uid(p1,p2) or generate_uid(p2,p1), whichever the spring was
    // originally created with), but at_uid() only looks up the exact key
    // it's given -- shared by every place that needs to find a
    // possibly-already-existing spring regardless of which order it was
    // created in.
    topology::Spring * _find_spring(topology::SpringCollection & collection, const topology::Particle & p1,
                                    const topology::Particle & p2) const;

    // Reads the current equilibrium of the spring already existing between
    // `p1`/`p2` into `equilibrium` and returns true, or returns false if no
    // such spring exists yet (STRETCH must run before BEND for this to find
    // every real 1-2 bond, including cross-residue ones -- see buildSprings).
    bool _existing_equilibrium(topology::Topology & topology, const topology::Particle & p1,
                                const topology::Particle & p2, double & equilibrium) const;

    // Retunes the existing --rigidbody spring between `p1`/`p2` in place
    // (the expected case), or adds a new one if none exists yet (fallback).
    void _retune_or_add_spring(topology::Topology & topology, topology::Particle & p1, topology::Particle & p2,
                                double equilibrium, double stiffness, const char * kind) const;

    // Adds a new dihedral ghost spring between `p1`/`p2` in `collection`, or
    // combines it with an already-existing one for that exact pair (two
    // different Fourier-term groups on the same axis may legitimately
    // target the same real substituent pair -- see the .cpp for the
    // combination formula).
    void _add_or_combine_dihedral_spring(topology::SpringCollection & collection, topology::Particle & p1,
                                         topology::Particle & p2, double equilibrium, double stiffness) const;
};

} // namespace rigidbodygroup
} // namespace biospring

#endif // __IO_BONDEDFORCEFIELDREADER_H__
