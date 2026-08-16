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
// grouping label (see DihedralEntry), used so the springs it produces can be
// routed to the matching Topology collection and, from there, selectively
// applied at build time (see buildSprings's enableDihedral* parameters and
// -dihedralbackbone/-dihedralsidechain in pdb2spn-cli.cpp) AND selectively
// enabled at runtime (see SpringNetwork's dihedral.phi/psi/omega/chi .msp
// settings). PHI/PSI/OMEGA were a single BACKBONE value until each got its
// own runtime .msp toggle (still built together, under the one
// -dihedralbackbone flag -- see buildSprings -- only their runtime
// enable/disable is independent).
// The list of families is spn::SpringNetwork's own DihedralFamilyIndex --
// the index the springs are stored under, all the way from the collection
// this reader fills to the .nc variable they are written to -- so a family
// is declared once and this reader only adds its .bi.ff spelling (see
// DIHEDRAL_FAMILY_KEYWORDS in the .cpp).
using DihedralFamily = spn::SpringNetwork::DihedralFamilyIndex;

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
    // This spring's share of its axis's exact dihedral-energy correction
    // (ring construction artifact minus AMBER's own real DC -- see
    // scripts/generate_bonded_forcefield.py's calibrate_ring/
    // emit_ghost_ring), in kJ.mol-1. Summed by buildSprings across every
    // DIHEDRAL entry actually applied and subtracted only when *reporting*
    // dihedral energy (SpringNetwork::getDihedralEnergy) -- never affects
    // forces, a constant has no gradient.
    double dc_offset;
};

// One massless virtual-site ("GHOSTPARTICLE" line): a ghost particle that
// does not correspond to any real PDB atom, placed algebraically from 3
// real anchor atoms (atom_B, atom_C -- the dihedral axis -- and atom_ref,
// which fixes the azimuthal reference direction) and 3 calibrated
// parameters (r, theta_deg, delta_deg -- see spn::GhostParticle for the
// placement formula). A DihedralEntry's atom_ref/atom_rotant may name a
// GhostParticleEntry's own `name` instead of a real atom name, for a
// ghost-to-ghost dihedral spring -- resolved against the ghost particles
// already created for the same residue (see buildSprings).
struct GhostParticleEntry
{
    std::string resname, name, atom_B, atom_C, atom_ref;
    double r;
    double theta_deg;
    double delta_deg;
};

// Parses a .bi.ff file ("bonded interaction .force field" -- see the
// .nbi.ff/.bi.ff naming note in data/forcefield/*.nbi.ff): lines of
// "STRETCH <name> <resname> <atom1> <atom2> <r0_A> <k_kJ.mol-1.A-2>",
// "BEND <name> <resname> <atom1> <atom2> <atom3> <theta0_deg>
// <k_kJ.mol-1.rad-2>", "GHOSTPARTICLE <name> <resname> <atom_B> <atom_C>
// <atom_ref> <r_A> <theta_deg> <delta_deg>" (a massless virtual site, see
// GhostParticleEntry), or "DIHEDRAL <name> <resname>
// <PHI|PSI|OMEGA|SIDECHAIN|PLANARITY> <atom_ref> <atom_rotant> <d0_A>
// <k_kJ.mol-1.A-2> <dc_offset_kJ.mol-1>" (atom_ref/atom_rotant may each name
// either a real atom or a GHOSTPARTICLE defined earlier in the same
// residue; dc_offset is this spring's share of its axis's exact
// dihedral-energy correction, see DihedralEntry). Atom names may carry a
// "+"/"-" prefix (CHARMM-style: next/previous residue), same convention as
// .rbody.
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
// and valence angle (see buildSprings). STRETCH and BEND each zero that
// existing --rigidbody spring's stiffness (kept, not removed, so nonbonded
// exclusion is unaffected) and add the real restraint as a brand new spring
// in its own dedicated collection instead of retuning in place -- STRETCH's
// new spring still connects the same two real atoms; BEND's connects two
// new ghost particles anchored to the vertex (see BendEntry's own comment
// for why a plain real-atom 1-3 spring is not accurate enough). This keeps
// each family's energy separately reportable (see
// spn::SpringNetwork::Energies) instead of folding everything into one
// generic "spring energy" that only means something for a plain ENM/
// rigid-body network in the first place.
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
    // reader needs, already validated to leave the right hinge axes free.
    // STRETCH and BEND never retune that existing --rigidbody spring's
    // parameters to the real ones in place any more: they zero its
    // stiffness instead (via _retune_or_add_spring, always called with
    // stiffness=0.0 -- kept, not removed, so the pair stays a
    // spring-neighbour for nonbonded exclusion) while still setting its
    // equilibrium to the real AMBER r0/r13, since BEND's own
    // _existing_equilibrium lookup only ever searches topology.springs()
    // and would otherwise silently read back whichever value --rigidbody
    // happened to set. The actual physics is added as a brand new spring in
    // a dedicated collection (Topology::add_stretch_spring /
    // add_bend_spring): STRETCH's connects the same two real atoms; BEND's
    // connects two new ghost particles anchored to the vertex along each
    // real bond direction, rescaled to the real bond length (see
    // BendEntry's own comment for why a plain real-atom 1-3 spring leaks
    // bond-stretch into the angle restraint). If a pair named in the file
    // has no existing --rigidbody spring at all (not the expected case),
    // _retune_or_add_spring adds a zero-stiffness placeholder instead of
    // retuning one, purely so the equilibrium/exclusion bookkeeping above
    // still has somewhere to live -- the real spring is added the same way
    // either case. The one instance seen so far (Proline's "-C"-N-CD angle:
    // its ring closes onto N via CD instead of an amide H, and P_PHI/P_RING
    // in ProteinAtomRigidGroups.rbody don't already connect "-C" to CD)
    // carries real physical meaning once BEND supplies the actual AMBER
    // angle constant -- it's what realistically restricts Proline's phi
    // (ring pucker allows some variation, unlike a free hinge). Deliberately
    // NOT added directly to ProteinAtomRigidGroups.rbody instead: with the
    // uniform --stiffness there and no real angle constant, it would freeze
    // phi completely rather than restrict it realistically -- worse than
    // leaving phi free, since it would assert a single arbitrary value with
    // artificial confidence instead of just not modelling the restriction.
    //
    // `translation` mirrors RigidBodyBuilder's own naming-translation table:
    // when -grp/--grp has renamed particles (e.g. CA -> ACA for Ala), the
    // .bi.ff file still uses original atom names, so this lets them be
    // resolved anyway. Must be an all-atom identity mapping (one atom per
    // rule, like amber.grp), same requirement as RigidBodyBuilder.
    //
    // DIHEDRAL entries are handled differently from STRETCH/BEND: a ghost
    // spring never corresponds to a real chemical bond, so it is always a
    // new addition (into Topology::dihedral_springs(entry.family)), never a
    // retune of an existing --rigidbody spring.
    //
    // `enableStretch`/`enableBend`/`enableDihedralBackbone`/
    // `enableDihedralSidechain` independently select which categories of
    // the .bi.ff file are actually applied (see -stretching/-bending/
    // -dihedralbackbone/-dihedralsidechain in pdb2spn-cli.cpp, where
    // -dihedral is a pure convenience alias that sets both dihedral flags
    // at once -- this reader only ever sees the two resolved booleans) --
    // none are applied unless explicitly requested. This is a build-time
    // decision: a family not requested here never gets a spring created at
    // all, so it costs nothing (no particle, no NetCDF entry) -- unrelated
    // to SpringNetwork's own dihedral.phi/psi/omega/chi/bending .msp
    // settings, a runtime debug on/off for whichever families WERE built
    // (see Configuration.hpp's own comment). Combining either dihedral
    // flag alone (without stretch/bend) gives the
    // "intermediate model" (see doc/BondedForceFieldSprings.md, Section 3.2
    // vs 3.1): bonds/angles stay at --rigidbody's uniform value, real
    // dihedral wells are layered on top anyway.
    //
    // `enableDihedralPlanarity` gates the PLANARITY entries (improper
    // dihedrals -- aromatic-ring/His hub planarity, generated since the
    // aromatic rings' rigid cliques were split into per-vertex hinge
    // groups; see -dihedralplanarity in pdb2spn-cli.cpp, folded into the
    // -dihedral convenience alias like the other dihedral families).
    void buildSprings(topology::Topology & topology, const reduce::ReduceRuleContainer * translation,
                      bool enableStretch, bool enableBend, bool enableDihedralBackbone,
                      bool enableDihedralSidechain, bool enableDihedralPlanarity) const;

    // Upper bound on how many ghost particles buildSprings will create for
    // `topology` (a residue-name match only, not a full anchor resolution
    // -- some may still be skipped for an unresolved anchor, e.g. a
    // chain-terminus residue, so this can overcount, never undercount).
    // MUST be called (and the result reserved via
    // Topology::reserve_particles) before any spring exists anywhere in
    // `topology` -- see Topology::reserve_particles's own comment for why:
    // ghost particle creation grows the particle vector, and any
    // reallocation silently invalidates every Spring's Particle&
    // reference created so far (that includes --rigidbody's springs,
    // built before -bondedinteraction ever runs -- see pdb2spn-cli.cpp's
    // call site, placed before RigidBodyBuilder for exactly this reason).
    size_t countExpectedGhostParticles(const topology::Topology & topology) const;

  protected:
    std::vector<StretchEntry> _stretch;
    std::vector<BendEntry> _bend;
    std::vector<DihedralEntry> _dihedral;
    std::vector<GhostParticleEntry> _ghostparticles;

    using ResidueParticleIndices = std::vector<size_t>;

    void _parse_line(const std::string & line, size_t line_id);

    std::vector<ResidueParticleIndices> _group_particles_by_residue(const topology::Topology & topology) const;

    // Creates the topology::Particle ghosts described by every
    // GhostParticleEntry matching `resname`, resolving their 3 anchors the
    // same way STRETCH/BEND/DIHEDRAL do (via _resolve_atom, so the same
    // +/- cross-residue convention applies). Each newly-created ghost's
    // index is appended to `residues[index]` so later DIHEDRAL entries in
    // the same residue can resolve a ghost particle's name exactly like a
    // real atom's, via the ordinary _resolve_atom lookup -- no separate
    // resolution path is needed. residues is intentionally non-const:
    // creating a particle changes topology.number_of_particles(), which is
    // exactly why this mutates the local grouping instead of the
    // once-computed groups staying accurate on their own.
    // Returns the number of ghost particles actually created (for
    // buildSprings's summary log -- see the .cpp).
    unsigned _create_ghost_particles(topology::Topology & topology, std::vector<ResidueParticleIndices> & residues,
                                     size_t index, const std::string & resname,
                                     const reduce::ReduceRuleContainer * translation) const;

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

    // Sets equilibrium/stiffness on the existing --rigidbody spring between
    // `p1`/`p2` in place (the expected case), or adds a new one if none
    // exists yet (fallback) -- generic in `stiffness`, but STRETCH/BEND both
    // only ever call this with stiffness=0.0 now, purely to zero out the old
    // spring's contribution while keeping its equilibrium at the real AMBER
    // value (see buildSprings's own comment for why).
    void _retune_or_add_spring(topology::Topology & topology, topology::Particle & p1, topology::Particle & p2,
                                double equilibrium, double stiffness, const char * kind) const;

    // Adds a new dihedral ghost spring between `p1`/`p2` in `collection`, or
    // combines it with an already-existing one for that exact pair (two
    // different Fourier-term groups on the same axis may legitimately
    // target the same real substituent pair -- see the .cpp for the
    // combination formula). `dc_offset` is this entry's own share of its
    // axis's dihedral-energy correction (see DihedralEntry): set directly
    // on a new spring, summed into an existing one's on combination (energy
    // contributions are additive, see topology::Spring::_dc_offset).
    void _add_or_combine_dihedral_spring(topology::SpringCollection & collection, topology::Particle & p1,
                                         topology::Particle & p2, double equilibrium, double stiffness,
                                         double dc_offset) const;
};

} // namespace rigidbodygroup
} // namespace biospring

#endif // __IO_BONDEDFORCEFIELDREADER_H__
