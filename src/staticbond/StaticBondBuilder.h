#ifndef __STATIC_BOND_BUILDER_H__
#define __STATIC_BOND_BUILDER_H__

#include <map>
#include <string>
#include <utility>

#include "topology/Topology.hpp"

namespace biospring
{
namespace staticbond
{

// Hydrogen bonds and disulfide bridges as ORDINARY SPRINGS, found once in the
// input structure and never re-paired.
//
// This is a deliberate step back from a dynamic hydrogen-bond term. Such a
// term can form and break bonds during the simulation, but it does not
// reproduce folding either, and it costs a mechanism that sits outside the
// spring model the rest of BioSpring is built on. A spring calibrated to a
// hydrogen bond's force constant delivers the same service -- holding a
// structure that is already folded -- inside the model that is already there.
//
// What is given up, and it should be said plainly: a spring never breaks. It
// holds a bond that a real molecule would open, and it cannot form one that
// was not in the input. These options are for keeping a known structure
// together, not for studying association.

// Force constant of a hydrogen-bond spring, in kJ.mol-1.A-2.
//
// Matched to the curvature of the Morse well the dynamic model used, at its
// minimum: V = De[1 - exp(-a(r-re))]^2 has V''(re) = 2.De.a^2, and BioSpring's
// spring energy is 0.5.k.dr^2, so k = 2.De.a^2 with De = 16.7 kJ.mol-1 and
// a = 1.34 A-1. The well DEPTH has no harmonic equivalent -- a spring is
// infinitely deep -- so the curvature is what can be matched, and it is what
// sets the force for the small displacements a held structure actually makes.
extern const float HYDROGEN_BOND_STIFFNESS;

// Donor/acceptor heavy-atom distance below which a hydrogen bond is taken to
// exist, in Angstrom. The usual crystallographic criterion (HBPLUS, DSSP).
extern const float HYDROGEN_BOND_CUTOFF;

// Force constant and search radius of a disulfide spring, in kJ.mol-1.A-2 and
// Angstrom. AMBER ff14SB's S-S bond is r0 = 2.038 A with k = 1389.1
// kJ.mol-1.A-2 in the 0.5.k.dr^2 convention, which is BioSpring's own, so the
// value transfers with no conversion. The search radius is generous next to
// 2.038 A: a real bridge is unambiguous, and anything near 2.5 A that is not
// one would have to be two sulfurs almost touching.
extern const float DISULFIDE_STIFFNESS;
extern const float DISULFIDE_CUTOFF;

// Residues i and j closer than this in sequence are not considered for a
// hydrogen bond. An i/i+1 donor-acceptor contact is a consequence of the
// backbone's own geometry, not a hydrogen bond, and springing it would just
// stiffen the chain.
extern const int HYDROGEN_BOND_MIN_SEQUENCE_GAP;

// (residue name, particle name) -> (donor capacity, acceptor capacity).
using DonorAcceptorTable = std::map<std::pair<std::string, std::string>, std::pair<int, int>>;

// Reads a .hbond table, keeping only what a static spring needs: which
// particles can donate and which can accept. The antecedent columns a dynamic
// term uses to orient the bond are parsed and ignored -- a spring pulls along
// the line between two atoms and has no direction of its own.
//
// Throws std::runtime_error if the file cannot be opened.
DonorAcceptorTable readDonorAcceptorTable(const std::string & path);

// Adds one spring per donor/acceptor pair within `cutoff`. A pair that already
// carries a spring is skipped rather than doubled, so these options compose
// with -cutoff, with -rigidbody, and with a network that has neither.
//
// Returns the number of springs added.
size_t addHydrogenBondSprings(topology::Topology & topology, const DonorAcceptorTable & table, float cutoff,
                              float stiffness);

// Adds one spring per pair of cysteine sulfurs within `cutoff`.
//
// A bridge that already carries a spring -- a PDB CONECT record makes one, at
// the global -stiffness -- is RECALIBRATED rather than skipped: it is the same
// bond, and leaving it at a mesh stiffness of several thousand would model the
// bridge far too rigidly without saying so.
//
// The sulfur is found by name -- a particle whose name ends in "SG" in a
// residue whose name starts with "CY" -- which covers both the raw PDB name
// and the type a .grp gives it (amber.grp names it CSG). Element names would
// be the honest test, but no reader in BioSpring populates them.
//
// Returns the number of springs added.
size_t addDisulfideSprings(topology::Topology & topology, float cutoff, float stiffness);

} // namespace staticbond
} // namespace biospring

#endif // __STATIC_BOND_BUILDER_H__
