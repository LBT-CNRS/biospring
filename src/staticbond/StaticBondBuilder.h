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

// Hydrogen bonds and disulfide bridges as ORDINARY SPRINGS, taking their force
// constant from chemistry.
//
// WHICH BONDS EXIST IS THE USER'S INPUT, NOT THIS CODE'S GUESS. The bonds come
// from the topology -- a PDB's CONECT records, anything pdb2spn already turned
// into a spring -- and all this does is give each declared bond the force
// constant its chemistry calls for, instead of the uniform -stiffness that
// pdb2spn applies to every CONECT alike. Nothing is searched for, nothing is
// invented, and a structure that declares no bond gets no spring.
//
// An earlier version scanned every donor/acceptor pair within 3.5 A instead.
// It found bonds nobody asked for: on 072's kinase, 122 of 263 sprung pairs
// sat between 3.2 and 3.5 A with a median angular weight of 0.01 -- distance
// coincidences, not hydrogen bonds -- and they dragged the median equilibrium
// length to 3.15 A against a real bond's 2.90. No cutoff or angular filter
// fixes that honestly, because "is this a bond" is the user's question to
// answer, not a threshold's.
//
// What a spring cannot do, and it should be said plainly: it never breaks, and
// it never forms. These options hold a declared structure together; they are
// not for studying association or dissociation.

// Force constant of a hydrogen-bond spring, in kJ.mol-1.A-2.
//
// Matched to the curvature of a Morse well at its minimum:
// V = De[1 - exp(-a(r-re))]^2 has V''(re) = 2.De.a^2, and BioSpring's spring
// energy is 0.5.k.dr^2, so k = 2.De.a^2 with De = 16.7 kJ.mol-1 and
// a = 1.34 A-1. The well DEPTH has no harmonic equivalent -- a spring is
// infinitely deep -- so the curvature is what can be matched, and it is what
// sets the restoring force for the small displacements a held structure makes.
// Measured against the Morse: +11 % at 0.05 A, +22 % at 0.1 A, +49 % at 0.2 A.
// Stiffer rather than softer, which is the right direction for something whose
// job is to hold.
extern const float HYDROGEN_BOND_STIFFNESS;

// Force constant of a disulfide spring, in kJ.mol-1.A-2. AMBER ff14SB's S-S
// bond is r0 = 2.038 A with k = 1389.1 in the 0.5.k.dr^2 convention, which is
// BioSpring's own, so the value transfers with no conversion.
extern const float DISULFIDE_STIFFNESS;

// What a .hbond rule says about one particle: how many hydrogens it can donate
// and how many lone pairs it can accept. That is all a declared bond needs --
// enough to tell a hydrogen bond from the covalent links a CONECT also carries.
struct DonorAcceptorRole
{
    int donorCapacity = 0;
    int acceptorCapacity = 0;
};

// (residue name, particle name) -> role.
using DonorAcceptorTable = std::map<std::pair<std::string, std::string>, DonorAcceptorRole>;

// Reads a .hbond table, keeping only the two capacities. The antecedent
// columns orient a dynamic bond and mean nothing to a spring, which pulls
// along the line between two atoms; they are parsed past so one table serves
// both mechanisms.
//
// Throws std::runtime_error if the file cannot be opened.
DonorAcceptorTable readDonorAcceptorTable(const std::string & path);

// Retunes every DECLARED bond whose two ends are a donor and an acceptor to
// the hydrogen-bond force constant. Equilibrium is left untouched: it is the
// declared bond's own length, so retuning cannot deform the input structure.
//
// Returns the number of springs retuned.
size_t retuneHydrogenBondSprings(topology::Topology & topology, const DonorAcceptorTable & table, float stiffness);

// Retunes every DECLARED bond between two cysteine sulfurs to the disulfide
// force constant. Without this a bridge declared by CONECT carries the global
// -stiffness -- 8000 in the rigid-body examples, 5.8x too rigid -- silently.
//
// The sulfur is found by name: a particle whose name ends in "SG" in a residue
// whose name starts with "CY", which covers both the raw PDB name and the type
// a .grp gives it (amber.grp names it CSG). Element names would be the honest
// test, but no reader in BioSpring populates them.
//
// Returns the number of springs retuned.
size_t retuneDisulfideSprings(topology::Topology & topology, float stiffness);

} // namespace staticbond
} // namespace biospring

#endif // __STATIC_BOND_BUILDER_H__
