#ifndef __HYDROGENBOND_ENERGY_HPP__
#define __HYDROGENBOND_ENERGY_HPP__

#include <cmath>

#include "../constants.hpp"

namespace biospring
{
namespace forcefield
{

// Morse potential between a donor and an acceptor heavy atom (no explicit
// hydrogen: distance only, no donor-H...acceptor angle). Parameters derived
// from the hydrogen-bond stretching mode reported in vibrational
// spectroscopy literature (~100-200 cm-1 for a moderate N-H...O=C bond,
// reduced mass of a N...O pair), not fit ad hoc:
//   De ~ 16.7 kJ.mol-1 (4 kcal.mol-1, typical backbone H-bond well depth)
//   re ~ 2.9 A (standard N...O distance for a backbone H-bond)
//   a  ~ 1.34 A-1 (from k = 2*De*a^2 with k ~ 60 kJ.mol-1.A-2, itself from
//        the ~150 cm-1 stretching frequency above)
// giving a practical range (<1% of De left) of ~7 A -- this sets the
// --hbond grid's cell size, not the long-range cutoff used for Coulomb.

/// @param distance Donor-acceptor heavy-atom distance, in Angstrom (A).
/// @param wellDepth Morse well depth De, in kJ.mol-1 (already molar).
/// @param equilibrium Equilibrium distance re, in Angstrom (A).
/// @param width Morse width parameter a, in A-1.
/// @return Hydrogen-bond energy, in kJ.mol-1 (0 at distance -> infinity, -De
///     at distance == equilibrium).
inline float hydrogen_bond_energy(float distance, float wellDepth, float equilibrium, float width)
{
    const float u = std::exp(-width * (distance - equilibrium));
    const float one_minus_u = 1.0f - u;
    return wellDepth * one_minus_u * one_minus_u - wellDepth;
}

/// @return dV/dr, i.e. the force module in the "self -> neighbor" convention
///     used throughout (positive = attractive), in kJ.mol-1.A-1. Multiply by
///     GLOBAL_SPRING_FORCE_CONVERT to use as a mechanical force module: the
///     well depth is already molar (kJ.mol-1), same convention as spring
///     stiffness and steric epsilon.
inline float hydrogen_bond_force_module(float distance, float wellDepth, float equilibrium, float width)
{
    const float u = std::exp(-width * (distance - equilibrium));
    float force_module = 2.0f * width * wellDepth * u * (1.0f - u);
    force_module *= GLOBAL_SPRING_FORCE_CONVERT;
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __HYDROGENBOND_ENERGY_HPP__
