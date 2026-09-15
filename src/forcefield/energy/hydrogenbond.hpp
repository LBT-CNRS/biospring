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

// ---- Angular weight -----------------------------------------------------
//
// A hydrogen bond is directional, and a distance-only Morse is not: on a
// B-DNA duplex it happily pairs two stacked same-strand groups sitting 2.75 A
// apart, closer to each other than to their real Watson-Crick partners at
// 2.90 and 2.94 A. Without an explicit hydrogen the direction still exists,
// carried by the donor's ANTECEDENT: the heavy atom it hangs off. The
// antecedent->donor vector stands in for where the H (or the lone pair)
// points, and how far the acceptor sits off that direction is what weights
// the well.
//
// w = max(cos(theta), 0)^2, theta the angle at the donor between
// (donor - antecedent) and (acceptor - donor). Measured over 48 real
// Watson-Crick bonds it runs 0.036 to 0.48, median 0.25; the stacked pair
// above gets 0.003. It does not need to separate them on its own -- donor
// and acceptor capacities let the parasitic pair form alongside the real one
// instead of starving it -- it only has to make it weigh nothing, and 0.003
// against 0.25 does.
//
// Clamped at zero rather than squared through: an acceptor BEHIND the donor
// (theta > 90 deg) is not a weak hydrogen bond, it is not one at all, and
// cos^2 alone would revive it.

/// @param cos_theta Cosine of the angle at the donor.
/// @return Dimensionless weight in [0, 1].
inline float hydrogen_bond_angular_factor(float cos_theta)
{
    return cos_theta > 0.0f ? cos_theta * cos_theta : 0.0f;
}

/// d(weight)/d(cos theta), for the force.
inline float hydrogen_bond_angular_derivative(float cos_theta)
{
    return cos_theta > 0.0f ? 2.0f * cos_theta : 0.0f;
}

} // namespace forcefield
} // namespace biospring

#endif // __HYDROGENBOND_ENERGY_HPP__
