#ifndef __STACKING_ENERGY_HPP__
#define __STACKING_ENERGY_HPP__

#include <cmath>

#include "../constants.hpp"
#include "../shared/stacking_shared.h"

namespace biospring
{
namespace forcefield
{

// Base stacking after HiRE-RNA v3. The functional form and why it cannot
// double count the hydrogen bond are documented once, in
// forcefield/shared/stacking_shared.h, which this file calls straight into.
//
// WHERE THE DEFAULTS COME FROM. Measured on 8PKL, the wild-type Escherichia
// coli 70S this project already uses, over its 4757 complete nucleobases:
// taking the 2317 chain-consecutive pairs whose angular product is above 0.7 --
// that is, the ones genuinely stacked rather than merely adjacent -- the
// frame-centre distance is 3.92 A median (IQR 3.72-4.12) with 2*var = 0.85 A^2.
// The same measurement is what separates a stack from everything else: over
// every base pair within 9 A, the angular product is 0.752 median along the
// chain against 0.313 for all other pairs, and the discrimination comes almost
// entirely from the two alignment factors (0.911/0.894 against 0.679/0.662)
// rather than from (ni.nj)^2, which barely moves (0.965 against 0.924) because
// in a folded RNA most nearby bases are roughly parallel whether they stack or
// not.
//
// The well depth is NOT measured here. Stacking free energies run about 4 to 12
// kJ.mol-1 per step experimentally, and HiRE-RNA fits its own eps against decoy
// ensembles; 8.0 kJ.mol-1 is a starting value in that range and it is a setting,
// not a constant.

/// e^(-(r-r_st)^2/sigma), the distance factor.
/// @param distance    Base frame centre-to-centre distance, in A.
/// @param equilibrium r_st, in A.
/// @param width       sigma, in A^2.
inline float stacking_distance_factor(float distance, float equilibrium, float width)
{
    return biospring_stacking_distance_factor(distance, equilibrium, width);
}

/// d(distance factor)/d(distance), in A-1.
inline float stacking_distance_derivative(float distance, float equilibrium, float width)
{
    return biospring_stacking_distance_derivative(distance, equilibrium, width);
}

/// A(u) = u^2 (2 - u^2), u the cosine between a base normal and the axis
/// joining the two bases: 1 face to face, 0 edge to edge.
inline float stacking_alignment_factor(float u) { return biospring_stacking_alignment_factor(u); }

/// dA/du = 4u(1 - u^2).
inline float stacking_alignment_derivative(float u) { return biospring_stacking_alignment_derivative(u); }

} // namespace forcefield
} // namespace biospring

#endif // __STACKING_ENERGY_HPP__
