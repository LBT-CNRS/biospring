#ifndef __STERIC_ENERGY_HPP__
#define __STERIC_ENERGY_HPP__

#include "../CombinationRules.hpp"
#include "../constants.hpp"
#include "../shared/steric_shared.h"

#include <cmath>

namespace biospring
{
namespace forcefield
{

// ======================================================================================
// Linear steric potential.
//
// radius_i, radius_j: particle radii, in Angstrom (A).
// distance: distance between the two particles, in Angstrom (A).

// Penalty stiffness for particle overlap, in kJ.mol-1.A-2 (molar
// convention, same unit as spring stiffness). Shared by energy and force so
// the reported energy is consistent with the force actually integrated
// (they used to diverge: 100 vs 1.0).
static const float STERIC_LINEAR_STIFFNESS = 1.0;

/// @return Steric energy, in kJ.mol-1 (0 when particles do not overlap).
inline float steric_energy_linear(float radius_i, float radius_j, float distance)
{
    float equilibrium = radius_i + radius_j;
    float distancevar = (distance - equilibrium);

    if (distancevar > 0)
        return 0.0;

    return 0.5 * STERIC_LINEAR_STIFFNESS * distancevar * distancevar;
}

/// @return Steric force module, in Da.A.fs-2 (see GLOBAL_SPRING_FORCE_CONVERT).
/// The arithmetic lives in ../shared/steric_shared.h, which the OpenCL kernel
/// compiles too. The constants stay here and are passed down.
inline float steric_force_module_linear(float radius_i, float radius_j, float distance)
{
    return biospring_steric_force_module_linear(radius_i, radius_j, distance, STERIC_LINEAR_STIFFNESS,
                                                static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT));
}

// ======================================================================================
// Amber 12-6 Lennard-Jones potential.
//
// radius_i, radius_j: particle radii (sigma), in Angstrom (A).
// epsilon_i, epsilon_j: particle well depths, in kJ.mol-1.
// distance: distance between the two particles, in Angstrom (A).
// Energies below are in kJ.mol-1, force modules in Da.A.fs-2 (converted via
// GLOBAL_SPRING_FORCE_CONVERT, see constants.hpp).

inline float steric_energy_amber(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    if (distance < MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon_ij = combination_rules::lorentz_berthelot::epsilon(epsilon_i, epsilon_j);
    float radius_ij = combination_rules::good_hope::radius(radius_i, radius_j);

    float repulsive = epsilon_ij * pow(radius_ij / distance, 12.0f);
    float attractive = -epsilon_ij * 2.0 * pow(radius_ij / distance, 6.0f);

    return repulsive + attractive;
}

inline float steric_force_module_amber(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    return biospring_steric_force_module_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                               static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF),
                                               static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT));
}

// ======================================================================================
// Lewitt 8-6 Lennard-Jones potential.
// Same units as the Amber 12-6 potential above (radius/distance in A,
// epsilon in kJ.mol-1, energy in kJ.mol-1, force module in Da.A.fs-2).

inline float steric_energy_lewitt(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    if (distance < MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon_ij = combination_rules::lorentz_berthelot::epsilon(epsilon_i, epsilon_j);
    float radius_ij = combination_rules::good_hope::radius(radius_i, radius_j);

    float repulsive = epsilon_ij * 3.0 * pow(radius_ij / distance, 8.0f);
    float attractive = -epsilon_ij * 4.0 * pow(radius_ij / distance, 6.0f);

    return repulsive + attractive;
}

inline float steric_force_module_lewitt(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    return biospring_steric_force_module_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF),
                                                static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT));
}

// ======================================================================================
// Zacharias 8-6 Lennard-Jones potential.
// Same units as the Amber 12-6 potential above (radius/distance in A,
// epsilon in kJ.mol-1, energy in kJ.mol-1, force module in Da.A.fs-2).

inline float steric_energy_zacharias(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    if (distance < MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon_ij = combination_rules::zacharias::epsilon(epsilon_i, epsilon_j);
    float radius_ij = combination_rules::zacharias::radius(radius_i, radius_j);

    float repulsive = epsilon_ij * pow(radius_ij / distance, 8.0f);
    float attractive = -epsilon_ij * pow(radius_ij / distance, 6.0f);

    return repulsive + attractive;
}

inline float steric_force_module_zacharias(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    return biospring_steric_force_module_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                   static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF),
                                                   static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT));
}

} // namespace forcefield
} // namespace biospring

#endif // __STERIC_ENERGY_HPP__