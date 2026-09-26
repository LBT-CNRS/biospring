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

// How each law COMBINES TWO RADII by default: the rule of the parameter set the
// law was published with, because that is the file it will be handed unless the
// user says otherwise. steric.radiusrule in the .msp overrides any of them; the
// table of what every shipped .ff actually stores is in ../shared/steric_shared.h.
//
//   AMBER 12-6      amber*.ff hold R* = rmin/2, so two of them ADD
//   Levitt 8-6      CAonlyLewitt.ff holds a CA-CA contact distance already, so
//                   two of them are MEANED -- summing gives a 10-13 A minimum
//                   where a real protein's closest non-bonded CA pair is 3.7 A
//   Zacharias 8-6   its own product rule, on epsilon as well as the radius
//   linear          two radii touch: a sum, and it always was one
static const int STERIC_RADIUS_RULE_AMBER = BIOSPRING_RADIUS_SUM;
static const int STERIC_RADIUS_RULE_LEWITT = BIOSPRING_RADIUS_GEOMETRIC_MEAN;
static const int STERIC_RADIUS_RULE_ZACHARIAS = BIOSPRING_RADIUS_PRODUCT;
static const int STERIC_RADIUS_RULE_LINEAR = BIOSPRING_RADIUS_SUM;

/// @return Steric energy, in kJ.mol-1 (0 when particles do not overlap).
inline float steric_energy_linear(float radius_i, float radius_j, float distance)
{
    return biospring_steric_energy_linear(radius_i, radius_j, distance, STERIC_LINEAR_STIFFNESS);
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

inline float steric_energy_amber(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                 float distance, int radiusrule = STERIC_RADIUS_RULE_AMBER)
{
    return biospring_steric_energy_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                        static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF), radiusrule);
}

inline float steric_force_module_amber(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                       float distance, int radiusrule = STERIC_RADIUS_RULE_AMBER)
{
    return biospring_steric_force_module_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                               static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF),
                                               static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT), radiusrule);
}

// ======================================================================================
// Lewitt 8-6 Lennard-Jones potential.
// Same units as the Amber 12-6 potential above (radius/distance in A,
// epsilon in kJ.mol-1, energy in kJ.mol-1, force module in Da.A.fs-2).

inline float steric_energy_lewitt(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                  float distance, int radiusrule = STERIC_RADIUS_RULE_LEWITT)
{
    return biospring_steric_energy_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                         static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF), radiusrule);
}

inline float steric_force_module_lewitt(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                        float distance, int radiusrule = STERIC_RADIUS_RULE_LEWITT)
{
    return biospring_steric_force_module_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF),
                                                static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT), radiusrule);
}

// ======================================================================================
// Zacharias 8-6 Lennard-Jones potential.
// Same units as the Amber 12-6 potential above (radius/distance in A,
// epsilon in kJ.mol-1, energy in kJ.mol-1, force module in Da.A.fs-2).

inline float steric_energy_zacharias(float radius_i, float radius_j, float epsilon_i, float epsilon_j, float distance)
{
    return biospring_steric_energy_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                            static_cast<float>(MINIMAL_DISTANCE_VDW_CUTOFF));
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