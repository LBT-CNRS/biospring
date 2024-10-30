#ifndef __STERIC_ENERGY_HPP__
#define __STERIC_ENERGY_HPP__

#include "../CombinationRules.hpp"
#include "../constants.hpp"

#include <cmath>

namespace biospring
{
namespace forcefield
{

// ======================================================================================
// Linear steric potential.

inline float steric_energy_linear(float radius_i, float radius_j, float distance)
{
    float equilibrium = radius_i + radius_j;
    float distancevar = (distance - equilibrium);

    if (distancevar > 0)
        return 0.0;

    float stiffness = 100;
    return 0.5 * (stiffness)*distancevar * distancevar;
}

inline float steric_force_module_linear(float radius_i, float radius_j, float distance)
{
    float equilibrium = radius_i + radius_j;
    float distancevar = (distance - equilibrium);

    if (distancevar > 0)
        return 0.0;

    float stiffness = 1.0;

    float force_module = -stiffness * fabs(distancevar);
    return force_module * GLOBAL_SPRING_FORCE_CONVERT;
}

// ======================================================================================
// Amber 12-6 Lennard-Jones potential.

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
    if (distance < MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon_ij = combination_rules::lorentz_berthelot::epsilon(epsilon_i, epsilon_j);
    float radius_ij = combination_rules::good_hope::radius(radius_i, radius_j);

    float repulsive = -epsilon_ij * 12.0 * (pow(radius_ij, 12.0f) / pow(distance, 13.0f));
    float attractive = epsilon_ij * 2.0 * 6.0 * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;
    return force_module * GLOBAL_SPRING_FORCE_CONVERT;
}

// ======================================================================================
// Lewitt 8-6 Lennard-Jones potential.

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
    if (distance < MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon_ij = combination_rules::lorentz_berthelot::epsilon(epsilon_i, epsilon_j);
    float radius_ij = combination_rules::good_hope::radius(radius_i, radius_j);

    float repulsive = -epsilon_ij * 3.0 * 8.0 * (pow(radius_ij, 8.0f) / pow(distance, 9.0f));
    float attractive = epsilon_ij * 4.0 * 6.0 * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;

    return force_module * GLOBAL_SPRING_FORCE_CONVERT;
}

// ======================================================================================
// Zacharias 8-6 Lennard-Jones potential.

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
    if (distance < MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon_ij = combination_rules::zacharias::epsilon(epsilon_i, epsilon_j);
    float radius_ij = combination_rules::zacharias::radius(radius_i, radius_j);

    float repulsive = -epsilon_ij * 8.0 * (pow(radius_ij, 8.0f) / pow(distance, 9.0f));
    float attractive = epsilon_ij * 6.0 * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;

    return force_module * GLOBAL_SPRING_FORCE_CONVERT;
}

} // namespace forcefield
} // namespace biospring

#endif // __STERIC_ENERGY_HPP__