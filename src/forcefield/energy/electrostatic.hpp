#ifndef __ELECTROSTATIC_ENERGY_HPP__
#define __ELECTROSTATIC_ENERGY_HPP__

#include "../constants.hpp"

namespace biospring
{
namespace forcefield
{

inline float electrostatic_energy(float charge1, float charge2, float distance, float dielectric)
{
    if (distance < MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF)
        return 0.0;

    // Coulomb's contant.
    double k = COULOMB_CONSTANT;

    // Convert charges to Coulomb.
    double q1 = charge1 * ELECTRONCHARGE_TO_COULOMB;
    double q2 = charge2 * ELECTRONCHARGE_TO_COULOMB;

    // Convert distance to meter.
    distance *= ANGSTROM_TO_METER;

    double energy = k * (q1 * q2) / (dielectric * distance);
    energy = energy * AVOGADRO_NUMBER; // in J.mol-1
    energy = energy * 1.0E-3;          // in kJ.mol-1
    return static_cast<float>(energy);
}

inline float electrostatic_force_module(float charge1, float charge2, float distance, float dielectric)
{
    if (distance < MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF)
        return 0.0;
    float force_module = -(charge1 * charge2) / (4.0 * PI * dielectric * distance * distance);
    force_module *= GLOBAL_ELECTROSTATIC_FORCE_CONVERT; // en N ou en kg.m.s-2;
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __ELECTROSTATIC_ENERGY_HPP__