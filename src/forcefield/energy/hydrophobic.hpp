#ifndef __HYDROPHOBIC_ENERGY_HPP__
#define __HYDROPHOBIC_ENERGY_HPP__

#include "../constants.hpp"
#include <cmath>

namespace biospring
{
namespace forcefield
{

inline float hydrophobic_energy(float hydrophobicity1, float hydrophobicity2, float distance)
{
    double energy = 0.0;
    energy = -(hydrophobicity1 * hydrophobicity2) * exp(-distance);
    energy = energy * AVOGADRO_NUMBER; // J/mol
    energy = energy * 1.0E-3;          // kJ/mol
    return energy;
}

inline float hydrophobic_force_module(float hydrophobicity1, float hydrophobicity2, float distance)
{
    float force_module = (hydrophobicity1 * hydrophobicity2) * exp(-distance);
    force_module *= GLOBAL_SPRING_FORCE_CONVERT;
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __HYDROPHOBIC_ENERGY_HPP__