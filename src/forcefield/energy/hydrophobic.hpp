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
    double force_module = 0.0;
    force_module = (hydrophobicity1 * hydrophobicity2) * exp(-distance);
    force_module = force_module * AVOGADRO_NUMBER * 1.0E-3; // J/mol
    force_module = force_module * 1.0E-3;                   // kJ/mol
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __HYDROPHOBIC_ENERGY_HPP__