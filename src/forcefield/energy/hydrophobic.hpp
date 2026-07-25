#ifndef __HYDROPHOBIC_ENERGY_HPP__
#define __HYDROPHOBIC_ENERGY_HPP__

#include "../constants.hpp"
#include <cmath>

namespace biospring
{
namespace forcefield
{

/// @param hydrophobicity1, hydrophobicity2 Particle hydrophobicity/transfer
///     scale, in kJ.mol-1 (see NetCDFWriter's "hydrophobicityscale" units).
/// @param distance Distance between the two particles, in Angstrom (A).
/// @return Pseudo-hydrophobicity energy, in kJ.mol-1. Note: this is an
///     empirical (not first-principles) potential, so hydrophobicity1 *
///     hydrophobicity2 is not literally an energy despite each factor being
///     in kJ.mol-1; the Avogadro/kJ scaling below only mirrors the reporting
///     convention used by the other energy terms in this module.
inline float hydrophobic_energy(float hydrophobicity1, float hydrophobicity2, float distance)
{
    double energy = 0.0;
    energy = -(hydrophobicity1 * hydrophobicity2) * exp(-distance);
    energy = energy * AVOGADRO_NUMBER; // J/mol
    energy = energy * 1.0E-3;          // kJ/mol
    return energy;
}

/// @return Pseudo-hydrophobicity force module, in Da.A.fs-2 (see
///     GLOBAL_SPRING_FORCE_CONVERT). hydrophobicity1/2 are treated as an
///     already-molar (kJ.mol-1) quantity here, same convention as spring
///     stiffness, so no separate Avogadro scaling is applied (unlike
///     hydrophobic_energy above).
inline float hydrophobic_force_module(float hydrophobicity1, float hydrophobicity2, float distance)
{
    float force_module = (hydrophobicity1 * hydrophobicity2) * exp(-distance);
    force_module *= GLOBAL_SPRING_FORCE_CONVERT;
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __HYDROPHOBIC_ENERGY_HPP__