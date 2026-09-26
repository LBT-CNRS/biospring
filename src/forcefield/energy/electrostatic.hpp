#ifndef __ELECTROSTATIC_ENERGY_HPP__
#define __ELECTROSTATIC_ENERGY_HPP__

#include "../constants.hpp"
#include "../shared/electrostatic_shared.h"

namespace biospring
{
namespace forcefield
{

/// @param charge1, charge2 Particle charges, in elementary charge units (e).
/// @param distance Distance between the two particles, in Angstrom (A).
/// @param dielectric Relative dielectric constant (dimensionless).
/// @return Coulomb energy, in kJ.mol-1 (real per-pair energy scaled to the
///     molar convention via AVOGADRO_NUMBER, see constants.hpp).
inline float electrostatic_energy(float charge1, float charge2, float distance, float dielectric)
{
    return biospring_electrostatic_energy(charge1, charge2, distance, dielectric,
                                         static_cast<float>(MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF),
                                         static_cast<float>(GLOBAL_ELECTROSTATIC_ENERGY_CONVERT));
}

/// @return Coulomb force module, in Da.A.fs-2 (see
///     GLOBAL_ELECTROSTATIC_FORCE_CONVERT: charges here are real per-particle
///     values, so unlike electrostatic_energy no Avogadro scaling applies).
///
/// The arithmetic lives in ../shared/electrostatic_shared.h, which the OpenCL
/// kernel compiles too, so the two backends cannot drift apart on it. The
/// constants stay here and are passed down: there is one definition of each.
inline float electrostatic_force_module(float charge1, float charge2, float distance, float dielectric)
{
    // charge1, charge2 in e; distance in A (not converted to meter here,
    // GLOBAL_ELECTROSTATIC_FORCE_CONVERT accounts for both units at once).
    return biospring_electrostatic_force_module(charge1, charge2, distance, dielectric,
                                                static_cast<float>(MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF),
                                                static_cast<float>(4.0 * PI),
                                                static_cast<float>(GLOBAL_ELECTROSTATIC_FORCE_CONVERT));
}

} // namespace forcefield
} // namespace biospring

#endif // __ELECTROSTATIC_ENERGY_HPP__