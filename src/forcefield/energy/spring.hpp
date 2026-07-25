#ifndef __SPRING_ENERGY_HPP__
#define __SPRING_ENERGY_HPP__

#include "../constants.hpp"

namespace biospring
{
namespace forcefield
{

/// @param distance Distance between the two particles, in Angstrom (A).
/// @param stiffness Spring constant, in kJ.mol-1.A-2 (molar convention;
///     typical literature value for Calpha elastic network models is
///     ~0.6 kcal.mol-1.A-2, i.e. ~2.5 kJ.mol-1.A-2).
/// @param equilibrium Rest length of the spring, in Angstrom (A).
/// @return Spring energy, in kJ.mol-1.
inline float spring_energy(float distance, float stiffness, float equilibrium)
{
    float distancevar = (distance - equilibrium);
    return 0.5 * stiffness * distancevar * distancevar;
}

/// @return Spring force module, in Da.A.fs-2 (see GLOBAL_SPRING_FORCE_CONVERT).
inline float spring_force_module(float distance, float stiffness, float equilibrium)
{
    float force_module = stiffness * (distance - equilibrium);
    force_module *= GLOBAL_SPRING_FORCE_CONVERT;
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __SPRING_ENERGY_HPP__