#ifndef __SPRING_ENERGY_HPP__
#define __SPRING_ENERGY_HPP__

#include "../constants.hpp"
#include "../shared/spring_shared.h"

namespace biospring
{
namespace forcefield
{

// The arithmetic lives in ../shared/spring_shared.h, which the OpenCL kernel
// source is built from as well. What stays here is this project's C++ naming
// and the unit convention; the formula itself is not restated.

/// @param distance Distance between the two particles, in Angstrom (A).
/// @param stiffness Spring constant, in kJ.mol-1.A-2 (molar convention;
///     typical literature value for Calpha elastic network models is
///     ~0.6 kcal.mol-1.A-2, i.e. ~2.5 kJ.mol-1.A-2).
/// @param equilibrium Rest length of the spring, in Angstrom (A).
/// @return Spring energy, in kJ.mol-1.
inline float spring_energy(float distance, float stiffness, float equilibrium)
{
    return biospring_spring_energy(distance, stiffness, equilibrium);
}

/// @return Spring force module, in Da.A.fs-2 (see GLOBAL_SPRING_FORCE_CONVERT).
inline float spring_force_module(float distance, float stiffness, float equilibrium)
{
    return biospring_spring_force_module(distance, stiffness, equilibrium,
                                         static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT));
}

} // namespace forcefield
} // namespace biospring

#endif // __SPRING_ENERGY_HPP__
