#ifndef __SPRING_ENERGY_HPP__
#define __SPRING_ENERGY_HPP__

#include "../constants.hpp"

namespace biospring
{
namespace forcefield
{

inline float spring_energy(float distance, float stiffness, float equilibrium)
{
    float distancevar = (distance - equilibrium);
    return 0.5 * stiffness * distancevar * distancevar;
}

inline float spring_force_module(float distance, float stiffness, float equilibrium)
{
    float force_module = stiffness * (distance - equilibrium);
    force_module *= GLOBAL_SPRING_FORCE_CONVERT;
    return force_module;
}

} // namespace forcefield
} // namespace biospring

#endif // __SPRING_ENERGY_HPP__