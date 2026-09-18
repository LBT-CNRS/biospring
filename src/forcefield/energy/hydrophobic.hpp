#ifndef __HYDROPHOBIC_ENERGY_HPP__
#define __HYDROPHOBIC_ENERGY_HPP__

#include "../constants.hpp"
#include "../shared/hydrophobic_shared.h"
#include <cmath>

namespace biospring
{
namespace forcefield
{

/// @param hydrophobicity1, hydrophobicity2 Per-particle hydrophobicity, from
///     the force field's seventh column (the .nc variable "hydrophobicity",
///     NOT "hydrophobicityscale", which is IMPALA's transfer energy). Their
///     product is in kJ.mol-1.A-1.
/// @param distance Distance between the two particles, in Angstrom (A).
/// @param decaylength The decay length of the attraction, in A.
/// @return Hydrophobic energy of the pair, in kJ.mol-1.
///
/// The arithmetic lives in ../shared/hydrophobic_shared.h next to the force,
/// which the OpenCL kernel compiles too, because the two are one function:
/// this is the integral of that gradient. They used to be written apart and
/// disagreed by N_A * 1e-3 -- the energy scaled as if the product were joules
/// per molecule while the force scaled as if it were already molar.
inline float hydrophobic_energy(float hydrophobicity1, float hydrophobicity2, float distance,
                                float decaylength)
{
    return biospring_hydrophobic_energy(hydrophobicity1, hydrophobicity2, distance, decaylength);
}

/// @return Hydrophobic force module, in Da.A.fs-2 (see
///     GLOBAL_SPRING_FORCE_CONVERT). The constant stays here and is passed
///     down, so the shared header carries no include.
inline float hydrophobic_force_module(float hydrophobicity1, float hydrophobicity2, float distance,
                                      float decaylength)
{
    return biospring_hydrophobic_force_module(hydrophobicity1, hydrophobicity2, distance, decaylength,
                                              static_cast<float>(GLOBAL_SPRING_FORCE_CONVERT));
}

} // namespace forcefield
} // namespace biospring

#endif // __HYDROPHOBIC_ENERGY_HPP__
