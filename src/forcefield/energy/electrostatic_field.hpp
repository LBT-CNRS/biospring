#ifndef __ELECTROSTATIC_FIELD_ENERGY_HPP__
#define __ELECTROSTATIC_FIELD_ENERGY_HPP__

#include "../constants.hpp"

namespace biospring
{
namespace forcefield
{

/// @param potential Electrostatic potential grid value at the particle's
///     position, in K/e (Kelvin per elementary charge; see BOLTZMANJPERK
///     below to convert to real energy).
/// @param charge Particle charge, in elementary charge units (e).
/// @return Electrostatic field energy, in kJ.mol-1.
inline float electrostatic_field_energy(float potential, float charge)
{
    float energy = potential * charge; // k.K.e-1.e ou en k.K
    energy = energy * BOLTZMANJPERK;   // k.K soit en J.K-1.K soit en J
    energy = energy * JOULE_TO_KJOULE; // kJ
    energy = energy * AVOGADRO_NUMBER; // kJ.mol-1
    return energy;
}

} // namespace forcefield
} // namespace biospring

#endif // __ELECTROSTATIC_FIELD_ENERGY_HPP__